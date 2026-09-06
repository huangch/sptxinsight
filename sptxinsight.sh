#!/usr/bin/env bash
# sptxinsight.sh - unified environment-aware runner for the sptxinsight CLI.
#
# Manages BOTH backends (native + docker) without depending on the legacy
# sptxinsight-docker-run.sh script, which has moved to bak_old_scripts/ and
# is kept only for reference.
#
# Subcommands:
#   run        --runner {native,docker} -d DIR [--gpu ID|all] [--tmpdir DIR]
#               [--no-pull] [--dry-run] [SPTXINSIGHT_ARGS ...]
#   status                                                  # current effective config
#   doctor     [--runner {native,docker}]                   # preflight
#   where                                                   # absolute path of this script
#   -h | --help
#   --version
#
#   `-d DIR` (a.k.a. `--data-dir`) is REQUIRED when --runner docker is in
#   effect: it names the host directory that gets bind-mounted to /workspace
#   inside the container. SPTXINSIGHT_DATA_DIR supplies the same value
#   non-interactively.
#
# Why this script uses --runner instead of -b/--backend (mirrors the
# same principle from wsinsight.sh: the wrapper's runner-selector flag must
# be in a DIFFERENT vocabulary from the inner CLI's own --backend):
#   sptxinsight's CLI itself has a global --backend flag (anndata, zarr,
#   spatialdata) that names the I/O backend. The wrapper's job is to select
#   where the CLI RUNS (native conda env vs docker container); to avoid
#   overloading sptxinsight's --backend and to keep the two flag spaces
#   disjoint, the wrapper uses --runner for native|docker. NOTE that
#   wsinsight.sh (the same-family wrapper) DOES use '-b' / '--backend' for
#   its own runner selector, because the wsinsight CLI does NOT have a
#   global --backend flag. The two wrappers diverge here ON PURPOSE:
#   never copy '-b' into this file without first removing --backend from
#   the sptxinsight CLI itself.
#
# Param-parsing rule:
#   Everything before the first sptxinsight subcommand name (run, ingest,
#   verify, niche, ...) is consumed by THIS script (env control: --runner,
#   -d, --gpu, --tmpdir, --no-pull, --dry-run). From (and including) the first
#   sptxinsight subcommand name, every token is passed through verbatim.
#
#   Boundary discovery:
#     1. Absorb the script's own flags and subcommands (status / doctor /
#        where take priority over sptxinsight subcommands of the same name).
#     2. Honor an explicit `--` delimiter.
#     3. The first position arg is checked against the list of known
#        sptxinsight subcommand names; if it matches, passthrough begins.
#     4. Backward-scan fallback: if the first position arg is not a known
#        subcommand, scan the remainder of argv for the first known one.
#     5. If no known sptxinsight subcommand is found anywhere, die with the
#        list of known commands.
#
#   Unknown script flags (-X that the wrapper does not recognize): in default
#   (lenient) mode, warn and treat as sptxinsight's; in SPTXINSIGHT_STRICT=1
#   mode, die.
#
# Exit code 0 on success, 1 on user error (bad args, no sptxinsight cmd),
# 2 on infrastructure failure (docker daemon down, sptxinsight not on PATH
# AND the schema --commands-only fetch failed).

set -euo pipefail

PROG="$(basename "$0")"
VERSION="1.0.0"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

IMAGE_ID="${SPTXINSIGHT_IMAGE:-huangchtw/sptxinsight:latest}"
HF_CACHE_VOLUME="${SPTXINSIGHT_HF_CACHE_VOLUME:-sptxinsight-hf-cache}"
COMMANDS_CACHE="$HOME/.cache/sptxinsight/commands.txt"
COMMANDS_CACHE_TTL_SECONDS="${SPTXINSIGHT_COMMANDS_TTL_SECONDS:-86400}"

# `-d` / `--data-dir` is REQUIRED when --runner docker is in effect.  It names
# the host directory that gets bind-mounted to /workspace inside the container;
# mirrors wsinsight.sh and clawpyter.sh. SPTXINSIGHT_DATA_DIR is the matching
# env var (kept for CI / non-interactive use; the flag takes precedence).
SCRIPT_DATA_DIR="${SPTXINSIGHT_DATA_DIR:-}"

# LAST-RESORT builtin subcommand list - updated as sptxinsight evolves.
# Used only when neither the cache nor a live `sptxinsight schema --commands-only`
# works. The bundled schema in sptxinsight/cli is the source of truth.
_SPT_BUILTIN_CMDS=(
    run ingest verify export schema niche niche-profile
    annotate marker-init marker-rerank
    hplot hplot-finalize cci agg
)

# ---------------------------------------------------------------------------
# usage
# ---------------------------------------------------------------------------
print_usage() {
    cat <<EOF
$PROG $VERSION - run sptxinsight with one of two backends (native | docker)

Usage:
  $PROG run     --runner {native,docker} -d DIR | --data-dir DIR
                [--gpu ID|all] [--tmpdir DIR] [--no-pull] [--dry-run] [SPTXINSIGHT_ARGS ...]
  $PROG status
  $PROG doctor  [--runner {native,docker}]
  $PROG where
  $PROG -h | --help
  $PROG --version

Backends (selected via --runner):
  native    Invoke the sptxinsight CLI on the host inside the activated env.
  docker    Run sptxinsight inside the $IMAGE_ID container (auto-pull; persist HF cache volume).
            Requires -d DIR (the host dir bind-mounted to /workspace inside the
            container). SPTXINSIGHT_DATA_DIR provides the same value non-interactively.

Environment overrides:
  SPTXINSIGHT_RUNNER              Default runner when --runner is not given (native | docker)
  SPTXINSIGHT_IMAGE                Override the docker image tag
  SPTXINSIGHT_HF_CACHE_VOLUME      Override the persistent HF model cache volume name
  SPTXINSIGHT_COMMANDS_TTL_SECONDS TTL for cached subcommand list (default 86400)
  SPTXINSIGHT_DATA_DIR             Default -d value when the flag is not given (used by docker)
  SPTXINSIGHT_STRICT=1             Die on unknown -X flags instead of warning + passthrough

Decision rule for argv parsing:
  Everything between the script's own flags and the first sptxinsight subcommand
  name (run/ingest/verify/...) is consumed by this script. From (and including)
  the first sptxinsight subcommand name, every remaining argument is passed
  verbatim to sptxinsight. Use -- to force passthrough explicitly.

Note: sptxinsight's CLI itself has a global --backend flag (anndata|zarr|spatialdata).
The wrapper's --runner is intentionally separate, to avoid overloading that flag.

Examples:
  $PROG run --samples ./samples --results ./results --base-type tumor --target-type lymphocyte       # native (default)
  $PROG --runner docker -d /workspace/project run --samples ./samples --results ./results
  $PROG --runner docker -d /workspace/project --gpu 0 --tmpdir /scratch --no-pull run --samples ./samples
  $PROG --runner docker -d /workspace/project --dry-run run --samples ./samples       # print resolved docker command, do not run
  $PROG doctor
  $PROG where
EOF
}

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
_log()  { printf '[%s] %s\n' "$PROG" "$*" >&2; }
_warn() { _log "WARNING: $*"; }
_die()  { _log "ERROR: $*"; exit "${2:-1}"; }

# ---------------------------------------------------------------------------
# sptxinsight subcommand discovery
# ---------------------------------------------------------------------------
_get_sptxinsight_command_list() {
    # 1. fresh cache wins
    # 2. live `sptxinsight schema --commands-only` writes the cache
    # 3. builtin static fallback

    if [[ -f "$COMMANDS_CACHE" ]]; then
        local cache_age
        cache_age=$(( $(date +%s) - $(stat -c %Y "$COMMANDS_CACHE" 2>/dev/null || echo 0) ))
        if [[ $cache_age -lt $COMMANDS_CACHE_TTL_SECONDS ]]; then
            cat "$COMMANDS_CACHE"
            return 0
        fi
    fi

    if command -v sptxinsight >/dev/null 2>&1; then
        local out
        if out="$(sptxinsight schema --commands-only 2>/dev/null)" \
           && [[ -n "$out" ]] \
           && command -v python3 >/dev/null 2>&1; then
            local extracted
            if extracted="$(printf '%s' "$out" | python3 -c "
import json, sys
try:
    d = json.load(sys.stdin)
    cmds = d.get('commands') or []
    if isinstance(cmds, dict): cmds = list(cmds.keys())
    for c in cmds: print(c)
except Exception:
    sys.exit(1)
")" && [[ -n "$extracted" ]]; then
                mkdir -p "$(dirname "$COMMANDS_CACHE")"
                printf '%s\n' "$extracted" > "$COMMANDS_CACHE"
                printf '%s\n' "$extracted"
                return 0
            fi
        fi
    fi

    # LAST RESORT
    _warn "could not refresh subcommand cache; using builtin list (last known good)"
    printf '%s\n' "${_SPT_BUILTIN_CMDS[@]}"
}

_is_sptxinsight_cmd() {
    local needle="$1" cmd
    for cmd in "${SPT_CMDS[@]}"; do
        [[ "$cmd" == "$needle" ]] && return 0
    done
    return 1
}

# Load subcommand list (once per invocation)
SPT_CMDS=()
while IFS= read -r c; do
    [[ -n "$c" ]] && SPT_CMDS+=("$c")
done < <(_get_sptxinsight_command_list)
[[ ${#SPT_CMDS[@]} -gt 0 ]] || _die "could not determine sptxinsight subcommand list (cache + fallback both empty)"

# ---------------------------------------------------------------------------
# Phase 1: parse argv (does NOT depend on runner choice)
# ---------------------------------------------------------------------------
SCRIPT_RUNNER=""        # native | docker
SCRIPT_GPU=""
SCRIPT_TMPDIR=""
SCRIPT_NO_PULL=0
SCRIPT_DRY_RUN=0
EXTRA_CMD=""           # status | doctor | where

while [[ $# -gt 0 ]]; do
    case "$1" in
        --runner)
            [[ -n "${2:-}" ]] || _die "--runner requires a value"
            SCRIPT_RUNNER="$2"; shift 2 ;;
        --runner=*)
            SCRIPT_RUNNER="${1#*=}"; shift ;;
        --gpu)
            [[ -n "${2:-}" ]] || _die "--gpu requires a value"
            SCRIPT_GPU="$2"; shift 2 ;;
        --gpu=*)
            SCRIPT_GPU="${1#*=}"; shift ;;
        --tmpdir)
            [[ -n "${2:-}" ]] || _die "--tmpdir requires a value"
            SCRIPT_TMPDIR="$2"; shift 2 ;;
        --tmpdir=*)
            SCRIPT_TMPDIR="${1#*=}"; shift ;;
        -d|--data-dir)
            [[ -n "${2:-}" ]] || _die "-d / --data-dir requires a value"
            SCRIPT_DATA_DIR="$2"; shift 2 ;;
        -d=*|--data-dir=*)
            SCRIPT_DATA_DIR="${1#*=}"; shift ;;
        --no-pull)
            SCRIPT_NO_PULL=1; shift ;;
        --dry-run|--dryrun)
            SCRIPT_DRY_RUN=1; shift ;;
        -h|--help)
            print_usage; exit 0 ;;
        --version)
            echo "$PROG $VERSION"; exit 0 ;;
        --)
            shift
            break ;;                       # explicit delimiter
        status|doctor|where)
            EXTRA_CMD="$1"; shift
            break ;;                       # script's own subcommands (priority)
        -b)
            # sptxinsight's CLI has its own global --backend flag, so this
            # wrapper does NOT claim '-b'. A user typing '-b docker' out of
            # habit from wsinsight.sh lands here. Surface a clear, named
            # hint instead of silently forwarding '-b docker ... schema' to
            # sptxinsight's --backend parser (which would 500 on the
            # remainder). The next token is the runner name; consume both.
            if [[ "${SPTXINSIGHT_STRICT:-0}" == "1" ]]; then
                _die "-b is reserved by sptxinsight's --backend flag; this wrapper's runner selector is --runner (see header)"
            fi
            [[ -n "${2:-}" ]] || _die "-b requires a value (e.g. -b docker); NOTE: -b is wsinsight.sh's runner shortcut. This wrapper uses --runner because sptxinsight's CLI has its own --backend flag."
            _warn "-b is the wsinsight.sh runner shortcut; on this wrapper use --runner $2 (sptxinsight's CLI also accepts --backend, which is why -b cannot be claimed here)."
            SCRIPT_RUNNER="$2"
            shift 2
            ;;
        -*)
            if [[ "${SPTXINSIGHT_STRICT:-0}" == "1" ]]; then
                _die "unknown option: $1 (SPTXINSIGHT_STRICT=1; known sptxinsight subcommands: ${SPT_CMDS[*]})"
            else
                _warn "unrecognized script flag: $1 -- treating as sptxinsight's (use SPTXINSIGHT_STRICT=1 to enforce)"
                break                       # let sptxinsight decide
            fi
            ;;
        *)
            # Position arg. Is it a known sptxinsight subcommand?
            if _is_sptxinsight_cmd "$1"; then
                break                       # passthrough starts here
            fi
            # Backward-scan fallback
            matched=0; found_idx=0; i=1
            for arg in "$@"; do
                if _is_sptxinsight_cmd "$arg"; then
                    matched=1; found_idx=$i; break
                fi
                i=$((i + 1))
            done
            if [[ $matched -eq 1 ]]; then
                shift $((found_idx - 1))
                break
            fi
            _die "no sptxinsight subcommand found in: $* (known: ${SPT_CMDS[*]})"
            ;;
    esac
done
SPTXINSIGHT_ARGS=("$@")

# ---------------------------------------------------------------------------
# Phase 2: dispatch on EXTRA_CMD (runner choice happens HERE)
# ---------------------------------------------------------------------------

if [[ -z "$SCRIPT_RUNNER" ]]; then
    SCRIPT_RUNNER="${SPTXINSIGHT_RUNNER:-native}"
fi
case "$SCRIPT_RUNNER" in
    native|docker) ;;
    *) _die "unknown runner: '$SCRIPT_RUNNER' (use 'native' or 'docker')" ;;
esac

# Map --gpu X to docker --gpus X with sane defaults. The result is APPENDED
# to the array named by $1 so that `--gpus all` stays as two separate argv
# tokens (`--gpus`, `all`), which docker requires; a single quoted string
# `--gpus all` is rejected as an unknown flag.
#
# Caller usage:
#   local -a gpus=()
#   resolve_docker_gpus_flag gpus
#   docker "${gpus[@]}" ...
resolve_docker_gpus_flag() {
    local -n _out="$1"
    case "$SCRIPT_GPU" in
        all|"")                 _out+=(--gpus all) ;;
        device=*|"capabilities"=*)
                                  _out+=(--gpus "$SCRIPT_GPU") ;;
        *)                       _out+=(--gpus "device=$SCRIPT_GPU") ;;
    esac
}

cmd_status() {
    cat <<EOF
PROG         : $PROG ($SCRIPT_DIR/$PROG)
VERSION      : $VERSION
Runner       : $SCRIPT_RUNNER$( [[ $SCRIPT_DRY_RUN -eq 1 ]] && echo " (dry-run)" )
Data dir(-d) : ${SCRIPT_DATA_DIR:-<unset - required only for --runner docker>}
GPU          : ${SCRIPT_GPU:-<default>$( [[ $SCRIPT_RUNNER == docker ]] && echo " (all)" )}
Tmpdir       : ${SCRIPT_TMPDIR:-<unchanged>}
No-pull      : $SCRIPT_NO_PULL
Image        : $IMAGE_ID
HF cache vol : $HF_CACHE_VOLUME
sptxinsight  : $(command -v sptxinsight 2>/dev/null || echo "<not on PATH>")
docker       : $(command -v docker 2>/dev/null || echo "<not on PATH>")
Subcommands  : ${#SPT_CMDS[@]} known; cache: $COMMANDS_CACHE
Pass-through : ${#SPTXINSIGHT_ARGS[@]} arg(s)${SPTXINSIGHT_ARGS[*]:+: ${SPTXINSIGHT_ARGS[*]}}
EOF
}

cmd_doctor() {
    local target="${DOC_RUNNER:-$SCRIPT_RUNNER}"
    local rc=0
    printf 'doctor (runner=%s)\n' "$target"
    case "$target" in
        native)
            local spti
            if spti="$(command -v sptxinsight)"; then
                printf '  [OK]  sptxinsight on PATH: %s\n' "$spti"
                if "$spti" --version >/dev/null 2>&1; then
                    printf '  [OK]  sptxinsight --version runs\n'
                else
                    printf '  [WARN] sptxinsight --version failed (env may be incomplete)\n'
                fi
            else
                printf '  [FAIL] sptxinsight not on PATH (activate the sptxinsight conda env, or use --runner docker)\n'
                rc=2
            fi
            if command -v nvidia-smi >/dev/null 2>&1; then
                local gpu_count
                gpu_count=$(nvidia-smi -L 2>/dev/null | wc -l)
                printf '  [OK]  nvidia-smi reports %s GPU(s)\n' "$gpu_count"
            else
                printf '  [INFO] nvidia-smi not on PATH (CPU-only ingest still works)\n'
            fi
            ;;
        docker)
            if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
                printf '  [OK]  docker daemon reachable\n'
                if docker image inspect "$IMAGE_ID" >/dev/null 2>&1; then
                    printf '  [OK]  image present locally: %s\n' "$IMAGE_ID"
                else
                    printf '  [INFO] image not local; will pull on first run: %s\n' "$IMAGE_ID"
                fi
                if docker volume inspect "$HF_CACHE_VOLUME" >/dev/null 2>&1; then
                    printf '  [OK]  HF cache volume exists: %s\n' "$HF_CACHE_VOLUME"
                else
                    printf '  [INFO] HF cache volume does not exist; will be auto-created: %s\n' "$HF_CACHE_VOLUME"
                fi
                if command -v nvidia-smi >/dev/null 2>&1; then
                    local gpu_count
                    gpu_count=$(nvidia-smi -L 2>/dev/null | wc -l)
                    printf '  [OK]  nvidia-smi reports %s GPU(s)\n' "$gpu_count"
                else
                    printf '  [WARN] nvidia-smi missing; GPU passthrough will not work\n'
                    rc=1
                fi
            else
                printf '  [FAIL] docker not reachable\n'
                rc=2
            fi
            ;;
    esac
    return $rc
}

cmd_where() {
    echo "$SCRIPT_DIR/$PROG"
}

# Build the resolved docker run args (printed by --dry-run). Takes the host
# data dir as $1 so this function doesn't depend on main()'s scope.
build_docker_command() {
    local data_dir="$1"
    local -a parts=(docker run --rm -it)
    resolve_docker_gpus_flag parts
    parts+=(--shm-size=32g --init)
    [[ -n "${HOST_UID:-}"   ]] && parts+=(-e HOST_UID)
    [[ -n "${HOST_GID:-}"   ]] && parts+=(-e HOST_GID)
    [[ -n "$SCRIPT_TMPDIR"  ]] && parts+=(-e TMPDIR="$SCRIPT_TMPDIR")
    parts+=(-v "$data_dir":/workspace -v "$HF_CACHE_VOLUME":/app/hf-cache)
    parts+=("$IMAGE_ID")
    if [[ ${#SPTXINSIGHT_ARGS[@]} -gt 0 ]]; then
        parts+=(sptxinsight "${SPTXINSIGHT_ARGS[@]}")
    fi
    # Shell-quote each token with printf %q. Tokens with embedded spaces
    # (e.g. `--gpus all` -> `--gpus` and `all`) end up as separate quoted
    # tokens, which is exactly what bash will see when re-parsing.
    local i
    for i in "${!parts[@]}"; do
        printf '%q' "${parts[$i]}"
        if [[ $i -lt $((${#parts[@]} - 1)) ]]; then printf ' '; fi
    done
}

case "$EXTRA_CMD" in
    where)   cmd_where; exit 0 ;;
    status)  cmd_status; exit 0 ;;
    doctor)  cmd_doctor; exit $? ;;
    "")      : ;;                       # fall through to run-dispatch
    *)       _die "unknown subcommand: $EXTRA_CMD" ;;
esac

# ---------------------------------------------------------------------------
# Phase 3: execute
# ---------------------------------------------------------------------------

main() {
    # No-arg invocation (or only --help / --version handled at the parser entry;
    # if we reach main with no sptxinsight subcommand and no script subcommand,
    # show usage).
    if [[ ${#SPTXINSIGHT_ARGS[@]} -eq 0 && -z "$EXTRA_CMD" ]]; then
        print_usage
        exit 0
    fi
    if [[ ${#SPTXINSIGHT_ARGS[@]} -eq 0 ]]; then
        if [[ "$SCRIPT_RUNNER" == "native" ]]; then
            _die "no sptxinsight subcommand supplied. Try: $PROG run --samples ./samples --results ./results --base-type tumor"
        fi
        # docker with no args -> interactive shell. Allowed by convention.
    fi

    case "$SCRIPT_RUNNER" in
        native)
            if [[ $SCRIPT_DRY_RUN -eq 1 ]]; then
                local spti_bin_dry
                if [[ -n "${SPTXINSIGHT_BIN:-}" ]]; then
                    spti_bin_dry="$SPTXINSIGHT_BIN"
                elif command -v sptxinsight >/dev/null 2>&1; then
                    spti_bin_dry="$(command -v sptxinsight)"
                else
                    spti_bin_dry="/opt/anaconda3/envs/sptxinsight/bin/sptxinsight"
                fi
                printf '+ SPTXINSIGHT_EXPERIMENTAL=%q %q %q\n' \
                    "${SPTXINSIGHT_EXPERIMENTAL:-1}" \
                    "$spti_bin_dry" "${SPTXINSIGHT_ARGS[*]}"
                exit 0
            fi
            local spti_bin
            if [[ -n "${SPTXINSIGHT_BIN:-}" ]]; then
                spti_bin="$SPTXINSIGHT_BIN"
            elif command -v sptxinsight >/dev/null 2>&1; then
                spti_bin="$(command -v sptxinsight)"
            else
                spti_bin="/opt/anaconda3/envs/sptxinsight/bin/sptxinsight"
            fi
            if [[ ! -x "$spti_bin" ]]; then
                _die "sptxinsight interpreter not found at '$spti_bin'. Activate the sptxinsight conda env or set SPTXINSIGHT_BIN."
            fi
            exec env -u CONDA_PREFIX -u CONDA_DEFAULT_ENV \
                PATH=/usr/bin:/bin \
                SPTXINSIGHT_EXPERIMENTAL="${SPTXINSIGHT_EXPERIMENTAL:-1}" \
                "$spti_bin" "${SPTXINSIGHT_ARGS[@]}"
            ;;
        docker)
            # Determine data dir: required for docker. The flag (-d / --data-dir)
            # takes precedence; SPTXINSIGHT_DATA_DIR is the env-var form.
            local DATA_DIR="$SCRIPT_DATA_DIR"
            if [[ -z "$DATA_DIR" ]]; then
                _die "docker runner needs a data dir. Pass -d DIR (or --data-dir DIR), or set SPTXINSIGHT_DATA_DIR=/path. The host dir is bind-mounted to /workspace inside the container; SPTXINSIGHT_ARGS sees /workspace as cwd."
            fi
            if [[ ! -d "$DATA_DIR" ]]; then
                _die "--data-dir '$DATA_DIR' does not exist (or is not a directory). Pass -d DIR pointing at an existing path."
            fi
            if [[ $SCRIPT_DRY_RUN -eq 1 ]]; then
                echo "+ $(build_docker_command "$DATA_DIR")"
                exit 0
            fi
            if [[ $SCRIPT_NO_PULL -eq 0 ]]; then
                docker pull "$IMAGE_ID" >/dev/null 2>&1 \
                    || _warn "docker pull failed; using local image (if any)"
            fi
            local -a docker_args=( run --rm -it )
            resolve_docker_gpus_flag docker_args
            docker_args+=( --shm-size=32g --init )
            [[ -n "${HOST_UID:-}"   ]] && docker_args+=( -e HOST_UID )
            [[ -n "${HOST_GID:-}"   ]] && docker_args+=( -e HOST_GID )
            [[ -n "$SCRIPT_TMPDIR"  ]] && docker_args+=( -e TMPDIR="$SCRIPT_TMPDIR" )
            docker_args+=( -v "$DATA_DIR":/workspace -v "$HF_CACHE_VOLUME":/app/hf-cache )
            docker_args+=( "$IMAGE_ID" )
            if [[ ${#SPTXINSIGHT_ARGS[@]} -gt 0 ]]; then
                docker_args+=( sptxinsight "${SPTXINSIGHT_ARGS[@]}" )
            fi
            exec docker "${docker_args[@]}"
            ;;
    esac
}

main "$@"
