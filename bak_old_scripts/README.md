# bak_old_scripts/

Legacy sptxinsight runner scripts that have been **superseded** by [`sptxinsight.sh`](../sptxinsight.sh) at the repo root.

## What moved here

| File | Date moved | Reason |
|---|---|---|
| `sptxinsight-docker-run.sh` | 2026-09-05 | Folded into `sptxinsight.sh --runner docker`. The new wrapper supersedes this script and fixes its `IMAGE_ID=sptxinsight:latest` (now `huangchtw/sptxinsight:latest`). |

## Why kept

- Source-of-truth for the historical docker invocation pattern (mount structure, env flags, GPU passthrough shape).
- The legacy script's `IMAGE_ID=sptxinsight:latest` references an **un-published** local-only image; the new wrapper uses **`huangchtw/sptxinsight:latest`** (matching Docker Hub).
- Useful if you need to reproduce exactly what the old wrapper did; **not** the recommended entry point going forward.

## Migration

Replace
```
bash sptxinsight-docker-run.sh /path/to/data --gpu 0 sptxinsight run ...
```
with
```
export SPTXINSIGHT_DATA_DIR=/path/to/data
./sptxinsight.sh --runner docker --gpu 0 run ...
```

Anything you used to pass through the legacy script after the data dir is now passed through `./sptxinsight.sh` after the first sptxinsight subcommand name (the wrapper auto-detects the boundary).
