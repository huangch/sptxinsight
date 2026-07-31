"""KurtoRank integration helpers used by sptxinsight commands."""

from __future__ import annotations

from dataclasses import dataclass
import multiprocessing
import os
from pathlib import Path
from typing import Sequence

import anndata as ad
import click

from ..uri_path import _LIST_SCHEMES
from ..uri_path import URIPath


@dataclass(frozen=True)
class AnnotateConfig:
    """Shared annotate settings independent of the input mode."""

    tissue_type: str
    spatial_key: str = "spatial"
    cell_type_key: str = "cell_type"
    markers_csv: Path | None = None
    common_only: bool = True
    normal_only: bool = False
    use_graphclust: bool = True
    chosen_leiden_res: float = 0.5
    min_genes: int = 10
    min_cells: int = 5
    lower_percentile: float = 1.0
    upper_percentile: float = 99.0
    n_top_genes: int = 3500
    n_perm: int = 1000
    n_jobs: int = 32
    allow_cuda_parallel: bool = False
    methods: tuple[str, ...] = ()
    generate_plots: bool = True
    overwrite: bool = False
    regenerate_plots: bool = False
    plot_format: str = "png"
    plot_dpi: int = 600
    use_top_k_markers: int | None = None
    verbose: bool = False


def _mode_from_uri(input_path: URIPath) -> str:
    if input_path.scheme in _LIST_SCHEMES:
        return "anndata"
    p = Path(input_path.local_path())
    if p.is_file() and p.suffix.lower() in {".h5ad", ".zarr"}:
        return "anndata"
    if p.is_dir():
        matrix_candidates = ("cell_by_feature_matrix.h5", "cell_feature_matrix.h5")
        if any((p / s).exists() for s in matrix_candidates) and (p / "cells.csv.gz").exists():
            return "xenium"
        if any(child.suffix.lower() in {".h5ad", ".zarr"} for child in p.iterdir()):
            return "anndata"
    return "unsupported"


def _resolve_xenium_dir(input_path: URIPath) -> Path:
    """Return the exact Xenium directory supplied by the caller."""
    return Path(input_path.local_path())


def _collect_anndata_samples(input_path: URIPath) -> list[URIPath]:
    """Collect AnnData sample URIs from file, directory, or manifest input."""
    if input_path.scheme in _LIST_SCHEMES:
        return list(input_path.iterdir())

    p = Path(input_path.local_path())
    if p.is_file():
        return [input_path]

    return [
        input_path / child.name
        for child in p.iterdir()
        if child.suffix.lower() in {".h5ad", ".zarr"}
    ]


def _load_anndata_sample(sample_path: URIPath) -> ad.AnnData:
    """Load one AnnData sample from a URIPath."""
    sample_local = Path(sample_path.local_path())
    if sample_local.suffix.lower() == ".h5ad":
        return ad.read_h5ad(sample_local)
    return ad.read_zarr(sample_local)


def _configure_runtime() -> None:
    """Match KurtoRank CLI runtime safeguards for process and BLAS behavior."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
    try:
        multiprocessing.set_start_method("spawn", force=True)
    except (RuntimeError, ValueError):
        pass
    try:
        import torch

        torch.multiprocessing.set_sharing_strategy("file_system")
    except Exception:
        pass


def _run_kurtorank_annotate_xenium(
    *,
    xenium_dir: Path,
    output_dir: Path,
    config: AnnotateConfig,
) -> None:
    """Execute KurtoRank annotate as an in-process Click command."""
    try:
        from kurtorank.annotate.main import annotate_cmd
    except Exception as exc:  # pragma: no cover - import-time environment issue
        raise click.ClickException(
            "Failed to import KurtoRank annotate pipeline. "
            "Ensure kurtorank and its scientific dependencies are installed. "
            f"Import error: {exc}"
        ) from None

    _configure_runtime()

    args: list[str] = [
        "--xenium-dir",
        str(xenium_dir),
        "--tissue-type",
        str(config.tissue_type),
        "--output-dir",
        str(output_dir),
        "--chosen-leiden-res",
        str(config.chosen_leiden_res),
        "--min-genes",
        str(config.min_genes),
        "--min-cells",
        str(config.min_cells),
        "--lower-percentile",
        str(config.lower_percentile),
        "--upper-percentile",
        str(config.upper_percentile),
        "--n-top-genes",
        str(config.n_top_genes),
        "--n-perm",
        str(config.n_perm),
        "--n-jobs",
        str(config.n_jobs),
        "--plot-format",
        str(config.plot_format),
        "--plot-dpi",
        str(config.plot_dpi),
    ]

    if config.markers_csv is not None:
        args.extend(["--markers-csv", str(config.markers_csv)])

    args.append("--common-only" if config.common_only else "--no-common-only")
    args.append("--normal-only" if config.normal_only else "--include-cancer")
    args.append("--use-graphclust" if config.use_graphclust else "--use-leiden")
    args.append(
        "--allow-cuda-parallel" if config.allow_cuda_parallel else "--no-allow-cuda-parallel"
    )
    args.append("--generate-plots" if config.generate_plots else "--no-generate-plots")
    args.append("--overwrite" if config.overwrite else "--no-overwrite")
    args.append(
        "--regenerate-plots" if config.regenerate_plots else "--no-regenerate-plots"
    )
    args.append("--verbose" if config.verbose else "--quiet")

    for method in config.methods:
        args.extend(["--method", str(method)])

    if config.use_top_k_markers is not None:
        args.extend(["--use-top-k-markers", str(config.use_top_k_markers)])

    try:
        annotate_cmd.main(args=args, prog_name="kurtorank annotate", standalone_mode=False)
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(f"KurtoRank annotate failed: {exc}") from None


def _run_kurtorank_annotate_anndata(
    *,
    adata: ad.AnnData,
    output_dir: Path,
    config: AnnotateConfig,
) -> ad.AnnData:
    """Execute KurtoRank internals on an in-memory AnnData object."""
    try:
        import kurtorank.annotate.main as _kmain
        from kurtorank.annotate.main import _default_markers_csv
        from kurtorank.annotate.main import export_qust_csvs
        from kurtorank.annotate.main import resolve_method_switch
        from kurtorank.annotate.main import run_kurtorank
        from kurtorank.annotate.main import run_leiden_scan
        from kurtorank.annotate.main import setup_logging
        from kurtorank.annotate.main import visualize_annotation
    except Exception as exc:  # pragma: no cover - import-time environment issue
        raise click.ClickException(
            "Failed to import KurtoRank annotate internals for AnnData mode. "
            "Ensure kurtorank and its scientific dependencies are installed. "
            f"Import error: {exc}"
        ) from None

    _configure_runtime()
    setup_logging(config.verbose)
    _kmain.PLOT_FORMAT = str(config.plot_format).lower()
    _kmain.PLOT_DPI = int(config.plot_dpi)

    output_dir.mkdir(parents=True, exist_ok=True)

    if config.spatial_key not in adata.obsm:
        raise click.ClickException(
            f"AnnData mode requires adata.obsm[{config.spatial_key!r}] coordinates."
        )
    if config.spatial_key != "spatial":
        adata.obsm["spatial"] = adata.obsm[config.spatial_key]

    use_graphclust = config.use_graphclust and "graphclust" in adata.obs

    method_switch = resolve_method_switch(config.methods)
    marker_path = config.markers_csv if config.markers_csv is not None else _default_markers_csv()

    primary_cluster = run_leiden_scan(
        adata,
        use_graphclust=use_graphclust,
        chosen_leiden_res=config.chosen_leiden_res,
        generate_plots=config.generate_plots,
        out_dir=output_dir,
    )
    adata, _ = run_kurtorank(
        adata,
        primary_cluster=primary_cluster,
        markers_csv=marker_path,
        tissue_type=config.tissue_type,
        common_only=config.common_only,
        normal_only=config.normal_only,
        method_switch=method_switch,
        n_perm=config.n_perm,
        n_jobs=config.n_jobs,
        generate_plots=config.generate_plots,
        out_dir=output_dir,
        use_top_k_markers=config.use_top_k_markers,
        allow_cuda_parallel=config.allow_cuda_parallel,
    )
    adata.uns["primary_cluster_key"] = primary_cluster

    if config.generate_plots:
        visualize_annotation(
            adata,
            primary_cluster=primary_cluster,
            generate_plots=config.generate_plots,
            out_dir=output_dir,
            method_switch=method_switch,
        )

    export_qust_csvs(adata, xenium_dir=output_dir, out_dir=output_dir)

    adata.obs[config.cell_type_key] = adata.obs["cell_subtype"].astype(str)
    adata.write_h5ad(output_dir / "annotated.h5ad")
    return adata


def run_kurtorank_annotate_input(
    *,
    input_path: URIPath,
    output_dir: Path,
    config: AnnotateConfig,
) -> None:
    """Dispatch KurtoRank annotate for Xenium or AnnData inputs."""
    mode = _mode_from_uri(input_path)
    if mode == "unsupported":
        raise click.ClickException(
            "Unsupported input for sptxinsight annotate. Provide a local Xenium "
            "directory containing Xenium output files directly, or an AnnData "
            "file / directory."
        )

    if mode == "xenium":
        xenium_dir = _resolve_xenium_dir(input_path)
        _run_kurtorank_annotate_xenium(xenium_dir=xenium_dir, output_dir=output_dir, config=config)
        return

    samples = _collect_anndata_samples(input_path)
    if not samples:
        raise click.ClickException(f"No .h5ad/.zarr samples found under {input_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    for idx, sample_path in enumerate(samples):
        sample_local = Path(sample_path.local_path())
        adata = _load_anndata_sample(sample_path)
        sample_out = output_dir if len(samples) == 1 else output_dir / sample_path.sample_id
        sample_out.mkdir(parents=True, exist_ok=True)
        _run_kurtorank_annotate_anndata(adata=adata, output_dir=sample_out, config=config)
        click.echo(f"[{idx + 1}/{len(samples)}] Annotated {sample_local.name} -> {sample_out}")
