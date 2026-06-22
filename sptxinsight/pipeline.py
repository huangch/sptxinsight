"""End-to-end H-Plot pipeline for spatial-transcriptomics samples.

Bridges AnnData samples to the vendored WSInsight H-Plot engine:

1. read each sample (``io.read_sample``),
2. write the per-sample CSV contract (``adapt.anndata_to_contract``),
3. run the WSI-free H-Plot generation + finalize.

Because spatial coordinates are already in microns, every sample is registered
with ``slide_mpp_lookup={slide_id: 1.0}`` and ``wsi_dir=None`` so no whole-slide
image is ever opened.
"""

from __future__ import annotations

import logging
from typing import List
from typing import Sequence

from .adapt import anndata_to_contract
from .insightlib.hplot_generation import hplot_finalize
from .insightlib.hplot_generation import hplot_generation
from .io import read_sample
from .uri_path import URIPath

_logger = logging.getLogger(__name__)


def adapt_samples(
    sample_uris: Sequence["str | URIPath"],
    results_dir: "str | URIPath",
    *,
    cell_type_key: str = "cell_type",
    spatial_key: str = "spatial",
    genes: Sequence[str] = (),
    expression_matrix: str = "X",
) -> tuple[list[URIPath], dict[str, float], list[str]]:
    """Read + adapt samples into the H-Plot CSV contract.

    Returns ``(slide_paths, slide_mpp_lookup, vocab)`` where ``slide_paths`` are
    name-only ``URIPath`` stems (never opened) and ``vocab`` is the union of all
    sanitized cell-type names across samples. ``genes`` (with
    ``expression_matrix``) are additionally written as ``expr_<gene>`` columns
    for gene-based base/target H-Plots.
    """
    results_dir = (
        results_dir if isinstance(results_dir, URIPath) else URIPath(str(results_dir))
    )
    slide_paths: list[URIPath] = []
    mpp_lookup: dict[str, float] = {}
    vocab: set[str] = set()

    for uri in sample_uris:
        upath = uri if isinstance(uri, URIPath) else URIPath(str(uri))
        slide_id = upath.sample_id
        _logger.info("Ingesting sample %s", slide_id)
        adata = read_sample(upath)
        types = anndata_to_contract(
            adata,
            slide_id,
            results_dir,
            cell_type_key=cell_type_key,
            spatial_key=spatial_key,
            genes=genes,
            expression_matrix=expression_matrix,
        )
        vocab.update(types)
        # Name-only stem: hplot_generation uses it solely to align CSV outputs.
        slide_paths.append(URIPath(f"{slide_id}.h5ad"))
        mpp_lookup[slide_id] = 1.0  # microns ⇒ mpp = 1

    return slide_paths, mpp_lookup, sorted(vocab)


def run_hplot(
    sample_uris: Sequence["str | URIPath"],
    results_dir: "str | URIPath",
    base_types: Sequence[str] | None,
    target_types: Sequence[str] | None,
    *,
    cell_type_key: str = "cell_type",
    spatial_key: str = "spatial",
    base_by: str = "celltype",
    target_by: str = "celltype",
    expression_matrix: str = "X",
    base_gene_threshold: float = 0.0,
    max_neighbor_distance_um: float = 25.0,
    hplot_k: int = 2,
    hplot_N: int = 8,
    hplot_R: float = 0.5,
    hplot_range_max: int | None = None,
    hplot_range_min: int | None = None,
    samples_with_valid_range_only: bool = False,
    num_workers: int = 8,
    overwrite: bool = False,
) -> List[str]:
    """Adapt samples and compute aggregated H-Plot outputs.

    Writes ``<results_dir>/hplot-outputs.csv`` via the vendored engine and
    returns the list of slide ids that FAILED H-Plot generation (empty on full
    success), mirroring ``hplot_generation``'s contract.

    ``base_by`` / ``target_by`` select whether ``base_types`` / ``target_types``
    name cell types (``"celltype"``) or genes (``"gene"``). In gene mode the
    listed names are gene symbols pulled from ``expression_matrix``; the target
    H-Plot value becomes mean expression per layer, and a gene-based base region
    is defined by mean expression exceeding ``base_gene_threshold``.
    """
    results_dir = (
        results_dir if isinstance(results_dir, URIPath) else URIPath(str(results_dir))
    )
    results_dir.mkdir(parents=True, exist_ok=True)

    base_types = list(base_types or [])
    target_types = list(target_types or [])
    genes: list[str] = []
    if base_by == "gene":
        genes.extend(base_types)
    if target_by == "gene":
        genes.extend(target_types)
    genes = list(dict.fromkeys(genes))  # de-dup, preserve order

    slide_paths, mpp_lookup, vocab = adapt_samples(
        sample_uris,
        results_dir,
        cell_type_key=cell_type_key,
        spatial_key=spatial_key,
        genes=genes,
        expression_matrix=expression_matrix,
    )
    _logger.info("Cell-type vocabulary: %s", vocab)
    if genes:
        _logger.info("Gene columns written from %s: %s", expression_matrix, genes)

    failed = hplot_generation(
        wsi_dir=None,
        slide_paths=slide_paths,
        results_dir=results_dir,
        base_type_list=list(base_types) if base_types else None,
        target_type_list=list(target_types) if target_types else None,
        base_by=base_by,
        target_by=target_by,
        base_gene_threshold=base_gene_threshold,
        max_neighbor_distance_um=max_neighbor_distance_um,
        hplot_k=hplot_k,
        hplot_N=hplot_N,
        hplot_R=hplot_R,
        hplot_range_max=hplot_range_max,
        hplot_range_min=hplot_range_min,
        hplot_samples_with_valid_range_only=samples_with_valid_range_only,
        num_workers=num_workers,
        slide_mpp_lookup=mpp_lookup,
        overwrite=overwrite,
    )
    # hplot_generation already wrote a hplot-outputs.csv; rebuild the canonical
    # gap-filled table from the per-sample CSVs (cheap + idempotent).
    hplot_finalize(results_dir, overwrite=True)
    if failed:
        _logger.warning(
            "H-Plot generation failed for %d sample(s): %s", len(failed), failed
        )
    return failed
