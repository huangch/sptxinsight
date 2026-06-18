"""sptxinsight: spatial-transcriptomics analysis sibling of WSInsight.

Reuses WSInsight's H-Plot spatial-heterogeneity engine (vendored under
``sptxinsight.insightlib``) but ingests AnnData/spatial-omics samples instead of
whole-slide images.
"""

from __future__ import annotations

try:
    from ._version import __version__
except Exception:  # pragma: no cover - version file is generated at build time
    __version__ = "0.0.0"

__all__ = ["__version__"]
