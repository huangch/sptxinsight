"""Compatibility shim for the WSI hooks used by the vendored ``hplot_generation``.

WSInsight's ``hplot_generation`` was written for whole-slide images and pulls the
micron-per-pixel spacing from the slide. In sptxinsight there is no slide: spatial
coordinates are already in microns, so the equivalent ``mpp`` is ``1.0``.

These two functions only execute on the WSI code path (``wsi_dir`` not ``None`` or
no ``slide_mpp_lookup`` entry). The sptxinsight drivers always pass
``wsi_dir=None`` and ``slide_mpp_lookup={slide_id: 1.0}``, so neither is reached in
normal operation; they exist so the vendored module imports cleanly and degrades
with a clear error if the WSI path is ever taken by mistake.
"""

from __future__ import annotations

from ..uri_path import URIPath


def get_avg_mpp(wsi_path: URIPath) -> float:
    """Return microns-per-pixel for a slide.

    sptxinsight has no slide to read. Callers must supply spacing through
    ``slide_mpp_lookup`` (microns ⇒ ``1.0``); reaching this means that lookup was
    missing for a sample.
    """
    raise RuntimeError(
        "sptxinsight has no whole-slide image to read mpp from. "
        "Pass slide_mpp_lookup={slide_id: 1.0} so spatial coordinates are treated "
        f"as microns (offending path: {wsi_path})."
    )


def _validate_wsi_directory(wsi_dir: URIPath) -> None:
    """No-op: sptxinsight never passes a real ``--wsi-dir``."""
    return None
