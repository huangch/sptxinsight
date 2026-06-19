#!/usr/bin/env python
"""Comprehensive HPV-stratified CCI border H-Plot screen.

Renders EVERY panel-valid ligand-receptor pair in BOTH directions
(in = receiver projection, out = sender projection) as an HPV+ vs Rest border
H-Plot, and writes a ranked summary (smallest per-layer Mann-Whitney p first)
so the most-separated pairs can be triaged quickly.

Outputs (under <results>/cci_hplot_screen/):
  cci_<LIG>_<REC>_<in|out>_hpv.png   per-panel figures (40 total)
  cci_screen_in.png / cci_screen_out.png   combined grids
  cci_screen_ranking.csv             ranked table for screening
"""
import json
import os
import glob
import math

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hplot.core import HPlot

R = "/workspace/wsinsight/wsinsight-model-development/data/results/head_neck/cme_run_gex"
CELLS = sorted(glob.glob(f"{R}/hplot-outputs-csv/cells/*.csv"))
OUT = f"{R}/cci_hplot_screen"
os.makedirs(OUT, exist_ok=True)
LAYER = "signed_distance_to_border"
LMIN, LMAX = -8, 20
PALETTE = {"HPV+": "#d62728", "Rest": "#1f77b4"}
ORDER = ["HPV+", "Rest"]
PYLIM = (1e-4, 1.0)  # fixed p-axis so panels are comparable across the screen

with open("/tmp/hn_hpv_map.json") as fp:
    HPV = json.load(fp)

def hpv_arm(stem):
    p16 = HPV.get(stem.split("_")[0], {}).get("p16", "N/A")
    return "HPV+" if p16 == "Positive" else "Rest"

# ---- assemble one long table: (sample, layer) mean per cci _in/_out column ----
recs = []
for path in CELLS:
    stem = os.path.basename(path)[:-4]
    df = pd.read_csv(path)
    df = df[(df[LAYER] >= LMIN) & (df[LAYER] <= LMAX)]
    cci_cols = [
        c for c in df.columns
        if c.startswith("cci_") and (c.endswith("_in_mean") or c.endswith("_out_mean"))
    ]
    agg = df.groupby(LAYER)[cci_cols].mean().reset_index()
    agg["sample"] = stem
    agg["hpv"] = hpv_arm(stem)
    recs.append(agg)
long = pd.concat(recs, ignore_index=True)
pairs = sorted(
    c[len("cci_"):-len("_in_mean")]
    for c in long.columns
    if c.endswith("_in_mean")
)
print(f"samples: {long['sample'].nunique()} | pairs: {len(pairs)} | "
      f"HPV+={long[long.hpv=='HPV+']['sample'].nunique()} "
      f"Rest={long[long.hpv=='Rest']['sample'].nunique()}")

def fit_one(col):
    sub = long[[LAYER, col, "sample", "hpv"]].dropna()
    h = HPlot()
    h.fit(sub, targets=col, layer=LAYER, group="hpv",
          color_map=PALETTE, legend_order=ORDER,
          pvalue=True, pvalue_test="mannwhitney", pvalue_correction="fdr_bh")
    return h, sub

def label_for(lig, rec, direction):
    if direction == "in":
        return f"{lig}\u2192{rec}\n(on {rec}\u207a receiver cells)"
    return f"{lig}\u2192{rec}\n(on {lig}\u207a sender cells)"

def summarize(h, col, direction):
    pv = h.layer_pvalues_.dropna(subset=["p_value"])
    if len(pv) == 0:
        return dict(min_p=np.nan, min_p_adj=np.nan, layer_at_min=np.nan,
                    higher_arm="", n_sig_layers=0)
    i = pv["p_value"].idxmin()
    lyr = int(round(pv.loc[i, "layer"]))
    g = long[long[LAYER] == lyr].groupby("hpv")[col].mean()
    higher = "HPV+" if g.get("HPV+", 0.0) >= g.get("Rest", 0.0) else "Rest"
    return dict(
        min_p=float(pv.loc[i, "p_value"]),
        min_p_adj=float(pv.loc[i, "p_adj"]),
        layer_at_min=lyr,
        higher_arm=higher,
        n_sig_layers=int((pv["p_value"] < 0.05).sum()),
    )

rows = []
for direction in ("in", "out"):
    ncol = 4
    nrow = math.ceil(len(pairs) / ncol)
    figG, axesG = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 4.4 * nrow), squeeze=False)
    flat = axesG.ravel()
    for k, stem in enumerate(pairs):
        lig, rec = stem.split("_", 1)
        col = f"cci_{stem}_{direction}_mean"
        h, _ = fit_one(col)
        s = summarize(h, col, direction)
        s.update(pair=stem, ligand=lig, receptor=rec, direction=direction)
        rows.append(s)
        lab = label_for(lig, rec, direction)
        # per-panel figure
        figP, axP = plt.subplots(figsize=(6.2, 4.4))
        h.plot(ax=axP, display_base_type="tumor", display_target_type=lab,
               value_kind="interaction", pvalue_show=True,
               pvalue_use_adjusted=False, pvalue_threshold=0.05, pvalue_ylim=PYLIM)
        figP.tight_layout()
        figP.savefig(f"{OUT}/cci_{stem}_{direction}_hpv.png", dpi=150, bbox_inches="tight")
        plt.close(figP)
        # grid panel
        h.plot(ax=flat[k], display_base_type="tumor", display_target_type=lab,
               value_kind="interaction", pvalue_show=True,
               pvalue_use_adjusted=False, pvalue_threshold=0.05, pvalue_ylim=PYLIM)
    for ax in flat[len(pairs):]:
        ax.axis("off")
    side = "receiver" if direction == "in" else "sender"
    figG.suptitle(f"CCI {direction.upper()} ({side}) vs tumor border \u2014 HPV+ vs Rest",
                  fontsize=16, y=1.004)
    figG.tight_layout()
    figG.savefig(f"{OUT}/cci_screen_{direction}.png", dpi=200, bbox_inches="tight")
    plt.close(figG)

rank = pd.DataFrame(rows)[
    ["pair", "direction", "ligand", "receptor",
     "min_p", "min_p_adj", "layer_at_min", "higher_arm", "n_sig_layers"]
].sort_values(["min_p", "n_sig_layers"], ascending=[True, False]).reset_index(drop=True)
rank.to_csv(f"{OUT}/cci_screen_ranking.csv", index=False)
print(f"\nwrote {len(pairs)*2} panels + 2 grids + ranking to {OUT}\n")
print("=== top 15 most-separated (pair x direction) ===")
print(rank.head(15).to_string(index=False))
