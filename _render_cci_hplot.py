#!/usr/bin/env python
"""Render HPV-stratified CCI border H-Plots for the head_neck cohort.

Reads the per-cell hplot cells CSVs (which carry signed_distance_to_border plus
all cci_* score columns), aggregates each sample to per-(sample, layer) means,
and draws HPV+ vs Rest H-Plots per ligand-receptor pair with the hplot package
(value_kind='interaction').
"""
import json
import os
import glob

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"  # keep text as <text> in SVG
import matplotlib.pyplot as plt
import pandas as pd

from hplot.core import HPlot

R = "/workspace/wsinsight/wsinsight-model-development/data/results/head_neck/cme_run_gex"
CELLS = sorted(glob.glob(f"{R}/hplot-outputs-csv/cells/*.csv"))
OUT = f"{R}/cci_hplot_hpv_panels"
os.makedirs(OUT, exist_ok=True)
LAYER = "signed_distance_to_border"
LMIN, LMAX = -8, 20

with open("/tmp/hn_hpv_map.json") as fp:
    HPV = json.load(fp)

def hpv_arm(stem):
    gsm = stem.split("_")[0]
    p16 = HPV.get(gsm, {}).get("p16", "N/A")
    return "HPV+" if p16 == "Positive" else "Rest"

# HPV-associated genes present in the head_neck panel: the HPV+ "immune-hot"
# checkpoint / cytotoxic / chemokine signature, plus proliferation & the
# HPV- skewed oncogenes (EGFR/MYC/SOX2). Any L-R pair touching one of these
# is rendered.
HPV_GENES = {
    # checkpoints / costimulation
    "CD274", "PDCD1", "CTLA4", "CD80", "CD86", "CD28",
    "LAG3", "HAVCR2", "TIGIT",
    # cytotoxic T cells
    "CD8A", "GZMB",
    # IFN-gamma / lymphoid-recruitment chemokines
    "CXCL9", "CXCL10", "CCL19", "CCR7", "CCL5",
    # proliferation surrogates (p16 biology)
    "MKI67", "PCNA",
    # oncogenes (HPV- skewed)
    "EGFR", "MYC", "SOX2",
    # broader immune context (costimulation, B/plasma, lymphocyte homing,
    # atypical chemokine receptors) -- not HPV-status markers but immune-related
    "CD70", "CD27", "TNFRSF13B", "TNFRSF17",
    "SELL", "ACKR1", "CXCL6",
}
PALETTE = {"HPV+": "#d62728", "Rest": "#1f77b4"}
ORDER = ["HPV+", "Rest"]

# ---- assemble one long table: (sample, layer) mean per cci column ----
recs = []
for path in CELLS:
    stem = os.path.basename(path)[:-4]
    df = pd.read_csv(path)
    df = df[(df[LAYER] >= LMIN) & (df[LAYER] <= LMAX)]
    cci_cols = [c for c in df.columns if c.startswith("cci_") and c.endswith("_in_mean")]
    agg = df.groupby(LAYER)[cci_cols].mean().reset_index()
    agg["sample"] = stem
    agg["hpv"] = hpv_arm(stem)
    recs.append(agg)
long = pd.concat(recs, ignore_index=True)
n_pos = long[long.hpv == "HPV+"]["sample"].nunique()
n_rest = long[long.hpv == "Rest"]["sample"].nunique()
print(f"samples: HPV+={n_pos}  Rest={n_rest}")

# derive available (ligand, receptor) pairs from the _in_mean columns and keep
# those that touch an HPV-associated gene.
_avail = sorted(
    c[len("cci_"):-len("_in_mean")]
    for c in long.columns
    if c.startswith("cci_") and c.endswith("_in_mean")
)
PAIRS = []
for stem in _avail:
    lig, rec = stem.split("_", 1)
    if lig in HPV_GENES or rec in HPV_GENES:
        PAIRS.append((lig, rec))
print(f"HPV-related pairs ({len(PAIRS)}): " + ", ".join(f"{a}->{b}" for a, b in PAIRS))

def render(ax, lig, rec):
    col = f"cci_{lig}_{rec}_in_mean"
    sub = long[[LAYER, col, "sample", "hpv"]].dropna()
    h = HPlot()
    h.fit(sub, targets=col, layer=LAYER, group="hpv",
          color_map=PALETTE, legend_order=ORDER,
          pvalue=True, pvalue_test="mannwhitney", pvalue_correction="fdr_bh")
    h.plot(ax=ax, display_base_type="tumor",
           display_target_type=f"{lig}\u2192{rec}\n(on {rec}\u207a receiver cells)",
           value_kind="interaction",
           pvalue_show=True, pvalue_use_adjusted=False, pvalue_threshold=0.05,
           pvalue_ylim=(1e-4, 1.0))
    return col

# ---- per-pair PNG + SVG ----
for lig, rec in PAIRS:
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    render(ax, lig, rec)
    fig.tight_layout()
    stem = f"{OUT}/cci_{lig}_{rec}_in_mean_hpv"
    fig.savefig(stem + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(stem + ".svg", bbox_inches="tight")
    plt.close(fig)

# ---- combined grid (4 columns) ----
import math
ncol = 4
nrow = math.ceil(len(PAIRS) / ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 4.4 * nrow), squeeze=False)
flat = axes.ravel()
for ax, (lig, rec) in zip(flat, PAIRS):
    render(ax, lig, rec)
for ax in flat[len(PAIRS):]:
    ax.axis("off")
fig.suptitle("CCI (ligand\u2192receptor, receiver) vs tumor border \u2014 HPV+ vs Rest",
             fontsize=15, y=1.02)
fig.tight_layout()
fig.savefig(f"{R}/cci_hplot_hpv.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{R}/cci_hplot_hpv.svg", bbox_inches="tight")
plt.close(fig)
print("wrote", OUT, "and cci_hplot_hpv.png/.svg")
