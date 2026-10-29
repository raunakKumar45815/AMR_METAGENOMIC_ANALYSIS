#!/usr/bin/env python3
"""
Bubble plots of the top MGE genes across STPs, months and seasons, using the
per-Gbp normalised matrix. Same structure and styling as the original script,
with four changes:

  1. input is MGE_gene_matrix_per_Gbp.tsv (Gene, Function, Origin, samples...)
  2. the third column (Origin) is ignored
  3. values stay on the MGE hits per Gbp scale - no relative-abundance step
  4. plots and the sheets behind them go to a separate output directory
"""
import os

import matplotlib
matplotlib.use("Agg")          # no display needed on a server
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

# -------- USER INPUTS --------
MATRIX_FILE = "results_1/MGE_gene_matrix_per_Gbp.tsv"
META_FILE   = "/data/Students/Raunak/RF-AMR-24-STP/MGE_RESULTS/new_results/MasterMetadataFinal.csv"
OUT_DIR     = "bubble_plots"
TOP_N       = 50
VALUE_LABEL = "MGE hits per Gbp"


# How values from several samples are combined into one bubble.
# "sum" reproduces the original script. Note that groups hold different numbers
# of samples (12 per STP, ~11 per month), so sums partly reflect group size;
# "mean" makes bubbles comparable across groups. Switch here if you prefer that.
AGG = "sum"

# -------- STYLING --------
FONT_SCALE   = 1.30          # multiplies every font size below
BASE_FONT    = 12            # tick labels
LABEL_FONT   = 15            # axis labels
TITLE_FONT   = 17            # plot title
LEGEND_FONT  = 12            # legend entries
BUBBLE_SIZES = (150, 1200)   # smallest and largest bubble area
BOLD         = "bold"        # set to "normal" to undo the bolding

# Figure sizes, one per panel. With 50 gene rows the canvas needs to be taller
# than it is wide, or the larger bubbles collide vertically.
FIG_STP    = (13, 18)
FIG_MONTH  = (15, 18)
FIG_SEASON = (10, 18)
# -----------------------------

plt.rcParams.update({
    "font.size":        BASE_FONT * FONT_SCALE,
    "font.weight":      BOLD,
    "axes.labelsize":   LABEL_FONT * FONT_SCALE,
    "axes.labelweight": BOLD,
    "axes.titlesize":   TITLE_FONT * FONT_SCALE,
    "axes.titleweight": BOLD,
    "xtick.labelsize":  BASE_FONT * FONT_SCALE,
    "ytick.labelsize":  BASE_FONT * FONT_SCALE,
    "legend.fontsize":  LEGEND_FONT * FONT_SCALE,
    "legend.title_fontsize": LEGEND_FONT * FONT_SCALE,
    "axes.linewidth":   1.4,
    "xtick.major.width": 1.4,
    "ytick.major.width": 1.4,
    "pdf.fonttype":     42,   # keep text editable in Illustrator/Inkscape
    "ps.fonttype":      42,
})

os.makedirs(OUT_DIR, exist_ok=True)

df = pd.read_csv(MATRIX_FILE, sep="\t")
meta = pd.read_csv(META_FILE, sep=None, engine="python")   # handles csv or tsv

GENE = df.columns[0]
CAT = df.columns[1]
IGNORED = df.columns[2]                                    # Origin, dropped
df = df.drop(columns=[IGNORED])
print(f"ignoring third column: {IGNORED}")

samples = [c for c in df.columns if c not in [GENE, CAT]]

# remove NA keyword
df = df[~df[GENE].astype(str).str.contains("NA:Keyword", case=False, na=False)].copy()

# merge exact duplicates (rows that differed only by the ignored Origin column)
df = df.groupby([GENE, CAT], as_index=False)[samples].sum()

# dominant category
df["_tot"] = df[samples].sum(axis=1)
df = (df.sort_values("_tot", ascending=False)
        .drop_duplicates(subset=GENE)
        .drop(columns="_tot"))

# values are already normalised (hits per Gbp) - no relative-abundance conversion
df["mean"] = df[samples].mean(axis=1)
top = df.sort_values("mean", ascending=False).head(TOP_N).copy()

# category order
cat_order = (top.groupby(CAT)[samples].sum().sum(axis=1)
                .sort_values(ascending=False).index.tolist())
palette = dict(zip(cat_order, sns.color_palette("tab20", len(cat_order))))

long = top.melt(id_vars=[GENE, CAT], value_vars=samples,
                var_name="Sample_Date_ID", value_name="abundance")

maps = {
    "stp":        dict(zip(meta["Sample_Date_ID"], meta["stp"])),
    "Month_Name": dict(zip(meta["Sample_Date_ID"], meta["Month_Name"])),
    "season":     dict(zip(meta["Sample_Date_ID"], meta["season"])),
}

# which matrix columns have no metadata row - they would silently become NaN
unmapped = sorted(set(samples) - set(meta["Sample_Date_ID"].astype(str)))
if unmapped:
    print(f"WARNING: {len(unmapped)} sample(s) absent from the metadata and excluded "
          f"from every plot:")
    for s in unmapped[:10]:
        print(f"  - {s}")
    if len(unmapped) > 10:
        print(f"  ... and {len(unmapped) - 10} more")


def bubble(df_long, x, out, title, order=None, size=(14, 12)):
    d = df_long.copy()
    d[x] = d["Sample_Date_ID"].map(maps[x])
    d = d[d[x].notna()]                                    # drop unmapped samples
    n_per_group = d.groupby(x, observed=True)["Sample_Date_ID"].nunique()
    d = d.groupby([GENE, CAT, x], as_index=False)["abundance"].agg(AGG)
    d = d[d["abundance"] > 0]
    if order is not None:
        d[x] = pd.Categorical(d[x], categories=order, ordered=True)
    cats = (d.groupby(CAT, observed=True)["abundance"].sum()
              .sort_values(ascending=False).index.tolist())
    gene_order = []
    for c in cats:
        tmp = (d[d[CAT] == c].groupby(GENE, observed=True)["abundance"].sum()
                 .sort_values(ascending=False))
        gene_order.extend(tmp.index.tolist())
    gene_order = list(dict.fromkeys(gene_order))[::-1]
    d[GENE] = pd.Categorical(d[GENE], categories=gene_order, ordered=True)

    fig, ax = plt.subplots(figsize=size)
    sns.scatterplot(data=d, x=x, y=GENE, hue=CAT, size="abundance",
                    sizes=BUBBLE_SIZES, palette=palette, edgecolor="black",
                    linewidth=.5, legend="brief", ax=ax)
    h, l = ax.get_legend_handles_labels()
    ax.get_legend().remove()
    split = l.index("abundance")
    leg1 = ax.legend(h[1:split], l[1:split], title="Function",
                     bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
                     labelspacing=1.0, borderpad=0.6)
    ax.add_artist(leg1)
    leg2 = ax.legend(h[split + 1:], l[split + 1:], title=VALUE_LABEL,
                     bbox_to_anchor=(1.02, .45), loc="upper left", frameon=False,
                     labelspacing=1.4, borderpad=0.6)
    for leg in (leg1, leg2):
        leg.get_title().set_fontweight(BOLD)
        for txt in leg.get_texts():
            txt.set_fontweight(BOLD)
    ax.set_title(title, fontweight=BOLD, fontsize=TITLE_FONT * FONT_SCALE)
    ax.set_xlabel("")
    ax.set_ylabel("Gene", fontweight=BOLD, fontsize=LABEL_FONT * FONT_SCALE)
    plt.xticks(rotation=45, ha="right")
    plt.setp(ax.get_xticklabels(), fontweight=BOLD)
    plt.setp(ax.get_yticklabels(), fontweight=BOLD, fontstyle="italic")
    ax.grid(axis="y", color="0.9", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    plt.tight_layout()
    # both legends sit outside the axes, so they must be named here or the tight
    # bounding box crops them mid-label
    plt.savefig(os.path.join(OUT_DIR, out), dpi=600, bbox_inches="tight",
                bbox_extra_artists=(leg1, leg2))
    plt.close()

    # ---- the sheet behind this plot
    stem = os.path.splitext(out)[0]
    d_out = d.copy()
    d_out[GENE] = d_out[GENE].astype(str)
    d_out[x] = d_out[x].astype(str)
    d_out = d_out.rename(columns={"abundance": f"{VALUE_LABEL} ({AGG})"})
    d_out.to_csv(os.path.join(OUT_DIR, f"{stem}_data_long.tsv"), sep="\t", index=False)

    wide = (d_out.pivot_table(index=[GENE, CAT], columns=x,
                              values=f"{VALUE_LABEL} ({AGG})", fill_value=0)
                 .reindex(gene_order[::-1], level=0))
    if order is not None:
        wide = wide.reindex(columns=[c for c in order if c in wide.columns])
    wide.to_csv(os.path.join(OUT_DIR, f"{stem}_data_wide.tsv"), sep="\t")

    n_per_group.rename("n_samples").to_csv(
        os.path.join(OUT_DIR, f"{stem}_samples_per_group.tsv"), sep="\t")

    print(f"  {out}: {d[GENE].nunique()} genes x {d[x].nunique()} groups, "
          f"{len(d)} bubbles")


bubble(long, "stp", "Top50_MGE_STP_Bubble.pdf", "Top 50 MGE genes across STPs",
       size=FIG_STP)
bubble(long, "Month_Name", "Top50_MGE_Month_Bubble.pdf", "Top 50 MGE genes across Months",
       ["January", "February", "March", "April", "May", "June", "July", "August",
        "September", "October", "November", "December"], FIG_MONTH)
bubble(long, "season", "Top50_MGE_Season_Bubble.pdf", "Top 50 MGE genes across Seasons",
       ["Winter", "Summer", "Monsoon", "Post-Monsoon"], FIG_SEASON)

# ---- the top-N table, per sample, on the per-Gbp scale
top.to_csv(os.path.join(OUT_DIR, "Top50_MGE_per_Gbp.tsv"), sep="\t", index=False)
long.to_csv(os.path.join(OUT_DIR, "Top50_MGE_per_Gbp_long.tsv"), sep="\t", index=False)

print(f"\nscale: {VALUE_LABEL}; groups combined by {AGG}")
print(f"Done - everything in {os.path.abspath(OUT_DIR)}")
