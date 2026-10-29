#!/usr/bin/env python3

from __future__ import annotations

import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# -------- USER INPUTS --------
MATRIX_FILE = "AMR-Matrix.tsv"
META_FILE   = "metadata.csv"
OUT_DIR     = "heatmaps"
TOP_N       = 50
VALUE_LABEL = "ARG hits per Gbp"

DROP_CLASSES = ["unclassified"]   # classes left out entirely; [] keeps all
RANK_BY      = "mean"      # "mean" | "prevalence" | "max"
ROW_ORDER    = "abundance" # "abundance"  most abundant at the top, least at the
                           #               bottom, ignoring class
                           # "class"      group rows by ARG_class, ranked within
LOG_SCALE    = True        # log10(value + 1); linear if False

SHOW_SAMPLE_LABELS = True  # the 132 tiny x labels; False leaves only STP blocks
SHORTEN_SAMPLE_LABELS = True   # PR-255-29-Dec-2024 -> PR-255-Dec-2024 (display only)

# Force specific x-axis labels. Keys may be the full name or the shortened one;
# this relabels the tick only, it does not move or alter any data, and the
# exported sheets keep the original column names.
RENAME_SAMPLES = {
    "PR-298-26-Aug-2024": "PR-293-Jul-2024",
    "PR-298-Aug-2024":    "PR-293-Jul-2024",
}
FIG_W, FIG_H = 26, 15
SAMPLE_FONT  = 5.5
GENE_FONT    = 9
GROUP_FONT   = 11
# -----------------------------

def short_sample(name: str) -> str:
    """Drop the day from a sample label: PR-255-29-Dec-2024 -> PR-255-Dec-2024.
    Anything that does not match the pattern is left exactly as it is."""
    m = re.match(r"^([A-Za-z]+-\d+)-\d{1,2}-([A-Za-z]{3,}-\d{4})$", str(name).strip())
    return f"{m.group(1)}-{m.group(2)}" if m else str(name)


CLASS_PAL = ["#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a", "#b15928",
             "#a6cee3", "#fb9a99", "#b2df8a", "#fdbf6f", "#cab2d6", "#8dd3c7",
             "#bebada", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f", "#80b1d3"]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # keep_default_na=False so a literal "NA" gene name stays a string
    mat = pd.read_csv(MATRIX_FILE, sep="\t", keep_default_na=False)
    gene_col, class_col = mat.columns[0], mat.columns[1]
    samples = [c for c in mat.columns if c not in (gene_col, class_col)]
    vals = mat[samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat = pd.concat([mat[[gene_col, class_col]], vals], axis=1)
    print(f"matrix: {len(mat)} genes x {len(samples)} samples")

    if DROP_CLASSES:
        drop = {c.strip().lower() for c in DROP_CLASSES}
        m = mat[class_col].astype(str).str.strip().str.lower().isin(drop)
        if m.any():
            print(f"dropping {int(m.sum())} row(s) in {sorted(drop)}")
            mat = mat[~m]

    # one row per gene: keep the class carrying the most signal
    mat["_tot"] = mat[samples].sum(axis=1)
    before = len(mat)
    mat = (mat.sort_values("_tot", ascending=False)
              .drop_duplicates(subset=gene_col)
              .drop(columns="_tot"))
    if len(mat) < before:
        print(f"collapsed {before - len(mat)} gene(s) listed under several classes")

    score = {"mean": mat[samples].mean(axis=1),
             "max": mat[samples].max(axis=1),
             "prevalence": (mat[samples] > 0).mean(axis=1)}[RANK_BY]
    mat = mat.assign(_score=score).sort_values("_score", ascending=False)
    top = mat.head(TOP_N).copy()
    print(f"top {TOP_N} by {RANK_BY}, highest first:")
    for _, r in top.head(3).iterrows():
        print(f"    {r[gene_col]:<28} {r['_score']:.3f}")
    print("    ...")
    for _, r in top.tail(2).iterrows():
        print(f"    {r[gene_col]:<28} {r['_score']:.3f}")

    if ROW_ORDER == "class":
        class_rank = (top.groupby(class_col)["_score"].sum()
                         .sort_values(ascending=False))
        top["_crank"] = top[class_col].map({c: i for i, c in
                                            enumerate(class_rank.index)})
        top = top.sort_values(["_crank", "_score"], ascending=[True, False])
    top = top.drop(columns=[c for c in ("_score", "_crank") if c in top.columns])

    # ---- sample order
    # Only Sample_Date_ID and stp are read from the metadata. Order within an
    # STP comes from the date embedded in Sample_Date_ID, so no Day/Month
    # columns are needed. Samples with no metadata row, or with a blank stp,
    # are dropped rather than collected into an extra block.
    group_blocks, ordered = None, samples
    if META_FILE and os.path.exists(META_FILE):
        meta = pd.read_csv(META_FILE, sep=None, engine="python", dtype=str)
        meta.columns = [c.strip() for c in meta.columns]
        need = {"Sample_Date_ID", "stp"}
        if not need <= set(meta.columns):
            sys.exit(f"metadata needs columns {sorted(need)}; found "
                     f"{list(meta.columns)}")

        meta = meta[["Sample_Date_ID", "stp"]].copy()
        # fillna before astype: pandas 3.0 keeps NA through astype(str) rather
        # than turning it into the string "nan", so a blank stp would survive
        for c in ("Sample_Date_ID", "stp"):
            meta[c] = meta[c].fillna("").astype(str).str.strip()
        blank = meta["stp"].isin(["", "nan", "NaN", "None", "NA", "-"])
        if blank.any():
            print(f"{int(blank.sum())} metadata row(s) have a blank stp and are "
                  f"ignored: {list(meta.loc[blank, 'Sample_Date_ID'])[:5]}")
            meta = meta[~blank]
        meta = meta.drop_duplicates(subset="Sample_Date_ID", keep="first")

        lookup = meta.set_index("Sample_Date_ID")["stp"]
        present = [s for s in samples if s in lookup.index]
        missing = [s for s in samples if s not in lookup.index]
        if missing:
            print(f"DROPPED {len(missing)} sample(s) with no metadata row: "
                  f"{missing[:8]}{' ...' if len(missing) > 8 else ''}")
        if not present:
            sys.exit("no matrix column matches Sample_Date_ID in the metadata")

        key = pd.DataFrame({"sample": present,
                            "stp": lookup.reindex(present).values})
        key["date"] = pd.to_datetime(
            key["sample"].str.extract(r"(\d{1,2}-[A-Za-z]{3}-\d{4})$",
                                      expand=False),
            format="%d-%b-%Y", errors="coerce")
        undated = int(key["date"].isna().sum())
        if undated:
            print(f"{undated} sample name(s) carry no readable date; placed last "
                  f"within their STP")
        key = key.sort_values(["stp", "date", "sample"],
                              na_position="last").reset_index(drop=True)
        ordered = list(key["sample"])

        blocks, start = [], 0
        stps = list(key["stp"])
        for i in range(1, len(stps) + 1):
            if i == len(stps) or stps[i] != stps[start]:
                blocks.append((stps[start], start, i))
                start = i
        group_blocks = blocks
        print(f"{len(ordered)} sample(s) in {len(blocks)} STP block(s): "
              f"{[b[0] for b in blocks]}")
    else:
        print("no metadata file found; leaving sample order as is")

    data = top[ordered].to_numpy(dtype=float)
    plot = np.log10(data + 1) if LOG_SCALE else data
    cbar_label = (f"log10({VALUE_LABEL} + 1)" if LOG_SCALE else VALUE_LABEL)

    # ---- figure
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    # the class strip sits immediately right of the cells: on the left it would
    # be separated from them by the gene labels and impossible to read against
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 0.02, 0.16], wspace=0.012)
    ax = fig.add_subplot(gs[0, 0])
    ax_cls = fig.add_subplot(gs[0, 1])
    ax_leg = fig.add_subplot(gs[0, 2]); ax_leg.axis("off")

    im = ax.imshow(plot, aspect="auto", cmap="viridis", interpolation="nearest")

    ax.set_xticks(np.arange(len(ordered)))
    if SHOW_SAMPLE_LABELS:
        xlabels = []
        renamed = []
        for s in ordered:
            lab = short_sample(s) if SHORTEN_SAMPLE_LABELS else s
            if s in RENAME_SAMPLES:
                lab = RENAME_SAMPLES[s]
            elif lab in RENAME_SAMPLES:
                lab = RENAME_SAMPLES[lab]
            if lab != (short_sample(s) if SHORTEN_SAMPLE_LABELS else s):
                renamed.append((s, lab))
            xlabels.append(lab)
        for orig, lab in renamed:
            print(f"relabelled on the x axis: {orig} -> {lab}")
        dup = {l for l in xlabels if xlabels.count(l) > 1}
        if dup:
            print(f"WARNING: {len(dup)} x-axis label(s) now appear more than once, "
                  f"so two columns look identical: {sorted(dup)}")
        ax.set_xticklabels(xlabels, rotation=90, fontsize=SAMPLE_FONT,
                           fontweight="bold")
    else:
        ax.set_xticklabels([])
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top[gene_col], fontsize=GENE_FONT, fontweight="bold")
    ax.tick_params(axis="both", length=2, pad=3)
    ax.set_ylabel("ARG", fontsize=13, fontweight="bold", labelpad=10)

    # faint cell separators, as in your reference
    ax.set_xticks(np.arange(-0.5, len(ordered), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(top), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.25)
    ax.tick_params(which="minor", length=0)

    # STP blocks: white separators plus one label per block
    if group_blocks:
        for name, s, e in group_blocks[:-1]:
            ax.axvline(e - 0.5, color="white", linewidth=2.2)
        for name, s, e in group_blocks:
            ax.text((s + e - 1) / 2, -1.2, name, ha="center", va="bottom",
                    fontsize=GROUP_FONT, fontweight="bold", rotation=0,
                    clip_on=False)

    # ARG_class strip down the left
    classes = list(dict.fromkeys(top[class_col]))
    cls_col = {c: CLASS_PAL[i % len(CLASS_PAL)] for i, c in enumerate(classes)}
    idx = np.array([[classes.index(c)] for c in top[class_col]])
    ax_cls.imshow(idx, aspect="auto", interpolation="nearest",
                  cmap=ListedColormap([cls_col[c] for c in classes]))
    ax_cls.set_xticks([]); ax_cls.set_yticks([])
    for sp in ax_cls.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.6)

    handles = [Patch(color=cls_col[c], label=c) for c in classes]
    ax_leg.legend(handles=handles, title="ARG class", loc="upper left",
                  frameon=False, fontsize=9,
                  title_fontproperties=dict(weight="bold", size=11),
                  prop=dict(weight="bold", size=9))

    cb = fig.colorbar(im, ax=ax_leg, fraction=0.35, pad=0.02, aspect=28)
    cb.set_label(cbar_label, fontsize=11, fontweight="bold")
    cb.ax.tick_params(labelsize=9)
    for t in cb.ax.get_yticklabels():
        t.set_fontweight("bold")

    title = f"Top {TOP_N} ARGs - sample wise"
    fig.suptitle(title, fontsize=18, fontweight="bold", y=0.965)

    out_pdf = os.path.join(OUT_DIR, f"Top{TOP_N}_ARG_heatmap.pdf")
    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ---- the sheet behind the plot
    sheet = top[[gene_col, class_col] + ordered]
    sheet.to_csv(os.path.join(OUT_DIR, f"Top{TOP_N}_ARG_heatmap_data.tsv"),
                 sep="\t", index=False)
    plotted = pd.DataFrame(plot, index=top[gene_col], columns=ordered)
    plotted.round(4).to_csv(
        os.path.join(OUT_DIR, f"Top{TOP_N}_ARG_heatmap_plotted_values.tsv"), sep="\t")

    print(f"\n{len(top)} genes x {len(ordered)} samples, "
          f"{'log10' if LOG_SCALE else 'linear'} scale")
    print(f"value range: {plot.min():.3f} - {plot.max():.3f}")
    print(f"zero cells: {(data == 0).sum()} of {data.size} "
          f"({100 * (data == 0).sum() / data.size:.0f}%)")
    print(f"wrote {out_pdf} (+ .png and 2 sheets) to {os.path.abspath(OUT_DIR)}")


if __name__ == "__main__":
    main()
