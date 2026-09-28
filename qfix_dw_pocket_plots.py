#!/usr/bin/env python3
"""Whole-ligand pocket D/W contact-occupancy figure, binder vs nonbinder, qfix cohort.

No existing figure covers the plain whole-ligand "dw_pocket" ML feature
family (D_<resid>/W_<resid>) as a binder-vs-nonbinder comparison --
analysis/core_vs_tail/ and analysis/core_vs_tail_ngs_observed/ only ever
plotted the core-restricted and tail-restricted D/W/R scores separately
(see qfix_core_tail_delta_plots.py's docstring), never the whole-ligand
version this ML feature group actually uses, and predate the qfix
pipeline besides.

Panels: the 8 D_<resid>/W_<resid> columns with the lowest qfix Mann-Whitney
p-value, out of 36 in the family. Box+jitter style and cohort construction
reuse qfix_gate_latch_stability_plots.py's build_qfix_cohort() / panel().

Usage:
    python qfix_dw_pocket_plots.py \
        --qfix_table qfix_500ns_feat_table.csv --target_seq_list seq_ids_qfix_all95.txt
"""
import argparse
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from qfix_gate_latch_stability_plots import build_qfix_cohort, panel
from ml_utils import build_feature_group_cols

OUT_DIR = "analysis/qfix"
N_PANELS = 8


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qfix_table", default="qfix_500ns_feat_table.csv")
    p.add_argument("--target_seq_list", default="seq_ids_qfix_all95.txt")
    args = p.parse_args()

    df_qfix = build_qfix_cohort(args.qfix_table, args.target_seq_list)
    n_binder = int((df_qfix["Group"] == "binder").sum())
    n_nonbinder = int((df_qfix["Group"] == "nonbinder").sum())
    print(f"Cohort: {len(df_qfix)} sequences, {n_binder} binders, {n_nonbinder} nonbinders")

    cols = build_feature_group_cols()["dw_pocket"]
    ranked = []
    for c in cols:
        b = df_qfix.loc[df_qfix["Group"] == "binder", c].dropna().values
        n = df_qfix.loc[df_qfix["Group"] == "nonbinder", c].dropna().values
        if len(b) < 3 or len(n) < 3:
            continue
        _, pval = mannwhitneyu(b, n, alternative="two-sided")
        ranked.append((c, pval))
    ranked.sort(key=lambda r: r[1])
    top_cols = [c for c, _ in ranked[:N_PANELS]]
    print(f"Top {N_PANELS} of {len(ranked)} scoreable columns by p-value: {top_cols}")

    fig, axes = plt.subplots(2, 4, figsize=(15, 8), dpi=150, constrained_layout=True)
    for ax, col in zip(axes.flat, top_cols):
        panel(ax, df_qfix, col, "contact occupancy (fraction of frames)")
    fig.suptitle("Pocket-residue contact with the whole ligand, binder vs nonbinder (qfix,\n"
                 f"95 ngs_observed sequences, {n_binder} binders / {n_nonbinder} nonbinders)\n"
                 "(D = direct heavy-atom contact; W = water-mediated contact; both are "
                 "whole-ligand, not core/tail-restricted)",
                 fontsize=10.5, fontweight="bold")
    out = f"{OUT_DIR}/dw_pocket_binder_vs_nonbinder_ngs_observed_95samples_qfix.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
