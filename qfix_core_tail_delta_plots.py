#!/usr/bin/env python3
"""Core-vs-tail contact-preference (delta_D/delta_W) figure, qfix cohort.

No existing figure covers the actual "core_tail_delta" ML feature family
(delta = D_core - D_tail / W_core - W_tail, a residue's *preference* for
the ligand's rigid steroid core vs. its flexible C20-C24 carboxylate tail
-- see the core_tail_delta method/reasoning summary from this session).
analysis/core_vs_tail/ and analysis/core_vs_tail_ngs_observed/ test core
and tail *separately* ("is there a Binder vs FP difference in this region
alone", by design never combined -- see core_vs_tail_regions.py's
docstring), which is a different question from whether a residue's core-
vs-tail *preference* differs between binders and nonbinders, and predates
the qfix pipeline besides.

Panels: the 8 delta_D/delta_W columns with the lowest qfix Mann-Whitney
p-value (all p<0.01), out of 36 in the family. Box+jitter style and cohort
construction reuse qfix_gate_latch_stability_plots.py's build_qfix_cohort()
/ box_jitter() / panel().

Usage:
    python qfix_core_tail_delta_plots.py \
        --qfix_table qfix_500ns_feat_table.csv --target_seq_list seq_ids_qfix_all95.txt
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from qfix_gate_latch_stability_plots import build_qfix_cohort, box_jitter, panel
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

    cols = build_feature_group_cols()["core_tail_delta"]
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
        panel(ax, df_qfix, col, "core - tail contact occupancy")
        ax.axhline(0, color="gray", linestyle="--", lw=1)
    fig.suptitle("Residues shift which part of the ligand they grip, binder vs nonbinder "
                 "(qfix,\n"
                 f"95 ngs_observed sequences, {n_binder} binders / {n_nonbinder} nonbinders)\n"
                 "(positive = prefers the rigid steroid core; negative = prefers the "
                 "C20-C24 carboxylate tail; dashed line = no preference)",
                 fontsize=10.5, fontweight="bold")
    out = f"{OUT_DIR}/core_tail_delta_binder_vs_nonbinder_ngs_observed_95samples_qfix.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
