#!/usr/bin/env python3
"""Binder vs nonbinder pocket-hydration figure, qfix 95-sequence cohort.

No existing figure covered this (repo-wide search found no hydration_pocket
plots), so this builds one from the same qfix-corrected dataframe
qfix_full_performance_eval.py evaluates the RF on. Box+jitter style mirrors
rmsd_to_ref/rmsd_to_ref_significance.py's box_jitter(), same as
qfix_gate_latch_stability_plots.py (whose build_qfix_cohort() this reuses).

Two panels: hydration_count_pocket_4A_mean (the whole-trajectory average,
most directly interpretable as "how wet is the pocket") and
hydration_count_pocket_4A_early20_mean (the single most SHAP-important
feature of any family, see earlier session analysis).

Usage:
    python qfix_pocket_hydration_plots.py \
        --qfix_table qfix_500ns_feat_table.csv --target_seq_list seq_ids_qfix_all95.txt
"""
import argparse
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from qfix_gate_latch_stability_plots import build_qfix_cohort, box_jitter

OUT_DIR = "analysis/qfix"

PANELS = [
    ("hydration_count_pocket_4A_mean", "Whole-trajectory mean pocket hydration"),
    ("hydration_count_pocket_4A_early20_mean", "Early-run (first 20%) mean pocket hydration"),
]


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

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), dpi=150, constrained_layout=True)
    for ax, (col, title) in zip(axes, PANELS):
        box_jitter(ax, df_qfix, col)
        b = df_qfix.loc[df_qfix["Group"] == "binder", col].dropna().values
        n = df_qfix.loc[df_qfix["Group"] == "nonbinder", col].dropna().values
        _, pval = mannwhitneyu(b, n, alternative="two-sided")
        ax.set_title(f"{title}\nMann-Whitney p={pval:.2g}", fontsize=10)
        ax.set_ylabel("Water-oxygen count within 4 Å of pocket", fontsize=9)
        ax.grid(True, axis="y", alpha=0.4)
        ax.set_axisbelow(True)

    fig.suptitle(f"Binders keep a drier pocket than nonbinders (qfix, 95 ngs_observed sequences,\n"
                 f"{n_binder} binders / {n_nonbinder} nonbinders)", fontsize=11, fontweight="bold")

    out_path = f"{OUT_DIR}/pocket_hydration_binder_vs_nonbinder_ngs_observed_95samples_qfix.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved -> {out_path}")


if __name__ == "__main__":
    main()
