#!/usr/bin/env python3
"""Whole-ligand contact-type composition figure, binder vs nonbinder, qfix cohort.

Companion to analysis/core_vs_tail_ngs_observed/contact_features_core_binder_vs_nonbinder.png
and _tail_binder_vs_nonbinder.png -- those two existing figures (pre-qfix,
dated Jul 15) showed the FULL contact-type composition (total contacts,
hydrophobic fraction/count, polar count, positively- and negatively-
charged count) split by ligand region (core/tail), and found strong
significance in several categories -- notably NOT in n_neg_charged, the
only category that made it into the ML "contact_type" feature group.

This script is the qfix equivalent, but WHOLE-LIGAND only, not split by
core/tail region: LIG_contacts/contact_features_all_40_500ns_qfix.csv (the
qfix reparameterization's contact-type output) was only ever computed for
the whole ligand -- there is no
contact_features_all_40_500ns_qfix_core.csv / _tail.csv, so a qfix
core/tail-split version of this figure isn't possible without first
generating that data. Flagged here rather than silently reusing the
pre-qfix core/tail split.

Usage:
    python qfix_contact_type_plots.py \
        --contact_table LIG_contacts/contact_features_all_40_500ns_qfix.csv \
        --target_seq_list seq_ids_qfix_all95.txt
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from qfix_gate_latch_stability_plots import build_qfix_cohort, box_jitter

OUT_DIR = "analysis/qfix"

PANELS = [
    ("mean_n_total", "Total contacts (mean)"),
    ("mean_frac_hydrophobic", "Fraction hydrophobic (mean)"),
    ("mean_n_hydrophobic", "Hydrophobic contacts (mean)"),
    ("mean_n_polar", "Polar contacts (mean)"),
    ("mean_n_pos_charged", "Positively-charged contacts (mean)"),
    ("mean_n_neg_charged", "Negatively-charged contacts (mean)"),
]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--contact_table", default="LIG_contacts/contact_features_all_40_500ns_qfix.csv")
    p.add_argument("--target_seq_list", default="seq_ids_qfix_all95.txt")
    args = p.parse_args()

    # Reuse build_qfix_cohort() purely for its authoritative Group labels
    # (binder/nonbinder) on this exact 95-sequence cohort; the D/W/RMSD/
    # RMSF feature columns it swaps in aren't used here.
    df_qfix = build_qfix_cohort("qfix_500ns_feat_table.csv", args.target_seq_list)
    labels = df_qfix[["name", "Group"]].rename(columns={"name": "seq_id"})

    contact = pd.read_csv(args.contact_table)
    df = labels.merge(contact, on="seq_id", how="inner")
    n_binder = int((df["Group"] == "binder").sum())
    n_nonbinder = int((df["Group"] == "nonbinder").sum())
    print(f"Cohort: {len(df)} sequences, {n_binder} binders, {n_nonbinder} nonbinders")

    fig, axes = plt.subplots(2, 3, figsize=(11.5, 8), dpi=150, constrained_layout=True)
    for ax, (col, title) in zip(axes.flat, PANELS):
        box_jitter(ax, df, col)
        b = df.loc[df["Group"] == "binder", col].dropna().values
        n = df.loc[df["Group"] == "nonbinder", col].dropna().values
        _, pval = mannwhitneyu(b, n, alternative="two-sided")
        ax.set_title(f"{title}\np={pval:.2g}", fontsize=10)
        ax.set_ylabel(col, fontsize=8)
        ax.grid(True, axis="y", alpha=0.4)
        ax.set_axisbelow(True)
    fig.suptitle("Whole-ligand contact-type composition, binder vs nonbinder (qfix,\n"
                 f"95 ngs_observed sequences, {n_binder} binders / {n_nonbinder} nonbinders)\n"
                 "(whole-ligand only -- qfix contact-type data was never split into "
                 "core/tail regions, unlike the pre-qfix comparison figure)",
                 fontsize=10.5, fontweight="bold")
    out = f"{OUT_DIR}/contact_type_binder_vs_nonbinder_ngs_observed_95samples_qfix.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
