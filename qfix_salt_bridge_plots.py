#!/usr/bin/env python3
"""Salt-bridge occupancy figure, binder vs nonbinder, qfix cohort.

No existing figure covers the "salt_bridge" ML feature family at all (a
repo-wide search found none, qfix or otherwise). Only 3 columns exist in
this family (see ml_utils.build_feature_group_cols), so all 3 are shown,
unlike the other qfix_*_plots.py scripts which select a top-N subset out
of a much larger family.

Box+jitter style and cohort construction reuse
qfix_gate_latch_stability_plots.py's build_qfix_cohort() / panel().

Usage:
    python qfix_salt_bridge_plots.py \
        --qfix_table qfix_500ns_feat_table.csv --target_seq_list seq_ids_qfix_all95.txt
"""
import argparse
import matplotlib.pyplot as plt

from qfix_gate_latch_stability_plots import build_qfix_cohort, panel
from ml_utils import build_feature_group_cols

OUT_DIR = "analysis/qfix"


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

    cols = build_feature_group_cols()["salt_bridge"]
    print(f"salt_bridge family: {cols}")

    fig, axes = plt.subplots(1, 3, figsize=(11, 4.5), dpi=150, constrained_layout=True)
    for ax, col in zip(axes, cols):
        panel(ax, df_qfix, col, "occupancy (%)")
    fig.suptitle("Salt-bridge anchoring of the ligand carboxylate, binder vs nonbinder "
                 "(qfix,\n"
                 f"95 ngs_observed sequences, {n_binder} binders / {n_nonbinder} nonbinders)",
                 fontsize=11, fontweight="bold")
    out = f"{OUT_DIR}/salt_bridge_binder_vs_nonbinder_ngs_observed_95samples_qfix.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
