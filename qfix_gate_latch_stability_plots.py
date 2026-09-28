#!/usr/bin/env python3
"""Gate/Latch RMSD-to-reference and RMSF figures, binder vs nonbinder, qfix cohort.

Companion to pocket_hydration_binder_vs_nonbinder_..._qfix.png: same 95-
sequence ngs_observed cohort, same qfix-corrected feature values the RF
model trains on (unlike analysis/rmsd_to_ref_ngs_observed/'s existing RMSD
figure, which predates the qfix pipeline and uses standard-parameterization
values). Box+jitter style mirrors rmsd_to_ref_significance.py's
box_jitter().

qfix_500ns_feat_table.csv only carries Gate/Latch RMSD-to-reference (not
Lb7a5/Whole/Recoil -- those extra regions were only ever computed in the
standalone rmsd_to_ref_significance.py analysis, never part of the RF's
own feature set, so there's nothing to re-run under qfix for them).

Finding worth knowing about if you re-run this: under qfix, only the Gate
shows a significant RMSD-to-reference separation (p=5.6e-06 mean, p=0.0064
SD) -- the Latch does not (p=0.96 mean, p=0.73 SD). This differs from
analysis/rmsd_to_ref_ngs_observed/'s pre-qfix figure, which found both Gate
(p=0.0014) and Latch (p=0.0037) significant under standard parameterization.
Flagged in Figure 1's title rather than silently reusing the old "gate/latch"
framing; worth reconciling (real effect of the ligand correction vs. an
apples-to-oranges difference in reference/alignment/window between the two
analyses) before treating either as final.

Usage:
    python qfix_gate_latch_stability_plots.py \
        --qfix_table qfix_500ns_feat_table.csv --target_seq_list seq_ids_qfix_all95.txt
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from model_swap_eval import load_baseline_df, load_target_seq_ids

GROUP_COLOR = {"binder": "#648FFF", "nonbinder": "#4D4D4D"}
GROUP_LABEL = {"binder": "Binder", "nonbinder": "Nonbinder"}

OUT_DIR = "analysis/qfix"


def box_jitter(ax, df, col, rng=None):
    """Draws a box-and-jitter plot for one column, Binder vs Nonbinder.

    Args:
        ax: Matplotlib Axes to draw on.
        df: DataFrame containing `col` and a "Group" column ("binder"/"nonbinder").
        col (str): Column to plot.
        rng: numpy Generator for jitter; a fixed-seed default is used if None.
    """
    rng = rng or np.random.default_rng(42)
    xpos = 0
    for group in ("binder", "nonbinder"):
        vals = df.loc[df["Group"] == group, col].dropna().values
        color = GROUP_COLOR[group]
        ax.boxplot(vals, positions=[xpos], widths=0.6, patch_artist=True,
                   medianprops=dict(color="black", linewidth=1.5),
                   boxprops=dict(facecolor=color, alpha=0.5),
                   whiskerprops=dict(color=color), capprops=dict(color=color),
                   flierprops=dict(marker="", linestyle="none"))
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(xpos + jitter, vals, color=color, s=16, alpha=0.85, zorder=3)
        xpos += 1
    ax.set_xticks([0, 1])
    ax.set_xticklabels([GROUP_LABEL["binder"], GROUP_LABEL["nonbinder"]])


def panel(ax, df, col, ylabel):
    """Draws one box_jitter panel plus a Mann-Whitney U p-value in the title.

    Args:
        ax: Matplotlib Axes to draw on.
        df: DataFrame containing `col` and a "Group" column.
        col (str): Column to plot.
        ylabel (str): Y-axis label.

    Returns:
        float: The Mann-Whitney two-sided p-value (Binder vs Nonbinder).
    """
    box_jitter(ax, df, col)
    b = df.loc[df["Group"] == "binder", col].dropna().values
    n = df.loc[df["Group"] == "nonbinder", col].dropna().values
    _, p = mannwhitneyu(b, n, alternative="two-sided")
    ax.set_title(f"{col}\np={p:.2g}", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(True, axis="y", alpha=0.4)
    ax.set_axisbelow(True)
    return p


def build_qfix_cohort(qfix_table, target_seq_list):
    """Builds the fully qfix-corrected cohort dataframe (see qfix_ml_performance_plots.py)."""
    target_seq_ids = load_target_seq_ids(target_seq_list)
    df, feature_group_cols = load_baseline_df()
    feature_cols = [c for cols in feature_group_cols.values() for c in cols]
    qfix = pd.read_csv(qfix_table).set_index("name")
    df_qfix = df.copy()
    for seq_id in target_seq_ids:
        row_idx = df_qfix.index[df_qfix["name"] == seq_id]
        df_qfix.loc[row_idx, feature_cols] = qfix.loc[seq_id, feature_cols].values
    return df_qfix


LATCH_RMSF_COLS = ["Latch (r114-118) mean (A)", "Latch (r114-118) SD (A)"]


def add_latch_rmsf(df_qfix):
    """Merges per-region Latch RMSF into df_qfix by Sequence name.

    Latch RMSF was never part of the ML feature table's "rmsf" group (that
    group only carries Gate mean/SD + 4 individual pocket residues -- see
    ml_utils.build_feature_group_cols), so it isn't in qfix_500ns_feat_table
    or reachable via build_qfix_cohort(). It does exist as a standalone
    per-region summary, computed the same qfix-vs-standard way as
    everything else in this cohort:
    analysis/rmsf_ca_per_seq_summary_500ns_qfix.csv covers 89 of the 95
    target sequences (the qfix reparameterization batch that had completed
    when it was generated); analysis/rmsf_ca_per_seq_summary_500ns.csv (no
    qfix suffix) covers all 95 under standard parameterization. Same
    override pattern as build_qfix_cohort(): start from standard values for
    every sequence, then overwrite with qfix values wherever qfix has them
    -- mirrors qfix_full_performance_eval.py's own "keeps standard values"
    fallback note for cohort sequences outside a --target_seq_list.

    Args:
        df_qfix: DataFrame from build_qfix_cohort(), with a "name" column.

    Returns:
        DataFrame: df_qfix with LATCH_RMSF_COLS added.
    """
    std = pd.read_csv("analysis/rmsf_ca_per_seq_summary_500ns.csv").set_index("Sequence")
    qfix89 = pd.read_csv("analysis/rmsf_ca_per_seq_summary_500ns_qfix.csv").set_index("Sequence")

    out = df_qfix.copy()
    for col in LATCH_RMSF_COLS:
        out[col] = out["name"].map(std[col])
    n_overridden = 0
    for seq_id in qfix89.index:
        row_idx = out.index[out["name"] == seq_id]
        if len(row_idx) == 0:
            continue
        out.loc[row_idx, LATCH_RMSF_COLS] = qfix89.loc[seq_id, LATCH_RMSF_COLS].values
        n_overridden += 1
    print(f"Latch RMSF: {n_overridden} sequences with qfix values, "
          f"{len(out) - n_overridden} kept standard values")
    return out


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

    # ── Figure 1: Gate/Latch RMSD-to-reference (mean + SD), qfix ─────────────
    rmsd_cols = ["Gate RMSD mean (A)", "Latch RMSD mean (A)",
                 "Gate RMSD SD (A)", "Latch RMSD SD (A)"]
    fig1, axes1 = plt.subplots(2, 2, figsize=(8, 8), dpi=150, constrained_layout=True)
    for ax, col in zip(axes1.flat, rmsd_cols):
        panel(ax, df_qfix, col, "RMSD to reference (Å)")
    fig1.suptitle("Binders' GATE loop stays closer to the reference pose -- LATCH does not "
                  "separate (qfix,\n"
                  f"95 ngs_observed sequences, {n_binder} binders / {n_nonbinder} nonbinders)",
                  fontsize=11, fontweight="bold")
    out1 = f"{OUT_DIR}/gate_latch_rmsd_binder_vs_nonbinder_ngs_observed_95samples_qfix.png"
    fig1.savefig(out1, dpi=300, bbox_inches="tight")
    print(f"Saved -> {out1}")

    # ── Figure 2: RMSF, qfix (Gate + Latch mean/SD, + 4 named pocket residues) ─
    df_rmsf = add_latch_rmsf(df_qfix)
    rmsf_cols = ["Gate (r84-90) mean (A)", "Latch (r114-118) mean (A)",
                 "Gate (r84-90) SD (A)", "Latch (r114-118) SD (A)",
                 "Y23 RMSF (A)", "R79 RMSF (A)", "I110 RMSF (A)", "G163 RMSF (A)"]
    fig2, axes2 = plt.subplots(2, 4, figsize=(15, 8), dpi=150, constrained_layout=True)
    for ax, col in zip(axes2.flat, rmsf_cols):
        panel(ax, df_rmsf, col, "RMSF (Å)")
    fig2.suptitle("Gate loop fluctuates less in binders; Latch does not separate (qfix,\n"
                  f"95 ngs_observed sequences, {n_binder} binders / {n_nonbinder} nonbinders)\n"
                  "(Y23/R79/I110/G163 are individual pocket residues, not gate/latch-region "
                  "aggregates; Latch RMSF is a standalone per-region summary, not part of the "
                  "ML \"rmsf\" feature group)",
                  fontsize=10.5, fontweight="bold")
    out2 = f"{OUT_DIR}/gate_rmsf_binder_vs_nonbinder_ngs_observed_95samples_qfix.png"
    fig2.savefig(out2, dpi=300, bbox_inches="tight")
    print(f"Saved -> {out2}")


if __name__ == "__main__":
    main()
