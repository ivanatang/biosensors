#!/usr/bin/env python3
"""Generates a SHAP summary plot for the qfix-feature RF classifier.

Fits the same scale -> GroupAwareSelector -> RandomForest pipeline used
throughout this project on the full qfix-corrected 95-sequence cohort
(disclosed as a full-data fit for interpretation, same caveat as the MDI/
permutation-importance figures), then explains it with a SHAP
TreeExplainer on the post-selection feature space (the RF never sees the
pre-selection columns, so SHAP must be computed in the same space MDI is).

Deliberately run in an isolated virtualenv (not the `biosensors` conda
env): shap>=0.51 hard-requires numpy>=2, which would silently break this
project's numpy==1.26.4/MDAnalysis==2.9.0 ABI pin if installed there (see
CLAUDE.md). This script only needs pandas/numpy/sklearn/scipy/shap/
matplotlib -- none of it touches MDAnalysis -- so it has no reason to
share an environment with the trajectory-analysis tooling.

Usage:
    python qfix_shap_summary.py --qfix_table qfix_500ns_feat_table.csv \
        --target_seq_list seq_ids_qfix_all95.txt
"""
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap

from model_swap_eval import load_baseline_df, load_target_seq_ids
from ml_utils import make_pipeline, RANDOM_STATE


def build_qfix_cohort(qfix_table, target_seq_list):
    """Builds the fully qfix-corrected cohort dataframe (see qfix_ml_performance_plots.py)."""
    target_seq_ids = load_target_seq_ids(target_seq_list)
    df, feature_group_cols = load_baseline_df()
    feature_cols = [c for cols in feature_group_cols.values() for c in cols]
    missing = set(target_seq_ids) - set(df["name"])
    if missing:
        raise ValueError(f"Target sequences not found in baseline cohort: {missing}")
    qfix = pd.read_csv(qfix_table).set_index("name")
    missing_cols = set(feature_cols) - set(qfix.columns)
    if missing_cols:
        raise ValueError(f"{qfix_table} is missing feature columns: {missing_cols}")

    df_qfix = df.copy()
    for seq_id in target_seq_ids:
        row_idx = df_qfix.index[df_qfix["name"] == seq_id]
        df_qfix.loc[row_idx, feature_cols] = qfix.loc[seq_id, feature_cols].values
    return df_qfix, feature_cols, feature_group_cols


def main():
    """Fits the full-data RF pipeline, computes SHAP values, and saves the summary plot."""
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qfix_table", default="qfix_500ns_feat_table.csv")
    p.add_argument("--target_seq_list", default="seq_ids_qfix_all95.txt")
    p.add_argument("--out_dir", default="analysis/ML")
    p.add_argument("--seq_source_tag", default="ngs_observed")
    p.add_argument("--suffix", default="_qfix")
    p.add_argument("--max_display", type=int, default=20)
    p.add_argument("--top_n_report", type=int, default=5)
    args = p.parse_args()

    df, feature_cols, feature_group_cols = build_qfix_cohort(args.qfix_table, args.target_seq_list)
    X = df[feature_cols].values
    y = df["Label"].values
    print(f"Qfix-corrected cohort: {len(y)} samples (Binders: {int(y.sum())} | "
          f"Nonbinders: {int((y == 0).sum())})")

    pipe = make_pipeline(feature_group_cols, feature_cols)
    pipe.fit(X, y)
    rf = pipe.named_steps["rf"]
    selected_idx = pipe.named_steps["select"].selected_indices_
    selected_cols = [feature_cols[i] for i in selected_idx]
    print(f"Full-data fit selected {len(selected_cols)} features")

    # SHAP must see exactly what the RF sees: scaled, then post-selection columns.
    X_scaled = pipe.named_steps["scale"].transform(X)
    X_selected = pipe.named_steps["select"].transform(X_scaled)

    explainer = shap.TreeExplainer(rf)
    sv = explainer.shap_values(X_selected)
    if isinstance(sv, list):
        sv_binder = sv[1]
    elif np.asarray(sv).ndim == 3:
        sv_binder = np.asarray(sv)[:, :, 1]
    else:
        sv_binder = np.asarray(sv)

    mean_abs_shap = np.abs(sv_binder).mean(axis=0)
    rank_df = pd.DataFrame({"feature": selected_cols, "mean_abs_shap": mean_abs_shap}) \
        .sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    print(f"\nTop {args.top_n_report} features by mean |SHAP value| (impact on P(Binder)):")
    print(rank_df.head(args.top_n_report).to_string(index=False))

    tag = f"{args.seq_source_tag}_{len(y)}samples{args.suffix}"
    os.makedirs(args.out_dir, exist_ok=True)

    shap.summary_plot(sv_binder, X_selected, feature_names=selected_cols,
                       max_display=args.max_display, show=False)
    fig = plt.gcf()
    fig.suptitle(f"SHAP summary -- qfix features -- {len(y)} samples", fontsize=11, fontweight="bold", y=1.02)
    out_path = os.path.join(args.out_dir, f"filtered_DW_RF_shap_summary_{tag}.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"\nSaved -> {out_path}")

    rank_csv = os.path.join(args.out_dir, f"filtered_DW_RF_shap_ranking_{tag}.csv")
    rank_df.to_csv(rank_csv, index=False)
    print(f"Saved -> {rank_csv}")


if __name__ == "__main__":
    main()
