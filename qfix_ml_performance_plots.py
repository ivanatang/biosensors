#!/usr/bin/env python3
"""Generates RF performance-visualization plots for the qfix-feature cohort.

Mirrors ML_classification.ipynb's binary-classification figures (CV metric
summary, confusion matrix, feature importances, ROC curve, permutation
test, and the precision-at-top-N cell), applied to the fully qfix-corrected
95-sequence cohort instead of the standard-parameterization one that
notebook trains on. The permutation-test and ROC-curve figures are computed
in that notebook but never saved to disk (only shown inline) -- they're
included here since they're directly relevant to "model performance."

Saved under analysis/ML/'s existing filt_RF_*_{SEQ_SOURCE}_{n}samples.png
naming convention, with a _qfix suffix so these sit alongside (without
overwriting) the existing standard-method plots for the same cohort size.

Usage:
    python qfix_ml_performance_plots.py --qfix_table qfix_500ns_feat_table.csv \
        --target_seq_list seq_ids_qfix_all95.txt
"""
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict, cross_validate
from sklearn.metrics import (confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc,
                              classification_report)
from sklearn.inspection import permutation_importance

from model_swap_eval import load_baseline_df, load_target_seq_ids
from ml_utils import make_pipeline, grouped_permutation_test, SCORING, N_SPLITS, RANDOM_STATE
from seq_utils import sequence_similarity_groups


def build_qfix_cohort(qfix_table, target_seq_list):
    """Builds the fully qfix-corrected cohort dataframe.

    Args:
        qfix_table (str): Path to the qfix feature table CSV.
        target_seq_list (str): seq_ids.txt-style list of sequences to swap
            qfix values in for.

    Returns:
        tuple: (df, feature_cols, feature_group_cols) -- the swapped-in
        DataFrame, flat feature column list, and group -> columns mapping.
    """
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
    """Runs the CV/permutation/full-data-fit pipeline and saves all figures."""
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qfix_table", default="qfix_500ns_feat_table.csv")
    p.add_argument("--target_seq_list", default="seq_ids_qfix_all95.txt")
    p.add_argument("--out_dir", default="analysis/ML")
    p.add_argument("--seq_source_tag", default="ngs_observed",
                    help="Matches ML_classification.ipynb's SEQ_SOURCE tag in the filename")
    p.add_argument("--suffix", default="_qfix", help="Filename suffix distinguishing these from the standard-method plots")
    p.add_argument("--n_permutations", type=int, default=100)
    args = p.parse_args()

    df, feature_cols, feature_group_cols = build_qfix_cohort(args.qfix_table, args.target_seq_list)
    X = df[feature_cols].values
    y = df["Label"].values
    seq_ids = df["name"].tolist()
    n_binders = int(y.sum())
    n_nonbinders = int((y == 0).sum())
    print(f"Qfix-corrected cohort: {len(y)} samples (Binders: {n_binders} | Nonbinders: {n_nonbinders})")

    sequences = df["Sequence"].astype(str).str.strip().tolist()
    seq_groups = sequence_similarity_groups(sequences, identity_threshold=0.95)
    cv = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    pipe = make_pipeline(feature_group_cols, feature_cols)

    cv_results = cross_validate(pipe, X, y, cv=cv, groups=seq_groups, scoring=SCORING,
                                 return_train_score=False)
    metrics = {k.replace("test_", ""): v for k, v in cv_results.items() if k.startswith("test_")}
    print(f"\n--- {N_SPLITS}-fold grouped stratified CV (qfix features) ---")
    for k, v in metrics.items():
        print(f"  {k:<20} {v.mean():.3f} +/- {v.std():.3f}")

    score, perm_scores, pval = grouped_permutation_test(
        pipe, X, y, seq_groups, cv, scoring="balanced_accuracy",
        n_permutations=args.n_permutations, random_state=RANDOM_STATE,
    )
    print(f"\nPermutation test | Observed balanced accuracy: {score:.3f} | p-value: {pval:.3f}")

    # Full-data fit for importance/interpretation only (disclosed as such,
    # not a generalization estimate -- see ROC curve title below).
    pipe.fit(X, y)
    rf = pipe.named_steps["rf"]
    selected_idx = pipe.named_steps["select"].selected_indices_
    selected_cols = [feature_cols[i] for i in selected_idx]
    print(f"\nFull-data fit selected {len(selected_cols)} features: {selected_cols}")
    y_pred = pipe.predict(X)
    y_prob = pipe.predict_proba(X)[:, 1]

    tag = f"{args.seq_source_tag}_{len(y)}samples{args.suffix}"
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    # ── Figure 1: CV metric summary ──────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(8, 4), dpi=150)
    metric_names = list(metrics.keys())
    means = [metrics[m].mean() for m in metric_names]
    sds = [metrics[m].std() for m in metric_names]
    x = np.arange(len(metric_names))
    ax1.bar(x, means, yerr=sds, capsize=5,
            color=plt.get_cmap("tab10")(np.linspace(0, 0.6, len(metric_names))),
            edgecolor="white")
    ax1.set_xticks(x)
    ax1.set_xticklabels(metric_names, rotation=20, ha="right", fontsize=9)
    ax1.set_ylabel("Score", fontsize=10)
    ax1.set_ylim(0, 1.15)
    ax1.axhline(0.5, color="gray", linestyle="--", lw=1, label="chance")
    ax1.set_title(f"{N_SPLITS}-fold grouped CV performance -- balanced RF (mean +/- SD)\n"
                  f"{len(y)} samples (Binders: {n_binders}  |  Nonbinders: {n_nonbinders})  --  qfix features",
                  fontsize=11, fontweight="bold")
    ax1.legend(fontsize=8)
    ax1.grid(True, axis="y", alpha=0.35)
    ax1.set_axisbelow(True)
    plt.tight_layout()
    fig1.savefig(os.path.join(out_dir, f"filt_RF_performance_metrics_{tag}.png"),
                 dpi=300, bbox_inches="tight")

    # ── Figure 2: permutation test ───────────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(7, 4), dpi=150)
    ax2.hist(perm_scores, bins=25, color="steelblue", alpha=0.7, edgecolor="white",
             label="Permuted scores")
    ax2.axvline(score, color="red", lw=2, label=f"Observed: {score:.3f}  (p={pval:.3f})")
    ax2.set_xlabel("Balanced accuracy", fontsize=10)
    ax2.set_ylabel("Count", fontsize=10)
    ax2.set_title("Permutation test -- qfix features", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.35)
    ax2.set_axisbelow(True)
    plt.tight_layout()
    fig2.savefig(os.path.join(out_dir, f"filt_RF_permutation_test_{tag}.png"),
                 dpi=300, bbox_inches="tight")

    # ── Figure 3: confusion matrix (grouped CV, out-of-fold predictions) ────
    y_pred_cv = cross_val_predict(pipe, X, y, cv=cv, groups=seq_groups)
    fig3, ax3 = plt.subplots(figsize=(7, 3), dpi=300)
    cm_mat = confusion_matrix(y, y_pred_cv)
    disp = ConfusionMatrixDisplay(cm_mat, display_labels=["Nonbinder", "Binder"])
    disp.plot(ax=ax3, colorbar=False, cmap="Blues")
    ax3.set_title(f"Confusion matrix ({N_SPLITS}-fold grouped CV, OOF predictions) -- qfix features",
                  fontsize=11, fontweight="bold")
    plt.tight_layout()
    fig3.savefig(os.path.join(out_dir, f"filt_RF_confusion_matrix_{tag}.png"),
                 dpi=300, bbox_inches="tight")

    # ── Figure 4: feature importances (MDI + permutation) ───────────────────
    mdi_imp = rf.feature_importances_
    perm_imp = permutation_importance(pipe, X, y, n_repeats=30, random_state=RANDOM_STATE, n_jobs=1)
    TOP_N = 20
    mdi_sorted_idx = np.argsort(mdi_imp)[-min(TOP_N, len(mdi_imp)):]
    mdi_names = [selected_cols[i] for i in mdi_sorted_idx]
    perm_sorted_idx = np.argsort(perm_imp.importances_mean)[-TOP_N:]
    perm_names = [feature_cols[i] for i in perm_sorted_idx]
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(14, 6), dpi=300, constrained_layout=True)
    colors_a = [prop_cycle[i % len(prop_cycle)] for i in range(len(mdi_sorted_idx))]
    ax4a.barh(range(len(mdi_sorted_idx)), mdi_imp[mdi_sorted_idx], color=colors_a, alpha=0.85)
    ax4a.set_yticks(range(len(mdi_sorted_idx)))
    ax4a.set_yticklabels(mdi_names, fontsize=9)
    ax4a.set_xlabel("Mean decrease in impurity", fontsize=10)
    ax4a.set_title(f"MDI -- top {len(mdi_sorted_idx)} of {len(selected_cols)} selected features (qfix)",
                   fontsize=11, fontweight="bold")
    ax4a.grid(True, axis="x", alpha=0.35)
    ax4a.set_axisbelow(True)

    colors_b = [prop_cycle[i % len(prop_cycle)] for i in range(len(perm_sorted_idx))]
    ax4b.barh(range(len(perm_sorted_idx)), perm_imp.importances_mean[perm_sorted_idx],
              xerr=perm_imp.importances_std[perm_sorted_idx], color=colors_b, alpha=0.85)
    ax4b.set_yticks(range(len(perm_sorted_idx)))
    ax4b.set_yticklabels(perm_names, fontsize=9)
    ax4b.set_xlabel("Mean accuracy decrease (permutation)", fontsize=10)
    ax4b.set_title(f"Permutation importance -- top {TOP_N} of all {len(feature_cols)} candidates (qfix)",
                   fontsize=11, fontweight="bold")
    ax4b.grid(True, axis="x", alpha=0.35)
    ax4b.set_axisbelow(True)
    fig4.savefig(os.path.join(out_dir, f"filt_RF_feat_importance_{tag}.png"),
                 dpi=300, bbox_inches="tight")

    # ── Figure 5: ROC curve ───────────────────────────────────────────────
    fpr, tpr, _ = roc_curve(y, y_prob)
    roc_auc_val = auc(fpr, tpr)
    fig5, ax5 = plt.subplots(figsize=(5, 5), dpi=150)
    ax5.plot(fpr, tpr, lw=2, color="steelblue", label=f"AUC = {roc_auc_val:.3f}")
    ax5.plot([0, 1], [0, 1], color="gray", linestyle="--", lw=1)
    ax5.set_xlabel("False positive rate", fontsize=10)
    ax5.set_ylabel("True positive rate", fontsize=10)
    ax5.set_title("ROC curve (full-data fit) -- qfix features", fontsize=11, fontweight="bold")
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.35)
    ax5.set_axisbelow(True)
    plt.tight_layout()
    fig5.savefig(os.path.join(out_dir, f"filt_RF_roc_curve_{tag}.png"),
                 dpi=300, bbox_inches="tight")

    # ── Classification report + incorrect predictions (printed, not plotted) ─
    print(f"\n--- Classification report ({N_SPLITS}-fold grouped CV, OOF predictions) ---")
    print(classification_report(y, y_pred_cv, target_names=["Nonbinder", "Binder"], zero_division=0))

    wrong_mask = y_pred_cv != y
    wrong_df = pd.DataFrame({
        "seq_id": [seq_ids[i] for i in np.where(wrong_mask)[0]],
        "true_label": ["Binder" if y[i] == 1 else "Nonbinder" for i in np.where(wrong_mask)[0]],
        "predicted_label": ["Binder" if y_pred_cv[i] == 1 else "Nonbinder" for i in np.where(wrong_mask)[0]],
        "error_type": ["False Negative" if y[i] == 1 else "False Positive" for i in np.where(wrong_mask)[0]],
    }).sort_values("error_type").reset_index(drop=True)
    print(f"\n--- Incorrectly predicted sequences ({len(wrong_df)} total) ---")
    print(wrong_df.to_string(index=False))

    # ── Precision-at-top-N ────────────────────────────────────────────────
    y_prob_cv = cross_val_predict(pipe, X, y, cv=cv, groups=seq_groups, method="predict_proba")[:, 1]
    order = np.argsort(-y_prob_cv)
    y_sorted = y[order]
    Ns = np.arange(1, len(y) + 1)
    precision_at_n = np.cumsum(y_sorted) / Ns
    base_rate = y.mean()

    fig_prec, ax_prec = plt.subplots(figsize=(8, 5), dpi=150, constrained_layout=True)
    ax_prec.plot(Ns, precision_at_n, color="#648FFF", lw=2, label="Precision@N (OOF, grouped CV)")
    ax_prec.axhline(base_rate, color="gray", linestyle="--", lw=1,
                     label=f"Base rate (random guessing) = {base_rate:.3f}")
    ax_prec.set_xlabel("N (top-N most-confident predicted binders)", fontsize=10)
    ax_prec.set_ylabel("Precision@N (fraction of top N that are true binders)", fontsize=10)
    ax_prec.set_title(f"Precision-at-top-N -- {N_SPLITS}-fold grouped CV, OOF predictions -- qfix features\n"
                       f"{len(y)} samples (Binders: {n_binders}  |  Nonbinders: {n_nonbinders})",
                       fontsize=11, fontweight="bold")
    ax_prec.set_ylim(0, 1.05)
    ax_prec.set_xlim(1, len(y))
    ax_prec.legend(fontsize=9)
    ax_prec.grid(True, alpha=0.4)
    ax_prec.set_axisbelow(True)
    fig_prec.savefig(os.path.join(out_dir, f"filt_RF_precision_at_topN_{tag}.png"),
                      dpi=300, bbox_inches="tight")

    practical_Ns = [n for n in [5, 10, 15, 20, 25, 30, 40, 50, 75, 100] if n <= len(y)]
    table_rows = [{"N": n, "true_binders_in_top_N": int(y_sorted[:n].sum()),
                   "precision_at_N": round(precision_at_n[n - 1], 3)} for n in practical_Ns]
    print("\n--- Precision at practical N (qfix features) ---")
    print(pd.DataFrame(table_rows).to_string(index=False))
    print(f"\nFor reference, base rate (random guessing) = {base_rate:.3f}")

    print(f"\nSaved 6 figures to {out_dir}/ with tag '{tag}':")
    for name in ("performance_metrics", "permutation_test", "confusion_matrix",
                 "feat_importance", "roc_curve", "precision_at_topN"):
        print(f"  filt_RF_{name}_{tag}.png")


if __name__ == "__main__":
    main()
