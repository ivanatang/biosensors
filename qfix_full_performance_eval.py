#!/usr/bin/env python3
"""Full CV performance report for the model on corrected (qfix) features.

Now that all 95 ngs_observed sequences have bond-order-fixed (qfix)
reparameterized feature values, this evaluates the RandomForest classifier
entirely on those corrected values (unlike model_swap_eval.py, which
compares standard vs. qfix side by side on identical folds -- here there is
no "baseline" run, since qfix is now the sole feature source for every row).

Reuses ml_utils.py's SCORING/make_pipeline/run_cv/grouped_permutation_test
so this mirrors ML_classification.ipynb's CV pipeline (cell 5) exactly,
plus a group-safe permutation test for significance against chance.

Usage:
    python qfix_full_performance_eval.py --qfix_table qfix_500ns_feat_table.csv \
        --target_seq_list seq_ids_qfix_all95.txt
"""
import argparse
import numpy as np
import pandas as pd

from model_swap_eval import load_baseline_df, load_target_seq_ids
from ml_utils import make_pipeline, run_cv, grouped_permutation_test, N_SPLITS, RANDOM_STATE


def main():
    """Builds the fully qfix-corrected cohort and reports CV + permutation results."""
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qfix_table", default="qfix_500ns_feat_table.csv")
    p.add_argument("--target_seq_list", default="seq_ids_qfix_all95.txt",
                    help="seq_ids.txt-style list of sequences to source qfix values for "
                         "(default: %(default)s)")
    p.add_argument("--n_permutations", type=int, default=100)
    p.add_argument("--permutation_scoring", default="balanced_accuracy")
    p.add_argument("--random_state", type=int, default=42)
    p.add_argument("--out", default="qfix_full_performance_results.csv")
    args = p.parse_args()

    target_seq_ids = load_target_seq_ids(args.target_seq_list)
    df, feature_group_cols = load_baseline_df()
    feature_cols = [c for cols in feature_group_cols.values() for c in cols]

    missing_targets = set(target_seq_ids) - set(df["name"])
    if missing_targets:
        raise ValueError(f"Target sequences not found in baseline cohort: {missing_targets}")
    missing_from_cohort = set(df["name"]) - set(target_seq_ids)
    if missing_from_cohort:
        print(f"NOTE: {len(missing_from_cohort)} baseline-cohort sequence(s) are NOT in "
              f"--target_seq_list, so they keep their standard (non-qfix) feature values: "
              f"{sorted(missing_from_cohort)}")

    qfix = pd.read_csv(args.qfix_table).set_index("name")
    missing_qfix_cols = set(feature_cols) - set(qfix.columns)
    if missing_qfix_cols:
        raise ValueError(f"{args.qfix_table} is missing feature columns: {missing_qfix_cols}")

    # Swap every target row's feature values to its qfix reparameterization,
    # keeping Sequence/Label (unaffected by the ligand fix) from the baseline df.
    df_qfix = df.copy()
    for seq_id in target_seq_ids:
        row_idx = df_qfix.index[df_qfix["name"] == seq_id]
        df_qfix.loc[row_idx, feature_cols] = qfix.loc[seq_id, feature_cols].values

    print(f"Fully qfix-corrected cohort: {len(df_qfix)} sequences "
          f"({len(target_seq_ids)} with qfix values swapped in)")

    metrics_df, cv_results, seq_groups = run_cv(df_qfix, feature_cols, feature_group_cols)
    print("\n7-fold grouped-stratified CV performance (qfix features):")
    print(metrics_df.to_string())

    y = df_qfix["Label"].values
    X = df_qfix[feature_cols].values
    from sklearn.model_selection import StratifiedGroupKFold
    cv = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    pipe = make_pipeline(feature_group_cols, feature_cols)

    print(f"\nRunning group-safe permutation test "
          f"({args.n_permutations} permutations, scoring={args.permutation_scoring!r})...")
    observed, perm_scores, pval = grouped_permutation_test(
        pipe, X, y, seq_groups, cv, scoring=args.permutation_scoring,
        n_permutations=args.n_permutations, random_state=args.random_state,
    )
    print(f"Observed {args.permutation_scoring}: {observed:.4f}")
    print(f"Permutation scores: mean={perm_scores.mean():.4f} SD={perm_scores.std():.4f} "
          f"max={perm_scores.max():.4f}")
    print(f"Permutation p-value: {pval:.4g}")

    out_rows = [{"metric": k, "mean": v.mean(), "sd": v.std()} for k, v in cv_results.items()]
    out_rows.append({"metric": f"permutation_{args.permutation_scoring}_observed", "mean": observed, "sd": np.nan})
    out_rows.append({"metric": f"permutation_{args.permutation_scoring}_null_mean", "mean": perm_scores.mean(), "sd": perm_scores.std()})
    out_rows.append({"metric": f"permutation_{args.permutation_scoring}_pvalue", "mean": pval, "sd": np.nan})
    pd.DataFrame(out_rows).to_csv(args.out, index=False)
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
