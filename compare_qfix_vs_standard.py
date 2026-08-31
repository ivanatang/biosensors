#!/usr/bin/env python3
"""
compare_qfix_vs_standard.py

Paired comparison of the 6 qfix pilot sequences' feature values (bond-order-
fixed reparameterization) against their original (standard) values in
feat_table_500ns.xlsx, per feature column. Since this is the SAME 6
sequences before/after (not two independent groups), this uses a paired
Wilcoxon signed-rank test rather than the Mann-Whitney/Cohen's d/rank-AUC
convention used elsewhere in this repo for independent Binder-vs-False
Positive comparisons (see agg_gate_latch_water_bridge.py etc.) -- that
convention doesn't apply to n=6 paired samples.

With n=6 pairs, the Wilcoxon test has limited power (the smallest possible
two-sided p-value is 1/32 = 0.03125), so this is exploratory: it tells you
which features moved the most and in which direction, not a definitive
significance claim.

Usage:
    python compare_qfix_vs_standard.py --qfix_table qfix_pilot_feat_table.csv
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from model_swap_eval import load_baseline_df, TARGET_SEQ_IDS


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qfix_table", default="qfix_pilot_feat_table.csv")
    p.add_argument("--out_long", default="qfix_vs_standard_deltas_long.csv",
                    help="Per-sequence per-feature delta table")
    p.add_argument("--out_summary", default="qfix_vs_standard_summary.csv",
                    help="Per-feature paired-comparison summary")
    args = p.parse_args()

    qfix = pd.read_csv(args.qfix_table).set_index("name")
    missing = set(TARGET_SEQ_IDS) - set(qfix.index)
    if missing:
        raise ValueError(f"{args.qfix_table} is missing sequences: {missing}")

    # Use the SAME fully-merged baseline df model_swap_eval.py trains on --
    # feat_table_500ns.xlsx alone only carries the dw_pocket/contact_type/
    # rmsf columns; gate-latch RMSD, salt bridge, core/tail delta, hydration,
    # and water bridge are merged in at runtime, same as the notebook does.
    std_full, feature_group_cols = load_baseline_df()
    std_full = std_full.set_index("name")
    feature_cols = [c for cols in feature_group_cols.values() for c in cols]
    missing_std = set(TARGET_SEQ_IDS) - set(std_full.index)
    if missing_std:
        raise ValueError(f"Baseline table is missing sequences: {missing_std}")
    missing_qfix_cols = set(feature_cols) - set(qfix.columns)
    if missing_qfix_cols:
        raise ValueError(f"{args.qfix_table} is missing feature columns: {missing_qfix_cols}")
    print(f"Comparing {len(feature_cols)} feature columns across {len(TARGET_SEQ_IDS)} paired sequences")

    std = std_full.loc[TARGET_SEQ_IDS, feature_cols]
    qfx = qfix.loc[TARGET_SEQ_IDS, feature_cols]

    # ── Long-format per-sequence deltas ──────────────────────────────────────
    long_rows = []
    for seq_id in TARGET_SEQ_IDS:
        for col in feature_cols:
            s, q = std.loc[seq_id, col], qfx.loc[seq_id, col]
            long_rows.append({"seq_id": seq_id, "feature": col,
                               "standard": s, "qfix": q, "delta": q - s})
    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(args.out_long, index=False)
    print(f"Saved per-sequence deltas -> {args.out_long}")

    # ── Per-feature paired summary ───────────────────────────────────────────
    summary_rows = []
    for col in feature_cols:
        s = std[col].values.astype(float)
        q = qfx[col].values.astype(float)
        delta = q - s

        mean_delta = float(np.mean(delta))
        median_delta = float(np.median(delta))
        # Percent change relative to the standard value's mean magnitude,
        # guarding against a near-zero denominator.
        denom = np.mean(np.abs(s))
        pct_change = float(mean_delta / denom * 100) if denom > 1e-9 else np.nan

        n_nonzero = int(np.sum(delta != 0))
        if n_nonzero >= 1:
            try:
                stat, pval = wilcoxon(delta)
            except ValueError:
                pval = np.nan
        else:
            pval = 1.0

        summary_rows.append({
            "feature": col, "n_nonzero_diffs": n_nonzero,
            "mean_standard": float(np.mean(s)), "mean_qfix": float(np.mean(q)),
            "mean_delta": mean_delta, "median_delta": median_delta,
            "pct_change": pct_change, "p": pval,
        })

    summary_df = pd.DataFrame(summary_rows)
    valid_p = summary_df["p"].notna()
    summary_df["q"] = np.nan
    summary_df.loc[valid_p, "q"] = bh_fdr(summary_df.loc[valid_p, "p"].values)
    summary_df = summary_df.sort_values("p", na_position="last")

    summary_df.to_csv(args.out_summary, index=False)
    print(f"Saved per-feature summary -> {args.out_summary}")

    n_sig = int((summary_df["q"] < 0.05).sum())
    print(f"\nFeatures with q < 0.05: {n_sig} of {len(feature_cols)}")
    print("\nTop 15 features by |mean_delta| as a fraction of the standard value:")
    top = summary_df.reindex(summary_df["pct_change"].abs().sort_values(ascending=False).index).head(15)
    print(top[["feature", "mean_standard", "mean_qfix", "mean_delta", "pct_change", "p", "q"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
