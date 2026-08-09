"""
agg_gate_latch_hbond.py
---------------------------------
Aggregates per-sequence gate/latch backbone-vs-side-chain real-H-bond
cross-reference results (from parse_gate_latch_hbond.py, run cohort-wide via
run_gate_latch_hbond_crossref_all.sh) into cohort-level averages, and tests
whether backbone/side-chain H-bond presence differs Binder vs False
Positive -- same significance-screening convention as
agg_gate_latch_water_bridge.py (Mann-Whitney U, Cohen's d, rank-AUC, BH-FDR).

Answers, across however many sequences have results:
  "On average, when the gate-latch-ligand network is actually forming
  (triple_bridge_active), what fraction of the time is a real (distance +
  angle) hydrogen bond present on the gate backbone vs. gate side chain,
  and on the latch backbone vs. latch side chain?"

Usage:
    python agg_gate_latch_hbond.py
    python agg_gate_latch_hbond.py --start-ns 40 --end-ns 500
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

out_dir      = "/projects/ivta1597/biosensors/water_analysis"
results_base = "/scratch/alpine/ivta1597/LCA_boltz_models"

GROUP_A, GROUP_B = "Binder", "False Positive"
PARTS = ["gate_backbone", "gate_sidechain", "latch_backbone", "latch_sidechain"]


def cohens_d(a, b):
    na, nb = len(a), len(b)
    pooled_var = ((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2)
    if pooled_var <= 0:
        return 0.0
    return (a.mean() - b.mean()) / np.sqrt(pooled_var)


def rank_auc(a, b):
    y = np.r_[np.ones(len(a)), np.zeros(len(b))]
    scores = np.r_[a, b]
    return roc_auc_score(y, scores)


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


def load_seq_type_map(seq_list_path):
    mapping = {}
    with open(seq_list_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            mapping[parts[0]] = parts[1] if len(parts) > 1 else "Unknown"
    return mapping


def compare_groups(df, col, group_a=GROUP_A, group_b=GROUP_B, min_n=5):
    a = df.loc[df["seq_type"] == group_a, col].dropna().values
    b = df.loc[df["seq_type"] == group_b, col].dropna().values
    if len(a) < min_n or len(b) < min_n:
        return None
    _, p = mannwhitneyu(a, b, alternative="two-sided")
    return dict(n_binder=len(a), n_fp=len(b), mean_binder=a.mean(), mean_fp=b.mean(),
                cohens_d=cohens_d(a, b), auc=rank_auc(a, b), p=p)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq_list', default="/projects/ivta1597/biosensors/seq_ids_orig.txt")
    parser.add_argument('--start-ns', type=float, default=40.0)
    parser.add_argument('--end-ns',   type=float, default=500.0)
    args = parser.parse_args()

    TAG = f"{int(args.start_ns)}_{int(args.end_ns)}ns"

    results_glob = os.path.join(
        results_base, "*", "*", "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd",
        f"*_gate_latch_hbond_crossref_{TAG}.csv"
    )
    long_out    = os.path.join(out_dir, f"gate_latch_hbond_crossref_all_{TAG}.csv")
    wide_out    = os.path.join(out_dir, f"gate_latch_hbond_crossref_wide_{TAG}.csv")
    cohort_out  = os.path.join(out_dir, f"gate_latch_hbond_cohort_avg_{TAG}.csv")
    stats_out   = os.path.join(out_dir, f"gate_latch_hbond_binder_vs_fp_stats_{TAG}.csv")

    print(f"Window      : {args.start_ns:.0f}-{args.end_ns:.0f} ns  (tag: {TAG})")
    print(f"Results glob: {results_glob}")

    files = sorted(glob.glob(results_glob))
    print(f"Found {len(files)} per-sequence crossref files")
    if not files:
        raise FileNotFoundError(
            f"No crossref CSVs found matching: {results_glob}\n"
            f"Run submit_gate_latch_water_bridge.sh ... true, then "
            f"submit_gate_latch_hbond.sh, then run_gate_latch_hbond_crossref_all.sh first."
        )

    # A per-sequence crossref CSV can be present but EMPTY (0 bytes) if
    # parse_gate_latch_hbond.py ran to completion without crashing but
    # found none of that sequence's 4 gmx hbond .xvg files -- e.g. its
    # gate_latch_hbond_gmx.sh job failed silently from the aggregator's
    # point of view (run_gate_latch_hbond_crossref_all.sh only reports a
    # PYTHON crash as "failed", not a python run that completed but wrote
    # nothing). Skip and report those instead of letting one bad sequence
    # kill the whole aggregation.
    dfs, empty_files = [], []
    for f in files:
        try:
            dfs.append(pd.read_csv(f))
        except pd.errors.EmptyDataError:
            empty_files.append(f)

    if empty_files:
        print(f"\nWARNING: {len(empty_files)} crossref CSV(s) were empty (likely missing "
              f"gmx hbond .xvg outputs for that sequence) -- skipped, not counted below:")
        for f in empty_files:
            print(f"  {f}")

    combined = pd.concat(dfs, ignore_index=True)

    seq_type_map = load_seq_type_map(args.seq_list)
    combined["seq_type"] = combined["seq_id"].map(seq_type_map)

    missing_type = combined.loc[combined["seq_type"].isna(), "seq_id"].unique()
    if len(missing_type):
        print(f"WARNING: {len(missing_type)} sequences not found in {args.seq_list}: "
              f"{list(missing_type)}")

    if os.path.exists(args.seq_list):
        with open(args.seq_list) as fh:
            all_ids = [l.split()[0] for l in fh if l.strip() and not l.startswith('#')]
        missing = set(all_ids) - set(combined["seq_id"])
        if missing:
            print(f"WARNING: {len(missing)} sequences missing hbond crossref results: {missing}")
        else:
            print("All sequences accounted for.")

    combined.to_csv(long_out, index=False)
    print(f"\nSaved combined long-format table: {long_out}  (shape={combined.shape})")

    # ── Pivot to one row per sequence: {part}_{condition}_frac columns ────
    combined["col"] = combined["group"] + "_" + combined["condition"] + "_frac"
    wide = combined.pivot_table(
        index=["seq_id", "seq_type"], columns="col", values="frac_hbond_present"
    ).reset_index()
    wide.columns.name = None
    wide.to_csv(wide_out, index=False)
    print(f"Saved per-sequence wide table: {wide_out}  (shape={wide.shape})")

    # ── Cohort-wide average, active-network frames only ────────────────────
    active = combined[combined["condition"] == "triple_bridge_active"]
    cohort_rows = []
    for part in PARTS:
        sub = active[active["group"] == part]["frac_hbond_present"].dropna()
        if len(sub) == 0:
            continue
        cohort_rows.append(dict(
            part=part, n_sequences=len(sub),
            mean_frac_hbond_present=round(float(sub.mean()), 4),
            median_frac_hbond_present=round(float(sub.median()), 4),
            std=round(float(sub.std()), 4),
        ))
        for st in sorted(combined["seq_type"].dropna().unique()):
            sub_st = active[(active["group"] == part) & (active["seq_type"] == st)]["frac_hbond_present"].dropna()
            if len(sub_st) == 0:
                continue
            cohort_rows.append(dict(
                part=f"{part} [{st}]", n_sequences=len(sub_st),
                mean_frac_hbond_present=round(float(sub_st.mean()), 4),
                median_frac_hbond_present=round(float(sub_st.median()), 4),
                std=round(float(sub_st.std()), 4),
            ))

    cohort_df = pd.DataFrame(cohort_rows)
    cohort_df.to_csv(cohort_out, index=False)
    print(f"\n-- Cohort-average: fraction of triple-bridge-active frames with a real H-bond --")
    print(cohort_df.to_string(index=False))
    print(f"\nSaved: {cohort_out}")

    # ── Binder vs False Positive stats on the active-frame fractions ───────
    stat_cols = [f"{p}_triple_bridge_active_frac" for p in PARTS]
    rows = []
    for col in stat_cols:
        if col not in wide.columns:
            continue
        result = compare_groups(wide, col)
        if result is None:
            continue
        result["feature"] = col
        rows.append(result)

    if rows:
        stats_df = pd.DataFrame(rows)
        stats_df["q"] = bh_fdr(stats_df["p"].values)
        stats_df = stats_df[["feature", "n_binder", "n_fp", "mean_binder", "mean_fp",
                              "cohens_d", "auc", "p", "q"]].sort_values("p")
        stats_df.to_csv(stats_out, index=False)
        print(f"\nGate/latch backbone-vs-side-chain H-bond presence -- Binder vs False Positive:")
        print(stats_df.to_string(index=False))
        print(f"\nSaved stats: {stats_out}")
    else:
        print("\nNot enough Binder/False Positive sequences with valid data to test.")


if __name__ == "__main__":
    main()
