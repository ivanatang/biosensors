"""Aggregates per-residue backbone/side-chain contact splits and tests
Binder vs False Positive separation.

Run on the login node after all residue_atom_split_worker.sh SLURM jobs
finish. Collects per-sequence *_res{RESSEQ}_atomsplit_{TAG}.csv files into
one table, then tests whether backbone_occupancy and sidechain_occupancy
separate Binder from False Positive, following this repo's established
significance-screening convention (core_vs_tail_regions.py /
rmsd_to_ref_significance.py): Mann-Whitney U, Cohen's d, rank-AUC, BH-FDR
across the two tests.

Usage:
    python agg_residue_atom_split.py                              # res 116, core, 40-500ns
    python agg_residue_atom_split.py --resseq 116 --ligand-region core
"""

import os
import argparse
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

out_dir = "/projects/ivta1597/biosensors/LIG_contacts"

GROUP_A, GROUP_B = "Binder", "False Positive"


def cohens_d(a, b):
    """Computes Cohen's d effect size between two samples.

    Args:
        a: First sample (array-like).
        b: Second sample (array-like).

    Returns:
        float: Standardized mean difference (a - b), pooled SD. 0.0 if
        pooled variance is non-positive.
    """
    na, nb = len(a), len(b)
    pooled_var = ((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2)
    if pooled_var <= 0:
        return 0.0
    return (a.mean() - b.mean()) / np.sqrt(pooled_var)


def rank_auc(a, b):
    """Computes the AUC of using the feature value to separate two groups.

    Args:
        a: Sample treated as the positive class (label 1).
        b: Sample treated as the negative class (label 0).

    Returns:
        float: ROC-AUC. 0.5 means no separation, >0.5 means `a` tends
        higher, <0.5 means `b` tends higher.
    """
    y = np.r_[np.ones(len(a)), np.zeros(len(b))]
    scores = np.r_[a, b]
    return roc_auc_score(y, scores)


def bh_fdr(pvals):
    """Applies a Benjamini-Hochberg FDR correction to a set of p-values.

    Args:
        pvals: Array-like of raw p-values.

    Returns:
        numpy.ndarray: FDR-adjusted q-values, same order as `pvals`.
    """
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
    """Loads seq_id -> seq_type from a tab-separated seq_ids file.

    Args:
        seq_list_path (str): Path to a seq_ids.txt-style file (columns:
            seq_id, seq_type, ...). Blank lines and lines starting with
            "#" are skipped.

    Returns:
        dict: Maps seq_id to seq_type ("Unknown" if the column is missing).
    """
    mapping = {}
    with open(seq_list_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            mapping[parts[0]] = parts[1] if len(parts) > 1 else "Unknown"
    return mapping


def load_seq_ids(seq_list_path):
    """Loads the seq_id column from a seq_ids.txt-style file.

    Args:
        seq_list_path (str): Path to the seq_ids file. Blank lines and
            lines starting with "#" are skipped.

    Returns:
        list[str]: seq_ids in file order.
    """
    with open(seq_list_path) as f:
        return [line.split()[0] for line in f
                if line.strip() and not line.startswith("#")]


def compare_groups(df, col, group_a=GROUP_A, group_b=GROUP_B, min_n=5):
    """Runs a Mann-Whitney U test between two groups for one column.

    Args:
        df: DataFrame containing `col` and a `seq_type` column.
        col: Name of the column to compare.
        group_a: `seq_type` value for the first group (default: GROUP_A).
        group_b: `seq_type` value for the second group (default: GROUP_B).
        min_n (int): Minimum non-NaN samples required per group; returns
            None below this (default: 5).

    Returns:
        dict | None: n, mean, Cohen's d, rank-AUC, and p-value for each
        group, or None if either group has fewer than `min_n` samples.
    """
    a = df.loc[df["seq_type"] == group_a, col].dropna().values
    b = df.loc[df["seq_type"] == group_b, col].dropna().values
    if len(a) < min_n or len(b) < min_n:
        return None
    _, p = mannwhitneyu(a, b, alternative="two-sided")
    return dict(n_binder=len(a), n_fp=len(b), mean_binder=a.mean(), mean_fp=b.mean(),
                cohens_d=cohens_d(a, b), auc=rank_auc(a, b), p=p)


def main():
    """Aggregates one residue's atom-split results and tests Binder vs FP."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq_list', default="/projects/ivta1597/biosensors/seq_ids_orig.txt")
    parser.add_argument('--resseq', type=int, default=116)
    parser.add_argument('--start-ns', type=float, default=40.0)
    parser.add_argument('--end-ns',   type=float, default=500.0)
    parser.add_argument('--ligand-region', choices=['whole', 'core', 'tail'], default='core')
    parser.add_argument('--suffix', default='',
                        help="Run-directory/output suffix, e.g. '_qfix', matching "
                             "whatever residue_atom_split_contact.py was run with "
                             "for these sequences (default: '', the standard "
                             "directory). NOTE: residue_atom_split_contact.py "
                             "itself does not currently have --suffix support, so "
                             "this only works once that script does too.")
    args = parser.parse_args()

    TAG        = f"{int(args.start_ns)}_{int(args.end_ns)}ns"
    REGION_TAG = "" if args.ligand_region == "whole" else f"_{args.ligand_region}"
    results_dir = os.path.join(out_dir, f"residue_atomsplit_results_{TAG}{REGION_TAG}{args.suffix}")
    combined_out = os.path.join(out_dir, f"residue{args.resseq}_atomsplit_all_{TAG}{REGION_TAG}{args.suffix}.csv")
    stats_out    = os.path.join(out_dir, f"residue{args.resseq}_atomsplit_binder_vs_fp_stats_{TAG}{REGION_TAG}{args.suffix}.csv")

    print(f"Residue     : {args.resseq}")
    print(f"Window      : {args.start_ns:.0f}-{args.end_ns:.0f} ns  (tag: {TAG})")
    print(f"Region      : {args.ligand_region}")
    print(f"Results dir : {results_dir}")

    # ── Build explicit per-sequence paths from seq_list (not a glob over the
    # whole shared results dir) -- a glob also picks up leftover output from
    # any sequence that ever had this step run, including ones no longer in
    # the current cohort (e.g. seq_ids_orig.txt's superseded 200-sequence
    # list), inflating the combined table with unlabeled rows. See the
    # equivalent fix in agg_gate_latch_water_bridge.py.
    seq_ids = load_seq_ids(args.seq_list)
    seq_type_map = load_seq_type_map(args.seq_list)

    rows    = []
    missing = []
    for seq_id in seq_ids:
        path = os.path.join(results_dir,
                             f"{seq_id}_res{args.resseq}_atomsplit_{TAG}{REGION_TAG}.csv")
        if not os.path.exists(path):
            print(f"  MISSING: {path}")
            missing.append(seq_id)
            continue
        df = pd.read_csv(path)
        df["seq_type"] = seq_type_map.get(seq_id, "Unknown")
        rows.append(df)

    print(f"Found {len(rows)} of {len(seq_ids)} summary files")
    if not rows:
        raise FileNotFoundError(
            f"No summary CSVs found for any sequence in {args.seq_list}.\n"
            f"Expected pattern: {results_dir}/<seq_id>_res{args.resseq}_atomsplit_{TAG}{REGION_TAG}.csv"
        )

    combined = pd.concat(rows, ignore_index=True)

    if missing:
        print(f"\nWARNING: {len(missing)} sequences missing atom-split results: {missing}")
    else:
        print("All sequences accounted for.")

    combined.to_csv(combined_out, index=False)
    print(f"\nSaved combined table: {combined_out}  (shape={combined.shape})")

    # ── Binder vs False Positive stats ────────────────────────────────────
    features = ["backbone_occupancy", "sidechain_occupancy",
                "backbone_mean_mindist_nm", "sidechain_mean_mindist_nm"]

    rows = []
    for feat in features:
        result = compare_groups(combined, feat)
        if result is None:
            continue
        result["feature"] = feat
        rows.append(result)

    if not rows:
        print("\nNot enough Binder/False Positive sequences with valid data to test.")
        return

    stats_df = pd.DataFrame(rows)
    stats_df["q"] = bh_fdr(stats_df["p"].values)
    stats_df = stats_df[["feature", "n_binder", "n_fp", "mean_binder", "mean_fp",
                          "cohens_d", "auc", "p", "q"]].sort_values("p")

    stats_df.to_csv(stats_out, index=False)
    print(f"\nResidue {args.resseq} backbone vs side-chain contact — Binder vs False Positive:")
    print(stats_df.to_string(index=False))
    print(f"\nSaved stats: {stats_out}")


if __name__ == "__main__":
    main()
