"""
agg_residue_atom_split.py
---------------------------
After all residue_atom_split_worker.sh SLURM jobs finish, run this on the
login node to (1) collect per-sequence *_res{RESSEQ}_atomsplit_{TAG}.csv
files into one table, and (2) test whether backbone_occupancy and
sidechain_occupancy separate Binder from False Positive, following this
repo's established significance-screening convention (core_vs_tail_regions.py /
rmsd_to_ref_significance.py): Mann-Whitney U, Cohen's d, rank-AUC, BH-FDR
across the two tests.

Usage:
    python agg_residue_atom_split.py                              # res 116, core, 40-500ns
    python agg_residue_atom_split.py --resseq 116 --ligand-region core
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

out_dir = "/projects/ivta1597/biosensors/LIG_contacts"

GROUP_A, GROUP_B = "Binder", "False Positive"


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

    # ── Load and combine per-sequence summaries ──────────────────────────
    summary_files = sorted(glob.glob(
        os.path.join(results_dir, f"*_res{args.resseq}_atomsplit_{TAG}{REGION_TAG}.csv")
    ))
    print(f"Found {len(summary_files)} summary files")
    if not summary_files:
        raise FileNotFoundError(
            f"No summary CSVs found in: {results_dir}\n"
            f"Expected pattern: *_res{args.resseq}_atomsplit_{TAG}{REGION_TAG}.csv"
        )

    combined = pd.concat([pd.read_csv(f) for f in summary_files], ignore_index=True)

    seq_type_map = load_seq_type_map(args.seq_list)
    combined["seq_type"] = combined["seq_id"].map(seq_type_map)

    missing_type = combined[combined["seq_type"].isna()]
    if len(missing_type):
        print(f"WARNING: {len(missing_type)} sequences not found in {args.seq_list}: "
              f"{missing_type['seq_id'].tolist()}")

    if os.path.exists(args.seq_list):
        with open(args.seq_list) as fh:
            all_ids = [l.split()[0] for l in fh if l.strip() and not l.startswith('#')]
        missing = set(all_ids) - set(combined["seq_id"])
        if missing:
            print(f"WARNING: {len(missing)} sequences missing atom-split results: {missing}")
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
