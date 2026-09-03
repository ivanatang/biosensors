"""Tests which RMSD-to-Boltz-reference regions separate Binder from False
Positive sequences.

Reads the per-sequence summary CSV produced by extract_gate_latch_rmsd_feats.py
(Gate, Latch, Lb7a5, Recoil, and Whole-protein Ca RMSD-to-reference) and, for
each region, tests its three recommended ML features (mean, wide-window
late mean, full-trajectory regression slope) against a Binder vs False
Positive split, following this repo's established significance-screening
convention (core_vs_tail_regions.py / water_density_analysis.py):
(1) Mann-Whitney U per region/feature, with a Cohen's d effect size and a
rank-AUC alongside the p-value; (2) a Benjamini-Hochberg FDR correction
across all region/feature tests together (5 regions x 3 features = 15
tests, so this matters far more than it did for the original 2-column
Gate/Latch-only Welch t-test in plot_gate_latch_rmsd_to_ref.ipynb, since a
handful of uncorrected "hits" is expected by chance alone across 15 tests);
and (3) plots, one figure per feature kind (mean / late window mean /
slope) with one panel per region, plus a significance summary and a plain
verdict table.

Input (produced by extract_gate_latch_rmsd_feats.py, which now runs on Alpine
against the PetaLibrary archive -- pull its output CSV back to this local
checkout, e.g. via dtn.rc.colorado.edu, before running this script):
    gate_latch_rmsd_to_ref_summary{TAG}.csv

Usage:
    python rmsd_to_ref_significance.py \
        --summary-csv /Users/ivanatang/Developer/biosensors/analysis/gate_latch_rmsd_to_ref_summary_500ns.csv \
        --out-dir     /Users/ivanatang/Developer/biosensors/analysis/rmsd_to_ref
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

# ── Constants ─────────────────────────────────────────────────────────────────
REPO_DIR = "/Users/ivanatang/Developer/biosensors"

GROUP_COLOR = {
    "Binder":         "#648FFF",
    "False Positive": "#DC267F",
    "Low Confidence": "#FE6100",
    "Fail Geometry":  "#FFB000",
}
GROUP_ORDER = ["Binder", "False Positive", "Low Confidence", "Fail Geometry"]

GROUP_A, GROUP_B = "Binder", "False Positive"   # the headline test

REGIONS = ["Gate", "Latch", "Lb7a5", "Recoil", "Whole"]

# Display label -> column-name template (region name is filled in below).
# These three are the ML features extract_gate_latch_rmsd_feats.py's own
# docstring recommends, since wide-window late mean / slope were the most
# robust binder/nonbinder discriminators found for gate/latch.
FEATURE_TEMPLATES = {
    "mean":      "{region} RMSD mean (A)",
    "late_mean": "{region} RMSD late{window} mean (A)",
    "slope":     "{region} slope (A/ns)",
}

# md_candidate_guide.csv's md_group values use different suffixes than the
# seq_id naming convention (pair_XXXX_binder / _nb / _low_pkt / _fail_gate)
# used everywhere else in the repo (mirrors core_vs_tail_regions.py).
MD_GROUP_SUFFIX = {
    "binder": "binder",
    "non_binder": "nb",
    "negative_low_pocket": "low_pkt",
    "negative_fail_gate": "fail_gate",
}


# ── Stats helpers (mirrors core_vs_tail_regions.py) ───────────────────────────
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


def compare_groups(df, col, group_a=GROUP_A, group_b=GROUP_B, min_n=5):
    """Runs a Mann-Whitney U test between two groups for one column.

    Args:
        df: DataFrame containing `col` and a `Group` column.
        col: Name of the column to compare.
        group_a: `Group` value for the first group (default: GROUP_A).
        group_b: `Group` value for the second group (default: GROUP_B).
        min_n (int): Minimum non-NaN samples required per group; returns
            None below this (default: 5).

    Returns:
        dict | None: n, mean, Cohen's d, rank-AUC, and p-value for each
        group, or None if either group has fewer than `min_n` samples.
    """
    a = df.loc[df["Group"] == group_a, col].dropna().values
    b = df.loc[df["Group"] == group_b, col].dropna().values
    if len(a) < min_n or len(b) < min_n:
        return None
    _, p = mannwhitneyu(a, b, alternative="two-sided")
    return dict(n_binder=len(a), n_fp=len(b), mean_binder=a.mean(), mean_fp=b.mean(),
                cohens_d=cohens_d(a, b), auc=rank_auc(a, b), p=p)


# ── Loading / optional sequencing-confirmed restriction ───────────────────────
def load_source_ids(guide_path, source):
    """Loads seq_ids from md_candidate_guide.csv matching one source value.

    Mirrors core_vs_tail_regions.py's load_source_ids.

    Args:
        guide_path (str): Path to md_candidate_guide.csv.
        source (str): Value to match in the `source` column (e.g.
            "ngs_observed" for sequencing-confirmed via Y2H/FACS sort-seq,
            vs. "designed_assumed").

    Returns:
        set: Sequence names (matching the summary CSV's "Sequence" column,
        e.g. pair_3098_binder) matching `source`.
    """
    guide = pd.read_csv(guide_path)
    matched = guide[guide["source"] == source].copy()
    matched["Sequence"] = matched.apply(
        lambda r: f"{r['pair_id']}_{MD_GROUP_SUFFIX.get(r['md_group'], r['md_group'])}",
        axis=1)
    return set(matched["Sequence"])


# ── Box + jitter helper (matches the repo's plot_Rg_sasa.py / core_vs_tail_regions.py style)
def box_jitter(ax, df, col, groups, rng=None):
    """Draws a box-and-jitter plot for one column, one box per group.

    Args:
        ax: Matplotlib Axes to draw on.
        df: DataFrame containing `col` and a `Group` column.
        col (str): Column to plot.
        groups (list[str]): `Group` values to plot, in order, each colored
            via GROUP_COLOR.
        rng: numpy Generator for jitter; a fixed-seed default is used if
            None.
    """
    rng = rng or np.random.default_rng(42)
    xpos = 0
    for group in groups:
        vals = df.loc[df["Group"] == group, col].dropna().values
        if len(vals) == 0:
            xpos += 1
            continue
        color = GROUP_COLOR[group]
        ax.boxplot(vals, positions=[xpos], widths=0.6, patch_artist=True,
                   medianprops=dict(color="black", linewidth=1.5),
                   boxprops=dict(facecolor=color, alpha=0.5),
                   whiskerprops=dict(color=color), capprops=dict(color=color),
                   flierprops=dict(marker="", linestyle="none"))
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(xpos + jitter, vals, color=color, s=14, alpha=0.8, zorder=3)
        xpos += 1


# ── Raw stats for all region x feature-kind combinations (no FDR yet) ─────────
def compute_all_stats(df, window):
    """Runs Mann-Whitney/Cohen's d/AUC for every (region, feature_kind)
    combination, with no FDR correction applied yet.

    FDR correction has to happen once, across every row this function
    returns together, not per feature-kind, since all of them compete for
    the same false-discovery budget.

    Args:
        df: Summary DataFrame with one column per region/feature-kind
            combination (see FEATURE_TEMPLATES) and a `Group` column.
        window (float): Wide-window size in ns, used to resolve the
            "late_mean" column names.

    Returns:
        pandas.DataFrame: One row per (feature_kind, region) with stats
        from compare_groups, uncorrected.
    """
    rows = []
    for feature_kind, template in FEATURE_TEMPLATES.items():
        for region in REGIONS:
            col = template.format(region=region, window=int(window))
            if col not in df.columns:
                continue
            res = compare_groups(df, col)
            if res is None:
                continue
            rows.append(dict(feature_kind=feature_kind, region=region, column=col, **res))
    return pd.DataFrame(rows)


# ── Per-feature-kind plot, using already globally-FDR-corrected stats ────────
def plot_feature_kind(df, stats_df, feature_kind, out_dir, sig_threshold):
    """Plots one feature kind across all regions, Binder vs False Positive.

    `stats_df` must already carry a 'q' column from the single global
    BH-FDR correction computed in main(), so stars here match the final
    verdict.

    Args:
        df: Summary DataFrame (see compute_all_stats).
        stats_df: FDR-corrected stats DataFrame (see compute_all_stats;
            must have a "q" column).
        feature_kind (str): Key of FEATURE_TEMPLATES to plot.
        out_dir (str): Directory to write the figure to.
        sig_threshold (float): FDR q-value significance threshold.

    Returns:
        dict | None: Summary (feature_kind, n tested/significant, top
        region) for the verdict table, or None if there wasn't enough data.
    """
    kind_df = stats_df[stats_df["feature_kind"] == feature_kind]
    if kind_df.empty:
        print(f"\n[{feature_kind}] Not enough Binder/False Positive data -- skipping.")
        return None

    n_sig = int((kind_df["q"] < sig_threshold).sum())
    print(f"[{feature_kind}] {n_sig}/{len(kind_df)} regions FDR-significant "
          f"({GROUP_A} vs {GROUP_B}, q<{sig_threshold})")

    n = len(kind_df)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                              dpi=300, constrained_layout=True, squeeze=False)
    for i, row in enumerate(kind_df.itertuples()):
        ax = axes[i // ncols][i % ncols]
        box_jitter(ax, df, row.column, GROUP_ORDER)
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right", fontsize=8)
        star = "  *SIGNIFICANT*" if row.q < sig_threshold else ""
        ax.set_title(f"{row.region}\np={row.p:.2g}, q={row.q:.2g}, d={row.cohens_d:.2f}{star}",
                     fontsize=9)
        ax.grid(True, alpha=0.4)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(f"RMSD-to-reference, {feature_kind.upper()} feature -- {GROUP_A} vs {GROUP_B}\n"
                 f"{n_sig}/{len(kind_df)} regions significant (FDR q<{sig_threshold}, "
                 f"corrected across all region/feature tests)",
                 fontsize=13)
    path = os.path.join(out_dir, f"rmsd_to_ref_{feature_kind}_binder_vs_nonbinder.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved {path}")

    top = kind_df.sort_values("p").iloc[0]
    return dict(feature_kind=feature_kind, n_regions_tested=len(kind_df),
                n_regions_significant=n_sig, top_region=top["region"],
                top_region_p=top["p"], top_region_q=top["q"], top_region_d=top["cohens_d"])


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    """Runs the RMSD-to-reference Binder vs False Positive analysis end to end."""
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summary-csv",
                        default=os.path.join(REPO_DIR, "analysis",
                                             "gate_latch_rmsd_to_ref_summary_500ns.csv"),
                        help="Per-sequence summary CSV from extract_gate_latch_rmsd_feats.py "
                             "(default: %(default)s)")
    parser.add_argument("--wide-window-ns", type=float, default=100.0,
                        help="Must match the --wide-window-ns used during extraction, so the "
                             "'late{window} mean' column name resolves correctly (default: 100.0)")
    parser.add_argument("--out-dir",
                        default=os.path.join(REPO_DIR, "analysis", "rmsd_to_ref"),
                        help="Output directory for stats CSVs and figures (default: %(default)s)")
    parser.add_argument("--structure-source", default="all",
                        choices=["ngs_observed", "designed_assumed", "all"],
                        help="Filter sequences by md_candidate_guide.csv's source column. "
                             "'all' (default) applies no filtering.")
    parser.add_argument("--structure-guide",
                        default=os.path.join(REPO_DIR, "md_candidate_guide.csv"),
                        help="Path to md_candidate_guide.csv (default: %(default)s)")
    parser.add_argument("--sig-threshold", type=float, default=0.05,
                        help="FDR q-value significance threshold (default: 0.05)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if not os.path.exists(args.summary_csv):
        raise FileNotFoundError(f"Summary CSV not found: {args.summary_csv}")
    df = pd.read_csv(args.summary_csv)
    print(f"Loaded {len(df)} sequences from {args.summary_csv}")

    if args.structure_source != "all":
        source_ids = load_source_ids(args.structure_guide, args.structure_source)
        before = len(df)
        df = df[df["Sequence"].isin(source_ids)].reset_index(drop=True)
        print(f"Restricting to source == {args.structure_source}: {before} -> {len(df)} sequences")

    n_binder = int((df["Group"] == GROUP_A).sum())
    n_fp = int((df["Group"] == GROUP_B).sum())
    print(f"{GROUP_A}: {n_binder}, {GROUP_B}: {n_fp}")

    full_stats = compute_all_stats(df, args.wide_window_ns)
    if full_stats.empty:
        raise RuntimeError("No region/feature columns found in the summary CSV -- check "
                            "--summary-csv and --wide-window-ns.")

    # ── ONE Benjamini-Hochberg FDR correction, across all region x feature
    # tests together (5 regions x 3 feature kinds = 15 tests here) -- all of
    # them compete for the same false-discovery budget ──
    full_stats["q"] = bh_fdr(full_stats["p"].values)
    full_stats = full_stats.sort_values("p")
    stats_path = os.path.join(args.out_dir, "rmsd_to_ref_significance_stats.csv")
    full_stats.to_csv(stats_path, index=False)
    print(f"\nSaved {stats_path}")
    print(full_stats.to_string(index=False))

    n_sig_overall = int((full_stats["q"] < args.sig_threshold).sum())
    print(f"\nOVERALL: {n_sig_overall}/{len(full_stats)} region/feature combinations "
          f"FDR-significant ({GROUP_A} vs {GROUP_B}, q<{args.sig_threshold})")

    verdicts = []
    for feature_kind in FEATURE_TEMPLATES:
        print("\n" + "=" * 70)
        print(f"RMSD-TO-REFERENCE -- feature kind: {feature_kind}")
        print("=" * 70)
        verdict = plot_feature_kind(df, full_stats, feature_kind, args.out_dir, args.sig_threshold)
        if verdict is not None:
            verdicts.append(verdict)

    # ── Verdict: a plain table, not a chart ──
    print("\n" + "=" * 70)
    print(f"VERDICT: {GROUP_A} vs {GROUP_B} difference, by feature kind")
    print("=" * 70)
    vdf = pd.DataFrame(verdicts)
    vdf.to_csv(os.path.join(args.out_dir, "rmsd_to_ref_verdict.csv"), index=False)
    print(vdf.to_string(index=False))

    print("\nDone. All outputs written to:", args.out_dir)


if __name__ == "__main__":
    main()
