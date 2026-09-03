"""Tests whether Binder vs nonbinder separation differs between the LCA
ligand's core (steroid ring system) and tail (C20-C24 carboxylate chain).

Primary comparison is Binder vs False Positive: both went through the same
computational pipeline, so the wet-lab outcome is the only difference,
making this the most direct binder-vs-nonbinder contrast available. Low
Confidence and Fail Geometry sequences appear in the plots for context only
(their pocket/geometry failed QC before this comparison even applies, and
there are only ~9-10 of each); they are not part of the headline test.

Core and tail are analyzed completely separately: every plot and stats
table answers "is there a Binder vs False Positive difference in this
region alone?" with no cross-region encoding (no shared panels or color
axis). The two regions are never compared to each other here; for that,
see compare_ligand_regions.py. They only appear together in the final
verdict table (region, n significant, top hit), which is plain text/CSV.

For both contact-type features and per-residue R-scores, per region, this
script: (1) runs Mann-Whitney U for every feature/residue with a Cohen's d
effect size and rank-AUC alongside the p-value, (2) applies a
Benjamini-Hochberg FDR correction, since testing many residues makes a
handful of uncorrected "hits" expected by chance, and (3) plots a
per-feature panel, a significance plot across all residues, and a
top-hit-residues panel, each scoped to that region only.

Inputs (produced by aggregate_r_scores.py / agg_contact_feats.py with
--ligand-region core|tail):
    r_scores_all_sequences_{TAG}_core.csv / _tail.csv
    contact_features_all_{TAG}_core.csv / _tail.csv

Usage:
    python core_vs_tail_regions.py \
        --r-scores-dir /projects/ivta1597/biosensors/water_analysis/agg_out \
        --contact-dir  /projects/ivta1597/biosensors/LIG_contacts \
        --seq-list     /projects/ivta1597/biosensors/seq_ids_orig.txt \
        --tag 40_500ns \
        --out-dir      /projects/ivta1597/biosensors/analysis/core_vs_tail
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

# ── Constants ─────────────────────────────────────────────────────────────────
GROUP_COLOR = {
    "Binder":         "#648FFF",
    "False Positive": "#DC267F",
    "Low Confidence": "#FE6100",
    "Fail Geometry":  "#FFB000",
}
GROUP_ORDER = ["Binder", "False Positive", "Low Confidence", "Fail Geometry"]

REGIONS = ["core", "tail"]

GROUP_A, GROUP_B = "Binder", "False Positive"   # the headline test

CONTACT_FEATURES = ["mean_frac_hydrophobic", "mean_n_total", "mean_n_hydrophobic",
                     "mean_n_polar", "mean_n_pos_charged", "mean_n_neg_charged"]


# ── Stats helpers ─────────────────────────────────────────────────────────────
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


# ── Loading ───────────────────────────────────────────────────────────────────
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


def load_r_scores(r_scores_dir, tag):
    """Loads per-region R-score CSVs written by aggregate_r_scores.py.

    Args:
        r_scores_dir (str): Directory containing
            r_scores_all_sequences_{tag}_{region}.csv files.
        tag (str): Run tag (e.g. "40_500ns").

    Returns:
        dict: Maps region ("core"/"tail") to its DataFrame. Missing files
        are skipped with a printed warning.
    """
    out = {}
    for region in REGIONS:
        path = os.path.join(r_scores_dir, f"r_scores_all_sequences_{tag}_{region}.csv")
        if os.path.exists(path):
            out[region] = pd.read_csv(path)
            print(f"Loaded {region} R-scores: {path}  ({len(out[region])} sequences)")
        else:
            print(f"  MISSING: {path}")
    return out


def load_dw_scores(dw_dir, tag):
    """Loads per-region D/W occupancy CSVs written by aggregate_r_scores.py.

    D and W are the two terms R = (D - W) * I is built from, kept here
    unmerged so they can be screened without R's ambiguity (R = 0 both
    when D and W are balanced and when neither ever contacts the ligand).

    Args:
        dw_dir (str): Directory containing
            dw_scores_all_sequences_{tag}_{region}.csv files (same run as
            load_r_scores).
        tag (str): Run tag (e.g. "40_500ns").

    Returns:
        dict: Maps region ("core"/"tail") to its DataFrame. Missing files
        are skipped with a printed warning.
    """
    out = {}
    for region in REGIONS:
        path = os.path.join(dw_dir, f"dw_scores_all_sequences_{tag}_{region}.csv")
        if os.path.exists(path):
            out[region] = pd.read_csv(path)
            print(f"Loaded {region} D/W scores: {path}  ({len(out[region])} sequences)")
        else:
            print(f"  MISSING: {path}")
    return out


def load_contact_features(contact_dir, tag, seq_type_map):
    """Loads per-region contact-feature CSVs written by agg_contact_feats.py.

    Args:
        contact_dir (str): Directory containing
            contact_features_all_{tag}_{region}.csv files.
        tag (str): Run tag (e.g. "40_500ns").
        seq_type_map (dict): seq_id -> seq_type, used to add a `seq_type`
            column (see load_seq_type_map).

    Returns:
        dict: Maps region ("core"/"tail") to its DataFrame. Missing files
        are skipped with a printed warning.
    """
    out = {}
    for region in REGIONS:
        path = os.path.join(contact_dir, f"contact_features_all_{tag}_{region}.csv")
        if os.path.exists(path):
            df = pd.read_csv(path)
            df["seq_type"] = df["seq_id"].map(seq_type_map)
            out[region] = df
            print(f"Loaded {region} contact features: {path}  ({len(df)} sequences)")
        else:
            print(f"  MISSING: {path}")
    return out


# md_candidate_guide.csv's md_group values use different suffixes than the
# seq_id naming convention (pair_XXXX_binder/_nb/_low_pkt/_fail_gate) used
# elsewhere in the repo.
MD_GROUP_SUFFIX = {
    "binder": "binder",
    "non_binder": "nb",
    "negative_low_pocket": "low_pkt",
    "negative_fail_gate": "fail_gate",
}


def load_source_ids(guide_path, source):
    """Loads seq_ids from md_candidate_guide.csv matching one source value.

    Args:
        guide_path (str): Path to md_candidate_guide.csv.
        source (str): Value to match in the `source` column (e.g.
            "ngs_observed" for sequencing-confirmed via Y2H/FACS sort-seq,
            vs. "designed_assumed").

    Returns:
        set: seq_ids (pair_id + mapped md_group suffix) matching `source`.
    """
    guide = pd.read_csv(guide_path)
    matched = guide[guide["source"] == source].copy()
    matched["seq_id"] = matched.apply(
        lambda r: f"{r['pair_id']}_{MD_GROUP_SUFFIX.get(r['md_group'], r['md_group'])}",
        axis=1)
    return set(matched["seq_id"])


def filter_to_source(data, source_ids, source, label):
    """Filters each region's DataFrame to a set of seq_ids, in place.

    Args:
        data (dict): Maps region to DataFrame (as returned by the load_*
            functions above). Modified and returned.
        source_ids (set): seq_ids to keep (see load_source_ids).
        source (str): Source name, used only in the printed summary.
        label (str): Data-kind label, used only in the printed summary.

    Returns:
        dict: `data`, with each region's DataFrame filtered to `source_ids`.
    """
    for region, df in data.items():
        before = len(df)
        data[region] = df[df["seq_id"].isin(source_ids)].reset_index(drop=True)
        print(f"  {label} {region}: {before} -> {len(data[region])} sequences "
              f"(source == {source})")
    return data


def resid_columns(df, prefix="R"):
    """Lists a DataFrame's per-residue columns, sorted by resSeq.

    Args:
        df: DataFrame with columns named "{prefix}_{resSeq}".
        prefix (str): Column-name prefix to match (default: "R").

    Returns:
        list[str]: Matching column names, sorted numerically by resSeq.
    """
    return sorted((c for c in df.columns if c.startswith(f"{prefix}_")),
                  key=lambda c: int(c.split("_")[1]))


def contactable_residues(df, prefix="R", min_n=5):
    """Finds residues with enough real ligand contact to be worth testing.

    R is NaN, and so silently dropped by compare_groups, exactly when a
    sequence never contacted that residue (D+W = 0). D and W instead record
    genuine zero occupancy in that case, so without this restriction a D/W
    screen would test all ~181 topology residues, including ones that never
    come near the ligand. This reproduces the same pocket restriction R
    gets for free from its NaN handling, so D/W are screened over the same
    residue population R was.

    Args:
        df: DataFrame with per-residue columns (see resid_columns).
        prefix (str): Column-name prefix identifying the metric (default:
            "R").
        min_n (int): Minimum Binder and False Positive sequences with a
            non-NaN score required to keep a residue (default: 5).

    Returns:
        list[int]: resSeq values for residues meeting the min_n threshold
        in both groups, sorted ascending.
    """
    keep = []
    for c in resid_columns(df, prefix):
        if compare_groups(df, c, min_n=min_n) is not None:
            keep.append(int(c.split("_")[1]))
    return sorted(keep)


# ── Box + jitter helper (matches the repo's plot_Rg_sasa.py style) ─────────────
def box_jitter(ax, df, col, groups, xpos_start=0, rng=None):
    """Draws a box-and-jitter plot for one column, one box per group.

    Args:
        ax: Matplotlib Axes to draw on.
        df: DataFrame containing `col` and a `seq_type` column.
        col (str): Column to plot.
        groups (list[str]): `seq_type` values to plot, in order, each
            colored via GROUP_COLOR.
        xpos_start (int): Starting x position (default: 0).
        rng: numpy Generator for jitter; a fixed-seed default is used if
            None.

    Returns:
        int: Next unused x position (xpos_start + len(groups)).
    """
    rng = rng or np.random.default_rng(42)
    xpos = xpos_start
    for group in groups:
        vals = df.loc[df["seq_type"] == group, col].dropna().values
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
    return xpos


# ── Contact-type features: Binder vs False Positive, ONE region at a time ─────
def contact_feature_region_analysis(df, region, out_dir):
    """Tests Binder vs False Positive contact-type composition, one region.

    Self-contained: nothing in this function's output references the other
    region, by design (see module docstring).

    Args:
        df: Contact-feature DataFrame for `region` (see
            load_contact_features).
        region (str): Region label, used in filenames/titles ("core" or
            "tail").
        out_dir (str): Directory to write the stats CSV and figure to.

    Returns:
        dict | None: Summary (region, n tested/significant, top feature)
        for the verdict table, or None if there wasn't enough data.
    """
    rows = []
    for feature in CONTACT_FEATURES:
        if feature not in df.columns:
            continue
        res = compare_groups(df, feature)
        if res is None:
            continue
        rows.append(dict(feature=feature, **res))
    if not rows:
        print(f"\n[Contact features - {region}] Not enough Binder/False Positive data — skipping.")
        return None
    stats_df = pd.DataFrame(rows)
    stats_df["q"] = bh_fdr(stats_df["p"].values)
    stats_path = os.path.join(out_dir, f"contact_features_{region}_binder_vs_fp_stats.csv")
    stats_df.sort_values("p").to_csv(stats_path, index=False)
    print(f"\nSaved {stats_path}")
    print(stats_df.sort_values("p").to_string(index=False))

    n_sig = int((stats_df["q"] < 0.05).sum())
    print(f"[Contact features - {region}] {n_sig}/{len(stats_df)} features "
          f"FDR-significant ({GROUP_A} vs {GROUP_B}, q<0.05)")

    # ── One figure, this region only: every feature, Binder vs False Positive ──
    features_present = stats_df["feature"].tolist()
    n = len(features_present)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                              dpi=300, constrained_layout=True, squeeze=False)
    for i, feature in enumerate(features_present):
        ax = axes[i // ncols][i % ncols]
        box_jitter(ax, df, feature, GROUP_ORDER)
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right", fontsize=8)
        row = stats_df[stats_df.feature == feature].iloc[0]
        star = "  *SIGNIFICANT*" if row["q"] < 0.05 else ""
        ax.set_title(f"{feature.replace('mean_', '')}\n"
                     f"p={row.p:.2g}, q={row.q:.2g}, d={row.cohens_d:.2f}{star}", fontsize=9)
        ax.grid(True, alpha=0.4)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(f"{region.upper()} REGION ONLY — contact-type composition by group\n"
                 f"{GROUP_A} vs {GROUP_B}: {n_sig}/{len(stats_df)} features significant (FDR q<0.05)",
                 fontsize=13)
    path = os.path.join(out_dir, f"contact_features_{region}_binder_vs_nonbinder.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved {path}")

    top = stats_df.sort_values("p").iloc[0]
    return dict(region=region, n_features_tested=len(stats_df), n_features_significant=n_sig,
                top_feature=top["feature"], top_feature_p=top["p"], top_feature_d=top["cohens_d"])


# ── Per-residue screens (R, D, or W): Binder vs False Positive, ONE region
# at a time ────────────────────────────────────────────────────────────────
def score_region_analysis(df, region, out_dir, prefix, label, top_n=6, restrict_to_resids=None):
    """Tests Binder vs False Positive per-residue scores, one region.

    Self-contained: nothing in this function's output references the other
    region, by design (see module docstring). Generalizes the original
    R-score-only screen so D and W (the two terms R = (D - W) * I is built
    from) can be screened the same way, without R's ambiguity where R = 0
    means either balanced direct/water-mediated contact or no contact at
    all. Each metric gets its own FDR family: all residues tested for
    *this* metric in *this* region, so D and W significance never inherits
    from R's or each other's p-values.

    Args:
        df: Score DataFrame with per-residue columns for `region` (see
            load_r_scores / load_dw_scores).
        region (str): Region label, used in filenames/titles ("core" or
            "tail").
        out_dir (str): Directory to write the stats CSV and figures to.
        prefix (str): Column-name prefix selecting the metric ("R", "D",
            or "W").
        label (str): Metric label used in titles/print statements (e.g.
            "R-score").
        top_n (int): Number of top-hit residues to plot individually
            (default: 6).
        restrict_to_resids (list[int] | None): If given, only test these
            resSeq values (used to screen D/W over the same pocket
            residues R was screened over; see contactable_residues).

    Returns:
        dict | None: Summary (region, metric, n tested/significant, top
        residue) for the verdict table, or None if there wasn't enough
        data.
    """
    stem = f"{prefix.lower()}_score"
    cols = resid_columns(df, prefix)
    if restrict_to_resids is not None:
        cols = [c for c in cols if int(c.split("_")[1]) in restrict_to_resids]
    rows = []
    for c in cols:
        res = compare_groups(df, c)
        if res is None:
            continue
        rows.append(dict(resSeq=int(c.split("_")[1]), **res))
    if not rows:
        print(f"\n[{label} - {region}] Not enough Binder/False Positive data — skipping.")
        return None
    stats_df = pd.DataFrame(rows)
    stats_df["q"] = bh_fdr(stats_df["p"].values)
    stats_path = os.path.join(out_dir, f"{stem}_{region}_binder_vs_fp_stats.csv")
    stats_df.sort_values("p").to_csv(stats_path, index=False)
    print(f"\nSaved {stats_path}")

    n_tested = len(stats_df)
    n_sig_raw = int((stats_df["p"] < 0.05).sum())
    n_sig_fdr = int((stats_df["q"] < 0.05).sum())
    print(f"[{label} - {region}] {n_tested} residues tested, {n_sig_raw} with p<0.05 "
          f"(uncorrected), {n_sig_fdr} with q<0.05 (FDR-corrected)")

    # ── Significance plot, this region only. Point color = which group is
    # higher (not region — there's only one region in this plot) ──
    fig, ax = plt.subplots(figsize=(max(10, stats_df["resSeq"].max() * 0.06), 5),
                            dpi=300, constrained_layout=True)
    colors = [GROUP_COLOR[GROUP_A] if d >= 0 else GROUP_COLOR[GROUP_B] for d in stats_df["cohens_d"]]
    ax.scatter(stats_df["resSeq"], -np.log10(stats_df["p"]), c=colors, s=25, alpha=0.85, zorder=2)
    sig = stats_df[stats_df["q"] < 0.05]
    if len(sig):
        sig_colors = [GROUP_COLOR[GROUP_A] if d >= 0 else GROUP_COLOR[GROUP_B] for d in sig["cohens_d"]]
        ax.scatter(sig["resSeq"], -np.log10(sig["p"]), facecolors="none", edgecolors=sig_colors,
                   s=110, linewidths=2, zorder=3)
        for _, row in sig.iterrows():
            ax.annotate(int(row["resSeq"]), (row["resSeq"], -np.log10(row["p"])),
                        textcoords="offset points", xytext=(0, 6), fontsize=7, ha="center")
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=1, label="p=0.05 (uncorrected)")
    ax.set_xlabel("Residue (resSeq)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(f"{region.upper()} REGION ONLY — per-residue {label}, {GROUP_A} vs {GROUP_B}\n"
                 f"{n_sig_fdr}/{n_tested} residues significant (FDR q<0.05); "
                 "open circle = FDR-significant")
    ax.grid(True, alpha=0.4)
    legend_handles = [
        mpatches.Patch(facecolor=GROUP_COLOR[GROUP_A], label=f"{GROUP_A} higher"),
        mpatches.Patch(facecolor=GROUP_COLOR[GROUP_B], label=f"{GROUP_B} higher"),
    ]
    ax.legend(handles=legend_handles, fontsize=8)
    path = os.path.join(out_dir, f"{stem}_{region}_significance.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved {path}")

    # ── Top-hit residues, this region only, one panel per residue ──
    top_hits = stats_df.sort_values("p").head(top_n)["resSeq"].tolist()
    print(f"[{label} - {region}] Top {len(top_hits)} residues by p-value: {top_hits}")
    if top_hits:
        n = len(top_hits)
        ncols = min(3, n)
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                                  dpi=300, constrained_layout=True, squeeze=False)
        for i, resSeq in enumerate(top_hits):
            ax = axes[i // ncols][i % ncols]
            col = f"{prefix}_{resSeq}"
            box_jitter(ax, df, col, GROUP_ORDER)
            ax.set_xticks(range(len(GROUP_ORDER)))
            ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right", fontsize=8)
            ax.axhline(0, color="grey", linewidth=0.6, alpha=0.5)
            row = stats_df[stats_df.resSeq == resSeq].iloc[0]
            ax.set_title(f"Residue {resSeq}\np={row.p:.2g}, q={row.q:.2g}, d={row.cohens_d:.2f}", fontsize=9)
            ax.grid(True, alpha=0.4)
        for j in range(n, nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.suptitle(f"{region.upper()} REGION ONLY — top {n} residues, "
                     f"{label} by group ({GROUP_A} vs {GROUP_B})", fontsize=13)
        path = os.path.join(out_dir, f"{stem}_{region}_top_residues.png")
        fig.savefig(path)
        plt.close(fig)
        print(f"Saved {path}")

    top = stats_df.sort_values("p").iloc[0]
    return dict(region=region, metric=label, n_residues_tested=n_tested, n_residues_significant=n_sig_fdr,
                top_residue=int(top["resSeq"]), top_residue_p=top["p"], top_residue_d=top["cohens_d"])


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    """Runs the core-vs-tail Binder vs False Positive analysis end to end."""
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--r-scores-dir", default="/projects/ivta1597/biosensors/water_analysis/agg_out")
    parser.add_argument("--contact-dir", default="/projects/ivta1597/biosensors/LIG_contacts")
    parser.add_argument("--seq-list", default="/projects/ivta1597/biosensors/seq_ids_orig.txt")
    parser.add_argument("--tag", default="40_500ns")
    parser.add_argument("--structure-source", default="all",
                        choices=["ngs_observed", "designed_assumed", "all"],
                        help="Filter sequences by md_candidate_guide.csv's source column. "
                             "'all' (default) applies no filtering.")
    parser.add_argument("--structure-guide",
                        default="/projects/ivta1597/biosensors/md_candidate_guide.csv",
                        help="Path to md_candidate_guide.csv (default: %(default)s)")
    parser.add_argument("--out-dir", default="/projects/ivta1597/biosensors/analysis/core_vs_tail")
    parser.add_argument("--top-n", type=int, default=6,
                        help="Number of top-hit residues to plot individually (default: 6)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    seq_type_map = load_seq_type_map(args.seq_list)

    print("=" * 70)
    print("Loading data (core/tail only)")
    print("=" * 70)
    r_data = load_r_scores(args.r_scores_dir, args.tag)
    dw_data = load_dw_scores(args.r_scores_dir, args.tag)
    contact_data = load_contact_features(args.contact_dir, args.tag, seq_type_map)

    if not r_data and not dw_data and not contact_data:
        raise FileNotFoundError(
            "No core/tail R-score, D/W-score, or contact-feature CSVs found. "
            "Check --r-scores-dir / --contact-dir / --tag."
        )

    if args.structure_source != "all":
        source_ids = load_source_ids(args.structure_guide, args.structure_source)
        print(f"\nRestricting to source == {args.structure_source}: "
              f"{len(source_ids)} in {args.structure_guide}")
        r_data = filter_to_source(r_data, source_ids, args.structure_source, "R-scores")
        dw_data = filter_to_source(dw_data, source_ids, args.structure_source, "D/W-scores")
        contact_data = filter_to_source(contact_data, source_ids, args.structure_source, "Contact feats")

    contact_summaries, r_summaries, d_summaries, w_summaries = [], [], [], []
    for region in REGIONS:
        if region in contact_data:
            print("\n" + "=" * 70)
            print(f"CONTACT FEATURES — {region.upper()} region only: {GROUP_A} vs {GROUP_B}")
            print("=" * 70)
            s = contact_feature_region_analysis(contact_data[region], region, args.out_dir)
            if s:
                contact_summaries.append(s)

    pocket_residues = {}
    for region in REGIONS:
        if region in r_data:
            print("\n" + "=" * 70)
            print(f"R-SCORE — {region.upper()} region only: {GROUP_A} vs {GROUP_B}, per residue")
            print("=" * 70)
            s = score_region_analysis(r_data[region], region, args.out_dir, "R", "R-score", top_n=args.top_n)
            if s:
                r_summaries.append(s)
            pocket_residues[region] = contactable_residues(r_data[region], "R")
            print(f"[pocket residues - {region}] {len(pocket_residues[region])} residues "
                  f"actually contact the ligand (>=5 Binder and >=5 FP sequences each); "
                  "D/W below are restricted to this same set.")

    # D and W are screened across the same pocket residues R was (see
    # contactable_residues), each with its own FDR family. This replaces an
    # earlier ad-hoc check that reused only the 5 residues R had already
    # flagged (optimistic FDR correction across just those 10 tests); this
    # screens the full pocket from scratch per metric, so a residue where D
    # and W cancel out in R but is real in D or W alone can still surface.
    for region in REGIONS:
        if region in dw_data:
            restrict = pocket_residues.get(region)
            print("\n" + "=" * 70)
            print(f"D OCCUPANCY — {region.upper()} region only: {GROUP_A} vs {GROUP_B}, per residue")
            print("=" * 70)
            s = score_region_analysis(dw_data[region], region, args.out_dir, "D", "D occupancy",
                                       top_n=args.top_n, restrict_to_resids=restrict)
            if s:
                d_summaries.append(s)

            print("\n" + "=" * 70)
            print(f"W OCCUPANCY — {region.upper()} region only: {GROUP_A} vs {GROUP_B}, per residue")
            print("=" * 70)
            s = score_region_analysis(dw_data[region], region, args.out_dir, "W", "W occupancy",
                                       top_n=args.top_n, restrict_to_resids=restrict)
            if s:
                w_summaries.append(s)

    # ── Verdict: a plain table, not a chart — each region's result stands on
    # its own above; this just puts the headline numbers side by side ──
    print("\n" + "=" * 70)
    print(f"VERDICT: {GROUP_A} vs {GROUP_B} difference, by region")
    print("=" * 70)
    if contact_summaries:
        cdf = pd.DataFrame(contact_summaries)
        cdf.to_csv(os.path.join(args.out_dir, "verdict_contact_features.csv"), index=False)
        print("\nContact-type features:")
        print(cdf.to_string(index=False))
    if r_summaries:
        rdf = pd.DataFrame(r_summaries)
        rdf.to_csv(os.path.join(args.out_dir, "verdict_r_scores.csv"), index=False)
        print("\nPer-residue R-scores:")
        print(rdf.to_string(index=False))
    if d_summaries:
        ddf = pd.DataFrame(d_summaries)
        ddf.to_csv(os.path.join(args.out_dir, "verdict_d_scores.csv"), index=False)
        print("\nPer-residue D occupancy:")
        print(ddf.to_string(index=False))
    if w_summaries:
        wdf = pd.DataFrame(w_summaries)
        wdf.to_csv(os.path.join(args.out_dir, "verdict_w_scores.csv"), index=False)
        print("\nPer-residue W occupancy:")
        print(wdf.to_string(index=False))

    print("\nDone. All outputs written to:", args.out_dir)


if __name__ == "__main__":
    main()
