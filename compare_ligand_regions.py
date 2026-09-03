"""Compares R-score and contact-type features across whole/core/tail LCA regions.

Tests the hypothesis (a collaborator's suggestion) that the steroid core and
carboxylate tail interact with the PYR1 pocket differently, and that
region-restricted metrics may separate Binders from nonbinders in ways the
whole-ligand analysis misses.

Inputs (produced by aggregate_r_scores.py / agg_contact_feats.py with
--ligand-region whole|core|tail):
    r_scores_all_sequences_{TAG}[_core|_tail].csv
    contact_features_all_{TAG}[_core|_tail].csv

A missing region is skipped with a warning rather than failing; all three
aren't required to get a comparison out of whichever regions are present.

Usage:
    python compare_ligand_regions.py \
        --r-scores-dir /projects/ivta1597/biosensors/water_analysis/agg_out \
        --contact-dir  /projects/ivta1597/biosensors/LIG_contacts \
        --seq-list     /projects/ivta1597/biosensors/seq_ids_orig.txt \
        --tag 40_500ns \
        --out-dir      /projects/ivta1597/biosensors/analysis/ligand_region

All p-values here are uncorrected for multiple comparisons (exploratory /
trend-finding, not a confirmatory test) — treat them as a ranking of
candidates worth a closer look, not a final result.
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import wilcoxon, mannwhitneyu

# ── Constants ─────────────────────────────────────────────────────────────────
GROUP_COLOR = {
    "Binder":         "#648FFF",
    "False Positive": "#DC267F",
    "Low Confidence": "#FE6100",
    "Fail Geometry":  "#FFB000",
}
GROUP_ORDER = ["Binder", "False Positive", "Low Confidence", "Fail Geometry"]

REGION_COLOR = {"whole": "#785EF0", "core": "#009E73", "tail": "#E69F00"}
REGION_ORDER = ["whole", "core", "tail"]

GROUP_A, GROUP_B = "Binder", "False Positive"   # primary separation test


# ── Loading ───────────────────────────────────────────────────────────────────
def region_suffix(region):
    """Maps a ligand region to its filename suffix.

    Args:
        region (str): "whole", "core", or "tail".

    Returns:
        str: "" for "whole", otherwise "_{region}".
    """
    return "" if region == "whole" else f"_{region}"


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
            seq_id = parts[0]
            seq_type = parts[1] if len(parts) > 1 else "Unknown"
            mapping[seq_id] = seq_type
    return mapping


def load_r_scores(r_scores_dir, tag):
    """Loads per-region R-score CSVs written by aggregate_r_scores.py.

    Args:
        r_scores_dir (str): Directory containing
            r_scores_all_sequences_{tag}[_region].csv files.
        tag (str): Run tag (e.g. "40_500ns").

    Returns:
        dict: Maps region ("whole"/"core"/"tail") to its DataFrame. Missing
        files are skipped with a printed warning.
    """
    out = {}
    for region in REGION_ORDER:
        path = os.path.join(r_scores_dir, f"r_scores_all_sequences_{tag}{region_suffix(region)}.csv")
        if os.path.exists(path):
            out[region] = pd.read_csv(path)
            print(f"Loaded {region:>5s} R-scores      : {path}  ({len(out[region])} sequences)")
        else:
            print(f"  (missing, skipping {region} R-scores) {path}")
    return out


def load_contact_features(contact_dir, tag, seq_type_map):
    """Loads per-region contact-feature CSVs written by agg_contact_feats.py.

    Args:
        contact_dir (str): Directory containing
            contact_features_all_{tag}[_region].csv files.
        tag (str): Run tag (e.g. "40_500ns").
        seq_type_map (dict): seq_id -> seq_type, used to add a `seq_type`
            column (see load_seq_type_map).

    Returns:
        dict: Maps region ("whole"/"core"/"tail") to its DataFrame. Missing
        files are skipped with a printed warning.
    """
    out = {}
    for region in REGION_ORDER:
        path = os.path.join(contact_dir, f"contact_features_all_{tag}{region_suffix(region)}.csv")
        if os.path.exists(path):
            df = pd.read_csv(path)
            df["seq_type"] = df["seq_id"].map(seq_type_map)
            n_unmapped = df["seq_type"].isna().sum()
            if n_unmapped:
                print(f"  WARNING: {n_unmapped} sequences in {region} contact features "
                      f"not found in --seq-list (seq_type will be NaN)")
            out[region] = df
            print(f"Loaded {region:>5s} contact feats : {path}  ({len(df)} sequences)")
        else:
            print(f"  (missing, skipping {region} contact features) {path}")
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


def resid_columns(df):
    """Lists a DataFrame's per-residue R-score columns, sorted by resSeq.

    Args:
        df: DataFrame with columns named "R_{resSeq}".

    Returns:
        list[str]: Matching column names, sorted numerically by resSeq.
    """
    return sorted((c for c in df.columns if c.startswith("R_")),
                  key=lambda c: int(c.split("_")[1]))


# ── R-score: per-residue region comparison ────────────────────────────────────
def r_score_region_comparison(r_data, out_dir, top_n=6):
    """Compares mean per-residue R-score between the core and tail regions.

    Plots mean R-score per residue for each available region, ranks
    residues by tail-minus-core delta, and runs a paired Wilcoxon test
    (core vs tail) per residue.

    Args:
        r_data (dict): Region -> R-score DataFrame (see load_r_scores).
            Requires both "core" and "tail".
        out_dir (str): Directory to write the plot and stats CSVs to.
        top_n (int): Number of most divergent residues to return from each
            end of the delta ranking (default: 6).

    Returns:
        list[int] | None: resSeq values of the top_n most core/tail-
        divergent residues (both directions), or None if "core" or "tail"
        is missing from `r_data`.
    """
    if "core" not in r_data or "tail" not in r_data:
        print("\n[R-score] Need both core and tail tables — skipping region comparison.")
        return None

    core_df, tail_df = r_data["core"], r_data["tail"]
    shared_ids = sorted(set(core_df["seq_id"]) & set(tail_df["seq_id"]))
    print(f"\n[R-score] {len(shared_ids)} sequences with both core and tail data")

    res_cols = sorted(set(resid_columns(core_df)) & set(resid_columns(tail_df)),
                       key=lambda c: int(c.split("_")[1]))
    core_idx = core_df.set_index("seq_id").loc[shared_ids, res_cols]
    tail_idx = tail_df.set_index("seq_id").loc[shared_ids, res_cols]

    mean_core = core_idx.mean(skipna=True)
    mean_tail = tail_idx.mean(skipna=True)
    delta = (mean_tail - mean_core).sort_values()

    # ── Plot: mean per-residue R-score, one line per region ──
    resnums = [int(c.split("_")[1]) for c in res_cols]
    fig, ax = plt.subplots(figsize=(max(10, len(res_cols) * 0.3), 4.5),
                            dpi=300, constrained_layout=True)
    ax.plot(resnums, mean_core.values, "o-", color=REGION_COLOR["core"],
            label="Core", markersize=3, linewidth=1)
    ax.plot(resnums, mean_tail.values, "o-", color=REGION_COLOR["tail"],
            label="Tail", markersize=3, linewidth=1)
    if "whole" in r_data:
        whole_df = r_data["whole"]
        whole_cols = [c for c in res_cols if c in whole_df.columns]
        whole_shared = [s for s in shared_ids if s in whole_df["seq_id"].values]
        mean_whole = whole_df.set_index("seq_id").loc[whole_shared, whole_cols].mean(skipna=True)
        ax.plot([int(c.split("_")[1]) for c in whole_cols], mean_whole.values, "--",
                color=REGION_COLOR["whole"], label="Whole", linewidth=1, alpha=0.7)
    ax.axhline(0, color="grey", linewidth=0.6, alpha=0.6)
    ax.set_xlabel("Residue (resSeq)")
    ax.set_ylabel("Mean R-score")
    ax.set_title(f"Mean per-residue R-score by ligand region (n={len(shared_ids)} sequences)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.4)
    path = os.path.join(out_dir, "r_score_by_region_per_residue.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved {path}")

    # ── Most core/tail-divergent residues ──
    div_table = pd.DataFrame({
        "resSeq": [int(c.split("_")[1]) for c in delta.index],
        "mean_R_core": mean_core[delta.index].values,
        "mean_R_tail": mean_tail[delta.index].values,
        "delta_tail_minus_core": delta.values,
    }).sort_values("delta_tail_minus_core")
    div_path = os.path.join(out_dir, "r_score_residue_deltas.csv")
    div_table.to_csv(div_path, index=False)
    print(f"Saved {div_path}  (all residues, sorted by tail-core delta)")

    top_div = pd.concat([div_table.head(top_n), div_table.tail(top_n)]).drop_duplicates()
    print(f"\nTop {len(top_div)} most core/tail-divergent residues:")
    print(top_div.to_string(index=False))

    # ── Paired Wilcoxon per residue: core R vs tail R ──
    stats_rows = []
    for c in res_cols:
        cvals, tvals = core_idx[c].values, tail_idx[c].values
        mask = ~np.isnan(cvals) & ~np.isnan(tvals)
        if mask.sum() < 8:
            continue
        try:
            _, p = wilcoxon(cvals[mask], tvals[mask])
        except ValueError:
            continue   # e.g. all paired differences are zero
        stats_rows.append(dict(resSeq=int(c.split("_")[1]), n_paired=int(mask.sum()),
                                mean_core=np.nanmean(cvals), mean_tail=np.nanmean(tvals),
                                wilcoxon_p=p))
    stats_df = pd.DataFrame(stats_rows).sort_values("wilcoxon_p")
    stats_path = os.path.join(out_dir, "r_score_core_vs_tail_wilcoxon.csv")
    stats_df.to_csv(stats_path, index=False)
    n_sig = (stats_df["wilcoxon_p"] < 0.05).sum()
    print(f"Saved {stats_path}")
    print(f"{n_sig}/{len(stats_df)} residues: core R-score significantly differs from "
          f"tail R-score (paired Wilcoxon, p<0.05, uncorrected)")

    return top_div["resSeq"].tolist()


def r_score_group_breakdown(r_data, resnums, out_dir):
    """Box-plots R-score by group x region for the given residues.

    One figure per residue in `resnums`, intended for the most
    core/tail-divergent residues (see r_score_region_comparison).

    Args:
        r_data (dict): Region -> R-score DataFrame (see load_r_scores).
        resnums (list[int]): resSeq values to plot.
        out_dir (str): Directory to write each figure to.
    """
    regions = [r for r in REGION_ORDER if r in r_data]
    for resSeq in resnums:
        col = f"R_{resSeq}"
        fig, ax = plt.subplots(figsize=(8, 5), dpi=300, constrained_layout=True)
        rng = np.random.default_rng(42)
        xpos, xticks, xticklabels = 0, [], []
        for region in regions:
            df = r_data[region]
            if col not in df.columns:
                continue
            for group in GROUP_ORDER:
                vals = df.loc[df["seq_type"] == group, col].dropna().values
                if len(vals) == 0:
                    continue
                color = GROUP_COLOR[group]
                ax.boxplot(vals, positions=[xpos], widths=0.6, patch_artist=True,
                           medianprops=dict(color="black", linewidth=1.5),
                           boxprops=dict(facecolor=color, alpha=0.5),
                           whiskerprops=dict(color=color), capprops=dict(color=color),
                           flierprops=dict(marker="", linestyle="none"))
                jitter = rng.uniform(-0.15, 0.15, len(vals))
                ax.scatter(xpos + jitter, vals, color=color, s=14, alpha=0.8, zorder=3)
                xticks.append(xpos)
                xticklabels.append(f"{region}\n{group}")
                xpos += 1
            xpos += 1   # gap between regions
        if not xticks:
            plt.close(fig)
            continue
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=7)
        ax.axhline(0, color="grey", linewidth=0.6, alpha=0.5)
        ax.set_ylabel("R-score")
        ax.set_title(f"Residue {resSeq}: R-score by group and ligand region")
        ax.grid(True, alpha=0.4)
        ax.legend(handles=[mpatches.Patch(facecolor=GROUP_COLOR[g], label=g) for g in GROUP_ORDER],
                  loc="best", fontsize=7)
        path = os.path.join(out_dir, f"r_score_resid{resSeq}_group_by_region.png")
        fig.savefig(path)
        plt.close(fig)
        print(f"Saved {path}")


def r_score_group_separation(r_data, out_dir):
    """Tests whether Binder vs False Positive separates better by region.

    Runs Mann-Whitney U per residue within each available region and flags
    residues significant under core/tail but not under whole.

    Args:
        r_data (dict): Region -> R-score DataFrame (see load_r_scores).
        out_dir (str): Directory to write the stats CSVs to.
    """
    regions = [r for r in REGION_ORDER if r in r_data]
    if len(regions) < 2:
        return
    rows = []
    for region in regions:
        df = r_data[region]
        for c in resid_columns(df):
            a = df.loc[df["seq_type"] == GROUP_A, c].dropna().values
            b = df.loc[df["seq_type"] == GROUP_B, c].dropna().values
            if len(a) < 5 or len(b) < 5:
                continue
            try:
                _, p = mannwhitneyu(a, b, alternative="two-sided")
            except ValueError:
                continue
            rows.append(dict(region=region, resSeq=int(c.split("_")[1]),
                              n_binder=len(a), n_fp=len(b),
                              mean_binder=a.mean(), mean_fp=b.mean(), mannwhitney_p=p))
    if not rows:
        print("\n[R-score] Not enough Binder/False Positive data for group separation test.")
        return
    sep_df = pd.DataFrame(rows)
    sep_path = os.path.join(out_dir, "r_score_binder_vs_fp_by_region.csv")
    sep_df.to_csv(sep_path, index=False)

    pivot = sep_df.pivot(index="resSeq", columns="region", values="mannwhitney_p")
    pivot_path = os.path.join(out_dir, "r_score_binder_vs_fp_pvalue_by_region.csv")
    pivot.to_csv(pivot_path)
    print(f"\nSaved {sep_path}")
    print(f"Saved {pivot_path}")

    if "whole" in pivot.columns:
        for region in [r for r in regions if r != "whole"]:
            newly_sig = pivot[(pivot[region] < 0.05) & (pivot["whole"] >= 0.05)].dropna(subset=[region])
            print(f"\n[R-score] Residues significant for {GROUP_A} vs {GROUP_B} under "
                  f"'{region}' but NOT under 'whole' (p<0.05, uncorrected): {len(newly_sig)}")
            if len(newly_sig):
                print(newly_sig.to_string())


# ── Contact-type composition comparison ────────────────────────────────────────
def contact_feature_region_comparison(contact_data, out_dir):
    """Compares contact-type composition across regions.

    Plots mean feature values by region, runs a paired Wilcoxon test
    (core vs tail) per feature if both are present, and tests whether
    Binder vs False Positive separation (Mann-Whitney U) differs by
    region.

    Args:
        contact_data (dict): Region -> contact-feature DataFrame (see
            load_contact_features). Needs at least two regions.
        out_dir (str): Directory to write the plot and stats CSVs to.
    """
    regions = [r for r in REGION_ORDER if r in contact_data]
    if len(regions) < 2:
        print("\n[Contact features] Need at least two regions — skipping.")
        return

    features = ["mean_frac_hydrophobic", "mean_n_total", "mean_n_hydrophobic",
                "mean_n_polar", "mean_n_pos_charged", "mean_n_neg_charged"]
    features = [f for f in features if all(f in contact_data[r].columns for r in regions)]
    if not features:
        print("\n[Contact features] None of the expected feature columns found — skipping.")
        return

    # ── Grouped bar: mean feature value by region ──
    fig, ax = plt.subplots(figsize=(9, 5), dpi=300, constrained_layout=True)
    width = 0.8 / len(regions)
    x = np.arange(len(features))
    for i, region in enumerate(regions):
        df = contact_data[region]
        means = [df[f].mean() for f in features]
        sems  = [df[f].std(ddof=1) / np.sqrt(len(df)) for f in features]
        ax.bar(x + i * width, means, width, yerr=sems, capsize=3,
               label=region.capitalize(), color=REGION_COLOR[region], alpha=0.85)
    ax.set_xticks(x + width * (len(regions) - 1) / 2)
    ax.set_xticklabels([f.replace("mean_", "") for f in features], rotation=25, ha="right")
    ax.set_ylabel("Mean value (per-sequence average over trajectory)")
    ax.set_title("Contact-type composition by ligand region")
    ax.grid(True, alpha=0.4)
    ax.legend()
    path = os.path.join(out_dir, "contact_features_by_region.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"\nSaved {path}")

    # ── Paired Wilcoxon: core vs tail, per feature ──
    if "core" in contact_data and "tail" in contact_data:
        core_df = contact_data["core"].set_index("seq_id")
        tail_df = contact_data["tail"].set_index("seq_id")
        shared = sorted(set(core_df.index) & set(tail_df.index))
        rows = []
        for f in features:
            cvals = core_df.loc[shared, f].values
            tvals = tail_df.loc[shared, f].values
            try:
                _, p = wilcoxon(cvals, tvals)
            except ValueError:
                p = np.nan
            rows.append(dict(feature=f, mean_core=cvals.mean(), mean_tail=tvals.mean(), wilcoxon_p=p))
        stats_df = pd.DataFrame(rows).sort_values("wilcoxon_p")
        stats_path = os.path.join(out_dir, "contact_features_core_vs_tail_wilcoxon.csv")
        stats_df.to_csv(stats_path, index=False)
        print(f"Saved {stats_path}")
        print(stats_df.to_string(index=False))

    # ── Does region restriction separate Binder from False Positive better? ──
    rows = []
    for region in regions:
        df = contact_data[region]
        for f in features:
            a = df.loc[df["seq_type"] == GROUP_A, f].dropna().values
            b = df.loc[df["seq_type"] == GROUP_B, f].dropna().values
            if len(a) < 3 or len(b) < 3:
                continue
            _, p = mannwhitneyu(a, b, alternative="two-sided")
            rows.append(dict(region=region, feature=f, n_binder=len(a), n_fp=len(b),
                              mean_binder=a.mean(), mean_fp=b.mean(), mannwhitney_p=p))
    if not rows:
        print(f"\n[Contact features] Not enough {GROUP_A}/{GROUP_B} data for group separation test.")
        return
    sep_df = pd.DataFrame(rows)
    sep_path = os.path.join(out_dir, "contact_features_binder_vs_fp_by_region.csv")
    sep_df.to_csv(sep_path, index=False)
    print(f"\nSaved {sep_path}")

    pivot = sep_df.pivot(index="feature", columns="region", values="mannwhitney_p")
    pivot_path = os.path.join(out_dir, "contact_features_binder_vs_fp_pvalue_by_region.csv")
    pivot.to_csv(pivot_path)
    print(f"Saved {pivot_path}")
    print(f"\n{GROUP_A} vs {GROUP_B} separation (Mann-Whitney p-value) by region:")
    print(pivot.to_string())


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    """Runs the whole/core/tail region comparison end to end."""
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--r-scores-dir", default="/projects/ivta1597/biosensors/water_analysis/agg_out",
                        help="Directory containing r_scores_all_sequences_{TAG}[_region].csv")
    parser.add_argument("--contact-dir", default="/projects/ivta1597/biosensors/LIG_contacts",
                        help="Directory containing contact_features_all_{TAG}[_region].csv")
    parser.add_argument("--seq-list", default="/projects/ivta1597/biosensors/seq_ids_orig.txt",
                        help="seq_id -> seq_type mapping, used for contact-feature group labels "
                             "(R-score CSVs already carry seq_type)")
    parser.add_argument("--tag", default="40_500ns", help="Time-window tag (default: 40_500ns)")
    parser.add_argument("--out-dir", default="/projects/ivta1597/biosensors/analysis/ligand_region",
                        help="Directory for output CSVs/PNGs")
    parser.add_argument("--top-n", type=int, default=6,
                        help="Number of most core/tail-divergent residues to plot individually "
                             "(default: 6, i.e. top 6 + bottom 6 by delta)")
    parser.add_argument("--structure-source", default="all",
                        choices=["ngs_observed", "designed_assumed", "all"],
                        help="Filter sequences by md_candidate_guide.csv's source column. "
                             "'all' (default) applies no filtering.")
    parser.add_argument("--structure-guide",
                        default="/projects/ivta1597/biosensors/md_candidate_guide.csv",
                        help="Path to md_candidate_guide.csv (default: %(default)s)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    seq_type_map = load_seq_type_map(args.seq_list)

    print("=" * 70)
    print("Loading R-score tables")
    print("=" * 70)
    r_data = load_r_scores(args.r_scores_dir, args.tag)

    print("\n" + "=" * 70)
    print("Loading contact-type feature tables")
    print("=" * 70)
    contact_data = load_contact_features(args.contact_dir, args.tag, seq_type_map)

    if not r_data and not contact_data:
        raise FileNotFoundError(
            "No R-score or contact-feature CSVs found. Check --r-scores-dir / "
            "--contact-dir / --tag."
        )

    if args.structure_source != "all":
        source_ids = load_source_ids(args.structure_guide, args.structure_source)
        print(f"\nRestricting to source == {args.structure_source}: "
              f"{len(source_ids)} in {args.structure_guide}")
        r_data = filter_to_source(r_data, source_ids, args.structure_source, "R-scores")
        contact_data = filter_to_source(contact_data, source_ids, args.structure_source, "Contact feats")

    print("\n" + "=" * 70)
    print("R-SCORE: core vs tail per-residue comparison")
    print("=" * 70)
    top_resnums = r_score_region_comparison(r_data, args.out_dir, top_n=args.top_n)
    if top_resnums:
        r_score_group_breakdown(r_data, top_resnums, args.out_dir)
    r_score_group_separation(r_data, args.out_dir)

    print("\n" + "=" * 70)
    print("CONTACT FEATURES: region comparison")
    print("=" * 70)
    contact_feature_region_comparison(contact_data, args.out_dir)

    print("\nDone. All outputs written to:", args.out_dir)


if __name__ == "__main__":
    main()
