#!/usr/bin/env python3
"""
build_feat_table.py

Assembles a self-contained feat_table_{N}ns.xlsx from the same per-sequence
feature-family CSVs that ML_classification.ipynb's first cell merges into
feat_table_500ns.xlsx at runtime -- but pointed at a requested analysis
window's outputs, and written to one fully-merged sheet (all_feats_{N}ns)
rather than requiring the same runtime merge every time the table is used.

Row/label anchor (name, ID, Group, Label, Binder_score, Sequence) is always
read from feat_table_500ns.xlsx's all_feats_500ns sheet and never re-derived,
so the 100ns/250ns/500ns tables stay row-for-row comparable. For --end-ns 500,
dw_pocket/contact_type/rmsf are reused directly from that same anchor sheet
(they're already pre-baked into it); for any other window they are computed
fresh from window-tagged CSVs, since no pre-built base sheet exists for
windows other than 500ns.

Usage:
    python build_feat_table.py --end-ns 250 --out feat_table_250ns.xlsx
    python build_feat_table.py --end-ns 100 --out feat_table_100ns.xlsx
    python build_feat_table.py --end-ns 500 --out feat_table_500ns_regen.xlsx   # validation only

Expects the same repo-relative input layout as ML_classification.ipynb's
first cell (analysis/, water_analysis/, water_analysis/agg_out/,
water_spatial/, salt_bridge/, LIG_contacts/), with each family's per-window
outputs already generated and, where they were computed on Alpine, synced
back into the local repo checkout.
"""
import os
import argparse
import pandas as pd

ANCHOR_TABLE = "feat_table_500ns.xlsx"
ANCHOR_SHEET = "all_feats_500ns"
ANCHOR_COLS  = ["name", "ID", "Group", "Label", "Binder_score", "Sequence"]

# Mirrors ML_classification.ipynb cell 1 exactly.
FLIPPED_LIGAND_NAMES = [
    "pair_3064_binder", "bind_033_binder", "bind_109_binder",
    "pair_1708_low_pkt", "nonb_012_nb", "nonb_020_nb", "nonb_055_nb",
]

POCKET_RESIDUES = [58, 59, 62, 83, 87, 88, 89, 92, 110, 115, 116, 117, 120, 122, 159, 160, 163, 164]
RMSD_REGIONS = ["Gate", "Latch"]

MD_GROUP_SUFFIX = {
    "binder":              "_binder",
    "non_binder":          "_nb",
    "negative_low_pocket": "_low_pkt",
    "negative_fail_gate":  "_fail_gate",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--end-ns", type=float, required=True,
                    help="End of analysis window in ns (100, 250, or 500)")
    p.add_argument("--out", required=True, help="Output .xlsx path")
    p.add_argument("--sheet-name", default=None,
                    help="Output sheet name (default: all_feats_{N}ns)")
    p.add_argument("--other-start-ns", type=float, default=40.0,
                    help="Start-ns used by every family except water_bridge "
                         "(default: 40, matches the existing 500ns table)")
    p.add_argument("--water-bridge-start-ns", type=float, default=0.0,
                    help="Start-ns used by the water_bridge family "
                         "(default: 0, matches the existing 500ns table)")
    return p.parse_args()


def report(df, cols, n_base, label, fillna=False):
    matched = df[cols[0]].notna().sum()
    print(f"  + {label:<22}: {matched}/{n_base} matched")
    if fillna:
        for c in cols:
            df[c] = df[c].fillna(0.0)


def main():
    args = parse_args()
    end_i      = int(args.end_ns)
    start_i    = int(args.other_start_ns)
    wb_start_i = int(args.water_bridge_start_ns)
    tag        = f"{start_i}_{end_i}ns"      # e.g. "40_250ns"
    wb_tag     = f"{wb_start_i}_{end_i}ns"   # e.g. "0_250ns"
    sheet_name = args.sheet_name or f"all_feats_{end_i}ns"
    is_500     = (end_i == 500 and start_i == 40)

    print(f"Building feat_table for window: {start_i}-{end_i}ns "
          f"(water_bridge: {wb_start_i}-{end_i}ns)")

    # ── Row/label anchor ──────────────────────────────────────────────────
    anchor_full = pd.read_excel(ANCHOR_TABLE, sheet_name=ANCHOR_SHEET)
    n_flipped = anchor_full["name"].isin(FLIPPED_LIGAND_NAMES).sum()
    anchor_full = anchor_full[~anchor_full["name"].isin(FLIPPED_LIGAND_NAMES)].reset_index(drop=True)
    n_base = len(anchor_full)
    print(f"Anchor cohort (from {ANCHOR_TABLE}::{ANCHOR_SHEET}): {n_base} sequences "
          f"({n_flipped} flipped-ligand QC exclusions removed)")

    df = anchor_full[ANCHOR_COLS].copy()

    # ── source (md_candidate_guide.csv) ─────────────────────────────────────
    mcg = pd.read_csv("md_candidate_guide.csv")
    mcg["name"] = mcg["pair_id"].astype(str) + mcg["md_group"].map(MD_GROUP_SUFFIX)
    df = df.merge(mcg[["name", "source"]], on="name", how="left")

    # ── dw_pocket / contact_type / rmsf ─────────────────────────────────────
    dw_cols = []
    for r in POCKET_RESIDUES:
        dw_cols += [f"D_{r}", f"W_{r}"]
    ct_type_cols = ["mean_n_neg_charged", "std_n_neg_charged", "occ_n_neg_charged_gt0"]
    rmsf_cols = ["Y23 RMSF (A)", "R79 RMSF (A)", "I110 RMSF (A)", "G163 RMSF (A)",
                 "Gate (r84-90) mean (A)", "Gate (r84-90) SD (A)"]

    if is_500:
        # Already pre-baked into the anchor sheet -- reuse directly.
        base_cols = dw_cols + ct_type_cols + rmsf_cols
        df = df.merge(anchor_full[["name"] + base_cols], on="name", how="left")
        report(df, dw_cols, n_base, "dw_pocket (base sheet)")
        report(df, ct_type_cols, n_base, "contact_type (base sheet)")
        report(df, rmsf_cols, n_base, "rmsf (base sheet)")
    else:
        # dw_pocket: whole-ligand D/W panel
        dw_path = f"water_analysis/agg_out/dw_scores_all_sequences_{tag}.csv"
        dw = pd.read_csv(dw_path)[["seq_id"] + dw_cols]
        df = df.merge(dw, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
        report(df, dw_cols, n_base, "dw_pocket")

        # contact_type
        ct_path = f"LIG_contacts/contact_features_all_{tag}.csv"
        ct = pd.read_csv(ct_path)[["seq_id"] + ct_type_cols]
        df = df.merge(ct, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
        report(df, ct_type_cols, n_base, "contact_type")

        # rmsf: single-residue + Ca region summary
        rmsf_single = pd.read_csv(f"analysis/rmsf_single_residues_per_seq_{end_i}ns.csv")
        rmsf_ca     = pd.read_csv(f"analysis/rmsf_ca_per_seq_summary_{end_i}ns.csv")
        rmsf = rmsf_single[["Sequence"] + [c for c in rmsf_single.columns if c in rmsf_cols]].merge(
            rmsf_ca[["Sequence"] + [c for c in rmsf_ca.columns if c in rmsf_cols]],
            on="Sequence", how="outer")
        df = df.merge(rmsf, left_on="name", right_on="Sequence", how="left").drop(columns=["Sequence"])
        report(df, rmsf_cols, n_base, "rmsf")

    # ── Gate-latch RMSD-to-reference ────────────────────────────────────────
    gl_cols = []
    for r in RMSD_REGIONS:
        gl_cols += [f"{r} RMSD mean (A)", f"{r} RMSD SD (A)", f"{r} drift100 (A)", f"{r} slope (A/ns)"]
    gl = pd.read_csv(f"analysis/gate_latch_rmsd_to_ref_summary_{end_i}ns.csv")
    gl = gl.rename(columns={"Sequence": "seq_id"})[["seq_id"] + gl_cols]
    df = df.merge(gl, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, gl_cols, n_base, "gate-latch RMSD-to-ref")

    # ── Salt bridges ─────────────────────────────────────────────────────────
    sb_cols = ["max_saltbridge_occupancy_pct", "n_saltbridges_gt50pct", "mean_top3_occupancy_pct"]
    sb_path = f"salt_bridge/saltbridge_features_all_seqs_{tag}.csv" if not is_500 \
        else "salt_bridge/saltbridge_features_all_seqs.csv"
    sb = pd.read_csv(sb_path)[["seq_id"] + sb_cols]
    df = df.merge(sb, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, sb_cols, n_base, "salt bridges", fillna=True)

    # ── Core vs. tail ligand-region D/W delta ───────────────────────────────
    core = pd.read_csv(f"water_analysis/agg_out/dw_scores_all_sequences_{tag}_core.csv")
    tail = pd.read_csv(f"water_analysis/agg_out/dw_scores_all_sequences_{tag}_tail.csv")
    dw_names = [f"D_{r}" for r in POCKET_RESIDUES] + [f"W_{r}" for r in POCKET_RESIDUES]
    core_sub = core[["seq_id"] + dw_names].rename(columns={c: f"{c}_core" for c in dw_names})
    tail_sub = tail[["seq_id"] + dw_names].rename(columns={c: f"{c}_tail" for c in dw_names})
    ct_delta = core_sub.merge(tail_sub, on="seq_id", how="inner")
    ct_cols = []
    for r in POCKET_RESIDUES:
        for letter in ("D", "W"):
            col = f"delta_{letter}_{r}"
            ct_delta[col] = ct_delta[f"{letter}_{r}_core"] - ct_delta[f"{letter}_{r}_tail"]
            ct_cols.append(col)
    df = df.merge(ct_delta[["seq_id"] + ct_cols], left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, ct_cols, n_base, "core/tail D/W delta", fillna=True)

    # ── Pocket-residue-referenced hydration ─────────────────────────────────
    hyd_cols = ["hydration_count_pocket_4A_mean", "hydration_count_pocket_4A_std",
                "hydration_count_pocket_4A_min", "hydration_count_pocket_4A_max",
                "hydration_count_pocket_4A_early20_mean", "hydration_count_pocket_4A_late20_mean",
                "hydration_count_pocket_4A_drift20", "hydration_count_pocket_4A_slope",
                "hydration_count_pocket_4A_slope_per_ns"]
    hyd_path = f"water_spatial/water_density_feats_pocket_{tag}.csv" if not is_500 \
        else "water_spatial/water_density_feats_pocket.csv"
    hyd = pd.read_csv(hyd_path)[["seq_id"] + hyd_cols]
    df = df.merge(hyd, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, hyd_cols, n_base, "pocket hydration (4A)")

    # ── Gate-latch-ligand water bridge ──────────────────────────────────────
    wb_cols = ["gate_bridge_occupancy", "triple_bridge_occupancy",
               "mean_n_triple_bridge_waters", "co_occurrence_occupancy",
               "mean_run_duration_ns"]
    wb_path = f"water_analysis/gate_latch_bridge_all_{wb_tag}_core.csv"
    wb = pd.read_csv(wb_path)[["seq_id"] + wb_cols]
    df = df.merge(wb, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, wb_cols, n_base, "gate-latch water bridge", fillna=True)

    print(f"\nFinal merged table: {df.shape[0]} sequences x {df.shape[1]} columns")

    feature_cols = dw_cols + ct_type_cols + rmsf_cols + gl_cols + sb_cols + ct_cols + hyd_cols + wb_cols
    print(f"Total candidate features: {len(feature_cols)}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with pd.ExcelWriter(args.out) as writer:
        df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"\nSaved -> {args.out}  (sheet: {sheet_name})")


if __name__ == "__main__":
    main()
