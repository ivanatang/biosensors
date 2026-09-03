#!/usr/bin/env python3
"""Builds a schema-matched feature table for the 6 qfix pilot sequences.

Covers bind_022/019/020_binder and nonb_006/008/009_nb, from their
bond-order-fixed reparameterization: same feature columns as
feat_table_500ns.xlsx, sourced from the _qfix-suffixed aggregated CSVs, for
comparison against those same 6 sequences' original (standard) feature
values and for a model-swap CV evaluation. Mirrors
build_new_ligand_feat_table.py's logic, with different sequences and
source paths.

Usage:
    python build_qfix_feat_table.py --seq_list seq_ids_qfix_pilot.txt \
        --out qfix_pilot_feat_table.csv
"""
import argparse
import pandas as pd

POCKET_RESIDUES = [58, 59, 62, 83, 87, 88, 89, 92, 110, 115, 116, 117, 120, 122, 159, 160, 163, 164]
RMSD_REGIONS = ["Gate", "Latch"]
TAG = "40_500ns"
WB_TAG = "0_500ns"


def parse_args():
    """Parses CLI args for the qfix pilot sequence list and output path.

    Returns:
        argparse.Namespace: Parsed arguments (seq_list, out).
    """
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seq_list", default="seq_ids_qfix_pilot.txt",
                    help="seq_ids.txt-style list of the qfix pilot sequences (default: %(default)s)")
    p.add_argument("--out", default="qfix_pilot_feat_table.csv", help="Output .csv path")
    return p.parse_args()


def report(df, cols, n_base, label):
    """Prints a match-rate summary for one merged feature family.

    Args:
        df: Merged feature table.
        cols (list[str]): Columns belonging to this feature family.
        n_base (int): Total row count, for the match-rate denominator.
        label (str): Feature-family name for the printed line.
    """
    matched = df[cols[0]].notna().sum()
    print(f"  + {label:<22}: {matched}/{n_base} matched")


def main():
    """Builds the qfix pilot feature table and writes it to CSV."""
    args = parse_args()

    rows = []
    with open(args.seq_list) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            name = parts[0]
            group_label = parts[1] if len(parts) > 1 else "Unknown"
            rows.append({"name": name, "Group": group_label,
                          "Label": 1 if group_label == "Binder" else 0})
    df = pd.DataFrame(rows)
    n_base = len(df)
    print(f"Qfix pilot cohort ({args.seq_list}): {n_base} sequences")

    dw_cols = []
    for r in POCKET_RESIDUES:
        dw_cols += [f"D_{r}", f"W_{r}"]
    ct_type_cols = ["mean_n_neg_charged", "std_n_neg_charged", "occ_n_neg_charged_gt0"]
    rmsf_cols = ["Y23 RMSF (A)", "R79 RMSF (A)", "I110 RMSF (A)", "G163 RMSF (A)",
                 "Gate (r84-90) mean (A)", "Gate (r84-90) SD (A)"]

    # ── dw_pocket: whole-ligand D/W panel ────────────────────────────────────
    dw = pd.read_csv(f"water_analysis/results_qfix_pilot/dw_scores_all_sequences_{TAG}.csv")[["seq_id"] + dw_cols]
    df = df.merge(dw, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, dw_cols, n_base, "dw_pocket")

    # ── contact_type ──────────────────────────────────────────────────────
    ct = pd.read_csv(f"LIG_contacts/contact_features_all_{TAG}_qfix.csv")[["seq_id"] + ct_type_cols]
    df = df.merge(ct, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, ct_type_cols, n_base, "contact_type")

    # ── rmsf: single-residue + Ca region summary ─────────────────────────────
    rmsf_single = pd.read_csv("analysis/rmsf_single_residues_per_seq_500ns_qfix.csv") \
        .rename(columns={"Sequence": "seq_id"})
    rmsf_ca = pd.read_csv("analysis/rmsf_ca_per_seq_summary_500ns_qfix.csv") \
        .rename(columns={"Sequence": "seq_id"})
    rmsf = rmsf_single[["seq_id"] + [c for c in rmsf_single.columns if c in rmsf_cols]].merge(
        rmsf_ca[["seq_id"] + [c for c in rmsf_ca.columns if c in rmsf_cols]],
        on="seq_id", how="outer")
    df = df.merge(rmsf, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, rmsf_cols, n_base, "rmsf")

    # ── Gate-latch RMSD-to-reference ─────────────────────────────────────────
    gl_cols = []
    for r in RMSD_REGIONS:
        gl_cols += [f"{r} RMSD mean (A)", f"{r} RMSD SD (A)", f"{r} drift100 (A)", f"{r} slope (A/ns)"]
    gl = pd.read_csv("analysis/gate_latch_rmsd_to_ref_summary_500ns_qfix.csv")
    gl = gl.rename(columns={"Sequence": "seq_id"})[["seq_id"] + gl_cols]
    df = df.merge(gl, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, gl_cols, n_base, "gate-latch RMSD-to-ref")

    # ── Salt bridges ──────────────────────────────────────────────────────
    sb_cols = ["max_saltbridge_occupancy_pct", "n_saltbridges_gt50pct", "mean_top3_occupancy_pct"]
    sb = pd.read_csv("salt_bridge/saltbridge_features_qfix_pilot.csv")[["seq_id"] + sb_cols]
    df = df.merge(sb, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, sb_cols, n_base, "salt bridges")
    df[sb_cols] = df[sb_cols].fillna(0.0)

    # ── Core vs. tail ligand-region D/W delta ────────────────────────────────
    core = pd.read_csv(f"water_analysis/results_qfix_pilot/dw_scores_all_sequences_{TAG}_core.csv")
    tail = pd.read_csv(f"water_analysis/results_qfix_pilot/dw_scores_all_sequences_{TAG}_tail.csv")
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
    report(df, ct_cols, n_base, "core/tail D/W delta")
    df[ct_cols] = df[ct_cols].fillna(0.0)

    # ── Pocket-residue-referenced hydration ──────────────────────────────────
    hyd_cols = ["hydration_count_pocket_4A_mean", "hydration_count_pocket_4A_std",
                "hydration_count_pocket_4A_min", "hydration_count_pocket_4A_max",
                "hydration_count_pocket_4A_early20_mean", "hydration_count_pocket_4A_late20_mean",
                "hydration_count_pocket_4A_drift20", "hydration_count_pocket_4A_slope",
                "hydration_count_pocket_4A_slope_per_ns"]
    hyd = pd.read_csv("water_spatial/water_density_feats_pocket_qfix_pilot.csv")[["seq_id"] + hyd_cols]
    df = df.merge(hyd, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, hyd_cols, n_base, "pocket hydration (4A)")

    # ── Gate-latch-ligand water bridge ───────────────────────────────────────
    wb_cols = ["gate_bridge_occupancy", "triple_bridge_occupancy",
               "mean_n_triple_bridge_waters", "co_occurrence_occupancy",
               "mean_run_duration_ns"]
    wb = pd.read_csv(f"water_analysis/gate_latch_bridge_all_{WB_TAG}_core_qfix.csv")[["seq_id"] + wb_cols]
    df = df.merge(wb, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    report(df, wb_cols, n_base, "gate-latch water bridge")
    df[wb_cols] = df[wb_cols].fillna(0.0)

    print(f"\nFinal qfix table: {df.shape[0]} sequences x {df.shape[1]} columns")

    feature_cols = dw_cols + ct_type_cols + rmsf_cols + gl_cols + sb_cols + ct_cols + hyd_cols + wb_cols
    print(f"Total feature columns (must match feat_table_500ns.xlsx's feature_cols exactly): {len(feature_cols)}")

    n_incomplete = df[feature_cols].isna().any(axis=1).sum()
    if n_incomplete:
        print(f"\nWARNING: {n_incomplete} sequence(s) have at least one missing feature -- "
              f"check the 'matched' counts above for which family is incomplete.")

    df.to_csv(args.out, index=False)
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
