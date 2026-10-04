#!/usr/bin/env python3
"""Builds the seq_ids list and arm metadata for the bile acid (BA) MD set.

"BA" (bile acid) follows the collaborator's BA_PYR1_XXXX design names; the
set covers LCA (LCAM in their files, the -1 monoanion) and LCA-3-sulfate.

Reads the collaborator's per-arm label table (orig_data/md_arm_set.csv; one
row per (design, ligand) Boltz-2 structure) and writes:

    seq_ids_ba.txt         standard 2-column seq_ids format (seq_id, label),
                           the single run list for both MD submission and
                           featurization. The MD submit scripts derive
                           prefix/id/dir_type by splitting the seq_id, so no
                           separate 4-column list is kept.
    ba_arm_metadata.csv    per-arm label, sequence, Boltz QC, and co-folding
                           scores, joined onto MD features in
                           build_ba_feat_table.py.

seq_id convention is <lig>_<id>_<suffix> (e.g. lca3s_0004_nb,
lcam_G0_binder). The ligand must be in the seq_id because 54 designs ship
with both ligands, and the MD scripts name run directories
<prefix>_<id>_<suffix>, so two arms of one design would otherwise collide.
Putting the ligand in the prefix (rather than a new seq_type) keeps the
existing binder/nb suffixes, so every hardcoded suffix -> subdirectory map in
the repo resolves BA runs without changes; only the base path differs.
That base comes from config_ba.yaml / the scripts' BASE override, not the
optional third seq_ids column: scripts treat that column as a full
per-sequence directory used as-is, and submit_multi_ID.sh skips any row that
has one.

Nonbinders get the "False Positive" label: every BA nonbinder arm was
selected as a co-folding false positive (its co-folding scores are balanced
against the binders'), matching that label's meaning in the LCA cohort.

Usage:
    python make_ba_seq_ids.py
"""
import argparse
import csv
import os
import re

REPO = os.path.dirname(os.path.abspath(__file__))

LIGAND_PREFIX = {"LCAM": "lcam", "LCA3S": "lca3s"}
CLASS_INFO = {
    # class -> (seq_ids label, seq_id suffix)
    "binder": ("Binder", "binder"),
    "nonbinder": ("False Positive", "nb"),
}

# Co-folding scores the collaborator balanced (or declared) per arm, kept for
# the co-folding-only baseline model in ba_ml_classification.py.
COFOLD_SCORES = [
    "plddt_ligand", "plddt_pocket", "binary_plddt_protein", "iptm",
    "binary_affinity_probability_binary", "binary_n_interface_unsatisfied",
    "n_valid_seeds", "geometry_score", "hbond_distance",
]
QC_COLS = [
    "binary_binding_mode", "needs_restrained_min", "pose_rmsd_spread",
    "seed_agreement", "flipped_fraction", "min_steric", "min_polar",
]


def parse_args() -> argparse.Namespace:
    """Parses CLI arguments."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--arm-csv", default=os.path.join(REPO, "orig_data", "md_arm_set.csv"))
    p.add_argument("--out-dir", default=REPO)
    return p.parse_args()


def design_short_id(sequence_id: str) -> str:
    """Converts a design name to the short id used in run names.

    Args:
        sequence_id: Design name from md_arm_set.csv, e.g. "BA_PYR1_0004",
            "POS_G0", or "POS_seq26".

    Returns:
        "0004", "G0", or "seq26" respectively.

    Raises:
        ValueError: The name matches neither the BA_PYR1_ nor POS_ pattern.
    """
    m = re.fullmatch(r"(?:BA_PYR1|POS)_(\w+)", sequence_id)
    if not m:
        raise ValueError(f"Unrecognized sequence_id: {sequence_id}")
    return m.group(1)


def main() -> None:
    """Writes the seq_ids list and the arm metadata table."""
    args = parse_args()
    with open(args.arm_csv, newline="") as f:
        rows = list(csv.DictReader(f))

    std_lines, meta_rows, seen = [], [], set()
    for r in rows:
        prefix = LIGAND_PREFIX[r["ligand"]]
        short_id = design_short_id(r["sequence_id"])
        label, suffix = CLASS_INFO[r["class"]]
        if (r["class"] == "binder") != (r["binds_this_ligand"] == "1"):
            raise ValueError(f"class/binds_this_ligand disagree for {r['sequence_id']} {r['ligand']}")
        seq_id = f"{prefix}_{short_id}_{suffix}"
        if seq_id in seen:
            raise ValueError(f"Duplicate seq_id: {seq_id}")
        seen.add(seq_id)

        std_lines.append(f"{seq_id}\t{label}")
        meta = {
            "seq_id": seq_id,
            "design": r["sequence_id"],
            "Ligand": r["ligand"],
            "Label": int(r["binds_this_ligand"]),
            "Group": label,
            "manifest_label": r["manifest_label"],
            "Sequence": r["protein_sequence"],
            "pocket_sequence": r["pocket_sequence"],
            "pdb": os.path.join("orig_data", "LCA_LCA3S", os.path.basename(r["pdb"])),
        }
        meta.update({c: r[c] for c in QC_COLS + COFOLD_SCORES})
        meta_rows.append(meta)

    with open(os.path.join(args.out_dir, "seq_ids_ba.txt"), "w") as f:
        f.write("\n".join(std_lines) + "\n")
    with open(os.path.join(args.out_dir, "ba_arm_metadata.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(meta_rows[0]))
        w.writeheader()
        w.writerows(meta_rows)

    n_bind = sum(m["Label"] for m in meta_rows)
    print(f"{len(meta_rows)} arms ({n_bind} binder, {len(meta_rows) - n_bind} nonbinder) "
          f"-> seq_ids_ba.txt, ba_arm_metadata.csv")


if __name__ == "__main__":
    main()
