"""
select_binding_pocket.py

Identifies the LCA binding site pocket from mdpocket exploration output for
each sequence in a seq_ids.txt file, by finding alpha spheres in the frequency
isosurface PDB within a distance cutoff of the ligand centroid.

Outputs a selected_pocket.pdb per sequence for use in mdpocket characterization:
    mdpocket --trajectory_file protein_only.xtc \
             --trajectory_format xtc \
             --selected_pocket selected_pocket.pdb \
             -f protein_only_ref.pdb \
             -o mdpocket_<seq_id>

Usage:
    python select_binding_pocket.py [seq_ids_orig.txt] [--cutoff 8.0] [--ligand_resname LIG]
                                     [--overwrite-existing]
"""

import os
import sys
import argparse
import numpy as np

# ── Configurable paths ────────────────────────────────────────────────────────
# medoid_PL.pdb is a raw MD-run file, never written to scratch by this
# pipeline, so it's looked up in the PetaLibrary archive only. The freq_iso
# PDB (mdpocket exploration output) is normally a scratch-only pipeline
# artifact, but for sequences processed a while ago it can itself have aged
# out of scratch's 90-day window -- if a prior pipeline run's outputs were
# captured in the same archive snapshot, fall back to archive for it too.
# selected_pocket.pdb (this script's own output) always writes to scratch.
# Mirrors the same fix already applied to water_analysis/R_score_calc.py and
# LIG_contacts/contact_type_analysis.py.
ARCHIVE_BASE = "/pl/active/shirts_archive/IvanaTang/biosensors"
BASE         = "/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
# ─────────────────────────────────────────────────────────────────────────────


def get_dir_type(seq_type):
    """Maps a seq_type label to its data subdirectory.

    Args:
        seq_type (str): One of "Binder", "False Positive", "Low
            Confidence", "Fail Geometry".

    Returns:
        str: Subdirectory name, or `seq_type` unchanged if unrecognized.
    """
    mapping = {
        "Binder":         "binders",
        "False Positive": "nonbinders",
        "Low Confidence": "neg_low_pkt",
        "Fail Geometry":  "neg_fail_gate",
    }
    return mapping.get(seq_type, seq_type)


def resolve_input_file(candidate_dirs, filename):
    """Finds a file across an ordered list of candidate directories.

    Directories don't necessarily hold every file for a sequence: some
    files landed directly in runrel/ while others for that same sequence
    live under runrel/500ns/, and scratch-only pipeline outputs older than
    90 days may only survive in the archive snapshot. Checks each candidate
    independently rather than resolving one shared directory. Matches
    pkt_vol_prep.sh's resolve_input_file.

    Args:
        candidate_dirs (list[str]): Directories to check, in order.
        filename (str): Filename to look for in each directory.

    Returns:
        str | None: Path to the first match found, or None.
    """
    for d in candidate_dirs:
        path = os.path.join(d, filename)
        if os.path.exists(path):
            return path
    return None


def parse_pdb_atoms(pdb_path, record_types=("ATOM", "HETATM")):
    """Parses atom coordinates from a PDB file.

    Args:
        pdb_path (str): Path to the PDB file.
        record_types (tuple[str]): Record types to include (default:
            ("ATOM", "HETATM")).

    Returns:
        tuple[numpy.ndarray, list[str]]: (n, 3) coordinate array (empty if
        no matching atoms) and the corresponding raw PDB lines.
    """
    coords = []
    lines  = []
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() in record_types:
                try:
                    coords.append([float(line[30:38]),
                                   float(line[38:46]),
                                   float(line[46:54])])
                    lines.append(line)
                except ValueError:
                    continue
    return np.array(coords) if coords else np.empty((0, 3)), lines


def get_ligand_coords(pdb_path, ligand_resname):
    """Extracts coordinates of all ligand atoms matching a residue name.

    Args:
        pdb_path (str): Path to the PDB file.
        ligand_resname (str): Residue name identifying the ligand.

    Returns:
        numpy.ndarray: (n, 3) coordinate array of matching atoms.

    Raises:
        ValueError: No atoms with `ligand_resname` found.
    """
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() in ("ATOM", "HETATM"):
                if line[17:20].strip() == ligand_resname:
                    try:
                        coords.append([float(line[30:38]),
                                       float(line[38:46]),
                                       float(line[46:54])])
                    except ValueError:
                        continue
    if not coords:
        raise ValueError(
            f"No atoms with resname '{ligand_resname}' found in {pdb_path}.\n"
            f"  Check with: grep HETATM {pdb_path} | awk '{{print $4}}' | sort -u"
        )
    return np.array(coords)


def process_sequence(seq_id, archive_flat_dir, run_dir, cutoff, ligand_resname, overwrite=False):
    """Selects the binding pocket for a single sequence and writes it out.

    Args:
        seq_id (str): Sequence identifier, used in log messages.
        archive_flat_dir (str): PetaLibrary archive directory holding
            medoid_PL.pdb (and, as a fallback, the freq_iso PDB).
        run_dir (str): scratch directory holding the freq_iso PDB and
            where selected_pocket.pdb is written.
        cutoff (float): Distance cutoff (A) from the ligand centroid.
        ligand_resname (str): Residue name identifying the ligand.
        overwrite (bool): Regenerate selected_pocket.pdb even if it already
            exists (default: False).

    Returns:
        str: "ok", "skip" (already done, or a required input is missing),
        or "fail" (an input was invalid, or no spheres passed the cutoff).
    """
    out_pdb = os.path.join(run_dir, "selected_pocket.pdb")

    archive_nested_dir = os.path.join(archive_flat_dir, "500ns")

    # ── Validate inputs ───────────────────────────────────────────────────────
    if not os.path.isdir(archive_flat_dir):
        print(f"  SKIP: archive directory not found: {archive_flat_dir}")
        return "skip"

    pl_pdb = resolve_input_file([archive_flat_dir, archive_nested_dir], "medoid_PL.pdb")
    if pl_pdb is None:
        print(f"  SKIP: medoid_PL.pdb not found in {archive_flat_dir} or {archive_nested_dir}")
        return "skip"

    freq_name = f"mdpocket_{seq_id}_freq_iso_0_5.pdb"
    freq_pdb  = resolve_input_file([run_dir, archive_flat_dir, archive_nested_dir], freq_name)
    if freq_pdb is None:
        print(f"  SKIP: {freq_name} not found in {run_dir}, {archive_flat_dir}, or {archive_nested_dir}")
        return "skip"

    # ── Skip if already done, unless overwriting ──────────────────────────────
    if os.path.exists(out_pdb) and not overwrite:
        print(f"  SKIP: selected_pocket.pdb already exists (use --overwrite-existing to regenerate)")
        return "skip"

    os.makedirs(run_dir, exist_ok=True)

    # ── Ligand centroid ───────────────────────────────────────────────────────
    try:
        lig_coords = get_ligand_coords(pl_pdb, ligand_resname)
    except ValueError as e:
        print(f"  FAIL: {e}")
        return "fail"

    centroid = lig_coords.mean(axis=0)
    print(f"  Ligand atoms    : {len(lig_coords)}")
    print(f"  Ligand centroid : ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f}) Å")

    # ── Load alpha spheres ────────────────────────────────────────────────────
    pocket_coords, pocket_lines = parse_pdb_atoms(freq_pdb)
    if len(pocket_coords) == 0:
        print(f"  FAIL: no atoms found in {freq_pdb}")
        return "fail"

    print(f"  Alpha spheres   : {len(pocket_coords)} total")

    # ── Distance filter ───────────────────────────────────────────────────────
    dists = np.linalg.norm(pocket_coords - centroid, axis=1)
    mask  = dists <= cutoff
    n_kept = mask.sum()

    print(f"  Cutoff          : {cutoff} Å  →  {n_kept} spheres kept")

    if n_kept == 0:
        nearest = dists.min()
        print(f"  WARN: nearest sphere is {nearest:.2f} Å away — "
              f"try --cutoff {int(np.ceil(nearest + 2))}")
        return "fail"

    # ── Write selected pocket ─────────────────────────────────────────────────
    with open(out_pdb, "w") as f:
        for line in np.array(pocket_lines)[mask]:
            f.write(line)
        f.write("END\n")

    print(f"  → {out_pdb}")
    return "ok"


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    """Selects the binding pocket for every sequence in a seq list."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list",        nargs="?",
                        default="/projects/ivta1597/biosensors/seq_ids_orig.txt",
                        help="Tab-separated seq_ids file (default: seq_ids_orig.txt)")
    parser.add_argument("--cutoff",        type=float, default=8.0,
                        help="Distance cutoff in Å from ligand centroid (default: 8.0)")
    parser.add_argument("--ligand_resname",default="LIG",
                        help="Residue name of the ligand in medoid_PL.pdb (default: LIG)")
    parser.add_argument("--overwrite-existing", action="store_true",
                        help="Regenerate selected_pocket.pdb even if it already exists, "
                             "instead of skipping.")
    parser.add_argument("--suffix", default="",
                        help="Run-directory suffix, e.g. '_qfix' for the "
                             "bond-order/charge-fix systems (default: '', the "
                             "standard production directory). Applies to every "
                             "sequence in seq_list -- run _qfix sequences via a "
                             "separate seq_list from standard ones.")
    args = parser.parse_args()

    global RUNREL
    RUNREL = f"{RUNREL}{args.suffix}"

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    counts = {"ok": 0, "skip": 0, "fail": 0}

    with open(args.seq_list) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            seq_id      = parts[0].strip()
            seq_type    = parts[1].strip() if len(parts) > 1 else ""
            custom_path = parts[2].strip() if len(parts) > 2 else ""

            print(f"\n{seq_id}  [{seq_type}]")

            if custom_path:
                archive_flat_dir = os.path.join(custom_path, RUNREL)
                run_dir          = os.path.join(custom_path, RUNREL)
            else:
                dir_type         = get_dir_type(seq_type)
                archive_flat_dir = os.path.join(ARCHIVE_BASE, dir_type, seq_id, RUNREL)
                run_dir          = os.path.join(BASE, dir_type, seq_id, RUNREL)

            result = process_sequence(seq_id, archive_flat_dir, run_dir, args.cutoff,
                                       args.ligand_resname, overwrite=args.overwrite_existing)
            counts[result] += 1

    print(f"\n{'='*35}")
    print(f" Processed : {counts['ok']}")
    print(f" Skipped   : {counts['skip']}")
    print(f" Failed    : {counts['fail']}")
    print(f"{'='*35}")


if __name__ == "__main__":
    main()
