"""Extracts a per-sequence scalar pocket water density from gmx spatial output.

Parses each sequence's water_density_{seq_id}.cube (a Gaussian cube file
produced by `gmx spatial`: a 3D spatial density function of water-oxygen
occupancy, in the water_spatial_prep.sh/water_spatial_run.sh fitted-
trajectory reference frame) and computes the mean voxel density within a
fixed radius of the ligand's mean position.

The plain `pocket_water_density_*` columns average over every voxel in
that flat 8 A sphere, including voxels simply inside the ligand or a
pocket-lining side chain, both of which physically exclude water
regardless of whether the sequence is a binder -- so a chunk of any
"depletion" there is anatomy, not a hydrophobic-exclusion signal, diluting
any real group difference. The `pocket_water_density_accessible_*`
columns restrict the same sphere to voxels not within a protein/ligand
heavy atom's van der Waals radius (mdtraj Element.radius, Bondi-based) for
most of a sampled set of frames, i.e. genuinely solvent-accessible space
near the ligand, not atom-occupied volume.

Ligand centroid is computed directly from this pipeline's own
fit_trim.xtc + fit_trim_ref.pdb (mean LIG heavy-atom position across the
fitted trajectory) rather than the separate RMSD pipeline's medoid_PL.pdb,
so everything needed to place the region of interest comes from this
pipeline's own output, without depending on two independently-computed
least-squares fits staying bit-identical.

Cube file format (Gaussian cube, as written by gmx spatial): lines 1-2 are
comments; line 3 is NATOMS and origin_x/y/z (bohr); lines 4-6 are the
NX/NY/NZ voxel counts and per-axis step vectors (bohr); the next NATOMS
lines are atom records (ignored here); the remainder is NX*NY*NZ density
values, whitespace-separated across however many values-per-line gmx
spatial writes (observed: one full Z-row per line, not the 6-per-line
convention some cube writers use). parse_cube() flattens all remaining
whitespace-separated tokens regardless of line breaks, so this only
matters if you're eyeballing the raw file.

1 bohr = 0.529177 Angstrom. Voxel step vectors are assumed axis-aligned
(gmx spatial writes an orthogonal grid); off-diagonal step components are
ignored. This parser (units, atom-record skip count, Z-fastest/Y/X reshape
order) was verified against a real GROMACS 2025.3 `gmx spatial` run
(bind_019_binder): the parsed density array's mean/min/max matched gmx
spatial's own reported "Raw data: average ..., min ..., max ..." line
exactly.

water_spatial_run.sh runs `gmx spatial -nodiv`, so voxel values here are
raw per-frame occupancy counts, not already normalized to bulk-water
density (`-div` turned out to be a boolean visualization-only flag, not a
numeric divisor; see water_spatial_run.sh). This script does the bulk-
water normalization below, in Angstrom^3 units for consistency with the
rest of the repo (everything else here uses Angstrom, not nm).

Output: water_spatial_feats.csv

Usage:
    python extract_water_spatial_feats.py [seq_ids_ngs_observed.txt]
                                           [--pocket-cutoff 8.0] [--out water_spatial_feats.csv]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import mdtraj as md
from scipy.spatial import cKDTree

BOHR_TO_ANG = 0.529177210903

# Bulk liquid water molecular number density, ~0.0334 waters/Angstrom^3
# (equivalent to the commonly-cited ~33.3 nm^-3 / ~1 g/cm^3 for water).
# Used only to express pocket_water_density_mean as a "fraction of bulk"
# for interpretability -- verify against this system's actual water model
# if that distinction matters (TIP3P vs SPC/E etc. differ slightly).
BULK_WATER_DENSITY_PER_ANG3 = 0.0334

# Generic heavy-atom van der Waals radius (Angstrom, ~carbon's Bondi radius)
# used only if mdtraj's own per-element radius is missing/zero for some atom.
FALLBACK_ATOM_RADIUS_ANG = 1.70

# solvent_accessible_mask() samples this many evenly-strided frames from the
# ~12,000-frame production trajectory rather than checking every frame --
# the protein backbone is fit-aligned and relatively static in this
# reference frame, so a coarse stride is a reasonable cost/accuracy
# tradeoff for identifying which voxels are typically atom-occupied.
DEFAULT_N_ACCESSIBILITY_FRAMES = 50

# A voxel counts as atom-occupied (and is excluded) if it falls within a
# protein/ligand heavy atom's van der Waals radius in at least this fraction
# of the sampled frames -- majority-vote rather than "ever occupied", so a
# voxel only transiently grazed by a fluctuating side chain still counts as
# solvent-accessible.
DEFAULT_OCCUPIED_FRAME_FRAC = 0.5

TYPE_SUBDIR = {
    "Binder":         "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry":  "neg_fail_gate",
}

BASE   = "/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"


def parse_cube(cube_path):
    """Parses a Gaussian cube file written by `gmx spatial`.

    Args:
        cube_path (str): Path to the .cube file.

    Returns:
        tuple: (origin_ang (3,), step_ang (3,), dims (3,) int, density
        (NX,NY,NZ)), all distances in Angstroms (see module docstring for
        the cube format and unit conversion).
    """
    with open(cube_path) as f:
        lines = f.readlines()

    header = lines[2].split()
    natoms = int(header[0])
    origin_bohr = np.array([float(x) for x in header[1:4]])

    dims = np.zeros(3, dtype=int)
    step_bohr = np.zeros(3)
    for axis in range(3):
        parts = lines[3 + axis].split()
        dims[axis] = int(parts[0])
        # Only the diagonal step component is used -- see module docstring
        # caveat about the axis-aligned-grid assumption.
        step_bohr[axis] = float(parts[1 + axis])

    data_start = 6 + natoms
    # np.array(..., dtype=float) on the pre-split token list avoids building a
    # list of individually-boxed Python float objects (float(x) per token) --
    # matters here since -nab 300 (needed to fix a separate gmx spatial
    # memory error) inflates the production-scale cube's padded bounding box
    # well beyond what the short local test window produced.
    tokens = []
    for line in lines[data_start:]:
        tokens.extend(line.split())
    density = np.array(tokens, dtype=np.float64).reshape(dims)  # Z fastest already matches this shape/order

    origin_ang = origin_bohr * BOHR_TO_ANG
    step_ang   = step_bohr * BOHR_TO_ANG
    return origin_ang, step_ang, dims, density


def ligand_centroid_ang(traj, lig_resname="LIG"):
    """Computes the ligand's mean heavy-atom position, in the cube grid's frame.

    Takes an already-loaded mdtraj Trajectory (see main()) rather than
    loading fit_trim.xtc itself, since solvent_accessible_mask() needs the
    same ~4 GB trajectory loaded too; loading it once and sharing it avoids
    doubling memory/IO cost per sequence.

    Args:
        traj: mdtraj Trajectory, fitted and production-windowed (the same
            reference frame the cube grid was generated in).
        lig_resname (str): Ligand residue name (default: "LIG").

    Returns:
        numpy.ndarray: (3,) mean ligand heavy-atom position, in Angstroms.

    Raises:
        ValueError: No residue named `lig_resname` is found.
    """
    top = traj.topology
    lig_atoms = [a.index for r in top.residues if r.name == lig_resname
                 for a in r.atoms if a.element.symbol != 'H']
    if not lig_atoms:
        raise ValueError(f"No residue named '{lig_resname}' found in topology")
    com_per_frame = traj.xyz[:, lig_atoms, :].mean(axis=1)   # (F, 3), nm
    return com_per_frame.mean(axis=0) * 10.0                  # -> Angstrom


def _local_index_bounds(dims, origin_ang, step_ang, centroid_ang, cutoff_ang):
    """Computes the local sub-box index range around a centroid.

    Shared by pocket_density() and solvent_accessible_mask() so both
    operate over the exact same sub-box.

    Args:
        dims (numpy.ndarray): (3,) voxel grid dimensions (NX, NY, NZ).
        origin_ang (numpy.ndarray): (3,) grid origin, in Angstroms.
        step_ang (numpy.ndarray): (3,) per-axis voxel step, in Angstroms.
        centroid_ang (numpy.ndarray): (3,) center of the sub-box, in
            Angstroms.
        cutoff_ang (float): Half-width of the sub-box, in Angstroms.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (lo_idx, hi_idx), each (3,)
        inclusive index bounds, clamped to the grid.
    """
    lo_idx = np.zeros(3, dtype=int)
    hi_idx = np.zeros(3, dtype=int)
    for axis in range(3):
        lo = int(np.floor((centroid_ang[axis] - cutoff_ang - origin_ang[axis]) / step_ang[axis]))
        hi = int(np.ceil((centroid_ang[axis] + cutoff_ang - origin_ang[axis]) / step_ang[axis]))
        lo_idx[axis] = max(lo, 0)
        hi_idx[axis] = min(hi, dims[axis] - 1)
    return lo_idx, hi_idx


def pocket_density(origin_ang, step_ang, density, centroid_ang, cutoff_ang,
                    accessible_mask=None):
    """Computes mean voxel occupancy within a radius of the ligand centroid.

    Restricts the meshgrid/distance computation to a small index sub-box
    around the centroid instead of the full (NX,NY,NZ) grid. The full-grid
    version OOM'd in production: -nab 300 (needed to fix a separate
    `gmx spatial` "item outside of allocated memory" error) pads the
    cube's bounding box well beyond the pocket region actually needed
    here, so a full-grid meshgrid could be tens of GB even though only a
    few thousand voxels near the ligand are ever used.

    Args:
        origin_ang (numpy.ndarray): (3,) grid origin, in Angstroms.
        step_ang (numpy.ndarray): (3,) per-axis voxel step, in Angstroms.
        density (numpy.ndarray): (NX,NY,NZ) raw per-frame occupancy counts.
        centroid_ang (numpy.ndarray): (3,) ligand centroid, in Angstroms.
        cutoff_ang (float): Sphere radius around the centroid, in
            Angstroms.
        accessible_mask (numpy.ndarray | None): Same-shape boolean mask
            from solvent_accessible_mask(); if given, voxels it marks
            atom-occupied are excluded from the mean in addition to the
            flat distance cutoff (see module docstring).

    Returns:
        tuple[float, int, float]: (mean_count, n_voxels, voxel_volume_ang3)
        -- mean raw per-frame occupancy count within the cutoff (NaN if no
        voxel qualifies), the number of voxels averaged, and the voxel
        volume in Angstrom^3 needed to convert the mean into a number
        density.
    """
    dims = density.shape
    lo_idx, hi_idx = _local_index_bounds(dims, origin_ang, step_ang, centroid_ang, cutoff_ang)

    voxel_volume_ang3 = float(step_ang[0] * step_ang[1] * step_ang[2])
    if np.any(hi_idx < lo_idx):
        return np.nan, 0, voxel_volume_ang3

    sub_density = density[lo_idx[0]:hi_idx[0] + 1,
                           lo_idx[1]:hi_idx[1] + 1,
                           lo_idx[2]:hi_idx[2] + 1]
    ix, iy, iz = np.meshgrid(np.arange(lo_idx[0], hi_idx[0] + 1),
                              np.arange(lo_idx[1], hi_idx[1] + 1),
                              np.arange(lo_idx[2], hi_idx[2] + 1),
                              indexing="ij")
    voxel_xyz = origin_ang + np.stack([ix, iy, iz], axis=-1) * step_ang
    dist = np.linalg.norm(voxel_xyz - centroid_ang, axis=-1)
    mask = dist <= cutoff_ang
    if accessible_mask is not None:
        mask = mask & accessible_mask
    n_voxels = int(mask.sum())
    if n_voxels == 0:
        return np.nan, 0, voxel_volume_ang3
    return float(sub_density[mask].mean()), n_voxels, voxel_volume_ang3


def solvent_accessible_mask(traj, dims, origin_ang, step_ang, centroid_ang, cutoff_ang,
                             n_sample_frames=DEFAULT_N_ACCESSIBILITY_FRAMES,
                             occupied_frame_frac=DEFAULT_OCCUPIED_FRAME_FRAC):
    """Builds a solvent-accessibility mask over pocket_density()'s sub-box.

    A voxel is "solvent-accessible" (True) if it is NOT within a
    protein/ligand heavy atom's van der Waals radius in at least
    `occupied_frame_frac` of `n_sample_frames` evenly-strided frames from
    `traj`. Only protein and ligand heavy atoms count as occluding; water/
    ions are exactly what pocket_density() is measuring, not something to
    exclude. Uses a nearest-atom-only check (per-frame KD-tree): for atoms
    this close in size (~1.5-1.8 A heavy-atom radii), a voxel within a
    farther atom's radius while just outside its nearest atom's radius is
    a rare edge case, not worth the cost of an exact multi-atom overlap
    test. Cheap because it only runs over the same small local sub-box
    pocket_density() restricts to (a few thousand voxels), not the full
    padded cube grid, avoiding the same OOM concern that motivated that
    restriction.

    Args:
        traj: mdtraj Trajectory.
        dims (numpy.ndarray): (3,) voxel grid dimensions (NX, NY, NZ).
        origin_ang (numpy.ndarray): (3,) grid origin, in Angstroms.
        step_ang (numpy.ndarray): (3,) per-axis voxel step, in Angstroms.
        centroid_ang (numpy.ndarray): (3,) ligand centroid, in Angstroms.
        cutoff_ang (float): Sphere radius around the centroid, in
            Angstroms.
        n_sample_frames (int): Number of evenly-strided frames to sample
            (default: DEFAULT_N_ACCESSIBILITY_FRAMES).
        occupied_frame_frac (float): Fraction of sampled frames a voxel
            must be atom-occupied in to be excluded (default:
            DEFAULT_OCCUPIED_FRAME_FRAC).

    Returns:
        numpy.ndarray: Boolean array, same shape as pocket_density()'s
        local sub-box (empty if the sub-box is empty).

    Raises:
        ValueError: No protein/LIG heavy atoms are found in the topology.
    """
    top = traj.topology
    atom_indices = [a.index for a in top.atoms
                    if a.element.symbol != 'H'
                    and (a.residue.is_protein or a.residue.name == "LIG")]
    if not atom_indices:
        raise ValueError("No protein/LIG heavy atoms found for exclusion mask")

    def _atom_radius_ang(atom):
        r = atom.element.radius
        return r * 10.0 if r and r > 0 else FALLBACK_ATOM_RADIUS_ANG

    atom_radii_ang = np.array([_atom_radius_ang(top.atom(i)) for i in atom_indices])

    lo_idx, hi_idx = _local_index_bounds(dims, origin_ang, step_ang, centroid_ang, cutoff_ang)
    if np.any(hi_idx < lo_idx):
        return np.zeros(0, dtype=bool)

    ix, iy, iz = np.meshgrid(np.arange(lo_idx[0], hi_idx[0] + 1),
                              np.arange(lo_idx[1], hi_idx[1] + 1),
                              np.arange(lo_idx[2], hi_idx[2] + 1),
                              indexing="ij")
    voxel_xyz = origin_ang + np.stack([ix, iy, iz], axis=-1) * step_ang
    sub_shape = voxel_xyz.shape[:3]
    voxel_xyz_flat = voxel_xyz.reshape(-1, 3)

    n_frames_total = traj.n_frames
    stride = max(1, n_frames_total // n_sample_frames)
    frame_idxs = np.arange(0, n_frames_total, stride)

    occupied_count = np.zeros(voxel_xyz_flat.shape[0], dtype=int)
    for fidx in frame_idxs:
        atom_xyz_ang = traj.xyz[fidx, atom_indices, :] * 10.0  # nm -> Angstrom
        tree = cKDTree(atom_xyz_ang)
        dist, nearest_idx = tree.query(voxel_xyz_flat, k=1)
        occupied_count += (dist <= atom_radii_ang[nearest_idx]).astype(int)

    occupied_frac = occupied_count / len(frame_idxs)
    accessible_flat = occupied_frac < occupied_frame_frac
    return accessible_flat.reshape(sub_shape)


def main():
    """Extracts pocket water density features for every sequence in a seq list."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list", nargs="?",
                        default="/projects/ivta1597/biosensors/seq_ids_ngs_observed.txt")
    parser.add_argument("--pocket-cutoff", type=float, default=8.0,
                        help="Radius in Angstrom around the ligand centroid to average "
                             "voxel density over (default: 8.0, matching "
                             "pkt_vol/select_binding_pocket.py's pocket-sphere cutoff).")
    parser.add_argument("--base", default=BASE)
    parser.add_argument("--output", default="water_spatial_feats.csv")
    parser.add_argument("--skip-solvent-accessible", action="store_true",
                        help="Skip the pocket_water_density_accessible_* columns "
                             "(cheaper, but loses the atom-occlusion-corrected metric).")
    parser.add_argument("--n-accessibility-frames", type=int,
                        default=DEFAULT_N_ACCESSIBILITY_FRAMES,
                        help=f"Frames sampled to build the solvent-accessible mask "
                             f"(default: {DEFAULT_N_ACCESSIBILITY_FRAMES}).")
    parser.add_argument("--occupied-frac-threshold", type=float,
                        default=DEFAULT_OCCUPIED_FRAME_FRAC,
                        help=f"A voxel is excluded if atom-occupied in at least this "
                             f"fraction of sampled frames (default: "
                             f"{DEFAULT_OCCUPIED_FRAME_FRAC}).")
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

    records = []
    missing = []

    with open(args.seq_list) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts       = line.split("\t")
            seq_id      = parts[0].strip()
            seq_type    = parts[1].strip() if len(parts) > 1 else ""
            custom_path = parts[2].strip() if len(parts) > 2 else ""

            if custom_path:
                run_dir = os.path.join(custom_path, RUNREL, "water_spatial")
            else:
                dir_type = TYPE_SUBDIR.get(seq_type, seq_type)
                run_dir  = os.path.join(args.base, dir_type, seq_id, RUNREL, "water_spatial")

            cube_path = os.path.join(run_dir, f"water_density_{seq_id}.cube")
            ref_pdb   = os.path.join(run_dir, "fit_trim_ref.pdb")
            fit_xtc   = os.path.join(run_dir, "fit_trim.xtc")

            if not (os.path.exists(cube_path) and os.path.exists(ref_pdb)
                    and os.path.exists(fit_xtc)):
                print(f"MISSING: {seq_id}  [{seq_type}]  ->  {run_dir}")
                missing.append(seq_id)
                continue

            try:
                origin_ang, step_ang, dims, density = parse_cube(cube_path)
                traj = md.load(fit_xtc, top=ref_pdb)
                centroid_ang = ligand_centroid_ang(traj)
                mean_count, n_voxels, voxel_vol_ang3 = pocket_density(
                    origin_ang, step_ang, density, centroid_ang, args.pocket_cutoff)

                # mean_count is a raw per-frame occupancy count per voxel
                # (gmx spatial -nodiv). Convert to a number density
                # (waters/Angstrom^3) by dividing by voxel volume, then
                # express as a fraction of bulk for interpretability.
                density_per_ang3 = (mean_count / voxel_vol_ang3
                                     if voxel_vol_ang3 > 0 else np.nan)
                frac_bulk = (density_per_ang3 / BULK_WATER_DENSITY_PER_ANG3
                             if not np.isnan(density_per_ang3) else np.nan)

                record = {
                    "seq_id":   seq_id,
                    "seq_type": seq_type,
                    "pocket_water_density_mean": mean_count,
                    "pocket_water_density_per_ang3": density_per_ang3,
                    "pocket_water_density_frac_bulk": frac_bulk,
                    "pocket_water_density_n_voxels":  n_voxels,
                    "grid_dims": f"{dims[0]}x{dims[1]}x{dims[2]}",
                }
                log_line = (f"OK: {seq_id}  [{seq_type}]  "
                            f"mean_count={mean_count:.4f}  frac_bulk={frac_bulk:.4f}  "
                            f"n_voxels={n_voxels}")

                if not args.skip_solvent_accessible:
                    acc_mask = solvent_accessible_mask(
                        traj, dims, origin_ang, step_ang, centroid_ang, args.pocket_cutoff,
                        n_sample_frames=args.n_accessibility_frames,
                        occupied_frame_frac=args.occupied_frac_threshold)
                    acc_mean, acc_n_voxels, _ = pocket_density(
                        origin_ang, step_ang, density, centroid_ang, args.pocket_cutoff,
                        accessible_mask=acc_mask)
                    acc_density_per_ang3 = (acc_mean / voxel_vol_ang3
                                             if voxel_vol_ang3 > 0 and not np.isnan(acc_mean)
                                             else np.nan)
                    acc_frac_bulk = (acc_density_per_ang3 / BULK_WATER_DENSITY_PER_ANG3
                                      if not np.isnan(acc_density_per_ang3) else np.nan)
                    record.update({
                        "pocket_water_density_accessible_mean": acc_mean,
                        "pocket_water_density_accessible_per_ang3": acc_density_per_ang3,
                        "pocket_water_density_accessible_frac_bulk": acc_frac_bulk,
                        "pocket_water_density_accessible_n_voxels": acc_n_voxels,
                        # what fraction of the flat sphere was actually water-
                        # accessible at all -- a diagnostic on its own, since a
                        # low value here is exactly what explained the
                        # unexpectedly-low plain pocket_water_density_mean.
                        "pocket_water_density_frac_solvent_accessible":
                            (acc_n_voxels / n_voxels if n_voxels else np.nan),
                    })
                    log_line += (f"  |  accessible_frac_bulk={acc_frac_bulk:.4f}  "
                                 f"accessible_n_voxels={acc_n_voxels}/{n_voxels}")

                records.append(record)
                print(log_line)
            except Exception as e:
                print(f"ERROR: {seq_id}  -  {e}")
                missing.append(seq_id)

    if not records:
        print("\nNo sequences loaded - nothing to write.")
        sys.exit(1)

    feat_df = pd.DataFrame(records)
    feat_df["pocket_water_density_missing"] = False
    feat_df.to_csv(args.output, index=False)

    print(f"\nFeatures written to: {args.output}")
    print(f"  Sequences processed : {len(records)}")
    print(f"  Sequences missing   : {len(missing)}")
    if missing:
        print(f"  Missing seq_ids     : {', '.join(missing)}")


if __name__ == "__main__":
    main()
