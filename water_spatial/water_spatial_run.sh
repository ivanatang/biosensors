#!/bin/bash

#SBATCH --job-name=water_spatial
#SBATCH --output=output_water_spatial_%j.out
#SBATCH --error=error_water_spatial_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# water_spatial_run.sh
#
# Builds a water-oxygen index group and runs `gmx spatial` on the sequence's
# fit_trim.xtc (produced by water_spatial_prep.sh) to generate a 3D spatial
# density function of water occupancy, output as a Gaussian .cube file.
#
# Requires water_spatial_prep.sh to have already produced
# <run_dir>/water_spatial/{protein_lig.ndx,fit_trim.xtc} for this sequence.
#
# Usage:
#   sbatch water_spatial_run.sh <seq_id> <dir_type> [bin_nm] [nab] [suffix]
#
# Arguments:
#   seq_id    - sequence identifier (e.g. pair_3059_binder)
#   dir_type  - directory group (binders | nonbinders | neg_low_pkt | neg_fail_gate)
#   bin_nm    - gmx spatial grid spacing in nm (default: 0.05 = 0.5 A)
#   suffix    - run-directory suffix, e.g. "_qfix" (default: ""), must match
#               what water_spatial_prep.sh was run with for this sequence
#   nab       - gmx spatial -nab ("number of ADDITIONAL BINS" per
#               `gmx spatial -h` -- padding added around the frame-0
#               bounding box is nab * bin_nm nm, NOT a fixed nm margin;
#               default here: 300, well above gmx's own default of 16).
#               VERIFIED against real bind_019_binder runs:
#               - gmx's own default -nab 16 failed even for a 21-frame/
#                 ~750ps test window ("item outside of the allocated
#                 memory"), the documented KNOWN ISSUES scenario in
#                 `gmx spatial -h`.
#               - -nab 100 (= 5 nm padding at bin=0.05) was sufficient for
#                 that short window but FAILED on the full ~12,000-frame,
#                 40-500ns production run -- a water molecule drifted
#                 ~0.02 nm past the padded box late in the trajectory,
#                 something the short test window had no time to expose.
#               If a job fails with "item outside of the allocated
#               memory", rerun with a higher nab (padding = nab * bin_nm).
#
# --time=08:00:00 below is an ESTIMATE, not a measurement: gmx spatial -h
# quotes ~30 min for a 50ns/32,000-atom trajectory; this system is a
# comparable ~29,265 atoms but the production window here is ~460ns (9.2x
# longer), so linear scaling suggests ~4.6h -- 8h leaves margin, but this
# has not been timed against a real full-length run.
#
# NOTE on -div: contrary to the original assumption here, `gmx spatial -h`
# shows -div is a BOOLEAN flag ("-[no]div", default yes), not a numeric
# bulk-density argument -- it auto-normalizes for VISUALIZATION only. This
# script uses -nodiv to get raw, physically-meaningful "counts per frame per
# voxel" instead; extract_water_spatial_feats.py does the actual bulk-water
# normalization in Python, where the divisor is an explicit, checkable
# constant rather than a black-box gmx internal computation. Verified
# against a real GROMACS 2025.3 run (bind_019_binder).
# ============================================================

set -euo pipefail

module purge
module load gcc
module load openmpi
module load anaconda
conda activate biosensors

# pandas/scipy (and other compiled-extension) imports need the envs own
# newer libstdc++, not the HPC images old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version GLIBCXX_3.4.29 not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

seq_id=$1
dir_type=$2
bin_nm=${3:-0.05}
nab=${4:-300}
suffix=${5:-}

# ── Configurable paths ────────────────────────────────────────────────────────
ARCHIVE_BASE="/pl/active/shirts_archive/IvanaTang/biosensors"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd${suffix}"
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
REF_TPR="prod_md_500ns.tpr"
# ─────────────────────────────────────────────────────────────────────────────

IN_FLAT_DIR="${ARCHIVE_BASE}/${dir_type}/${seq_id}/${RUNREL}"
RUN_DIR="${BASE}/${dir_type}/${seq_id}/${RUNREL}/water_spatial"

resolve_input_file() {
    local flat_dir="$1"
    local filename="$2"
    if [[ -f "${flat_dir}/${filename}" ]]; then
        echo "${flat_dir}/${filename}"
    elif [[ -f "${flat_dir}/500ns/${filename}" ]]; then
        echo "${flat_dir}/500ns/${filename}"
    else
        echo ""
    fi
}

TPR=$(resolve_input_file "$IN_FLAT_DIR" "$REF_TPR")
FIT_TRIM_XTC="${RUN_DIR}/fit_trim.xtc"
WATER_NDX="${RUN_DIR}/water_ow.ndx"
CUBE="${RUN_DIR}/water_density_${seq_id}.cube"

echo "============================================================"
echo "  gmx spatial water density"
echo "  seq_id      : $seq_id"
echo "  dir_type    : $dir_type"
echo "  bin (nm)    : $bin_nm"
echo "  nab         : $nab"
echo "  start       : $(date)"
echo "============================================================"

if [[ -z "$TPR" ]]; then
    echo "ERROR: $REF_TPR not found in $IN_FLAT_DIR or $IN_FLAT_DIR/500ns"
    exit 1
fi
if [[ ! -f "$FIT_TRIM_XTC" ]]; then
    echo "ERROR: $FIT_TRIM_XTC not found — run water_spatial_prep.sh for $seq_id first"
    exit 1
fi

# ── Step 1: build water-oxygen index group ──────────────────────────────────
# VERIFIED against a real .tpr (bind_019_binder, via the SOL .itp and .gro):
# residue name is "SOL" (as assumed) but the oxygen ATOM name is "O", not
# "OW" as originally assumed -- water in this system's topology uses generic
# per-element atom names (O, H, H), not the TIP3P-convention OW/HW1/HW2.
# `resname SOL and name OW` would have silently matched zero atoms.
#
# Also: `gmx select`'s "name = expression" naming syntax (e.g.
# 'water_ow = resname SOL and name O') fails outright with "Too few
# selections provided" in GROMACS 2025.3 -- confirmed empirically, not just
# assumed. Omit the name; gmx auto-generates one, and since it's the only
# group in this file it can be referenced by index (0) when piped to
# gmx spatial below.
echo ""
echo "── Step 1: build water-oxygen index group ──────────────────────────"
if [[ -f "$WATER_NDX" ]]; then
    echo "SKIP: $WATER_NDX already exists"
else
    "$GMX" select -s "$TPR" -on "$WATER_NDX" \
        -select 'resname SOL and name O'
    if [[ $? -ne 0 || ! -f "$WATER_NDX" ]]; then
        echo "ERROR: failed to build $WATER_NDX"
        exit 1
    fi
    echo "OK: $WATER_NDX"
fi

# ── Step 2: gmx spatial ──────────────────────────────────────────────────────
# `gmx spatial` has NO -o flag (contrary to the original assumption here) --
# it always writes a fixed "grid.cube" in the current working directory.
# Run it from inside RUN_DIR (already unique per sequence, so no collision
# risk across concurrent SLURM jobs) and rename the result afterward.
#
# It also prompts for TWO groups ("generate SDF" and "output coords, e.g.
# solute"), not one -- confirmed empirically; supplying only one answer left
# the second read uninitialized and produced nonsense huge bounding-box
# values ("item outside of allocated memory"). Both prompts get "0" (the
# only group in water_ow.ndx) since a separate small "solute" reference
# group isn't being generated here.
echo ""
echo "── Step 2: gmx spatial ─────────────────────────────────────────────"
if [[ -f "$CUBE" ]]; then
    echo "SKIP: $CUBE already exists"
else
    cd "$RUN_DIR"
    rm -f grid.cube
    printf '%s\n%s\n' "0" "0" | \
        "$GMX" spatial \
            -s "$TPR" -f "$FIT_TRIM_XTC" -n "$WATER_NDX" \
            -pbc -nodiv -bin "$bin_nm" -nab "$nab"

    if [[ ${PIPESTATUS[1]} -ne 0 || ! -f grid.cube ]]; then
        echo "ERROR: gmx spatial failed for $seq_id"
        exit 1
    fi
    mv grid.cube "$CUBE"
    echo "OK: $CUBE"
fi

echo ""
echo "============================================================"
echo "  Done: $seq_id"
echo "  End: $(date)"
echo "============================================================"
