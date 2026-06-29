#!/bin/bash
#SBATCH --job-name=rg_pocket
#SBATCH --output=logs/output_%A_%a.out
#SBATCH --error=logs/error_%A_%a.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-128%20

# Usage: sbatch --array=1-<N>%20 submit_Rg_pocket.sh [seq_ids.txt]
#   N = number of lines in seq_ids.txt  (wc -l seq_ids.txt)

set -uo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
SEQ_IDS="${1:-seq_ids.txt}"

# 27 consensus pocket positions (from config.yaml: medoid.pocket_resids)
POCKET_RESIDS="59 60 61 62 79 81 83 87 88 89 91 92 94 108 109 110 115 117 120 122 141 158 159 160 163 164 167"

declare -A TYPE_SUBDIR=(
    ["Binder"]="binders"
    ["False Positive"]="nonbinders"
    ["Low Confidence"]="neg_low_pkt"
    ["Fail Geometry"]="neg_fail_gate"
)

# ── Load modules ──────────────────────────────────────────────────────────────
module purge
module load gcc
module load openmpi

# ── Get this task's line from seq_ids.txt ─────────────────────────────────────
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SEQ_IDS")
if [[ -z "$line" ]]; then
    echo "No entry for array task ${SLURM_ARRAY_TASK_ID} — exiting"
    exit 0
fi

folder_name=$(echo "$line" | cut -f1)
label=$(echo "$line"       | cut -f2)
custom_base=$(echo "$line" | cut -f3)

echo "Task ${SLURM_ARRAY_TASK_ID}: $folder_name ($label)"

# ── Resolve run directory ─────────────────────────────────────────────────────
subdir="${TYPE_SUBDIR[$label]:-}"
if [[ -z "$subdir" ]]; then
    echo "ERROR: unknown label '$label' — exiting"
    exit 1
fi

if [[ -n "$custom_base" ]]; then
    rundir="${custom_base}/${RUNREL}"
else
    rundir="${BASE}/${subdir}/${folder_name}/${RUNREL}"
fi

echo "rundir: $rundir"

if [[ ! -d "$rundir" ]]; then
    echo "ERROR: directory not found — $rundir"
    exit 1
fi

# ── Check required inputs ─────────────────────────────────────────────────────
missing=""
for f in medoid_PL.pdb PL_only_40_500ns.xtc; do
    [[ ! -f "${rundir}/${f}" ]] && missing+=" $f"
done
if [[ -n "$missing" ]]; then
    echo "ERROR: missing inputs:$missing"
    exit 1
fi

cd "$rundir"

# ── Skip if already done ──────────────────────────────────────────────────────
if [[ -f "Rg_pocket.xvg" ]]; then
    echo "SKIP: Rg_pocket.xvg already exists"
    exit 0
fi

# ── Build pocket-residue index (all atoms of the 27 pocket positions) ─────────
echo "Building pocket residue index..."
"$GMX" select \
    -s medoid_PL.pdb \
    -on pocket_res.ndx \
    -select "protein and resid $POCKET_RESIDS" \
    2>> gyrate_pocket.log

if [[ ! -f "pocket_res.ndx" ]]; then
    echo "ERROR: failed to create pocket_res.ndx (check gyrate_pocket.log)"
    exit 1
fi

# ── Run gmx gyrate on pocket residues ─────────────────────────────────────────
echo "Running gmx gyrate on pocket residues..."
echo "0" | "$GMX" gyrate \
    -s medoid_PL.pdb \
    -f PL_only_40_500ns.xtc \
    -n pocket_res.ndx \
    -o Rg_pocket.xvg \
    2>> gyrate_pocket.log || true

if [[ -f "Rg_pocket.xvg" ]]; then
    echo "  → Rg_pocket.xvg done"
else
    echo "  → ERROR: Rg_pocket.xvg not created (check gyrate_pocket.log)"
    exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID} complete"
