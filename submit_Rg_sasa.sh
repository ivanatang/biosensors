#!/bin/bash
#SBATCH --job-name=rg_sasa
#SBATCH --output=logs/output_%A_%a.out
#SBATCH --error=logs/error_%A_%a.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-128%20

# Usage: sbatch --array=1-<N>%20 submit_rg_sasa.sh [seq_ids.txt]
#   N       = number of lines in seq_ids.txt  (wc -l seq_ids.txt)
#   %20     = max 20 jobs running simultaneously
#   seq_ids.txt can be passed as argument or defaults to ./seq_ids.txt

set -uo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
SEQ_IDS="${1:-seq_ids.txt}"

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
custom_base=$(echo "$line" | cut -f3)   # empty string if column absent

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
for f in medoid_PL.pdb PL_only_40_500ns.xtc index.ndx; do
    [[ ! -f "${rundir}/${f}" ]] && missing+=" $f"
done
if [[ -n "$missing" ]]; then
    echo "ERROR: missing inputs:$missing"
    exit 1
fi

cd "$rundir"

# ── Rg ────────────────────────────────────────────────────────────────────────
if [[ -f "Rg_PL.xvg" ]]; then
    echo "SKIP: Rg_PL.xvg already exists"
else
    echo "Running gmx gyrate..."
    echo "Protein_LIG" | "$GMX" gyrate \
        -s medoid_PL.pdb \
        -f PL_only_40_500ns.xtc \
        -n index.ndx \
        -o Rg_PL.xvg \
        2>> gyrate_rg.log || true

    if [[ -f "Rg_PL.xvg" ]]; then
        echo "  → Rg_PL.xvg done"
    else
        echo "  → ERROR: Rg_PL.xvg not created (check gyrate_rg.log)"
    fi
fi

# ── SASA ──────────────────────────────────────────────────────────────────────
if [[ -f "sasa_PL.xvg" ]]; then
    echo "SKIP: sasa_PL.xvg already exists"
else
    echo "Running gmx sasa..."
    echo "Protein_LIG" | "$GMX" sasa \
        -s medoid_PL.pdb \
        -f PL_only_40_500ns.xtc \
        -n index.ndx \
        -o sasa_PL.xvg \
        2>> sasa.log || true

    if [[ -f "sasa_PL.xvg" ]]; then
        echo "  → sasa_PL.xvg done"
    else
        echo "  → ERROR: sasa_PL.xvg not created (check sasa.log)"
    fi
fi

echo "Task ${SLURM_ARRAY_TASK_ID} complete"
