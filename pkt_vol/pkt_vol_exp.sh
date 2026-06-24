#!/usr/bin/env bash
# =============================================================================
# run_mdpocket_exploration.sh  —  SLURM worker script
#
# Runs mdpocket exploration for a single sequence.
# cd into run directory before calling mdpocket to avoid malloc crash
# from long path lengths (known fpocket bug).
#
# Called by submit_mdpocket_exploration.sh:
#   sbatch run_mdpocket_exploration.sh <seq_id> <dir_type>
# =============================================================================

#SBATCH --job-name=mdpocket_exp
#SBATCH --output=logs/mdpocket_exp_%j.out
#SBATCH --error=logs/mdpocket_exp_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export TMPDIR=$SLURM_SCRATCH
export SLURM_EXPORT_ENV=ALL

module purge
module load gcc
module load openmpi
module load anaconda
conda activate fpocket_env

# ── Configurable paths ────────────────────────────────────────────────────────
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
PROT_XTC="protein_only.xtc"
PROT_PDB="protein_only.pdb"
# ─────────────────────────────────────────────────────────────────────────────

SEQ_ID=$1
DIR_TYPE=$2

if [[ -z "$SEQ_ID" || -z "$DIR_TYPE" ]]; then
    echo "ERROR: usage: sbatch run_mdpocket_exploration.sh <seq_id> <dir_type>"
    exit 1
fi

RUN_DIR="${BASE}/${DIR_TYPE}/${SEQ_ID}/${RUNREL}"
FREQ_ISO="${RUN_DIR}/mdpocket_${SEQ_ID}_freq_iso_0_5.pdb"

echo "seq_id   : $SEQ_ID"
echo "dir_type : $DIR_TYPE"
echo "run_dir  : $RUN_DIR"
echo ""

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -d "$RUN_DIR" ]]; then
    echo "ERROR: run directory not found: $RUN_DIR"; exit 1
fi
if [[ ! -f "${RUN_DIR}/${PROT_XTC}" ]]; then
    echo "ERROR: protein_only.xtc not found in $RUN_DIR"; exit 1
fi
if [[ ! -f "${RUN_DIR}/${PROT_PDB}" ]]; then
    echo "ERROR: protein_only.pdb not found in $RUN_DIR"; exit 1
fi

# ── Skip if already done ──────────────────────────────────────────────────────
if [[ -f "$FREQ_ISO" ]]; then
    echo "SKIP: freq_iso_0_5.pdb already exists"
    exit 0
fi

# ── Run mdpocket exploration ──────────────────────────────────────────────────
# cd into run directory — avoids malloc crash from long path lengths (fpocket bug)
echo "Running mdpocket exploration..."
cd "$RUN_DIR"

mdpocket \
    --trajectory_file   "$PROT_XTC"        \
    --trajectory_format xtc                \
    -f                  "$PROT_PDB"        \
    -o                  "mdpocket_${SEQ_ID}"

EXIT_CODE=$?

if [[ $EXIT_CODE -eq 0 && -f "$FREQ_ISO" ]]; then
    echo ""
    echo "OK: $FREQ_ISO"
    ls -lh mdpocket_${SEQ_ID}_*.pdb mdpocket_${SEQ_ID}_*.dx 2>/dev/null
else
    echo ""
    echo "FAILED: mdpocket exited with code $EXIT_CODE"
    exit 1
fi
