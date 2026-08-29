#!/usr/bin/env bash
# =============================================================================
# run_mdpocket_exploration.sh  —  SLURM worker script
#
# Runs mdpocket exploration for a single sequence.
# cd into run directory before calling mdpocket to avoid malloc crash
# from long path lengths (known fpocket bug).
#
# Called by submit_mdpocket_exploration.sh:
#   sbatch run_mdpocket_exploration.sh <seq_id> <dir_type> [overwrite] [suffix]
#
# overwrite ("true"/"false", default false) forces mdpocket to rerun and
# removes existing mdpocket_<seq_id>_* outputs first, instead of skipping
# when freq_iso_0_5.pdb already has pocket atoms. Needed after a stale
# input gets fixed upstream (e.g. a PetaLibrary archive resync).
#
# suffix - run-directory suffix, e.g. "_qfix" (default: "")
# =============================================================================

#SBATCH --job-name=mdpocket_exp
#SBATCH --output=logs/mdpocket_exp_%j.out
#SBATCH --error=logs/mdpocket_exp_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
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
OVERWRITE=${3:-false}
SUFFIX=${4:-}
RUNREL="${RUNREL}${SUFFIX}"

if [[ -z "$SEQ_ID" || -z "$DIR_TYPE" ]]; then
    echo "ERROR: usage: sbatch run_mdpocket_exploration.sh <seq_id> <dir_type> [overwrite] [suffix]"
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

# ── Skip if already done (and non-empty), unless overwriting ──────────────────
# Existence alone isn't a reliable "done" signal -- an interrupted job (or a
# past bug) can leave a 0-atom freq_iso file behind, which would otherwise be
# skipped forever. Require at least one ATOM/HETATM record too.
if [[ "$OVERWRITE" == "true" ]]; then
    echo "OVERWRITE: removing existing mdpocket exploration outputs"
    rm -f "${RUN_DIR}"/mdpocket_${SEQ_ID}_*
elif [[ -f "$FREQ_ISO" ]] && grep -qE '^(ATOM|HETATM)' "$FREQ_ISO"; then
    echo "SKIP: freq_iso_0_5.pdb already exists and has pocket atoms"
    exit 0
elif [[ -f "$FREQ_ISO" ]]; then
    echo "RERUN: freq_iso_0_5.pdb exists but has no ATOM/HETATM records -- regenerating"
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
    if grep -qE '^(ATOM|HETATM)' "$FREQ_ISO"; then
        echo "OK: $FREQ_ISO"
    else
        echo "OK (mdpocket exit 0) but $FREQ_ISO has no ATOM/HETATM records --"
        echo "  no alpha sphere was present in >=50% of frames for this sequence."
        echo "  This can be a genuine result (unstable/collapsed pocket) rather than a bug."
    fi
    ls -lh mdpocket_${SEQ_ID}_*.pdb mdpocket_${SEQ_ID}_*.dx 2>/dev/null
else
    echo ""
    echo "FAILED: mdpocket exited with code $EXIT_CODE"
    exit 1
fi
