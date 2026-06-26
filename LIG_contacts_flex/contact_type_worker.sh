#!/bin/bash
# contact_type_worker.sh
# ----------------------
# SLURM worker script for a single sequence.
# SEQ_ID, START_NS, END_NS are passed in via --export by submit_contact_analysis.sh
#
# Usage (via submit_contact_analysis.sh):
#   sbatch --export=SEQ_ID=pair_3059_binder,START_NS=40,END_NS=500 contact_type_worker.sh
#   sbatch --export=SEQ_ID=pair_3059_binder,START_NS=40,END_NS=250 contact_type_worker.sh

#SBATCH --job-name=contact_type
#SBATCH --output=logs/contact_type_%j.out
#SBATCH --error=logs/contact_type_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ─────────────────────────────────────────────
# USER VARIABLES
# ─────────────────────────────────────────────
PYTHON_SCRIPT="/scratch/alpine/ivta1597/LCA_boltz_models/LIG_contacts_flex/contact_type_analysis.py"
CONDA_ENV="biosensors"
LOG_DIR="/scratch/alpine/ivta1597/LCA_boltz_models/LIG_contacts_flex/logs"

# ─────────────────────────────────────────────
set -euo pipefail
mkdir -p "${LOG_DIR}"

# ── Validate required env vars ────────────────────────────────────
if [[ -z "${SEQ_ID}" ]]; then
    echo "ERROR: SEQ_ID is not set. Submit via submit_contact_analysis.sh" >&2
    exit 1
fi

# Apply defaults for optional window args
START_NS="${START_NS:-40}"
END_NS="${END_NS:-500}"

echo "──────────────────────────────────────────"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Seq ID     : ${SEQ_ID}"
echo "Window     : ${START_NS}–${END_NS} ns"
echo "Node       : $(hostname)"
echo "Start time : $(date)"
echo "──────────────────────────────────────────"

module purge
module load anaconda
conda activate "${CONDA_ENV}"

python "${PYTHON_SCRIPT}" "${SEQ_ID}" \
    --start-ns "${START_NS}"          \
    --end-ns   "${END_NS}"

echo "Finished at: $(date)"
