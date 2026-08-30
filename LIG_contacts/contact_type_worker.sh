#!/bin/bash
# contact_type_worker.sh
# ----------------------
# SLURM worker script for a single sequence.
# SEQ_ID, START_NS, END_NS, LIGAND_REGION are passed in via --export by
# submit_contact_analysis.sh
#
# Usage (via submit_contact_analysis.sh):
#   sbatch --export=SEQ_ID=pair_3059_binder,START_NS=40,END_NS=500 contact_type_worker.sh
#   sbatch --export=SEQ_ID=pair_3059_binder,START_NS=40,END_NS=250 contact_type_worker.sh
#   sbatch --export=SEQ_ID=pair_3059_binder,START_NS=40,END_NS=500,LIGAND_REGION=core contact_type_worker.sh

#SBATCH --job-name=contact_type
#SBATCH --output=/projects/ivta1597/biosensors/LIG_contacts/logs/contact_type_%j.out
#SBATCH --error=/projects/ivta1597/biosensors/LIG_contacts/logs/contact_type_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ─────────────────────────────────────────────
# USER VARIABLES
# ─────────────────────────────────────────────
PYTHON_SCRIPT="/projects/ivta1597/biosensors/LIG_contacts/contact_type_analysis.py"
CONDA_ENV="biosensors"
LOG_DIR="/projects/ivta1597/biosensors/LIG_contacts/logs"

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
LIGAND_REGION="${LIGAND_REGION:-whole}"
SUFFIX="${SUFFIX:-}"

echo "──────────────────────────────────────────"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Seq ID     : ${SEQ_ID}"
echo "Window     : ${START_NS}–${END_NS} ns"
echo "Region     : ${LIGAND_REGION}"
echo "Suffix     : '${SUFFIX}'"
echo "Node       : $(hostname)"
echo "Start time : $(date)"
echo "──────────────────────────────────────────"

module purge
module load anaconda
conda activate "${CONDA_ENV}"

# pandas/scipy (and other compiled-extension) imports need the env's own
# newer libstdc++, not the HPC image's old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version `GLIBCXX_3.4.29' not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/${CONDA_ENV}/lib:$LD_LIBRARY_PATH"

python "${PYTHON_SCRIPT}" "${SEQ_ID}"    \
    --start-ns "${START_NS}"             \
    --end-ns   "${END_NS}"               \
    --ligand-region "${LIGAND_REGION}"   \
    --suffix "${SUFFIX}"

echo "Finished at: $(date)"
