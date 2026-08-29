#!/bin/bash

#SBATCH --job-name=rmsd_to_ref_extract
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# Single job, not an array: extract_gate_latch_rmsd_feats.py already loops
# over every sequence in one seq list in a single Python process.
#
# Usage:
#   sbatch submit_extract_rmsd_to_ref.sh [seq_list] [extra extract_gate_latch_rmsd_feats.py args...]
#
# Examples:
#   sbatch submit_extract_rmsd_to_ref.sh                              # defaults to seq_ids_orig.txt
#   sbatch submit_extract_rmsd_to_ref.sh /projects/ivta1597/biosensors/seq_ids_orig.txt --tag _500ns
#
# --time above is a rough estimate for 194 sequences; adjust after the
# first run if it under/over-shoots.
#
# NOTE: sbatch copies this script into a per-job spool directory
# (/var/spool/slurmd/...) before running it, so ${BASH_SOURCE[0]} does NOT
# point at this file's real location -- cd to the hardcoded repo path below
# instead of trying to derive it from the script's own path.
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

SEQ_LIST="${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}"
if [[ $# -gt 0 ]]; then
    shift
fi

# Resolve a relative SEQ_LIST against the directory `sbatch` was invoked
# from (SLURM sets SLURM_SUBMIT_DIR) *before* cd'ing into rmsd_to_ref/ below
# -- otherwise a relative path silently resolves against the wrong
# directory once we cd, since it's not re-resolved at open() time.
if [[ "$SEQ_LIST" != /* ]]; then
    SEQ_LIST="${SLURM_SUBMIT_DIR}/${SEQ_LIST}"
fi

cd /projects/ivta1597/biosensors/rmsd_to_ref

echo "Running extract_gate_latch_rmsd_feats.py against $SEQ_LIST"
python extract_gate_latch_rmsd_feats.py "$SEQ_LIST" "$@"
