#!/bin/bash
#SBATCH --job-name=salt_bridge
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

CONFIG="$1"
SEQ_ID="$2"
SEQ_TYPE="$3"
SCRIPT_DIR="$4"

export TMPDIR=$SLURM_SCRATCH
export SLURM_EXPORT_ENV=ALL

module purge
module load gcc
module load openmpi
module load anaconda
conda activate biosensors

python "${SCRIPT_DIR}/salt_bridge_analysis.py" \
       "$CONFIG" "$SEQ_ID" "$SEQ_TYPE"
