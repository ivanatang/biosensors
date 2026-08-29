#!/bin/bash
#SBATCH --job-name=salt_bridge
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

CONFIG="$1"
SEQ_ID="$2"
SEQ_TYPE="$3"
SCRIPT_DIR="$4"
START_NS="${5:-40}"
END_NS="${6:-500}"

export TMPDIR=$SLURM_SCRATCH
export SLURM_EXPORT_ENV=ALL

module purge
module load gcc
module load openmpi
module load anaconda
conda activate biosensors

# pandas/scipy (and other compiled-extension) imports need the envs own
# newer libstdc++, not the HPC images old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version GLIBCXX_3.4.29 not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

python "${SCRIPT_DIR}/salt_bridge_analysis.py" \
       "$CONFIG" "$SEQ_ID" "$SEQ_TYPE" \
       --start-ns "$START_NS" --end-ns "$END_NS"
