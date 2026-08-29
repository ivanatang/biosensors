#!/bin/bash

#SBATCH --job-name=extract_water_spatial_feats
#SBATCH --output=output_extract_water_spatial_feats_%j.out
#SBATCH --error=error_extract_water_spatial_feats_%j.err
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
# SLURM wrapper for extract_water_spatial_feats.py -- parsing a cube file
# loads its full (NX,NY,NZ) density array into memory, and once -nab was
# raised to 300 in water_spatial_run.sh (needed to fix a separate
# `gmx spatial` "item outside of allocated memory" error at production
# scale), that padded grid got large enough to OOM a login node. Run this
# on a compute allocation instead.
#
# Usage:
#   sbatch extract_water_spatial_feats_job.sh [seq_list] [--pocket-cutoff 8.0] [--output out.csv]
#
# Example:
#   sbatch extract_water_spatial_feats_job.sh seq_ids_test_one.txt --output water_spatial_feats_test.csv
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

echo "============================================================"
echo "  extract_water_spatial_feats.py"
echo "  args  : $*"
echo "  start : $(date)"
echo "============================================================"

python extract_water_spatial_feats.py "$@"

echo ""
echo "============================================================"
echo "  End: $(date)"
echo "============================================================"
