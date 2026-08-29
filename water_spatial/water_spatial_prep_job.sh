#!/bin/bash

#SBATCH --job-name=water_spatial_prep
#SBATCH --output=output_water_spatial_prep_%j.out
#SBATCH --error=error_water_spatial_prep_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# SLURM wrapper for water_spatial_prep.sh -- that script is a plain bash
# loop (no SBATCH header of its own, meant to be run directly per
# pkt_vol_prep.sh's convention), but the PBC-correct/fit/trim gmx trjconv
# calls are heavy enough over the full 40-500ns window to need a real
# compute allocation rather than running on a login node.
#
# Usage:
#   sbatch water_spatial_prep_job.sh [seq_list] [--overwrite-existing] [--start-ns N] [--end-ns N]
#
# Example:
#   sbatch water_spatial_prep_job.sh seq_ids_test_one.txt
#   sbatch water_spatial_prep_job.sh seq_ids_orig.txt   # full 193-sequence cohort
#
# water_spatial_prep.sh loops through the whole seq_list serially inside
# this one job (per-sequence should_run/skip logic means a timed-out job
# can just be resubmitted with the same seq_list to resume where it left
# off -- already-complete sequences are skipped, not redone). 24h walltime
# is a generous estimate: bind_019_binder's full prep+run+extract test
# took ~15 min total, of which prep is only one (not the slowest) stage.
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
echo "  water_spatial_prep.sh"
echo "  args  : $*"
echo "  start : $(date)"
echo "============================================================"

bash water_spatial_prep.sh "$@"

echo ""
echo "============================================================"
echo "  End: $(date)"
echo "============================================================"
