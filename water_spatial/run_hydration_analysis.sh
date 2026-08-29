#!/bin/bash

#SBATCH --job-name=hydration_calc
#SBATCH --output=output_hydration_%j.out
#SBATCH --error=error_hydration_%j.err
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
# Usage:
#   sbatch run_hydration_analysis.sh <seq_id> <seq_type> [reference_region] [start_ns] [end_ns] [cutoffs] [stride] [suffix]
#
# Arguments:
#   seq_id            - sequence identifier (e.g. pair_3059_binder)
#   seq_type          - directory group (binders | nonbinders | neg_low_pkt | neg_fail_gate)
#   reference_region  - ligand | pocket_residues (default: ligand)
#   start_ns          - start of analysis window in ns (default: 40)
#   end_ns            - end of analysis window in ns   (default: 500)
#   cutoffs           - comma-separated hydration-shell cutoffs in Angstrom (default: 3,4,5,6,8)
#   stride            - frame stride (default: 10)
#   suffix            - run-directory/output suffix, e.g. "_qfix" (default: "")
#
# Example:
#   sbatch run_hydration_analysis.sh pair_3059_binder binders
#   sbatch run_hydration_analysis.sh pair_3059_binder binders pocket_residues
# ============================================================

set -euo pipefail

# Absolute path, not relative -- SLURM copies this script into a per-job
# spool directory before executing it, and the job's CWD is wherever `sbatch`
# was invoked from (not this script's real location), so a relative script
# path silently fails with FileNotFoundError depending on the caller's CWD.
HYDRATION_CALC_SCRIPT="/projects/ivta1597/biosensors/water_spatial/hydration_calc.py"

module purge
module load anaconda
conda activate biosensors

# pandas/scipy (and other compiled-extension) imports need the envs own
# newer libstdc++, not the HPC images old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version GLIBCXX_3.4.29 not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

seq_id=$1
seq_type=$2
reference_region=${3:-ligand}
start_ns=${4:-40}
end_ns=${5:-500}
cutoffs=${6:-3,4,5,6,8}
stride=${7:-10}
suffix=${8:-}

echo "============================================================"
echo "  Hydration-shell water count"
echo "  seq_id   : $seq_id"
echo "  seq_type : $seq_type"
echo "  region   : ${reference_region}"
echo "  window   : ${start_ns}-${end_ns} ns"
echo "  cutoffs  : ${cutoffs} A"
echo "  stride   : ${stride}"
echo "  suffix   : '${suffix}'"
echo "  start    : $(date)"
echo "============================================================"

python "$HYDRATION_CALC_SCRIPT" \
    --seq_id   "$seq_id"   \
    --seq_type "$seq_type" \
    --reference-region "$reference_region" \
    --start-ns "$start_ns" \
    --end-ns   "$end_ns"   \
    --cutoffs  "$cutoffs"  \
    --stride   "$stride"   \
    --suffix   "$suffix"

echo ""
echo "============================================================"
echo "  Done: $seq_id"
echo "  End: $(date)"
echo "============================================================"
