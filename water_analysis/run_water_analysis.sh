#!/bin/bash

#SBATCH --job-name=water_contact
#SBATCH --output=output_water_%j.out
#SBATCH --error=error_water_%j.err
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
#   sbatch run_water_analysis.sh <seq_id> <seq_type> [start_ns] [end_ns] [ligand_region] [suffix]
#
# Arguments:
#   seq_id         - sequence identifier (e.g. pair_3059_binder)
#   seq_type       - directory group (binders | nonbinders | neg_low_pkt | neg_fail_gate)
#   start_ns       - start of analysis window in ns (default: 40)
#   end_ns         - end of analysis window in ns   (default: 500)
#   ligand_region  - whole | core | tail (default: whole)
#   suffix         - run-directory/output suffix, e.g. "_qfix" for the
#                    bond-order/charge-fix systems (default: "", standard)
#
# Examples:
#   sbatch run_water_analysis.sh pair_3059_binder binders                    # full 500 ns
#   sbatch run_water_analysis.sh pair_3059_binder binders 40 250             # 250 ns window
#   sbatch run_water_analysis.sh pair_3059_binder binders 40 300             # 300 ns window
#   sbatch run_water_analysis.sh pair_3059_binder binders 40 500 core        # steroid core only
#   sbatch run_water_analysis.sh pair_3059_binder binders 40 500 tail        # carboxylate tail only
#   sbatch run_water_analysis.sh bind_022_binder binders 40 500 whole _qfix  # _qfix system
#
# NOTE: Hbond_threshold.py and water_hbond_stability.py (steps 2-3) always
# select the whole ligand and read/write the unsuffixed water_contacts_{TAG}
# path, so they are only meaningful for ligand_region=whole. They are
# skipped automatically for core/tail runs.
# ============================================================

set -euo pipefail

# Absolute paths, not relative -- SLURM copies this script into a per-job
# spool directory before executing it, and the job's CWD is wherever `sbatch`
# was invoked from (not this script's real location), so relative script
# paths silently fail with FileNotFoundError depending on the caller's CWD.
R_SCORE_SCRIPT="/projects/ivta1597/biosensors/water_analysis/R_score_calc.py"
HBOND_THRESHOLD_SCRIPT="/projects/ivta1597/biosensors/water_analysis/Hbond_threshold.py"
HBOND_STABILITY_SCRIPT="/projects/ivta1597/biosensors/water_analysis/water_hbond_stability.py"

module purge
module load anaconda
conda activate biosensors

# pandas/scipy (and other compiled-extension) imports need the envs own
# newer libstdc++, not the HPC images old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version GLIBCXX_3.4.29 not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

seq_id=$1
seq_type=$2
start_ns=${3:-40}
end_ns=${4:-500}
ligand_region=${5:-whole}
suffix=${6:-}

echo "============================================================"
echo "  Water contact analysis"
echo "  seq_id   : $seq_id"
echo "  seq_type : $seq_type"
echo "  window   : ${start_ns}–${end_ns} ns"
echo "  region   : ${ligand_region}"
echo "  suffix   : '${suffix}'"
echo "  start    : $(date)"
echo "============================================================"

echo ""
echo "=== Step 1: R_score_calc.py ==="
python "$R_SCORE_SCRIPT" --seq_id $seq_id --seq_type $seq_type --start-ns $start_ns --end-ns $end_ns --ligand-region $ligand_region --suffix "$suffix"
if [ $? -ne 0 ]; then
    echo "ERROR: R_score_calc.py failed for $seq_id"
    exit 1
fi

if [ "$ligand_region" != "whole" ]; then
    echo ""
    echo "============================================================"
    echo "  Region-restricted run ($ligand_region): skipping"
    echo "  Hbond_threshold.py / water_hbond_stability.py (whole-ligand only)."
    echo "  End: $(date)"
    echo "============================================================"
    exit 0
fi

echo ""
echo "=== Step 2: Hbond_threshold.py ==="
python "$HBOND_THRESHOLD_SCRIPT" \
    --seq_id   $seq_id   \
    --seq_type $seq_type \
    --start-ns $start_ns \
    --end-ns   $end_ns \
    --suffix   "$suffix"
if [ $? -ne 0 ]; then
    echo "ERROR: Hbond_threshold.py failed for $seq_id"
    exit 1
fi

echo ""
echo "=== Step 3: water_hbond_stability.py ==="
python "$HBOND_STABILITY_SCRIPT" \
    --seq_id   $seq_id   \
    --seq_type $seq_type \
    --start-ns $start_ns \
    --end-ns   $end_ns \
    --suffix   "$suffix"
if [ $? -ne 0 ]; then
    echo "ERROR: water_hbond_stability.py failed for $seq_id"
    exit 1
fi

echo ""
echo "============================================================"
echo "  All steps completed for $seq_id"
echo "  Output directory: water_contacts_${start_ns}_${end_ns}ns/"
echo "  End: $(date)"
echo "============================================================"
