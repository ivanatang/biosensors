#!/bin/bash

#SBATCH --job-name=water_contact
#SBATCH --output=output_water_%j.out
#SBATCH --error=error_water_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# Usage:
#   sbatch run_water_analysis.sh <seq_id> <seq_type> [start_ns] [end_ns]
#
# Arguments:
#   seq_id    - sequence identifier (e.g. pair_3059_binder)
#   seq_type  - directory group (binders | nonbinders | neg_low_pkt | neg_fail_gate)
#   start_ns  - start of analysis window in ns (default: 40)
#   end_ns    - end of analysis window in ns   (default: 500)
#
# Examples:
#   sbatch run_water_analysis.sh pair_3059_binder binders           # full 500 ns
#   sbatch run_water_analysis.sh pair_3059_binder binders 40 250    # 250 ns window
#   sbatch run_water_analysis.sh pair_3059_binder binders 40 300    # 300 ns window
# ============================================================

set -euo pipefail

module purge
module load anaconda
conda activate biosensors

seq_id=$1
seq_type=$2
start_ns=${3:-40}
end_ns=${4:-500}

echo "============================================================"
echo "  Water contact analysis"
echo "  seq_id   : $seq_id"
echo "  seq_type : $seq_type"
echo "  window   : ${start_ns}–${end_ns} ns"
echo "  start    : $(date)"
echo "============================================================"

echo ""
echo "=== Step 1: R_score_calc.py ==="
python R_score_calc.py --seq_id   $seq_id   --seq_type $seq_type --start-ns $start_ns --end-ns   $end_ns
if [ $? -ne 0 ]; then
    echo "ERROR: R_score_calc.py failed for $seq_id"
    exit 1
fi

echo ""
echo "=== Step 2: Hbond_threshold.py ==="
python Hbond_threshold.py \
    --seq_id   $seq_id   \
    --seq_type $seq_type \
    --start-ns $start_ns \
    --end-ns   $end_ns
if [ $? -ne 0 ]; then
    echo "ERROR: Hbond_threshold.py failed for $seq_id"
    exit 1
fi

echo ""
echo "=== Step 3: water_hbond_stability.py ==="
python water_hbond_stability.py \
    --seq_id   $seq_id   \
    --seq_type $seq_type \
    --start-ns $start_ns \
    --end-ns   $end_ns
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
