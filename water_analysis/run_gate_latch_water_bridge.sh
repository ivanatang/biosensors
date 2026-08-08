#!/bin/bash

#SBATCH --job-name=gl_water_bridge
#SBATCH --output=output_glbridge_%j.out
#SBATCH --error=error_glbridge_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# Usage:
#   sbatch run_gate_latch_water_bridge.sh <seq_id> <seq_type> [start_ns] [end_ns] [ligand_region] [dump_frames]
#
# Arguments:
#   seq_id         - sequence identifier (e.g. pair_3059_binder)
#   seq_type       - directory group (binders | nonbinders | neg_low_pkt | neg_fail_gate)
#   start_ns       - start of analysis window in ns (default: 40)
#   end_ns         - end of analysis window in ns   (default: 500)
#   ligand_region  - whole | core | tail (default: core)
#   dump_frames    - "true" to also write the per-frame triple_bridge CSV
#                     needed by parse_gate_latch_hbond.py (default: false)
#
# Examples:
#   sbatch run_gate_latch_water_bridge.sh pair_3059_binder binders
#   sbatch run_gate_latch_water_bridge.sh pair_3059_binder binders 40 500 core
#   sbatch run_gate_latch_water_bridge.sh pair_3059_binder binders 40 500 core true
# ============================================================

PYTHON_SCRIPT="/projects/ivta1597/biosensors/water_analysis/gate_latch_water_bridge.py"

set -euo pipefail

module purge
module load anaconda
conda activate biosensors

seq_id=$1
seq_type=$2
start_ns=${3:-0}   # default 0 (not 40) so first_appearance_ns is meaningful
end_ns=${4:-500}
ligand_region=${5:-core}
dump_frames=${6:-false}

echo "============================================================"
echo "  Gate-latch-ligand water bridge analysis"
echo "  seq_id      : $seq_id"
echo "  seq_type    : $seq_type"
echo "  window      : ${start_ns}-${end_ns} ns"
echo "  region      : ${ligand_region}"
echo "  dump_frames : ${dump_frames}"
echo "  start       : $(date)"
echo "============================================================"

extra_args=()
[[ "$dump_frames" == "true" ]] && extra_args+=(--dump-frames)

python "$PYTHON_SCRIPT" \
    --seq_id "$seq_id" --seq_type "$seq_type" \
    --start-ns "$start_ns" --end-ns "$end_ns" \
    --ligand-region "$ligand_region" \
    "${extra_args[@]}"

echo "Finished at: $(date)"
