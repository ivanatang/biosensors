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
#   sbatch run_gate_latch_water_bridge.sh <seq_id> <seq_type> [start_ns] [end_ns] [ligand_region] [suffix]
#
# Arguments:
#   seq_id         - sequence identifier (e.g. pair_3059_binder)
#   seq_type       - directory group (binders | nonbinders | neg_low_pkt | neg_fail_gate)
#   start_ns       - start of analysis window in ns (default: 40)
#   end_ns         - end of analysis window in ns   (default: 500)
#   ligand_region  - whole | core | tail (default: core)
#   suffix         - run-directory/output suffix, e.g. "_qfix" (default: "")
#
# Examples:
#   sbatch run_gate_latch_water_bridge.sh pair_3059_binder binders
#   sbatch run_gate_latch_water_bridge.sh pair_3059_binder binders 40 500 core
#   sbatch run_gate_latch_water_bridge.sh bind_022_binder binders 40 500 whole _qfix
# ============================================================

PYTHON_SCRIPT="/projects/ivta1597/biosensors/water_analysis/gate_latch_water_bridge.py"

set -euo pipefail

module purge
module load anaconda
conda activate biosensors

# pandas/scipy (and other compiled-extension) imports need the envs own
# newer libstdc++, not the HPC images old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version GLIBCXX_3.4.29 not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

seq_id=$1
seq_type=$2
start_ns=${3:-0}   # default 0 (not 40) so first_appearance_ns is meaningful
end_ns=${4:-500}
ligand_region=${5:-core}
suffix=${6:-}

echo "============================================================"
echo "  Gate-latch-ligand water bridge analysis"
echo "  seq_id   : $seq_id"
echo "  seq_type : $seq_type"
echo "  window   : ${start_ns}-${end_ns} ns"
echo "  region   : ${ligand_region}"
echo "  suffix   : '${suffix}'"
echo "  start    : $(date)"
echo "============================================================"

python "$PYTHON_SCRIPT" \
    --seq_id "$seq_id" --seq_type "$seq_type" \
    --start-ns "$start_ns" --end-ns "$end_ns" \
    --ligand-region "$ligand_region" --suffix "$suffix"

echo "Finished at: $(date)"
