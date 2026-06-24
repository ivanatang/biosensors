#!/bin/bash
#SBATCH --job-name=gate-latch
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc3
#SBATCH --partition=amilan
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# --- Parse arguments ---
WORKDIR="${1}"

# --- Validate directory ---
if [[ ! -d "$WORKDIR" ]]; then
    echo "ERROR: Directory not found: $WORKDIR"
    exit 1
fi

cd "$WORKDIR" || { echo "ERROR: Cannot cd into $WORKDIR"; exit 1; }

# --- Check required files ---
for f in prod_md_500ns.tpr prod_md_500ns_fit.xtc prod_md_40_500ns.xtc index.ndx; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $WORKDIR/$f"
        exit 1
    fi
done

# Load GROMACS if on HPC (comment out if gmx is already in PATH)
# module load gromacs

set -euo pipefail   # exit on error, undefined var, or pipe failure

GRO="prod_md_500ns.gro"
GATE_RES=88    # P88
LATCH_RES=116  # L116

if [[ ! -f "$GRO" ]]; then
    echo "ERROR: GRO file not found: $WORKDIR/$GRO"
    exit 1
fi

CA_GATE=$(grep -E "^\s+${GATE_RES}[A-Z]+" "$GRO" | grep -E "\s+CA\s+" | awk '{print $3}')
CA_LATCH=$(grep -E "^\s+${LATCH_RES}[A-Z]+" "$GRO" | grep -E "\s+CA\s+" | awk '{print $3}')

if [[ -z "$CA_GATE" || -z "$CA_LATCH" ]]; then
    echo "ERROR: Could not extract atom indices for gate (res ${GATE_RES}) or latch (res ${LATCH_RES}) from $GRO"
    exit 1
fi

echo "  Gate  CA atom index (res ${GATE_RES}): $CA_GATE"
echo "  Latch CA atom index (res ${LATCH_RES}): $CA_LATCH"

echo "[ gate_latch_dist ]"  > gate_latch116.ndx
echo "$CA_GATE $CA_LATCH" >> gate_latch116.ndx

gmx distance \
    -f prod_md_500ns_fit.xtc \
    -n gate_latch116.ndx \
    -select '"gate_latch_dist"' \
    -oall gate_latch116_timeseries.xvg

