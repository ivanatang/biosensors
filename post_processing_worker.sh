#!/bin/bash
#SBATCH --job-name=postproc
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

# Usage: sbatch postproc_worker.sh <WORKDIR> [CONFIG]
#   WORKDIR  full path to the run subdirectory (prod_md_0p9_cutoff_3dt_64x1_16PME_642dd)
#   CONFIG   path to config.yaml (default: same directory as this script)

set -euo pipefail

WORKDIR="${1:?ERROR: WORKDIR argument required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CONFIG="${2:-${SCRIPT_DIR}/config.yaml}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.yaml not found: $CONFIG" >&2
    exit 1
fi

if [[ ! -d "$WORKDIR" ]]; then
    echo "ERROR: WORKDIR not found: $WORKDIR" >&2
    exit 1
fi

# ── HPC environment ─────────────────────────────────────────────────────────────
module purge
module load gcc
module load openmpi
module load anaconda
conda activate biosensors

# ── Parse config.yaml ───────────────────────────────────────────────────────────
eval "$(python3 << 'PYEOF'
import yaml, os
with open(os.environ['CONFIG']) as f:
    d = yaml.safe_load(f)
pp = d['postproc']
tr = d['trajectory']
o  = pp['outputs']
print('GMX='       + repr(pp['gmx']))
print('TPR='       + repr(tr['tpr']))
print('XTC='       + repr(tr['xtc']))
print('HINGE='     + repr(pp['hinge_resids']))
print('B='         + str(pp['equil_start_ps']))
print('E='         + str(pp['equil_end_ps']))
print('NDX='       + repr(o['index']))
print('PBC_XTC='   + repr(o['pbc_xtc']))
print('FIT_XTC='   + repr(o['fit_xtc']))
print('RMSD='      + repr(o['rmsd_ca']))
print('TRIM_XTC='  + repr(o['trimmed_xtc']))
PYEOF
)"

cd "$WORKDIR"
echo "── Post-processing: $WORKDIR ──"

# ── Step 1: make_ndx ────────────────────────────────────────────────────────────
if [[ ! -f "$NDX" ]]; then
    echo "[1/5] make_ndx → $NDX"
    "$GMX" make_ndx -f "$TPR" -o "$NDX" << EOF
r LIG & ! a H*
name 20 LIG_heavy
1 & ! a H*
name 21 Protein_heavy
8 & ! a H*
name 22 Protein_SC_heavy
1 | 13
name 23 Protein_LIG
ri ${HINGE} & a CA
name 24 CA_hinge
q
EOF
else
    echo "[1/5] SKIP make_ndx ($NDX exists)"
fi

# ── Step 2: trjconv — PBC correction (center on Protein_LIG, output System) ─────
if [[ ! -f "$PBC_XTC" ]]; then
    echo "[2/5] trjconv PBC → $PBC_XTC"
    echo "23 0" | "$GMX" trjconv \
        -s "$TPR" -f "$XTC" -n "$NDX" -o "$PBC_XTC" \
        -pbc mol -center
else
    echo "[2/5] SKIP PBC trjconv ($PBC_XTC exists)"
fi

# ── Step 3: trjconv — rot+trans fit (fit to Protein_LIG, output System) ─────────
if [[ ! -f "$FIT_XTC" ]]; then
    echo "[3/5] trjconv fit → $FIT_XTC"
    echo "23 0" | "$GMX" trjconv \
        -s "$TPR" -f "$PBC_XTC" -n "$NDX" -o "$FIT_XTC" \
        -fit rot+trans
else
    echo "[3/5] SKIP fit trjconv ($FIT_XTC exists)"
fi

# ── Step 4: RMSD on C-alpha ──────────────────────────────────────────────────────
if [[ ! -f "$RMSD" ]]; then
    echo "[4/5] rms (C-alpha) → $RMSD"
    echo "3 3" | "$GMX" rms \
        -s "$TPR" -f "$FIT_XTC" -o "$RMSD"
else
    echo "[4/5] SKIP rms ($RMSD exists)"
fi

# ── Step 5: trjconv — trim to equilibrated window ────────────────────────────────
if [[ ! -f "$TRIM_XTC" ]]; then
    echo "[5/5] trjconv trim (${B} – ${E} ps) → $TRIM_XTC"
    echo "0" | "$GMX" trjconv \
        -s "$TPR" -f "$FIT_XTC" -o "$TRIM_XTC" \
        -b "$B" -e "$E"
else
    echo "[5/5] SKIP trim ($TRIM_XTC exists)"
fi

echo "── Done: $(basename "$WORKDIR") ──"
