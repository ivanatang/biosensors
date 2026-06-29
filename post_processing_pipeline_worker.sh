#!/bin/bash
#SBATCH --job-name=pp_pipeline
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
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

# Usage: sbatch post_processing_pipeline_worker.sh <WORKDIR> [CONFIG] [END_TIME_PS]
#   WORKDIR     full path to the run subdirectory (prod_md_0p9_cutoff_3dt_64x1_16PME_642dd)
#   CONFIG      path to config.yaml (default: same directory as this script)
#   END_TIME_PS end of the analysis window in ps (default: equil_end_ps from config, 500000)
#               e.g. 250000 to analyse only the first 250 ns of production
#
# Phase 1 (Steps  1-5) always processes the full trajectory for QC/preprocessing.
# Phases 2-3 honour END_TIME_PS: the medoid is found within that window and all
# RMSD/RMSF/distance outputs are computed up to that time point.
# If you re-run with a different END_TIME_PS, delete the existing Phase 3 outputs
# first so the skip-if-exists checks do not reuse stale files.
#
# Runs the full post-processing pipeline for one sequence in three phases:
#   Phase 1 (Steps  1-5):  PBC correction, fitting, C-alpha RMSD, trajectory trim
#   Phase 2 (Step   6):    Find structural medoid of the production ensemble
#   Phase 3 (Steps 7-20):  Medoid-referenced RMSD, RMSF, gate-latch distance
#
# NOTE: Step 19 (RMSF for gate/latch/Lb7a5/recoil loop subgroups) requires
#       groups 30-33 to be present in index.ndx. Run run_gate_latch.sh first
#       to define these groups; the step is skipped with a warning if they
#       are absent.

set -euo pipefail

WORKDIR="${1:?ERROR: WORKDIR argument required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CONFIG="${2:-${SCRIPT_DIR}/config.yaml}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.yaml not found: $CONFIG" >&2; exit 1
fi
if [[ ! -d "$WORKDIR" ]]; then
    echo "ERROR: WORKDIR not found: $WORKDIR" >&2; exit 1
fi

FIND_MEDOID="${SCRIPT_DIR}/get_medoid.py"
if [[ ! -f "$FIND_MEDOID" ]]; then
    echo "ERROR: get_medoid.py not found: $FIND_MEDOID" >&2; exit 1
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
pp  = d['postproc']
tr  = d['trajectory']
med = d['medoid']
o   = pp['outputs']
print('GMX='        + repr(pp['gmx']))
print('TPR='        + repr(tr['tpr']))
print('XTC='        + repr(tr['xtc']))
print('HINGE='      + repr(pp['hinge_resids']))
print('GATE_RES='   + str(pp['gate_res']))
print('LATCH_RES='  + str(pp['latch_res']))
print('B='          + str(pp['equil_start_ps']))
print('E='          + str(pp['equil_end_ps']))
print('NDX='        + repr(o['index']))
print('PBC_XTC='    + repr(o['pbc_xtc']))
print('FIT_XTC='    + repr(o['fit_xtc']))
print('RMSD='       + repr(o['rmsd_ca']))
print('TRIM_XTC='   + repr(o['trimmed_xtc']))
print('MED_STRIDE=' + str(med['stride']))
print('MEDOID_TXT=' + repr(med['medoid_txt']))
PYEOF
)"

END_TIME_PS="${3:-$E}"
END_NS=$(( END_TIME_PS / 1000 ))

cd "$WORKDIR"
echo "════════════════════════════════════════════════════════════════"
echo "  Post-processing pipeline: $WORKDIR"
echo "  Analysis window: $(( B / 1000 ))–${END_NS} ns"
echo "  Start time: $(date)"
echo "════════════════════════════════════════════════════════════════"

# ──────────────────────────────────────────────────────────────────────────────
# PHASE 1: PBC correction, rotational/translational fitting, and trajectory trim
# ──────────────────────────────────────────────────────────────────────────────

# Step 1: Generate atom index groups
# Create named atom groups needed throughout the pipeline: ligand and protein
# heavy atoms, the combined Protein+LIG group used for centering/fitting, and
# the gate-loop CA selection used for hinge RMSF calculations.
if [[ ! -f "$NDX" ]]; then
    echo "[1/20] make_ndx → $NDX"
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
    echo "[1/20] SKIP make_ndx ($NDX exists)"
fi

# Step 2: PBC correction
# Re-image molecules into the primary unit cell and center the protein+ligand
# complex, removing periodic boundary condition artifacts before fitting.
if [[ ! -f "$PBC_XTC" ]]; then
    echo "[2/20] trjconv PBC → $PBC_XTC"
    echo "23 0" | "$GMX" trjconv \
        -s "$TPR" -f "$XTC" -n "$NDX" -o "$PBC_XTC" \
        -pbc mol -center
else
    echo "[2/20] SKIP PBC trjconv ($PBC_XTC exists)"
fi

# Step 3: Rotational and translational fitting
# Align each frame to the initial structure by minimizing the RMSD of
# protein+ligand, removing overall tumbling for cleaner structural analysis.
if [[ ! -f "$FIT_XTC" ]]; then
    echo "[3/20] trjconv fit → $FIT_XTC"
    echo "23 0" | "$GMX" trjconv \
        -s "$TPR" -f "$PBC_XTC" -n "$NDX" -o "$FIT_XTC" \
        -fit rot+trans
else
    echo "[3/20] SKIP fit trjconv ($FIT_XTC exists)"
fi

# Step 4: C-alpha RMSD vs initial frame
# Measure backbone drift over the full trajectory to verify equilibration
# and identify the production-phase start; output used as a QC plot.
if [[ ! -f "$RMSD" ]]; then
    echo "[4/20] rms (C-alpha) → $RMSD"
    echo "3 3" | "$GMX" rms \
        -s "$TPR" -f "$FIT_XTC" -o "$RMSD"
else
    echo "[4/20] SKIP rms ($RMSD exists)"
fi

# Step 5: Trim trajectory to equilibrated production window
# Discard the first ${B} ps (equilibration warm-up) and retain only the
# production-phase frames used in all downstream structural analyses.
if [[ ! -f "$TRIM_XTC" ]]; then
    echo "[5/20] trjconv trim (${B}–${E} ps) → $TRIM_XTC"
    echo "0" | "$GMX" trjconv \
        -s "$TPR" -f "$FIT_XTC" -o "$TRIM_XTC" \
        -b "$B" -e "$E"
else
    echo "[5/20] SKIP trim ($TRIM_XTC exists)"
fi

# ──────────────────────────────────────────────────────────────────────────────
# PHASE 2: Identify the structural medoid of the production ensemble
# ──────────────────────────────────────────────────────────────────────────────

# Step 6: Find the structural medoid frame
# Identify the single trajectory frame that minimizes the sum of squared
# pairwise C-alpha RMSD distances — the most representative structure of the
# ensemble. Used as the reference for all subsequent RMSD and RMSF calculations
# to avoid introducing crystal-structure or first-frame bias.
if [[ ! -f "$MEDOID_TXT" ]]; then
    echo "[6/20] get_medoid.py (stride=${MED_STRIDE}, end=${END_NS} ns) → $MEDOID_TXT"
    python3 "$FIND_MEDOID" \
        -s "$TPR" -f "$TRIM_XTC" \
        --select "protein and name CA" \
        --stride "$MED_STRIDE" \
        --end-ns "$END_NS" \
        --out "$MEDOID_TXT"
else
    echo "[6/20] SKIP get_medoid ($MEDOID_TXT exists)"
fi

DUMP_TIME=$(grep '^medoid_time_ps:' "$MEDOID_TXT" | awk '{print $2}')
if [[ -z "$DUMP_TIME" ]]; then
    echo "ERROR: could not parse medoid_time_ps from $MEDOID_TXT" >&2; exit 1
fi
echo "  Medoid time: ${DUMP_TIME} ps"

# ──────────────────────────────────────────────────────────────────────────────
# PHASE 3: Medoid-referenced structural analysis
# ──────────────────────────────────────────────────────────────────────────────

# Step 7: Extract medoid frame — Protein+Ligand
# Dump only protein and ligand atoms at the medoid time point; used as the
# reference topology for RMSF calculations to avoid solvent-inflated file size.
if [[ ! -f "medoid_PL.pdb" ]]; then
    echo "[7/20] Extract medoid (Protein+LIG, group 23) → medoid_PL.pdb"
    echo 23 | "$GMX" trjconv \
        -s "$TPR" -f "$FIT_XTC" -n "$NDX" \
        -o medoid_PL.pdb -dump "$DUMP_TIME"
else
    echo "[7/20] SKIP medoid_PL.pdb (exists)"
fi

# Step 8: Extract medoid frame — full system
# Dump the complete system (including solvent) at the medoid time point;
# used as the reference topology for distance and RMSD calculations.
if [[ ! -f "medoid_system.pdb" ]]; then
    echo "[8/20] Extract medoid (System, group 0) → medoid_system.pdb"
    echo 0 | "$GMX" trjconv \
        -s "$TPR" -f "$FIT_XTC" -n "$NDX" \
        -o medoid_system.pdb -dump "$DUMP_TIME"
else
    echo "[8/20] SKIP medoid_system.pdb (exists)"
fi

# Step 9: C-alpha RMSD vs medoid
# Quantify per-frame structural deviation from the representative medoid
# structure over the production trajectory; used for conformational clustering
# and to distinguish open/closed states.
if [[ ! -f "rmsd_CA_to_medoid.xvg" ]]; then
    echo "[9/20] rms (CA vs medoid, up to ${END_NS} ns) → rmsd_CA_to_medoid.xvg"
    echo "3 3" | "$GMX" rms \
        -s medoid_system.pdb -f "$TRIM_XTC" -n "$NDX" \
        -o rmsd_CA_to_medoid.xvg -e "$END_TIME_PS"
else
    echo "[9/20] SKIP rmsd_CA_to_medoid.xvg (exists)"
fi

# Step 10: Extract protein+ligand trajectory
# Strip solvent and ions from the trimmed trajectory, keeping only protein and
# ligand atoms to reduce file size and accelerate MDAnalysis-based analyses
# (R-score, water contacts, etc.).
if [[ ! -f "PL_only_40_500ns.xtc" ]]; then
    echo "[10/20] Extract Protein+LIG trajectory (up to ${END_NS} ns) → PL_only_40_500ns.xtc"
    echo 23 | "$GMX" trjconv \
        -s medoid_system.pdb -f "$TRIM_XTC" -n "$NDX" \
        -o PL_only_40_500ns.xtc -e "$END_TIME_PS"
else
    echo "[10/20] SKIP PL_only_40_500ns.xtc (exists)"
fi

# Step 11: Select pocket-proximal C-alpha atoms
# Identify CA atoms within 5 Å of the ligand heavy atoms in the medoid
# structure; these pocket-lining residues capture local flexibility around
# the binding site independent of global backbone motion.
if [[ ! -f "CA_near_LIG_5A.ndx" ]]; then
    echo "[11/20] Select CA within 5 Å of LIG_heavy → CA_near_LIG_5A.ndx"
    "$GMX" select \
        -s medoid_system.pdb -n "$NDX" \
        -select 'name CA and group "Protein" and within 0.5 of group "LIG_heavy"' \
        -on CA_near_LIG_5A.ndx
else
    echo "[11/20] SKIP CA_near_LIG_5A.ndx (exists)"
fi

# Step 12: Add pocket CA group to main index
# gmx select writes the group header as the raw selection string (long and
# unwieldy), so rename it to CA_near_LIG_5A, then patch that block into the
# main index.ndx in-place so group 25 is available to downstream gmx commands.
echo "[12/20] Patch CA_near_LIG_5A into $NDX"
sed -i 's/^\[.*\]/[ CA_near_LIG_5A ]/' CA_near_LIG_5A.ndx

export NDX
python3 - << 'PYEOF'
import re, sys, os
NDX_FILE   = os.environ["NDX"]
PATCH_FILE = "CA_near_LIG_5A.ndx"
TARGET     = "CA_near_LIG_5A"

with open(NDX_FILE) as f:
    content = f.read()
with open(PATCH_FILE) as f:
    new_block = f.read().strip() + "\n"

block_re = re.compile(r'(\[[ \t]*[^\]]+[ \t]*\].*?)(?=\[[ \t]*[^\]]+[ \t]*\]|\Z)', re.DOTALL)
blocks   = block_re.findall(content)

replaced   = False
out_blocks = []
for block in blocks:
    m = re.match(r'\[[ \t]*([^\]]+?)[ \t]*\]', block.strip())
    if m and m.group(1).strip() == TARGET:
        out_blocks.append(new_block)
        replaced = True
    else:
        out_blocks.append(block)

if not replaced:
    sys.stderr.write(f"WARNING: group [ {TARGET} ] not found in {NDX_FILE}; appending.\n")
    out_blocks.append("\n" + new_block)

with open(NDX_FILE, "w") as f:
    f.write("".join(out_blocks))
print(f"  {'Replaced' if replaced else 'Appended'} group [ {TARGET} ] in {NDX_FILE}.")
PYEOF

# Step 13: RMSD — pocket-proximal CA atoms (group 25)
# Track conformational stability of the binding-pocket lining residues
# relative to the medoid, revealing local flexibility independent of global
# backbone drift.
if [[ ! -f "rmsd_CA_near_LIG_5A.xvg" ]]; then
    echo "[13/20] RMSD CA_near_LIG_5A (group 25, up to ${END_NS} ns) → rmsd_CA_near_LIG_5A.xvg"
    echo "3 25" | "$GMX" rms \
        -s medoid_system.pdb -f "$TRIM_XTC" -n "$NDX" \
        -o rmsd_CA_near_LIG_5A.xvg -e "$END_TIME_PS"
else
    echo "[13/20] SKIP rmsd_CA_near_LIG_5A.xvg (exists)"
fi

# Step 14: RMSD — gate-loop hinge CA (group 24)
# Monitor the conformational dynamics of the gate loop (residues 84–90),
# the primary structural element that controls binding-pocket accessibility
# and distinguishes open vs. closed states.
if [[ ! -f "rmsd_CA_hinge.xvg" ]]; then
    echo "[14/20] RMSD CA_hinge (group 24, up to ${END_NS} ns) → rmsd_CA_hinge.xvg"
    echo "3 24" | "$GMX" rms \
        -s medoid_system.pdb -f "$TRIM_XTC" -n "$NDX" \
        -o rmsd_CA_hinge.xvg -e "$END_TIME_PS"
else
    echo "[14/20] SKIP rmsd_CA_hinge.xvg (exists)"
fi

# Step 15: RMSD — ligand heavy atoms (group 20)
# Quantify ligand pose stability within the pocket over the production
# trajectory, discriminating well-anchored bound states from loosely
# associated or partially dissociated conformations.
if [[ ! -f "rmsd_lig_heavy.xvg" ]]; then
    echo "[15/20] RMSD LIG_heavy (group 20, up to ${END_NS} ns) → rmsd_lig_heavy.xvg"
    echo "3 20" | "$GMX" rms \
        -s medoid_system.pdb -f "$TRIM_XTC" -n "$NDX" \
        -o rmsd_lig_heavy.xvg -e "$END_TIME_PS"
else
    echo "[15/20] SKIP rmsd_lig_heavy.xvg (exists)"
fi

# Step 16: Build gate–latch distance index
# Extract the atom numbers of P88 CA (gate loop) and L117 CA (latch) from the
# GRO file and write a two-atom index group; gate–latch separation is the key
# closure metric that separates binders from nonbinders.
GRO="prod_md_500ns.gro"
if [[ ! -f "$GRO" ]]; then
    echo "ERROR: $GRO not found in $WORKDIR" >&2; exit 1
fi
echo "[16/20] Build gate–latch distance index → gate_latch.ndx"
CA_GATE=$(grep -E "^\s+${GATE_RES}[A-Z]+" "$GRO" | grep -E "\s+CA\s+" | awk '{print $3}')
CA_LATCH=$(grep -E "^\s+${LATCH_RES}[A-Z]+" "$GRO" | grep -E "\s+CA\s+" | awk '{print $3}')
if [[ -z "$CA_GATE" || -z "$CA_LATCH" ]]; then
    echo "ERROR: could not find CA atom indices for gate res ${GATE_RES} or latch res ${LATCH_RES} in $GRO" >&2
    exit 1
fi
echo "  Gate  CA (res ${GATE_RES}): atom $CA_GATE"
echo "  Latch CA (res ${LATCH_RES}): atom $CA_LATCH"
printf '[ gate_latch_dist ]\n%s %s\n' "$CA_GATE" "$CA_LATCH" > gate_latch.ndx

# Step 17: Gate–latch distance timeseries
# Compute the time-resolved CA–CA distance between the gate and latch over
# the full fitted trajectory; the primary observable for tracking pocket
# closure dynamics and classifying binder vs. nonbinder behaviour.
if [[ ! -f "gate_latch_timeseries.xvg" ]]; then
    echo "[17/20] Gate–latch distance ($(( B / 1000 ))–${END_NS} ns) → gate_latch_timeseries.xvg"
    "$GMX" distance \
        -f "$FIT_XTC" -n gate_latch.ndx \
        -select '"gate_latch_dist"' \
        -oall gate_latch_timeseries.xvg \
        -b "$B" -e "$END_TIME_PS"
else
    echo "[17/20] SKIP gate_latch_timeseries.xvg (exists)"
fi

# Step 18: RMSF — full Protein+Ligand (group 23)
# Compute per-residue root-mean-square fluctuation referenced to the medoid
# structure, capturing average mobility across the production trajectory for
# every residue.
if [[ ! -f "rmsf_PL.xvg" ]]; then
    echo "[18/20] RMSF Protein+LIG (group 23) → rmsf_PL.xvg"
    echo 23 | "$GMX" rmsf \
        -s medoid_PL.pdb -f PL_only_40_500ns.xtc -n "$NDX" \
        -o rmsf_PL.xvg -res
else
    echo "[18/20] SKIP rmsf_PL.xvg (exists)"
fi

# Step 19: RMSF — loop subgroups (groups 30–33)
# Compute per-residue RMSF separately for gate, latch, Lb7a5, and recoil loop
# subgroups, isolating region-specific flexibility signatures used as ML
# features. Groups 30–33 must be present in index.ndx (created by
# run_gate_latch.sh); this step is skipped with a warning if they are absent.
echo "[19/20] RMSF loop subgroups (groups 30–33)"
for GROUP_NUM in 30 31 32 33; do
    case $GROUP_NUM in
        30) LABEL="ca_gate";   OUTFILE="rmsf_PL_ca_gate.xvg"   ;;
        31) LABEL="ca_latch";  OUTFILE="rmsf_PL_ca_latch.xvg"  ;;
        32) LABEL="ca_Lb7a5";  OUTFILE="rmsf_PL_ca_Lb7a5.xvg"  ;;
        33) LABEL="ca_recoil"; OUTFILE="rmsf_PL_ca_recoil.xvg" ;;
    esac
    if grep -qE "^\[[ \t]*${LABEL}[ \t]*\]" "$NDX" 2>/dev/null; then
        if [[ ! -f "$OUTFILE" ]]; then
            echo "  Group ${GROUP_NUM} (${LABEL}) → ${OUTFILE}"
            echo "$GROUP_NUM" | "$GMX" rmsf \
                -s medoid_PL.pdb -f PL_only_40_500ns.xtc -n "$NDX" \
                -o "$OUTFILE" -res
        else
            echo "  SKIP ${OUTFILE} (exists)"
        fi
    else
        echo "  WARNING: group '${LABEL}' (${GROUP_NUM}) not found in $NDX — skipping. Run run_gate_latch.sh first."
    fi
done

# Step 20: Remove GROMACS backup files
# Clean up the #filename# backup files that GROMACS creates when overwriting
# outputs, preventing unnecessary scratch storage consumption.
echo "[20/20] Removing GROMACS backup files (#*#)"
rm -f \#*\#

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Pipeline complete: $(basename "$(dirname "$WORKDIR")")/$(basename "$WORKDIR")"
echo "  End time: $(date)"
echo "  Output files:"
echo "    $NDX  (groups 20–25 + any loop subgroups)"
echo "    $RMSD"
echo "    $TRIM_XTC"
echo "    $MEDOID_TXT"
echo "    medoid_PL.pdb, medoid_system.pdb"
echo "    rmsd_CA_to_medoid.xvg, rmsd_CA_near_LIG_5A.xvg"
echo "    rmsd_CA_hinge.xvg, rmsd_lig_heavy.xvg"
echo "    PL_only_40_500ns.xtc"
echo "    gate_latch.ndx, gate_latch_timeseries.xvg"
echo "    rmsf_PL.xvg  (+ loop subgroup files if groups 30–33 present)"
echo "════════════════════════════════════════════════════════════════"
