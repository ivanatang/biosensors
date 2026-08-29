#!/bin/bash
#SBATCH --job-name=pp_pipeline
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
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

# Usage: sbatch post_processing_pipeline_worker.sh <WORKDIR> [CONFIG] [END_TIME_PS] [PHASE] [FORCE]
#   WORKDIR     full path to the run subdirectory (prod_md_0p9_cutoff_3dt_64x1_16PME_642dd)
#   CONFIG      path to config.yaml (default: same directory as this script)
#   END_TIME_PS end of the analysis window in ps (default: equil_end_ps from config, 500000)
#               e.g. 250000 to analyse only the first 250 ns of production
#   PHASE       which phase(s) to run: 1, 2, 3, or all (default: all)
#               Phase 2 and 3 require Phase 1 outputs to already exist.
#               Phase 3 requires Phase 2 outputs to already exist.
#   FORCE       true to overwrite existing outputs; false to skip them (default: false)
#
# Phases 2-3 honour END_TIME_PS: the medoid is found within that window and all
# RMSD/RMSF/distance outputs are written to WORKDIR/<END_NS>ns/.
# Re-running with a different END_TIME_PS produces a separate output directory
# so results from different windows never overwrite each other.
#
#   Phase 1 (Steps  1-5):  PBC correction, fitting, C-alpha RMSD, trajectory trim
#   Phase 2 (Step   6):    Find structural medoid of the production ensemble
#   Phase 3 (Steps 7-21):  Medoid-referenced RMSD, RMSF, gate-latch distance
#
# NOTE: Step 1 creates loop subgroups 25-28 (ca_gate, ca_latch, ca_Lb7a5,
#       ca_recoil) in index.ndx. Step 19 uses these for per-loop RMSF.

set -euo pipefail

WORKDIR="${1:?ERROR: WORKDIR argument required}"
# SLURM copies the script to a spool dir, so BASH_SOURCE[0] is not the repo.
# Use SLURM_SUBMIT_DIR when available (sbatch), fall back to BASH_SOURCE (local).
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export CONFIG="${2:-${SCRIPT_DIR}/config.yaml}"
PHASE="${4:-all}"
FORCE="${5:-false}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.yaml not found: $CONFIG" >&2; exit 1
fi
if [[ ! -d "$WORKDIR" ]]; then
    echo "ERROR: WORKDIR not found: $WORKDIR" >&2; exit 1
fi
if [[ "$PHASE" != "all" && "$PHASE" != "1" && "$PHASE" != "2" && "$PHASE" != "3" ]]; then
    echo "ERROR: PHASE must be 1, 2, 3, or all (got: '$PHASE')" >&2; exit 1
fi
if [[ "$FORCE" != "true" && "$FORCE" != "false" ]]; then
    echo "ERROR: FORCE must be true or false (got: '$FORCE')" >&2; exit 1
fi

# Returns 0 (run the step) if FORCE=true or the output file does not exist yet.
should_run() { [[ "$FORCE" == "true" || ! -f "$1" ]]; }

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

# pandas/scipy (and other compiled-extension) imports need the envs own
# newer libstdc++, not the HPC images old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version GLIBCXX_3.4.29 not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

# ── Parse config.yaml ───────────────────────────────────────────────────────────
eval "$(python3 << 'PYEOF'
import yaml, os
with open(os.environ['CONFIG']) as f:
    d = yaml.safe_load(f)
pp  = d['postproc']
tr  = d['trajectory']
med = d['medoid']
o   = pp['outputs']
print('GMX='            + repr(pp['gmx']))
print('TPR='            + repr(tr['tpr']))
print('XTC='            + repr(tr['xtc']))
print('HINGE='          + repr(pp['hinge_resids']))
print('GATE_RES='       + str(pp['gate_res']))
print('LATCH_RES='      + str(pp['latch_res']))
print('B='              + str(pp['equil_start_ps']))
print('E='              + str(pp['equil_end_ps']))
print('NDX='            + repr(o['index']))
print('PBC_XTC='        + repr(o['pbc_xtc']))
print('FIT_XTC='        + repr(o['fit_xtc']))
print('RMSD='           + repr(o['rmsd_ca']))
print('TRIM_XTC='       + repr(o['trimmed_xtc']))
print('MED_STRIDE='      + str(med['stride']))
print('MEDOID_TXT_FNAME=' + repr(med['medoid_txt']))
print('POCKET_RESIDS='  + repr(med['pocket_resids']))
PYEOF
)"

# ── Derive analysis window and output directory ──────────────────────────────────
END_TIME_PS="${3:-$E}"
END_NS=$(( END_TIME_PS / 1000 ))
START_NS=$(( B / 1000 ))

# All Phase 3 outputs live in a subdirectory named after the analysis end time.
# Phase 1 preprocessing outputs (index, pbc/fit/trim XTC, RMSD) remain in WORKDIR.
OUT_DIR="${WORKDIR}/${END_NS}ns"

# Phase 3 file paths — all inside OUT_DIR
MEDOID_TXT="${OUT_DIR}/${MEDOID_TXT_FNAME}"
MEDOID_PL="${OUT_DIR}/medoid_PL.pdb"
MEDOID_SYS="${OUT_DIR}/medoid_system.pdb"
PL_ONLY="${OUT_DIR}/PL_only_40_${END_NS}ns.xtc"
CA_NDX="${OUT_DIR}/CA_near_LIG_5A.ndx"
OUT_NDX="${OUT_DIR}/index.ndx"   # copy of NDX with CA_near_LIG_5A patched in
GL_NDX="${OUT_DIR}/gate_latch.ndx"

cd "$WORKDIR"
mkdir -p "$OUT_DIR"

echo "════════════════════════════════════════════════════════════════"
echo "  Post-processing pipeline: $WORKDIR"
echo "  Analysis window : ${START_NS}–${END_NS} ns"
echo "  Output directory: $OUT_DIR"
echo "  Start time      : $(date)"
echo "════════════════════════════════════════════════════════════════"

# ──────────────────────────────────────────────────────────────────────────────
# PHASE 1: PBC correction, rotational/translational fitting, and trajectory trim
# These outputs are shared across all analysis windows and stay in WORKDIR.
# ──────────────────────────────────────────────────────────────────────────────
if [[ "$PHASE" == "all" || "$PHASE" == "1" ]]; then

# Step 1: Generate atom index groups
# Create named atom groups needed throughout the pipeline: ligand and protein
# heavy atoms, the combined Protein+LIG group used for centering/fitting,
# the gate-loop CA selection used for hinge RMSF, and loop subgroups 25-28
# (ca_gate 84-90, ca_latch 114-118, ca_Lb7a5 148-155, ca_recoil 154-166).
if should_run "$NDX"; then
    echo "[1/21] make_ndx → $NDX"
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
ri 84-90 & a CA
name 25 ca_gate
ri 114-118 & a CA
name 26 ca_latch
ri 148-155 & a CA
name 27 ca_Lb7a5
ri 154-166 & a CA
name 28 ca_recoil
q
EOF
else
    echo "[1/21] SKIP make_ndx ($NDX exists)"
fi

# Step 2: PBC correction
# Re-image molecules into the primary unit cell and center the protein+ligand
# complex, removing periodic boundary condition artifacts before fitting.
if should_run "$PBC_XTC"; then
    echo "[2/21] trjconv PBC → $PBC_XTC"
    echo "23 0" | "$GMX" trjconv \
        -s "$TPR" -f "$XTC" -n "$NDX" -o "$PBC_XTC" \
        -pbc mol -center
else
    echo "[2/21] SKIP PBC trjconv ($PBC_XTC exists)"
fi

# Step 3: Rotational and translational fitting
# Align each frame to the initial structure by minimizing the RMSD of
# protein+ligand, removing overall tumbling for cleaner structural analysis.
if should_run "$FIT_XTC"; then
    echo "[3/21] trjconv fit → $FIT_XTC"
    echo "23 0" | "$GMX" trjconv \
        -s "$TPR" -f "$PBC_XTC" -n "$NDX" -o "$FIT_XTC" \
        -fit rot+trans
else
    echo "[3/21] SKIP fit trjconv ($FIT_XTC exists)"
fi

# Step 4: C-alpha RMSD vs initial frame
# Measure backbone drift over the full trajectory to verify equilibration
# and identify the production-phase start; output used as a QC plot.
if should_run "$RMSD"; then
    echo "[4/21] rms (C-alpha) → $RMSD"
    echo "3 3" | "$GMX" rms \
        -s "$TPR" -f "$FIT_XTC" -o "$RMSD"
else
    echo "[4/21] SKIP rms ($RMSD exists)"
fi

# Step 5: Trim trajectory to equilibrated production window
# Discard the first ${B} ps (equilibration warm-up) and retain only the
# production-phase frames used in all downstream structural analyses.
if should_run "$TRIM_XTC"; then
    echo "[5/21] trjconv trim (${B}–${E} ps) → $TRIM_XTC"
    echo "0" | "$GMX" trjconv \
        -s "$TPR" -f "$FIT_XTC" -o "$TRIM_XTC" \
        -b "$B" -e "$E"
else
    echo "[5/21] SKIP trim ($TRIM_XTC exists)"
fi

fi # end PHASE 1

# ──────────────────────────────────────────────────────────────────────────────
# PHASE 2: Identify the structural medoid of the production ensemble
# ──────────────────────────────────────────────────────────────────────────────
if [[ "$PHASE" == "all" || "$PHASE" == "2" || "$PHASE" == "3" ]]; then

# Step 6: Find the structural medoid frame
# Identify the single trajectory frame that minimizes the sum of squared
# pairwise C-alpha RMSD distances — the most representative structure of the
# ensemble. Used as the reference for all subsequent RMSD and RMSF calculations
# to avoid introducing crystal-structure or first-frame bias.
if should_run "$MEDOID_TXT"; then
    echo "[6/21] get_medoid.py (stride=${MED_STRIDE}, window=${START_NS}–${END_NS} ns) → $MEDOID_TXT"
    python3 "$FIND_MEDOID" \
        -s "$TPR" -f "$TRIM_XTC" \
        --select "protein and name CA" \
        --stride "$MED_STRIDE" \
        --end-ns "$END_NS" \
        --out "$MEDOID_TXT"
else
    echo "[6/21] SKIP get_medoid ($MEDOID_TXT exists)"
fi

fi # end PHASE 2

# Read medoid time — required for Phase 3 regardless of whether Phase 2 just ran
if [[ "$PHASE" == "all" || "$PHASE" == "2" || "$PHASE" == "3" ]]; then
    if [[ ! -f "$MEDOID_TXT" ]]; then
        echo "ERROR: $MEDOID_TXT not found — run Phase 2 first" >&2; exit 1
    fi
    DUMP_TIME=$(grep '^medoid_time_ps:' "$MEDOID_TXT" | awk '{print $2}')
    if [[ -z "$DUMP_TIME" ]]; then
        echo "ERROR: could not parse medoid_time_ps from $MEDOID_TXT" >&2; exit 1
    fi
    echo "  Medoid time: ${DUMP_TIME} ps"
fi

# ──────────────────────────────────────────────────────────────────────────────
# PHASE 3: Medoid-referenced structural analysis  →  all outputs in $OUT_DIR
# ──────────────────────────────────────────────────────────────────────────────
if [[ "$PHASE" == "all" || "$PHASE" == "3" ]]; then

# Step 7: Extract medoid frame — Protein+Ligand
# Dump only protein and ligand atoms at the medoid time point; used as the
# reference topology for RMSF calculations to avoid solvent-inflated file size.
if should_run "$MEDOID_PL"; then
    echo "[7/21] Extract medoid (Protein+LIG, group 23) → $MEDOID_PL"
    echo 23 | "$GMX" trjconv \
        -s "$TPR" -f "$FIT_XTC" -n "$NDX" \
        -o "$MEDOID_PL" -dump "$DUMP_TIME"
else
    echo "[7/21] SKIP $MEDOID_PL (exists)"
fi

# Step 8: Extract medoid frame — full system
# Dump the complete system (including solvent) at the medoid time point;
# used as the reference topology for distance and RMSD calculations.
if should_run "$MEDOID_SYS"; then
    echo "[8/21] Extract medoid (System, group 0) → $MEDOID_SYS"
    echo 0 | "$GMX" trjconv \
        -s "$TPR" -f "$FIT_XTC" -n "$NDX" \
        -o "$MEDOID_SYS" -dump "$DUMP_TIME"
else
    echo "[8/21] SKIP $MEDOID_SYS (exists)"
fi

# Step 9: C-alpha RMSD vs medoid
# Quantify per-frame structural deviation from the representative medoid
# structure over the production trajectory; used for conformational clustering
# and to distinguish open/closed states.
if should_run "${OUT_DIR}/rmsd_CA_to_medoid.xvg"; then
    echo "[9/21] rms (CA vs medoid, ${START_NS}–${END_NS} ns) → ${OUT_DIR}/rmsd_CA_to_medoid.xvg"
    echo "3 3" | "$GMX" rms \
        -s "$MEDOID_SYS" -f "$TRIM_XTC" -n "$NDX" \
        -o "${OUT_DIR}/rmsd_CA_to_medoid.xvg" -e "$END_TIME_PS"
else
    echo "[9/21] SKIP rmsd_CA_to_medoid.xvg (exists)"
fi

# Step 10: Extract protein+ligand trajectory
# Strip solvent and ions from the trimmed trajectory, keeping only protein and
# ligand atoms to reduce file size and accelerate MDAnalysis-based analyses
# (R-score, water contacts, etc.).
if should_run "$PL_ONLY"; then
    echo "[10/21] Extract Protein+LIG trajectory (${START_NS}–${END_NS} ns) → $PL_ONLY"
    echo 23 | "$GMX" trjconv \
        -s "$MEDOID_SYS" -f "$TRIM_XTC" -n "$NDX" \
        -o "$PL_ONLY" -e "$END_TIME_PS"
else
    echo "[10/21] SKIP $PL_ONLY (exists)"
fi

# Step 11: Select pocket-proximal C-alpha atoms
# Identify CA atoms within 5 Å of the ligand heavy atoms in the medoid
# structure; these pocket-lining residues capture local flexibility around
# the binding site independent of global backbone motion.
if should_run "$CA_NDX"; then
    echo "[11/21] Select CA within 5 Å of LIG_heavy → $CA_NDX"
    "$GMX" select \
        -s "$MEDOID_SYS" -n "$NDX" \
        -select 'name CA and group "Protein" and within 0.5 of group "LIG_heavy"' \
        -on "$CA_NDX"
else
    echo "[11/21] SKIP $CA_NDX (exists)"
fi

# Step 12: Add pocket CA group to the output index
# Copy the shared index.ndx into the output directory so each analysis window
# has its own independent index. Rename the gmx select header to CA_near_LIG_5A
# and patch it into the copy so group 25 is available to downstream gmx commands
# without modifying the original index.ndx in WORKDIR.
echo "[12/21] Copy index.ndx → $OUT_NDX and patch CA_near_LIG_5A"
cp "$NDX" "$OUT_NDX"
sed -i 's/^\[.*\]/[ CA_near_LIG_5A ]/' "$CA_NDX"

export OUT_NDX CA_NDX
python3 - << 'PYEOF'
import re, sys, os
NDX_FILE   = os.environ["OUT_NDX"]
PATCH_FILE = os.environ["CA_NDX"]
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
if should_run "${OUT_DIR}/rmsd_CA_near_LIG_5A.xvg"; then
    echo "[13/21] RMSD CA_near_LIG_5A (group 29, ${START_NS}–${END_NS} ns) → ${OUT_DIR}/rmsd_CA_near_LIG_5A.xvg"
    echo "3 29" | "$GMX" rms \
        -s "$MEDOID_SYS" -f "$TRIM_XTC" -n "$OUT_NDX" \
        -o "${OUT_DIR}/rmsd_CA_near_LIG_5A.xvg" -e "$END_TIME_PS"
else
    echo "[13/21] SKIP rmsd_CA_near_LIG_5A.xvg (exists)"
fi

# Step 14: RMSD — gate-loop hinge CA (group 24)
# Monitor the conformational dynamics of the gate loop (residues 84–90),
# the primary structural element that controls binding-pocket accessibility
# and distinguishes open vs. closed states.
if should_run "${OUT_DIR}/rmsd_CA_hinge.xvg"; then
    echo "[14/21] RMSD CA_hinge (group 24, ${START_NS}–${END_NS} ns) → ${OUT_DIR}/rmsd_CA_hinge.xvg"
    echo "3 24" | "$GMX" rms \
        -s "$MEDOID_SYS" -f "$TRIM_XTC" -n "$OUT_NDX" \
        -o "${OUT_DIR}/rmsd_CA_hinge.xvg" -e "$END_TIME_PS"
else
    echo "[14/21] SKIP rmsd_CA_hinge.xvg (exists)"
fi

# Step 15: RMSD — ligand heavy atoms (group 20)
# Quantify ligand pose stability within the pocket over the production
# trajectory, discriminating well-anchored bound states from loosely
# associated or partially dissociated conformations.
if should_run "${OUT_DIR}/rmsd_lig_heavy.xvg"; then
    echo "[15/21] RMSD LIG_heavy (group 20, ${START_NS}–${END_NS} ns) → ${OUT_DIR}/rmsd_lig_heavy.xvg"
    echo "3 20" | "$GMX" rms \
        -s "$MEDOID_SYS" -f "$TRIM_XTC" -n "$OUT_NDX" \
        -o "${OUT_DIR}/rmsd_lig_heavy.xvg" -e "$END_TIME_PS"
else
    echo "[15/21] SKIP rmsd_lig_heavy.xvg (exists)"
fi

# Step 16: Build gate–latch distance index
# Extract the atom numbers of P88 CA (gate loop) and L117 CA (latch) from the
# GRO file and write a two-atom index group; gate–latch separation is the key
# closure metric that separates binders from nonbinders.
GRO="prod_md_500ns.gro"
if [[ ! -f "$GRO" ]]; then
    echo "ERROR: $GRO not found in $WORKDIR" >&2; exit 1
fi
if should_run "$GL_NDX"; then
    echo "[16/21] Build gate–latch distance index → $GL_NDX"
    CA_GATE=$(grep -E "^\s+${GATE_RES}[A-Z]+" "$GRO" | grep -E "\s+CA\s+" | awk '{print $3}')
    CA_LATCH=$(grep -E "^\s+${LATCH_RES}[A-Z]+" "$GRO" | grep -E "\s+CA\s+" | awk '{print $3}')
    if [[ -z "$CA_GATE" || -z "$CA_LATCH" ]]; then
        echo "ERROR: could not find CA atom indices for gate res ${GATE_RES} or latch res ${LATCH_RES} in $GRO" >&2
        exit 1
    fi
    echo "  Gate  CA (res ${GATE_RES}): atom $CA_GATE"
    echo "  Latch CA (res ${LATCH_RES}): atom $CA_LATCH"
    printf '[ gate_latch_dist ]\n%s %s\n' "$CA_GATE" "$CA_LATCH" > "$GL_NDX"
else
    echo "[16/21] SKIP $GL_NDX (exists)"
fi

# Step 17: Gate–latch distance timeseries
# Compute the time-resolved CA–CA distance between the gate and latch; the
# primary observable for tracking pocket closure dynamics and classifying
# binder vs. nonbinder behaviour.
if should_run "${OUT_DIR}/gate_latch_timeseries.xvg"; then
    echo "[17/21] Gate–latch distance (${START_NS}–${END_NS} ns) → ${OUT_DIR}/gate_latch_timeseries.xvg"
    "$GMX" distance \
        -f "$FIT_XTC" -n "$GL_NDX" \
        -select '"gate_latch_dist"' \
        -oall "${OUT_DIR}/gate_latch_timeseries.xvg" \
        -b "$B" -e "$END_TIME_PS"
else
    echo "[17/21] SKIP gate_latch_timeseries.xvg (exists)"
fi

# Step 18: RMSF — full Protein+Ligand (group 23)
# Compute per-residue root-mean-square fluctuation referenced to the medoid
# structure, capturing average mobility across the production trajectory for
# every residue.
if should_run "${OUT_DIR}/rmsf_PL.xvg"; then
    echo "[18/21] RMSF Protein+LIG (group 23) → ${OUT_DIR}/rmsf_PL.xvg"
    echo 23 | "$GMX" rmsf \
        -s "$MEDOID_PL" -f "$PL_ONLY" -n "$OUT_NDX" \
        -o "${OUT_DIR}/rmsf_PL.xvg" -res
else
    echo "[18/21] SKIP rmsf_PL.xvg (exists)"
fi

# Step 19: RMSF — loop subgroups (groups 30–33)
# Compute per-residue RMSF separately for gate, latch, Lb7a5, and recoil loop
# subgroups, isolating region-specific flexibility signatures used as ML
# features. Groups 30–33 must be present in index.ndx (created by
# run_gate_latch.sh); this step is skipped with a warning if absent.
echo "[19/21] RMSF loop subgroups (groups 30–33)"
for GROUP_NUM in 25 26 27 28; do
    case $GROUP_NUM in
        25) LABEL="ca_gate";   OUTFILE="${OUT_DIR}/rmsf_PL_ca_gate.xvg"   ;;
        26) LABEL="ca_latch";  OUTFILE="${OUT_DIR}/rmsf_PL_ca_latch.xvg"  ;;
        27) LABEL="ca_Lb7a5";  OUTFILE="${OUT_DIR}/rmsf_PL_ca_Lb7a5.xvg"  ;;
        28) LABEL="ca_recoil"; OUTFILE="${OUT_DIR}/rmsf_PL_ca_recoil.xvg" ;;
    esac
    if should_run "$OUTFILE"; then
        echo "  Group ${GROUP_NUM} (${LABEL}) → ${OUTFILE}"
        echo "$GROUP_NUM" | "$GMX" rmsf \
            -s "$MEDOID_PL" -f "$PL_ONLY" -n "$OUT_NDX" \
            -o "$OUTFILE" -res
    else
        echo "  SKIP $OUTFILE (exists)"
    fi
done

# Step 20: Radius of gyration and solvent-accessible surface area
# Compute Rg and SASA for two regions using the per-window PL trajectory and
# medoid reference, mirroring the logic in compute_Rg_sasa.sh:
#   pocket — 27 consensus binding-site residues (Leonard et al. 2024), capturing
#            pocket opening/closing dynamics independent of overall backbone motion
#   whole  — full Protein+Ligand complex, characterising overall compactness
echo "[20/21] Rg and SASA"

POCKET_NDX="${OUT_DIR}/pocket_res.ndx"

# Build pocket residue index from the medoid structure (reuse if already present)
if should_run "$POCKET_NDX"; then
    echo "  Building pocket residue index → $POCKET_NDX"
    "$GMX" select \
        -s "$MEDOID_PL" \
        -on "$POCKET_NDX" \
        -select "resid ${POCKET_RESIDS}"
fi

# Pocket Rg
if should_run "${OUT_DIR}/Rg_pocket.xvg"; then
    echo "  gyrate (pocket) → ${OUT_DIR}/Rg_pocket.xvg"
    echo "0" | "$GMX" gyrate \
        -s "$MEDOID_PL" -f "$PL_ONLY" -n "$POCKET_NDX" \
        -o "${OUT_DIR}/Rg_pocket.xvg"
else
    echo "  SKIP Rg_pocket.xvg (exists)"
fi

# Pocket SASA
if should_run "${OUT_DIR}/sasa_pocket.xvg"; then
    echo "  sasa (pocket) → ${OUT_DIR}/sasa_pocket.xvg"
    echo "0" | "$GMX" sasa \
        -s "$MEDOID_PL" -f "$PL_ONLY" -n "$POCKET_NDX" \
        -o "${OUT_DIR}/sasa_pocket.xvg"
else
    echo "  SKIP sasa_pocket.xvg (exists)"
fi

# Whole Protein+Ligand Rg
if should_run "${OUT_DIR}/Rg_PL.xvg"; then
    echo "  gyrate (Protein+LIG) → ${OUT_DIR}/Rg_PL.xvg"
    echo "Protein_LIG" | "$GMX" gyrate \
        -s "$MEDOID_PL" -f "$PL_ONLY" -n "$OUT_NDX" \
        -o "${OUT_DIR}/Rg_PL.xvg"
else
    echo "  SKIP Rg_PL.xvg (exists)"
fi

# Whole Protein+Ligand SASA
if should_run "${OUT_DIR}/sasa_PL.xvg"; then
    echo "  sasa (Protein+LIG) → ${OUT_DIR}/sasa_PL.xvg"
    echo "Protein_LIG" | "$GMX" sasa \
        -s "$MEDOID_PL" -f "$PL_ONLY" -n "$OUT_NDX" \
        -o "${OUT_DIR}/sasa_PL.xvg"
else
    echo "  SKIP sasa_PL.xvg (exists)"
fi

# Step 21: Remove GROMACS backup files
# Clean up the #filename# backup files that GROMACS creates when overwriting
# outputs, preventing unnecessary scratch storage consumption.
echo "[21/21] Removing GROMACS backup files (#*#)"
rm -f \#*\# "${OUT_DIR}"/\#*\#

fi # end PHASE 3

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Pipeline complete (phase: $PHASE): $(basename "$(dirname "$WORKDIR")")/$(basename "$WORKDIR")"
echo "  End time: $(date)"
if [[ "$PHASE" == "all" || "$PHASE" == "1" ]]; then
echo ""
echo "  Phase 1 outputs (WORKDIR):"
echo "    $NDX  (groups 20–24)"
echo "    $RMSD"
echo "    $TRIM_XTC"
fi
if [[ "$PHASE" == "all" || "$PHASE" == "2" || "$PHASE" == "3" ]]; then
echo ""
echo "  Phase 2-3 outputs ($OUT_DIR/):"
echo "    $MEDOID_TXT_FNAME"
echo "    medoid_PL.pdb, medoid_system.pdb"
echo "    index.ndx  (copy of shared index with CA_near_LIG_5A as group 29)"
echo "    rmsd_CA_to_medoid.xvg, rmsd_CA_near_LIG_5A.xvg"
echo "    rmsd_CA_hinge.xvg, rmsd_lig_heavy.xvg"
echo "    $(basename "$PL_ONLY")"
echo "    gate_latch.ndx, gate_latch_timeseries.xvg"
echo "    rmsf_PL.xvg  (+ loop subgroup files if groups 30–33 present)"
echo "    Rg_pocket.xvg, sasa_pocket.xvg, pocket_res.ndx"
echo "    Rg_PL.xvg, sasa_PL.xvg"
fi
echo "════════════════════════════════════════════════════════════════"
