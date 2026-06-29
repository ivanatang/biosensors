#!/bin/bash
#SBATCH --job-name=medoid_analysis_250ns
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

# ============================================================
# Usage:
#   Bash:  bash run_medoid_analysis_250ns.sh <workdir> <dump_time> [end_time_ps]
#   SLURM: sbatch run_medoid_analysis_250ns.sh <workdir> <dump_time> [end_time_ps]
#
# Arguments:
#   workdir      - path to prod_md_0p9_cutoff_3dt_64x1_16PME_642dd directory
#   dump_time    - medoid frame time in ps from get_medoid.py
#                  (e.g. 187612.500)
#   end_time_ps  - end of analysis window in ps (default: 250000)
#
# All output files are written to <workdir>/half_time/
# No existing files are modified or overwritten.
# ============================================================

# --- Parse arguments ---
WORKDIR="${1}"
DUMP_TIME="${2}"
END_TIME_PS="${3:-250000}"

if [[ -z "$WORKDIR" || -z "$DUMP_TIME" ]]; then
    echo "ERROR: Missing required arguments."
    echo "Usage: bash $0 <workdir> <dump_time> [end_time_ps]"
    exit 1
fi

# Derive tag and output directory
END_NS=$(awk "BEGIN {printf \"%d\", ${END_TIME_PS}/1000}")
TAG="${END_NS}ns"
OUT="half_time"   # subdirectory for all _250ns output files

# Tagged output file paths (all inside $OUT/)
NDX_TAG="${OUT}/index_${TAG}.ndx"
CA_NDX_TAG="${OUT}/CA_near_LIG_5A_${TAG}.ndx"
GL_NDX_TAG="${OUT}/gate_latch_${TAG}.ndx"

# --- Validate directory ---
if [[ ! -d "$WORKDIR" ]]; then
    echo "ERROR: Directory not found: $WORKDIR"
    exit 1
fi

echo "============================================================"
echo "  Medoid Analysis Pipeline (windowed)"
echo "  Working directory : $WORKDIR"
echo "  Output directory  : $WORKDIR/$OUT/"
echo "  Dump time (ps)    : $DUMP_TIME"
echo "  Analysis window   : 40000–${END_TIME_PS} ps  (40–${END_NS} ns)"
echo "  Start time        : $(date)"
echo "============================================================"

cd "$WORKDIR" || { echo "ERROR: Cannot cd into $WORKDIR"; exit 1; }

# --- Check required files ---
for f in prod_md_500ns.tpr prod_md_500ns_fit.xtc prod_md_40_500ns.xtc index.ndx prod_md_500ns.gro; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $WORKDIR/$f"
        exit 1
    fi
done

# Load GROMACS if on HPC (comment out if gmx is already in PATH)
# module load gromacs

set -euo pipefail

# --- Create output directory ---
mkdir -p "${OUT}"
echo "Output directory: $WORKDIR/$OUT/"

# -------------------------------------------------------
# Step 1: RMSD of CA — windowed fitted trajectory
# -------------------------------------------------------
echo ""
echo "[1/15] Computing CA RMSD (prod_md_500ns_fit.xtc, 40–${END_NS} ns)..."
echo "3 3" | gmx rms \
    -s prod_md_500ns.tpr \
    -f prod_md_500ns_fit.xtc \
    -o ${OUT}/rmsd_CA_${TAG}.xvg \
    -b 40000 -e "${END_TIME_PS}"

# -------------------------------------------------------
# Step 2: Extract medoid frame — protein+ligand (group 23)
# -------------------------------------------------------
echo ""
echo "[2/15] Extracting medoid frame (Protein+Ligand) at t=${DUMP_TIME} ps..."
echo 23 | gmx trjconv \
    -s prod_md_500ns.tpr \
    -f prod_md_500ns_fit.xtc \
    -n index.ndx \
    -o ${OUT}/medoid_PL_${TAG}.pdb \
    -dump "${DUMP_TIME}"

# -------------------------------------------------------
# Step 3: Extract medoid frame — full system (group 0)
# -------------------------------------------------------
echo ""
echo "[3/15] Extracting medoid frame (System) at t=${DUMP_TIME} ps..."
echo 0 | gmx trjconv \
    -s prod_md_500ns.tpr \
    -f prod_md_500ns_fit.xtc \
    -n index.ndx \
    -o ${OUT}/medoid_system_${TAG}.pdb \
    -dump "${DUMP_TIME}"

# -------------------------------------------------------
# Step 4: RMSD of CA against windowed medoid
# -------------------------------------------------------
echo ""
echo "[4/15] Computing CA RMSD (prod_md_40_500ns.xtc -> medoid_system_${TAG}.pdb)..."
echo "3 3" | gmx rms \
    -s ${OUT}/medoid_system_${TAG}.pdb \
    -f prod_md_40_500ns.xtc \
    -o ${OUT}/rmsd_CA_to_medoid_${TAG}.xvg \
    -e "${END_TIME_PS}"

# -------------------------------------------------------
# Step 5: Extract windowed protein+ligand trajectory
# -------------------------------------------------------
echo ""
echo "[5/15] Extracting Protein+Ligand trajectory (40–${END_NS} ns)..."
echo 23 | gmx trjconv \
    -s ${OUT}/medoid_system_${TAG}.pdb \
    -f prod_md_40_500ns.xtc \
    -n index.ndx \
    -o ${OUT}/PL_only_40_${TAG}.xtc \
    -e "${END_TIME_PS}"

# -------------------------------------------------------
# Step 6: Select CA atoms within 0.5 nm of LIG_heavy
# -------------------------------------------------------
echo ""
echo "[6/15] Selecting CA atoms within 5 Å of LIG_heavy..."
gmx select \
    -s ${OUT}/medoid_system_${TAG}.pdb \
    -n index.ndx \
    -select 'name CA and group "Protein" and within 0.5 of group "LIG_heavy"' \
    -on "${CA_NDX_TAG}"

# -------------------------------------------------------
# Step 7: Copy index.ndx -> half_time/index_250ns.ndx,
# then patch the CA_near_LIG_5A block in the copy.
# The original index.ndx is never modified.
# -------------------------------------------------------
echo ""
echo "[7/15] Creating ${NDX_TAG} and patching CA_near_LIG_5A block..."
cp index.ndx "${NDX_TAG}"

sed -i 's/^\[.*\]/[ CA_near_LIG_5A ]/' "${CA_NDX_TAG}"

python3 - "${NDX_TAG}" "${CA_NDX_TAG}" <<'PYEOF'
import re, sys

NDX_FILE   = sys.argv[1]
PATCH_FILE = sys.argv[2]
TARGET     = "CA_near_LIG_5A"

with open(NDX_FILE) as f:
    content = f.read()

with open(PATCH_FILE) as f:
    new_block = f.read().strip() + "\n"

block_re = re.compile(r'(\[[ \t]*[^\]]+[ \t]*\].*?)(?=\[[ \t]*[^\]]+[ \t]*\]|\Z)', re.DOTALL)
blocks = block_re.findall(content)

replaced = False
out_blocks = []
for block in blocks:
    m = re.match(r'\[[ \t]*([^\]]+?)[ \t]*\]', block.strip())
    if m and m.group(1).strip() == TARGET:
        out_blocks.append(new_block)
        replaced = True
    else:
        out_blocks.append(block)

if not replaced:
    sys.stderr.write(f"WARNING: group [ {TARGET} ] not found in {NDX_FILE}. Appending instead.\n")
    out_blocks.append("\n" + new_block)

with open(NDX_FILE, "w") as f:
    f.write("".join(out_blocks))

status = "Replaced" if replaced else "Appended"
print(f"  {status} group [ {TARGET} ] in {NDX_FILE}.")
PYEOF

# -------------------------------------------------------
# Step 8: RMSD — CA near LIG (group 25)
# -------------------------------------------------------
echo ""
echo "[8/15] Computing RMSD for CA_near_LIG_5A (group 25)..."
echo "3 25" | gmx rms \
    -s ${OUT}/medoid_system_${TAG}.pdb \
    -f prod_md_40_500ns.xtc \
    -n "${NDX_TAG}" \
    -o ${OUT}/rmsd_CA_near_LIG_5A_${TAG}.xvg \
    -e "${END_TIME_PS}"

# -------------------------------------------------------
# Step 9: RMSD — CA hinge (group 24)
# -------------------------------------------------------
echo ""
echo "[9/15] Computing RMSD for CA hinge (group 24)..."
echo "3 24" | gmx rms \
    -s ${OUT}/medoid_system_${TAG}.pdb \
    -f prod_md_40_500ns.xtc \
    -n "${NDX_TAG}" \
    -o ${OUT}/rmsd_CA_hinge_${TAG}.xvg \
    -e "${END_TIME_PS}"

# -------------------------------------------------------
# Step 10: RMSD — ligand heavy atoms (group 20)
# -------------------------------------------------------
echo ""
echo "[10/15] Computing RMSD for ligand heavy atoms (group 20)..."
echo "3 20" | gmx rms \
    -s ${OUT}/medoid_system_${TAG}.pdb \
    -f prod_md_40_500ns.xtc \
    -n "${NDX_TAG}" \
    -o ${OUT}/rmsd_lig_heavy_${TAG}.xvg \
    -e "${END_TIME_PS}"

# -------------------------------------------------------
# Step 11: Gate–latch distance index (P88 CA — L117 CA)
# -------------------------------------------------------
echo ""
echo "[11/15] Building gate–latch distance index (P88 CA – L117 CA)..."

GRO="prod_md_500ns.gro"
GATE_RES=88
LATCH_RES=117

CA_GATE=$(grep -E "^\s+${GATE_RES}[A-Z]+"  "$GRO" | grep -E "\s+CA\s+" | awk '{print $3}')
CA_LATCH=$(grep -E "^\s+${LATCH_RES}[A-Z]+" "$GRO" | grep -E "\s+CA\s+" | awk '{print $3}')

if [[ -z "$CA_GATE" || -z "$CA_LATCH" ]]; then
    echo "ERROR: Could not extract CA atom indices for gate (res ${GATE_RES}) or latch (res ${LATCH_RES})"
    exit 1
fi

echo "  Gate  CA atom index (res ${GATE_RES}): $CA_GATE"
echo "  Latch CA atom index (res ${LATCH_RES}): $CA_LATCH"

echo "[ gate_latch_dist ]"  > "${GL_NDX_TAG}"
echo "$CA_GATE $CA_LATCH" >> "${GL_NDX_TAG}"

# -------------------------------------------------------
# Step 12: Gate–latch distance timeseries — windowed
# -------------------------------------------------------
echo ""
echo "[12/15] Computing gate–latch distance timeseries (40–${END_NS} ns)..."
gmx distance \
    -f prod_md_500ns_fit.xtc \
    -n "${GL_NDX_TAG}" \
    -select '"gate_latch_dist"' \
    -oall ${OUT}/gate_latch_timeseries_${TAG}.xvg \
    -b 40000 -e "${END_TIME_PS}"

# -------------------------------------------------------
# Step 13: RMSF — full protein+ligand (group 23)
# -------------------------------------------------------
echo ""
echo "[13/15] Computing per-residue RMSF (Protein+Ligand, 40–${END_NS} ns)..."
echo 23 | gmx rmsf \
    -s ${OUT}/medoid_PL_${TAG}.pdb \
    -f ${OUT}/PL_only_40_${TAG}.xtc \
    -n "${NDX_TAG}" \
    -o ${OUT}/rmsf_PL_${TAG}.xvg \
    -res

# -------------------------------------------------------
# Step 14: RMSF — loop subgroups (groups 30–33)
# -------------------------------------------------------
echo ""
echo "[14/15] Computing per-residue RMSF for loop subgroups (40–${END_NS} ns)..."

for GROUP_NUM in 30 31 32 33; do
    case $GROUP_NUM in
        30) LABEL="ca_gate"   ; OUTFILE="${OUT}/rmsf_PL_ca_gate_${TAG}.xvg"   ;;
        31) LABEL="ca_latch"  ; OUTFILE="${OUT}/rmsf_PL_ca_latch_${TAG}.xvg"  ;;
        32) LABEL="ca_Lb7a5"  ; OUTFILE="${OUT}/rmsf_PL_ca_Lb7a5_${TAG}.xvg"  ;;
        33) LABEL="ca_recoil" ; OUTFILE="${OUT}/rmsf_PL_ca_recoil_${TAG}.xvg" ;;
    esac
    echo "  Group ${GROUP_NUM} (${LABEL}) -> ${OUTFILE}"
    echo $GROUP_NUM | gmx rmsf \
        -s ${OUT}/medoid_PL_${TAG}.pdb \
        -f ${OUT}/PL_only_40_${TAG}.xtc \
        -n "${NDX_TAG}" \
        -o "$OUTFILE" \
        -res
done

# -------------------------------------------------------
# Step 15: Clean up GROMACS backup files from half_time/
# -------------------------------------------------------
echo ""
echo "[15/15] Removing GROMACS backup files (#*#)..."
rm -f "${OUT}"/\#*\#

echo ""
echo "============================================================"
echo "  All steps completed successfully."
echo "  End time: $(date)"
echo "  Window: 40–${END_NS} ns  |  Output: $WORKDIR/$OUT/"
echo "  Output files:"
echo "    $OUT/rmsd_CA_${TAG}.xvg"
echo "    $OUT/medoid_PL_${TAG}.pdb"
echo "    $OUT/medoid_system_${TAG}.pdb"
echo "    $OUT/rmsd_CA_to_medoid_${TAG}.xvg"
echo "    $OUT/PL_only_40_${TAG}.xtc"
echo "    $OUT/CA_near_LIG_5A_${TAG}.ndx"
echo "    $OUT/index_${TAG}.ndx"
echo "    $OUT/rmsd_CA_near_LIG_5A_${TAG}.xvg"
echo "    $OUT/rmsd_CA_hinge_${TAG}.xvg"
echo "    $OUT/rmsd_lig_heavy_${TAG}.xvg"
echo "    $OUT/gate_latch_${TAG}.ndx"
echo "    $OUT/gate_latch_timeseries_${TAG}.xvg"
echo "    $OUT/rmsf_PL_${TAG}.xvg"
echo "    $OUT/rmsf_PL_ca_gate_${TAG}.xvg"
echo "    $OUT/rmsf_PL_ca_latch_${TAG}.xvg"
echo "    $OUT/rmsf_PL_ca_Lb7a5_${TAG}.xvg"
echo "    $OUT/rmsf_PL_ca_recoil_${TAG}.xvg"
echo "============================================================"
