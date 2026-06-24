#!/bin/bash
#SBATCH --job-name=medoid_analysis
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

# ============================================================
# Usage:
#   Bash:  bash run_medoid_analysis.sh <workdir> <dump_time>
#   SLURM: sbatch run_medoid_analysis.sh <workdir> <dump_time>
#
# Arguments:
#   workdir   - path to directory containing simulation files
#   dump_time - time (ps) for the medoid frame (e.g. 40151.175)
#
# Example:
#   bash run_medoid_analysis.sh /scratch/alpine/ivta1597/LCA_boltz_models/seq01/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd 40151.175
# ============================================================

# --- Parse arguments ---
WORKDIR="${1}"
DUMP_TIME="${2}"

if [[ -z "$WORKDIR" || -z "$DUMP_TIME" ]]; then
    echo "ERROR: Missing required arguments."
    echo "Usage: bash $0 <workdir> <dump_time>"
    echo "  workdir   : directory containing prod_md_500ns.tpr, xtc, index.ndx"
    echo "  dump_time : medoid frame time in ps (e.g. 40151.175)"
    exit 1
fi

# --- Validate directory ---
if [[ ! -d "$WORKDIR" ]]; then
    echo "ERROR: Directory not found: $WORKDIR"
    exit 1
fi

echo "============================================================"
echo "  Medoid Analysis Pipeline"
echo "  Working directory : $WORKDIR"
echo "  Dump time (ps)    : $DUMP_TIME"
echo "  Start time        : $(date)"
echo "============================================================"

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

# -------------------------------------------------------
# Step 1: RMSD of CA against fitted trajectory
# -------------------------------------------------------
echo ""
echo "[1/15] Computing CA RMSD (prod_md_500ns_fit.xtc -> prod_md_500ns.tpr)..."
echo "3 3" | gmx rms \
    -s prod_md_500ns.tpr \
    -f prod_md_500ns_fit.xtc \
    -o rmsd_CA.xvg

# -------------------------------------------------------
# Step 2: Extract medoid frame — protein+ligand (group 23)
# -------------------------------------------------------
echo ""
echo "[2/15] Extracting medoid frame (Protein+Ligand) at t=${DUMP_TIME} ps..."
echo 23 | gmx trjconv \
    -s prod_md_500ns.tpr \
    -f prod_md_500ns_fit.xtc \
    -n index.ndx \
    -o medoid_PL.pdb \
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
    -o medoid_system.pdb \
    -dump "${DUMP_TIME}"

# -------------------------------------------------------
# Step 4: RMSD of CA against medoid structure
# -------------------------------------------------------
echo ""
echo "[4/15] Computing CA RMSD (prod_md_40_500ns.xtc -> medoid_system.pdb)..."
echo "3 3" | gmx rms \
    -s medoid_system.pdb \
    -f prod_md_40_500ns.xtc \
    -o rmsd_CA_to_medoid.xvg

# -------------------------------------------------------
# Step 5: Extract protein+ligand trajectory (group 23)
# -------------------------------------------------------
echo ""
echo "[5/15] Extracting Protein+Ligand trajectory from prod_md_40_500ns.xtc..."
echo 23 | gmx trjconv \
    -s medoid_system.pdb \
    -f prod_md_40_500ns.xtc \
    -n index.ndx \
    -o PL_only_40_500ns.xtc

# -------------------------------------------------------
# Step 6: Select CA atoms within 0.5 nm of LIG_heavy
# -------------------------------------------------------
echo ""
echo "[6/15] Selecting CA atoms within 5 Å of LIG_heavy..."
gmx select \
    -s medoid_system.pdb \
    -n index.ndx \
    -select 'name CA and group "Protein" and within 0.5 of group "LIG_heavy"' \
    -on CA_near_LIG_5A.ndx

# -------------------------------------------------------
# Step 7: Rename group header, then replace the existing
#   [ CA_near_LIG_5A ] block in index.ndx in-place.
#   gmx select names the group after the selection string
#   (long/unwieldy), so first fix the header with sed,
#   then use Python to find+replace the matching block.
# -------------------------------------------------------
echo ""
echo "[7/15] Renaming group to 'CA_near_LIG_5A' and replacing block in index.ndx..."
sed -i 's/^\[.*\]/[ CA_near_LIG_5A ]/' CA_near_LIG_5A.ndx

python3 - <<'PYEOF'
import re, sys

NDX_FILE   = "index.ndx"
PATCH_FILE = "CA_near_LIG_5A.ndx"
TARGET     = "CA_near_LIG_5A"

with open(NDX_FILE) as f:
    content = f.read()

with open(PATCH_FILE) as f:
    new_block = f.read().strip() + "\n"

# Split ndx into blocks; each block begins with a [ name ] header
# and runs up to (but not including) the next header or EOF.
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
    -s medoid_system.pdb \
    -f prod_md_40_500ns.xtc \
    -n index.ndx \
    -o rmsd_CA_near_LIG_5A.xvg

# -------------------------------------------------------
# Step 9: RMSD — CA hinge (group 24)
# -------------------------------------------------------
echo ""
echo "[9/15] Computing RMSD for CA hinge (group 24)..."
echo "3 24" | gmx rms \
    -s medoid_system.pdb \
    -f prod_md_40_500ns.xtc \
    -n index.ndx \
    -o rmsd_CA_hinge.xvg

# -------------------------------------------------------
# Step 10: RMSD — ligand heavy atoms (group 20)
# -------------------------------------------------------
echo ""
echo "[10/15] Computing RMSD for ligand heavy atoms (group 20)..."
echo "3 20" | gmx rms \
    -s medoid_system.pdb \
    -f prod_md_40_500ns.xtc \
    -n index.ndx \
    -o rmsd_lig_heavy.xvg

# -------------------------------------------------------
# Step 11: Gate–latch distance index (P88 CA — L117 CA)
# -------------------------------------------------------
echo ""
echo "[11/15] Building gate–latch distance index (P88 CA – L117 CA)..."

GRO="prod_md_500ns.gro"
GATE_RES=88    # P88
LATCH_RES=117  # L117

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

echo "[ gate_latch_dist ]"  > gate_latch.ndx
echo "$CA_GATE $CA_LATCH" >> gate_latch.ndx

# -------------------------------------------------------
# Step 12: Gate–latch distance timeseries
# -------------------------------------------------------
echo ""
echo "[12/15] Computing gate–latch distance timeseries..."
gmx distance \
    -f prod_md_500ns_fit.xtc \
    -n gate_latch.ndx \
    -select '"gate_latch_dist"' \
    -oall gate_latch_timeseries.xvg

# -------------------------------------------------------
# Step 13: RMSF — full protein+ligand (group 23)
# -------------------------------------------------------
echo ""
echo "[13/15] Computing per-residue RMSF (Protein+Ligand, group 23)..."
echo 23 | gmx rmsf \
    -s medoid_PL.pdb \
    -f PL_only_40_500ns.xtc \
    -n index.ndx \
    -o rmsf_PL.xvg \
    -res

# -------------------------------------------------------
# Step 14: RMSF — loop subgroups (groups 30–33)
# -------------------------------------------------------
echo ""
echo "[14/15] Computing per-residue RMSF for loop subgroups..."

for GROUP_NUM in 30 31 32 33; do
    case $GROUP_NUM in
        30) LABEL="ca_gate"   ; OUTFILE="rmsf_PL_ca_gate.xvg"   ;;
        31) LABEL="ca_latch"  ; OUTFILE="rmsf_PL_ca_latch.xvg"  ;;
        32) LABEL="ca_Lb7a5"  ; OUTFILE="rmsf_PL_ca_Lb7a5.xvg"  ;;
        33) LABEL="ca_recoil" ; OUTFILE="rmsf_PL_ca_recoil.xvg" ;;
    esac
    echo "  Group ${GROUP_NUM} (${LABEL}) -> ${OUTFILE}"
    echo $GROUP_NUM | gmx rmsf \
        -s medoid_PL.pdb \
        -f PL_only_40_500ns.xtc \
        -n index.ndx \
        -o "$OUTFILE" \
        -res
done

# -------------------------------------------------------
# Step 15: Clean up GROMACS backup files (#file#)
# -------------------------------------------------------
echo ""
echo "[15/15] Removing GROMACS backup files (#*#)..."
rm -f \#*\#

echo ""
echo "============================================================"
echo "  All steps completed successfully."
echo "  End time: $(date)"
echo "  Output files:"
echo "    rmsd_CA.xvg"
echo "    medoid_PL.pdb"
echo "    medoid_system.pdb"
echo "    rmsd_CA_to_medoid.xvg"
echo "    PL_only_40_500ns.xtc"
echo "    CA_near_LIG_5A.ndx         (standalone)"
echo "    index.ndx                  (CA_near_LIG_5A group replaced in-place)"
echo "    rmsd_CA_near_LIG_5A.xvg"
echo "    rmsd_CA_hinge.xvg"
echo "    rmsd_lig_heavy.xvg"
echo "    gate_latch.ndx"
echo "    gate_latch_timeseries.xvg"
echo "    rmsf_PL.xvg"
echo "    rmsf_PL_ca_gate.xvg"
echo "    rmsf_PL_ca_latch.xvg"
echo "    rmsf_PL_ca_Lb7a5.xvg"
echo "    rmsf_PL_ca_recoil.xvg"
echo "============================================================"
