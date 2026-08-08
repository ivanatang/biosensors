#!/bin/bash

#SBATCH --job-name=gl_hbond_gmx
#SBATCH --output=output_glhbond_%j.out
#SBATCH --error=error_glhbond_%j.err
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
# gate_latch_hbond_gmx.sh
# ------------------------------------------------------------
# Companion to gate_latch_water_bridge.py, which finds frames where a
# bridging water is simultaneously within 4 A (heavy-atom DISTANCE only) of
# ligand + gate + latch. That 4 A test doesn't distinguish a real hydrogen
# bond from a water that just happens to be nearby, and doesn't say whether
# a genuine H-bond (if any) lands on the residue's backbone or side chain.
#
# This script answers that with gmx hbond, which uses the actual geometric
# H-bond definition: heavy-atom donor...acceptor distance < R_CUT (default
# 0.35 nm) AND H-donor...acceptor angle < A_CUT (default 30 deg from
# linear) -- see `gmx hbond -h` on this cluster's build to confirm exact
# flag names/defaults before trusting these blindly; -r/-a matched the
# classic gmx hbond interface as of GROMACS <=2022 but were NOT verified
# against gromacs-2025.3 directly (no way to run gmx from this dev
# environment) -- if `gmx hbond -h` shows different flags, fix the
# GMX_HBOND_CMD template below before submitting at scale.
#
# Backbone/side-chain groups are built with `gmx select`, reusing the
# `group "Protein" and resid ...` syntax already confirmed to parse on this
# GROMACS 2025.3 build in compute_Rg_sasa.sh (the bare `protein` selection
# macro fails to parse on this build; `group "Protein"` works around it).
#
# Groups (residue numbers per CLAUDE.md: gate 84-90, latch 114-118):
#   gate_backbone   : gate residues,  atoms named N/CA/C/O
#   gate_sidechain  : gate residues,  every other atom (incl. its own H's --
#                     needed since side-chain N-H/O-H groups, e.g. His
#                     imidazole, Ser/Thr -OH, are real H-bond donors)
#   latch_backbone  : latch residues, atoms named N/CA/C/O
#   latch_sidechain : latch residues, every other atom
#   water_sol       : resname SOL
#
# For each of the 4 protein groups, gmx hbond is run against water_sol
# alone (each pair written to its OWN 2-group .ndx, always group 0=water,
# group 1=protein-part, so the interactive group prompt can always be
# answered "0\n1" -- avoids depending on GROMACS's own default-group
# numbering, which varies run to run and has already caused off-by-one
# problems elsewhere in this repo).
#
# -num output only (time vs. H-bond count meeting the real geometric
# criteria) -- NOT -hbm/-dist/-ang, to keep this a first-pass yes/no signal
# per residue-part per frame. -num .xvg files parse with this repo's
# standard convention (skip lines starting with # or @). Re-run with -hbm
# added later if per-specific-atom-pair detail turns out to be needed.
#
# gmx hbond has no per-frame concept of "which water" the way
# gate_latch_water_bridge.py does -- it just reports how many D-A pairs
# meeting the geometric criteria exist between the two groups each frame,
# pooled over ALL waters. Cross-reference against the triple-bridge frame
# mask (gate_latch_water_bridge.py --dump-frames) with
# parse_gate_latch_hbond.py to ask "in frames where the literal single-
# water triple bridge is active, is a real H-bond also present, and on
# which side of the residue" -- that join is what actually answers the
# backbone-vs-side-chain question; the raw -num output alone does not
# (it would also count ordinary surface hydration unrelated to the
# gate-latch-ligand network).
#
# Usage:
#   sbatch gate_latch_hbond_gmx.sh <seq_id> <seq_type> [start_ns] [end_ns]
#
# Example:
#   sbatch gate_latch_hbond_gmx.sh pair_3085_binder binders 40 500
# ============================================================

set -euo pipefail

module purge
module load gcc
module load openmpi

GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

seq_id=$1
seq_type=$2
start_ns=${3:-40}
end_ns=${4:-500}
start_ps=$(awk -v n="$start_ns" 'BEGIN{printf "%d", n*1000}')
end_ps=$(awk -v n="$end_ns"   'BEGIN{printf "%d", n*1000}')

R_CUT=0.35   # nm  -- heavy-atom donor...acceptor distance cutoff
A_CUT=30     # deg -- H-donor...acceptor angle cutoff (deviation from linear)

rundir="${BASE}/${seq_type}/${seq_id}/${RUNREL}"
xtc="${rundir}/prod_md_500ns.xtc"
gro="${rundir}/prod_md_500ns.gro"
tpr="${rundir}/prod_md_500ns.tpr"

# Prefer a .tpr if one exists in the run directory -- it carries explicit
# bonded topology, so gmx hbond doesn't have to guess H attachment from
# interatomic distances. Falls back to the .gro (same file the mdtraj-based
# scripts in this repo already use as their topology reference) if no .tpr
# is present. CONFIRM which is actually available in this run directory --
# it wasn't checked before writing this script since none of the existing
# Python analysis here needed bonded topology.
if [[ -f "$tpr" ]]; then
    struct="$tpr"
    echo "Using .tpr for topology (explicit bonds): $struct"
else
    struct="$gro"
    echo "No .tpr found at $tpr -- falling back to $struct"
    echo "  (gmx hbond will infer H attachment without explicit bond records;"
    echo "   double check H-bond counts look sane, e.g. against Hbond_threshold.py)"
fi

for f in "$xtc" "$struct"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

echo "============================================================"
echo "  Gate-latch backbone/side-chain H-bond analysis (gmx hbond)"
echo "  seq_id   : $seq_id"
echo "  seq_type : $seq_type"
echo "  window   : ${start_ns}-${end_ns} ns"
echo "  r_cut    : ${R_CUT} nm   a_cut: ${A_CUT} deg"
echo "  start    : $(date)"
echo "============================================================"

cd "$rundir"

BACKBONE_NAMES="(name N or name CA or name C or name O)"

build_group () {
    local out_ndx=$1 group_name=$2 selection=$3
    local raw_ndx="${out_ndx}.raw"
    "$GMX" select -s "$struct" -on "$raw_ndx" -select "$selection"
    # gmx select names the output group after the selection text itself,
    # which is long/messy -- rewrite the header to a clean, predictable name.
    awk -v name="$group_name" '
        /^\[/ { print "[ " name " ]"; next }
        { print }
    ' "$raw_ndx" > "$out_ndx"
    rm -f "$raw_ndx"
}

build_group gate_backbone.ndx   gate_backbone   "group \"Protein\" and resid 84 to 90   and ${BACKBONE_NAMES}"
build_group gate_sidechain.ndx  gate_sidechain  "group \"Protein\" and resid 84 to 90   and not ${BACKBONE_NAMES}"
build_group latch_backbone.ndx  latch_backbone  "group \"Protein\" and resid 114 to 118 and ${BACKBONE_NAMES}"
build_group latch_sidechain.ndx latch_sidechain "group \"Protein\" and resid 114 to 118 and not ${BACKBONE_NAMES}"
build_group water_sol.ndx       water_sol       "resname SOL"

for grp in gate_backbone gate_sidechain latch_backbone latch_sidechain; do
    n_atoms=$(grep -A1 "\[ ${grp} \]" "${grp}.ndx" | tail -n1 | wc -w)
    echo "  ${grp}: ${n_atoms} atoms"
    if [[ "$n_atoms" -eq 0 ]]; then
        echo "  ERROR: ${grp} selection matched 0 atoms -- check resid numbering / selection syntax"
        exit 1
    fi
done

run_hbond () {
    local grp=$1
    local pair_ndx="${grp}_pair.ndx"
    cat water_sol.ndx "${grp}.ndx" > "$pair_ndx"

    echo "── gmx hbond: water_sol vs ${grp} ──"
    printf "0\n1\n" | "$GMX" hbond \
        -s "$struct" -f "$xtc" -n "$pair_ndx" \
        -b "$start_ps" -e "$end_ps" \
        -r "$R_CUT" -a "$A_CUT" \
        -num "hbond_${grp}_${start_ns}_${end_ns}ns_num.xvg" \
        2> "hbond_${grp}_${start_ns}_${end_ns}ns.log"
}

for grp in gate_backbone gate_sidechain latch_backbone latch_sidechain; do
    run_hbond "$grp"
done

echo "Finished at: $(date)"
echo "Outputs written to: $rundir"
echo "  hbond_{gate,latch}_{backbone,sidechain}_${start_ns}_${end_ns}ns_num.xvg"
