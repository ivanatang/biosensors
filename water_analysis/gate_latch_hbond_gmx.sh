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
# 0.35 nm, -hbr) AND H-donor...acceptor angle < A_CUT (default 30 deg from
# linear, -hba).
#
# CONFIRMED against gromacs-2025.3's own `gmx hbond -h` (previously guessed
# and got it wrong twice -- see git history): this build's gmx hbond is the
# NEW implementation added in GROMACS 2024 ("If you need the old one, use
# gmx hbond-legacy" per its own help text), rewritten onto the same
# selection-language engine as `gmx select` (the "Too few selections
# provided" error from an earlier version of this script traced back to
# selectioncollection.cpp, the exact source file behind gmx select's own
# parser). There is NO interactive index-number prompt in this version --
# -r/-t take selection EXPRESSIONS directly (reference/target selections,
# NOT a distance cutoff -- -r means something different than it did in the
# pre-2024 tool most tutorials still describe), and the real geometric
# cutoff flags are -hbr/-hba, not -r/-a.
#
# Backbone/side-chain groups are built with `gmx select`, reusing the
# `group "Protein" and resid ...` syntax already confirmed to parse on this
# GROMACS 2025.3 build in compute_Rg_sasa.sh (the bare `protein` selection
# macro fails to parse on this build; `group "Protein"` works around it).
# The same `group "name"` syntax is reused directly in gmx hbond's -r/-t
# expressions below, referencing groups supplied via -n.
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
# For each of the 4 protein groups, gmx hbond -r/-t is run against
# water_sol alone (each pair's groups supplied via a small 2-group -n file,
# just to keep each invocation minimal/uncluttered -- no numbering games
# needed anymore since selections reference groups by NAME).
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
# Read the raw trajectory from PetaLibrary (persistent), NOT
# /scratch/alpine (auto-deletes after 90 days per CLAUDE.md) -- same split
# R_score_calc.py and gate_latch_water_bridge.py already use, for the same
# reason ("scratch auto-deletes ... older runs' xtc/gro are already
# gone"). An earlier version of this script read the trajectory from
# scratch too (copied from compute_Rg_sasa.sh, which only gets away with
# it because it operates on small already-extracted derived files, not the
# full raw trajectory) -- that silently worked for whichever sequences
# hadn't been purged yet and silently failed (missing .xtc) for the rest.
INPUT_BASE="/pl/active/shirts_archive/IvanaTang/biosensors"
OUTPUT_BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

seq_id=$1
seq_type=$2
start_ns=${3:-40}
end_ns=${4:-500}
start_ps=$(awk -v n="$start_ns" 'BEGIN{printf "%d", n*1000}')
end_ps=$(awk -v n="$end_ns"   'BEGIN{printf "%d", n*1000}')

R_CUT=0.35   # nm  -- heavy-atom donor...acceptor distance cutoff
A_CUT=30     # deg -- H-donor...acceptor angle cutoff (deviation from linear)

indir="${INPUT_BASE}/${seq_type}/${seq_id}/${RUNREL}"
outdir="${OUTPUT_BASE}/${seq_type}/${seq_id}/${RUNREL}"
mkdir -p "$outdir"

xtc="${indir}/prod_md_500ns.xtc"
gro="${indir}/prod_md_500ns.gro"
tpr="${indir}/prod_md_500ns.tpr"

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

cd "$outdir"

BACKBONE_NAMES="(name N or name CA or name C or name O)"

build_group () {
    local out_ndx=$1 group_name=$2 selection=$3
    # Must already end in .ndx -- gmx select silently APPENDS .ndx to -on
    # filenames that don't (e.g. "foo.ndx.raw" -> "foo.ndx.raw.ndx"), which
    # left the awk step below looking for a file that was never created.
    local raw_ndx="${out_ndx%.ndx}_raw.ndx"
    # gmx select writes its normal startup banner/progress text to stderr
    # unconditionally, success or failure -- redirect it to its own log
    # (same pattern as the gmx hbond calls below) so it doesn't pollute
    # SLURM's error_*.err and look like a failure when nothing broke.
    "$GMX" select -s "$struct" -on "$raw_ndx" -select "$selection" \
        2> "${group_name}_select.log"
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
    # Sum atom-index tokens across ALL lines belonging to this group, not
    # just the first one -- GROMACS wraps .ndx atom lists at ~15 per line,
    # so a plain "-A1" undercounts any group with more than ~15 atoms.
    n_atoms=$(awk -v name="$grp" '
        $0 == "[ " name " ]" { found=1; next }
        /^\[/ { found=0 }
        found { n += NF }
        END { print n+0 }
    ' "${grp}.ndx")
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

    # GROMACS 2024+ rewrote gmx hbond onto the same selection-language
    # engine as gmx select (confirmed: the "Too few selections provided"
    # error traced back to selectioncollection.cpp, the same source file
    # behind the group "..." syntax already proven to work for gmx select
    # in compute_Rg_sasa.sh). There is no interactive index-number prompt
    # in this version -- -r/-t take selection EXPRESSIONS directly, and -n
    # just supplies the named groups those expressions can reference via
    # group "name". -hbr/-hba are the real distance/angle cutoff flags
    # (NOT -r/-a, which don't mean what the classic pre-2024 gmx hbond's
    # docs/tutorials imply: -r here means "reference selection").
    echo "── gmx hbond: water_sol vs ${grp} ──"
    "$GMX" hbond \
        -s "$struct" -f "$xtc" -n "$pair_ndx" \
        -r 'group "water_sol"' -t "group \"${grp}\"" \
        -b "$start_ps" -e "$end_ps" \
        -hbr "$R_CUT" -hba "$A_CUT" \
        -num "hbond_${grp}_${start_ns}_${end_ns}ns_num.xvg" \
        -o "hbond_${grp}_${start_ns}_${end_ns}ns_pairs.ndx" \
        2> "hbond_${grp}_${start_ns}_${end_ns}ns.log"
}

for grp in gate_backbone gate_sidechain latch_backbone latch_sidechain; do
    run_hbond "$grp"
done

echo "Finished at: $(date)"
echo "Outputs written to: $outdir"
echo "  hbond_{gate,latch}_{backbone,sidechain}_${start_ns}_${end_ns}ns_num.xvg"
