#!/usr/bin/env bash
# =============================================================================
# water_spatial_prep.sh
#
# Prep stage for the gmx-spatial water-density pipeline: for each sequence in
# SEQ_LIST, produce a PBC-corrected, rotation/translation-fit, solvent-
# retained, production-windowed trajectory (fit_trim.xtc) ready for
# `gmx spatial`.
#
# This is deliberately self-contained rather than depending on
# post_processing_pipeline_worker.sh's Phase 1 output (prod_md_40_500ns.xtc):
# that file, and the index.ndx it depends on, are scratch-only pipeline
# artifacts (90-day TTL, not guaranteed to exist yet for every sequence).
# Reading raw inputs from the PetaLibrary archive instead — same as
# water_analysis/R_score_calc.py and pkt_vol/pkt_vol_prep.sh already do —
# avoids that missing/incomplete-file risk entirely.
#
# The fit/centering group ("Protein_LIG": Protein + LIG) and output group
# ("System": everything, including water/ions) mirror the *logic* of
# post_processing_pipeline_worker.sh's Steps 2/3/5, just built from scratch
# from the archived TPR instead of reusing that script's shared index.ndx.
#
# VERIFIED against a real local file set (bind_019_binder, GROMACS 2025.3):
# the index group MUST be built with `gmx make_ndx` (which preserves the
# auto-generated default groups: System, Protein, LIG, SOL, ...), not
# `gmx select -on` (which writes an index file containing ONLY the newly
# selected group). Answering "System" as the output group against a
# select-only index silently mismatches -- gmx falls back to reusing
# whatever single group IS in the file, so trjconv -pbc mol -center
# silently produced a protein+ligand-only trajectory with water already
# stripped, defeating the entire point of this pipeline, with no error.
# This matches why post_processing_pipeline_worker.sh's own Step 1 uses
# make_ndx (heredoc) rather than gmx select for exactly this group.
#
# Usage:
#   bash water_spatial_prep.sh [seq_ids_ngs_observed.txt] [--overwrite-existing]
#                               [--start-ns 40] [--end-ns 500] [--suffix _qfix]
#
# Input  (per sequence, archive only): <ARCHIVE_BASE>/<dir_type>/<seq_id>/<RUNREL>/prod_md_500ns.{tpr,xtc}
# Output (per sequence, scratch):      <BASE>/<dir_type>/<seq_id>/<RUNREL>/water_spatial/{protein_lig.ndx,pbc.xtc,fit_trim.xtc}
# =============================================================================

set -uo pipefail

# ── Configurable paths ────────────────────────────────────────────────────────
ARCHIVE_BASE="/pl/active/shirts_archive/IvanaTang/biosensors"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
REF_TPR="prod_md_500ns.tpr"
REF_XTC="prod_md_500ns.xtc"
OUT_SUBDIR="water_spatial"
FIT_GROUP="protein_lig"
OUTPUT_GROUP="System"
# ─────────────────────────────────────────────────────────────────────────────

OVERWRITE=false
START_NS=40
END_NS=500
SUFFIX=""
SEQ_LIST="/projects/ivta1597/biosensors/seq_ids_ngs_observed.txt"
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --overwrite-existing) OVERWRITE=true; shift ;;
        --start-ns)            START_NS="$2"; shift 2 ;;
        --end-ns)               END_NS="$2"; shift 2 ;;
        --suffix)               SUFFIX="$2"; shift 2 ;;
        *)                     POSITIONAL+=("$1"); shift ;;
    esac
done
[[ ${#POSITIONAL[@]} -gt 0 ]] && SEQ_LIST="${POSITIONAL[0]}"

# "" reproduces the original RUNREL exactly; e.g. "_qfix" points this at the
# bond-order/charge-fix systems' differently-named run directory.
RUNREL="${RUNREL}${SUFFIX}"

START_PS=$(( START_NS * 1000 ))
END_PS=$(( END_NS * 1000 ))

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"       ;;
        "False Positive") echo "nonbinders"    ;;
        "Low Confidence") echo "neg_low_pkt"   ;;
        "Fail Geometry")  echo "neg_fail_gate" ;;
        *)                echo "$1"            ;;
    esac
}

# Files land in either runrel/ or runrel/500ns/ independently of each other,
# same caveat noted throughout this repo -- check each independently.
resolve_input_file() {
    local flat_dir="$1"
    local filename="$2"
    if [[ -f "${flat_dir}/${filename}" ]]; then
        echo "${flat_dir}/${filename}"
    elif [[ -f "${flat_dir}/500ns/${filename}" ]]; then
        echo "${flat_dir}/500ns/${filename}"
    else
        echo ""
    fi
}

should_run() { [[ "$OVERWRITE" == "true" || ! -f "$1" ]]; }

total=0; failed=0; skipped=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue
    ((total++))

    if [[ -n "$custom_path" ]]; then
        IN_FLAT_DIR="${custom_path}/${RUNREL}"
        OUT_DIR="${custom_path}/${RUNREL}/${OUT_SUBDIR}"
    else
        dir_type=$(get_dir_type "$seq_type")
        IN_FLAT_DIR="${ARCHIVE_BASE}/${dir_type}/${seq_id}/${RUNREL}"
        OUT_DIR="${BASE}/${dir_type}/${seq_id}/${RUNREL}/${OUT_SUBDIR}"
    fi

    echo ""
    echo "========================================"
    echo " $seq_id  [$seq_type]"
    echo "========================================"

    TPR=$(resolve_input_file "$IN_FLAT_DIR" "$REF_TPR")
    XTC=$(resolve_input_file "$IN_FLAT_DIR" "$REF_XTC")

    if [[ -z "$TPR" || -z "$XTC" ]]; then
        echo "SKIP: raw TPR/XTC not found in $IN_FLAT_DIR or $IN_FLAT_DIR/500ns"
        ((skipped++)); continue
    fi

    mkdir -p "$OUT_DIR"
    NDX="${OUT_DIR}/protein_lig.ndx"
    PBC_XTC="${OUT_DIR}/pbc.xtc"
    FIT_TRIM_XTC="${OUT_DIR}/fit_trim.xtc"

    # ── Step 1: build the Protein_LIG index group (+ all default groups) ───
    echo ""
    echo "── Step 1: build Protein_LIG index group ──────────────────────────"
    if ! should_run "$NDX"; then
        echo "SKIP: $NDX already exists"
    else
        # "1 | 13" (Protein | LIG) mirrors post_processing_pipeline_worker.sh's
        # own group 23 exactly -- verified against a real .tpr that GROMACS's
        # auto-generated default group order is 0 System, 1 Protein, ...,
        # 13 LIG for this system's topology. gmx make_ndx (not gmx select)
        # is required here so System/Protein/LIG/SOL etc. all stay available
        # alongside the new custom group -- see comment above.
        "$GMX" make_ndx -f "$TPR" -o "$NDX" << EOF
1 | 13
name 20 ${FIT_GROUP}
q
EOF
        if [[ $? -ne 0 || ! -f "$NDX" ]]; then
            echo "FAILED: step 1 — $seq_id"; ((failed++)); continue
        fi
        echo "OK: $NDX"
    fi

    # ── Steps 2-3: PBC correction, then rotation/translation fit + trim ────
    # pbc.xtc is a throwaway intermediate (deleted below), so both steps are
    # gated on fit_trim.xtc's existence rather than pbc.xtc's — otherwise
    # pbc.xtc would look "missing" and get needlessly regenerated on every
    # rerun even after fit_trim.xtc is already done.
    echo ""
    echo "── Steps 2-3: PBC correct, fit + trim (${START_NS}-${END_NS} ns) ──"
    if ! should_run "$FIT_TRIM_XTC"; then
        echo "SKIP: $FIT_TRIM_XTC already exists"
    else
        printf '%s\n%s\n' "$FIT_GROUP" "$OUTPUT_GROUP" | \
            "$GMX" trjconv \
                -s "$TPR" -f "$XTC" -n "$NDX" -o "$PBC_XTC" \
                -pbc mol -center

        if [[ ${PIPESTATUS[1]} -ne 0 || ! -f "$PBC_XTC" ]]; then
            echo "FAILED: step 2 (PBC correction) — $seq_id"; ((failed++)); continue
        fi
        echo "OK: $PBC_XTC"

        printf '%s\n%s\n' "$FIT_GROUP" "$OUTPUT_GROUP" | \
            "$GMX" trjconv \
                -s "$TPR" -f "$PBC_XTC" -n "$NDX" -o "$FIT_TRIM_XTC" \
                -fit rot+trans -b "$START_PS" -e "$END_PS"

        if [[ ${PIPESTATUS[1]} -ne 0 || ! -f "$FIT_TRIM_XTC" ]]; then
            echo "FAILED: step 3 (fit + trim) — $seq_id"; ((failed++)); continue
        fi
        echo "OK: $FIT_TRIM_XTC"
    fi

    # Clean up the PBC-only intermediate — fit_trim.xtc is all downstream
    # steps need, and solvated trajectories are large.
    rm -f "$PBC_XTC"

    # ── Step 4: single-frame reference PDB (frame 0 of fit_trim.xtc) ────────
    # Gives extract_water_spatial_feats.py an mdtraj-loadable topology, and
    # lets it compute the ligand's mean position directly from THIS
    # pipeline's own fitted trajectory (rather than pulling medoid_PL.pdb
    # from the separate RMSD pipeline, which -- despite fitting to the same
    # embedded TPR reference structure and a matching Protein+LIG selection,
    # so in practice likely equivalent -- is an independently-computed fit
    # this script has no way to verify stays bit-identical to). Mirrors
    # pkt_vol_prep.sh's own `-dump 0` reference-PDB pattern.
    echo ""
    echo "── Step 4: extract reference PDB (frame 0) ─────────────────────────"
    REF_PDB="${OUT_DIR}/fit_trim_ref.pdb"
    if ! should_run "$REF_PDB"; then
        echo "SKIP: $REF_PDB already exists"
    else
        printf '%s\n' "$OUTPUT_GROUP" | \
            "$GMX" trjconv \
                -s "$TPR" -f "$FIT_TRIM_XTC" -n "$NDX" \
                -o "$REF_PDB" -dump 0

        if [[ ${PIPESTATUS[1]} -ne 0 || ! -f "$REF_PDB" ]]; then
            echo "FAILED: step 4 — $seq_id"; ((failed++)); continue
        fi
        echo "OK: $REF_PDB"
    fi

done < "$SEQ_LIST"

echo ""
echo "============================="
echo " Total    : $total"
echo " Skipped  : $skipped"
echo " Failed   : $failed"
echo " Complete : $((total - failed - skipped))"
echo "============================="
