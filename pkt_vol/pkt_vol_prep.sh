#!/usr/bin/env bash
# =============================================================================
# run_pocket_pipeline_manual.sh
#
# Steps 1-2 of pocket volume pipeline for all sequences in SEQ_LIST.
#   1) Strip ligand    → protein_only.xtc
#   2) Extract ref PDB → protein_only.pdb
#
# Reads the raw MD topology/trajectory (TPR, PL_only_40_500ns.xtc) from the
# PetaLibrary archive, not scratch -- scratch auto-deletes after 90 days and
# older runs' inputs are already gone. Outputs (protein_only.*) are still
# written to scratch, same as before. Mirrors the same fix already applied to
# water_analysis/R_score_calc.py and LIG_contacts/contact_type_analysis.py.
#
# After this script completes, run mdpocket exploration manually per sequence:
#   mdpocket --trajectory_file protein_only.xtc --trajectory_format xtc \
#            -f protein_only.pdb -o mdpocket_<seq_id>
#
# Usage:
#   bash run_pocket_pipeline_manual.sh [seq_ids_orig.txt] [--overwrite-existing] [--suffix _qfix]
#
# --overwrite-existing regenerates protein_only.xtc/.pdb even if they already
# exist, instead of skipping. Needed after a stale input gets fixed upstream
# (e.g. a PetaLibrary archive resync) -- without it, the "already exists"
# check would keep skipping sequences whose scratch outputs were built from
# the old, truncated/stale input.
# =============================================================================

# ── Configurable paths ────────────────────────────────────────────────────────
INPUT_BASE="/pl/active/shirts_archive/IvanaTang/biosensors"   # raw TPR/trajectory (archive)
OUTPUT_BASE="/scratch/alpine/ivta1597/LCA_boltz_models"       # pipeline outputs (scratch)
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
REF_TPR="prod_md_500ns.tpr"
PL_XTC="PL_only_40_500ns.xtc"
PROT_XTC="protein_only.xtc"
PROT_PDB="protein_only.pdb"
PROT_GROUP="Protein"
# ─────────────────────────────────────────────────────────────────────────────

OVERWRITE=false
SEQ_LIST="/projects/ivta1597/biosensors/seq_ids_orig.txt"
SUFFIX=""
next_is_suffix=false
for arg in "$@"; do
    if [[ "$next_is_suffix" == "true" ]]; then
        SUFFIX="$arg"
        next_is_suffix=false
        continue
    fi
    case "$arg" in
        --overwrite-existing) OVERWRITE=true ;;
        --suffix)              next_is_suffix=true ;;
        *)                    SEQ_LIST="$arg" ;;
    esac
done

# "" reproduces the original RUNREL exactly; e.g. "_qfix" points this at the
# bond-order/charge-fix systems' differently-named run directory.
RUNREL="${RUNREL}${SUFFIX}"

# Returns 0 (run/regenerate) if OVERWRITE=true or the output file doesn't
# exist yet. Mirrors post_processing_pipeline_worker.sh's should_run/FORCE.
should_run() { [[ "$OVERWRITE" == "true" || ! -f "$1" ]]; }

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

# Individual files land in either runrel/ or runrel/500ns/ independently of
# each other -- e.g. many sequences have the .tpr directly in runrel/ even
# though other files for that same sequence are under runrel/500ns/. So look
# each file up on its own rather than resolving one shared directory.
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

total=0; failed=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue
    ((total++))

    if [[ -n "$custom_path" ]]; then
        IN_FLAT_DIR="${custom_path}/${RUNREL}"
        OUT_RUN_DIR="${custom_path}/${RUNREL}"
    else
        dir_type=$(get_dir_type "$seq_type")
        IN_FLAT_DIR="${INPUT_BASE}/${dir_type}/${seq_id}/${RUNREL}"
        OUT_RUN_DIR="${OUTPUT_BASE}/${dir_type}/${seq_id}/${RUNREL}"
    fi

    TPR=$(resolve_input_file "$IN_FLAT_DIR" "$REF_TPR")
    IN_XTC=$(resolve_input_file "$IN_FLAT_DIR" "$PL_XTC")
    PROT_XTC_PATH="${OUT_RUN_DIR}/${PROT_XTC}"
    PROT_PDB_PATH="${OUT_RUN_DIR}/${PROT_PDB}"

    echo ""
    echo "========================================"
    echo " $seq_id  [$seq_type]"
    echo "========================================"

    # ── Validate inputs ───────────────────────────────────────────────────────
    abort=0
    [[ ! -d "$IN_FLAT_DIR" ]] && echo "ERROR: archive input directory not found: $IN_FLAT_DIR" && abort=1
    [[ -z "$TPR"    ]] && echo "ERROR: TPR not found in $IN_FLAT_DIR or $IN_FLAT_DIR/500ns"                && abort=1
    [[ -z "$IN_XTC" ]] && echo "ERROR: PL trajectory not found in $IN_FLAT_DIR or $IN_FLAT_DIR/500ns"      && abort=1
    if [[ $abort -eq 1 ]]; then ((failed++)); continue; fi

    mkdir -p "$OUT_RUN_DIR"

    # ── Step 1: strip ligand ──────────────────────────────────────────────────
    echo ""
    echo "── Step 1: strip ligand ─────────────────────────────────────────────"
    if ! should_run "$PROT_XTC_PATH"; then
        echo "SKIP: $PROT_XTC already exists (use --overwrite-existing to regenerate)"
    else
        printf '%s\n' "$PROT_GROUP" | \
            "$GMX" trjconv \
                -s "$TPR"           \
                -f "$IN_XTC"        \
                -o "$PROT_XTC_PATH"

        if [[ ${PIPESTATUS[1]} -ne 0 || ! -f "$PROT_XTC_PATH" ]]; then
            echo "FAILED: step 1 — $seq_id"; ((failed++)); continue
        fi
        echo "OK: $PROT_XTC_PATH"
    fi

    # ── Step 2: extract reference PDB ────────────────────────────────────────
    echo ""
    echo "── Step 2: extract reference PDB ───────────────────────────────────"
    if ! should_run "$PROT_PDB_PATH"; then
        echo "SKIP: $PROT_PDB already exists (use --overwrite-existing to regenerate)"
    else
        printf '%s\n' "$PROT_GROUP" | \
            "$GMX" trjconv \
                -s    "$TPR"           \
                -f    "$PROT_XTC_PATH" \
                -o    "$PROT_PDB_PATH" \
                -dump 0

        if [[ ${PIPESTATUS[1]} -ne 0 || ! -f "$PROT_PDB_PATH" ]]; then
            echo "FAILED: step 2 — $seq_id"; ((failed++)); continue
        fi
        echo "OK: $PROT_PDB_PATH"
    fi

done < "$SEQ_LIST"

echo ""
echo "============================="
echo " Total    : $total"
echo " Failed   : $failed"
echo " Complete : $((total - failed))"
echo "============================="
echo ""
echo "Next: run mdpocket exploration manually for each sequence:"
echo "  cd <run_dir>"
echo "  mdpocket --trajectory_file protein_only.xtc --trajectory_format xtc \\"
echo "           -f protein_only.pdb -o mdpocket_<seq_id>"
