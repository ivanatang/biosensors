#!/usr/bin/env bash
# =============================================================================
# run_pocket_pipeline_manual.sh
#
# Steps 1-2 of pocket volume pipeline for all sequences in SEQ_LIST.
#   1) Strip ligand    → protein_only.xtc
#   2) Extract ref PDB → protein_only.pdb
#
# After this script completes, run mdpocket exploration manually per sequence:
#   mdpocket --trajectory_file protein_only.xtc --trajectory_format xtc \
#            -f protein_only.pdb -o mdpocket_<seq_id>
#
# Usage:
#   bash run_pocket_pipeline_manual.sh [seq_ids.txt]
# =============================================================================

# ── Configurable paths ────────────────────────────────────────────────────────
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
REF_TPR="prod_md_500ns.tpr"
PL_XTC="PL_only_40_500ns.xtc"
PROT_XTC="protein_only.xtc"
PROT_PDB="protein_only.pdb"
PROT_GROUP="Protein"
# ─────────────────────────────────────────────────────────────────────────────

SEQ_LIST=${1:-seq_ids.txt}

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

total=0; failed=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue
    ((total++))

    if [[ -n "$custom_path" ]]; then
        RUN_DIR="${custom_path}/${RUNREL}"
    else
        dir_type=$(get_dir_type "$seq_type")
        RUN_DIR="${BASE}/${dir_type}/${seq_id}/${RUNREL}"
    fi

    TPR="${RUN_DIR}/${REF_TPR}"
    IN_XTC="${RUN_DIR}/${PL_XTC}"
    PROT_XTC_PATH="${RUN_DIR}/${PROT_XTC}"
    PROT_PDB_PATH="${RUN_DIR}/${PROT_PDB}"

    echo ""
    echo "========================================"
    echo " $seq_id  [$seq_type]"
    echo "========================================"

    # ── Validate inputs ───────────────────────────────────────────────────────
    abort=0
    [[ ! -d "$RUN_DIR" ]] && echo "ERROR: run directory not found: $RUN_DIR" && abort=1
    [[ ! -f "$TPR"     ]] && echo "ERROR: TPR not found: $TPR"               && abort=1
    [[ ! -f "$IN_XTC"  ]] && echo "ERROR: PL trajectory not found: $IN_XTC"  && abort=1
    if [[ $abort -eq 1 ]]; then ((failed++)); continue; fi

    # ── Step 1: strip ligand ──────────────────────────────────────────────────
    echo ""
    echo "── Step 1: strip ligand ─────────────────────────────────────────────"
    if [[ -f "$PROT_XTC_PATH" ]]; then
        echo "SKIP: $PROT_XTC already exists"
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
    if [[ -f "$PROT_PDB_PATH" ]]; then
        echo "SKIP: $PROT_PDB already exists"
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
