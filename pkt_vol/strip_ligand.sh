#!/usr/bin/env bash
# =============================================================================
# strip_ligand.sh
#
# For each sequence in SEQ_LIST, extract a protein-only trajectory from the
# existing protein+ligand production XTC using gmx trjconv.
#
# Usage:
#   bash strip_ligand_trjconv.sh [seq_ids.txt]
#
# Input  (per sequence): <BASE>/<dir_type>/<seq_id>/<RUNREL>/PL_only_40_500ns.xtc
# Output (per sequence): <BASE>/<dir_type>/<seq_id>/<RUNREL>/protein_only.xtc
# =============================================================================

# ── Configurable paths ────────────────────────────────────────────────────────
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
INPUT_XTC="PL_only_40_500ns.xtc"
REF_PDB="prod_md_500ns.tpr"
INDEX="index.ndx"
OUTPUT_XTC="protein_only.xtc"
CENTER_GROUP="Protein"
OUTPUT_GROUP="Protein"
# ─────────────────────────────────────────────────────────────────────────────

SEQ_LIST=${1:-seq_ids.txt}

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

# ── Map display seq_type → directory name used in the file system ─────────────
get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"       ;;
        "False Positive") echo "nonbinders"    ;;
        "Low Confidence") echo "neg_low_pkt"   ;;
        "Fail Geometry")  echo "neg_fail_gate" ;;
        *)                echo "$1"            ;;   # fallback: use as-is
    esac
}

processed=0
skipped=0
failed=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    # Skip empty lines and comments
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    # ── Resolve run directory ─────────────────────────────────────────────────
    if [[ -n "$custom_path" ]]; then
        RUN_DIR="${custom_path}/${RUNREL}"
        echo "INFO (custom path): $seq_id → $RUN_DIR"
    else
        dir_type=$(get_dir_type "$seq_type")
        RUN_DIR="${BASE}/${dir_type}/${seq_id}/${RUNREL}"
    fi

    # ── Validate inputs ───────────────────────────────────────────────────────
    if [[ ! -d "$RUN_DIR" ]]; then
        echo "SKIP (dir not found): $seq_id  [$RUN_DIR]"
        ((skipped++))
        continue
    fi

    IN_XTC="${RUN_DIR}/${INPUT_XTC}"
    REF="${RUN_DIR}/${REF_PDB}"
    NDX="${RUN_DIR}/${INDEX}"
    OUT_XTC="${RUN_DIR}/${OUTPUT_XTC}"

    if [[ ! -f "$IN_XTC" ]]; then
        echo "SKIP (no input XTC): $seq_id  [$IN_XTC]"
        ((skipped++))
        continue
    fi

    if [[ ! -f "$REF" ]]; then
        echo "SKIP (no ref PDB): $seq_id  [$REF]"
        ((skipped++))
        continue
    fi

    # ── Skip if output already exists ─────────────────────────────────────────
    if [[ -f "$OUT_XTC" ]]; then
        echo "SKIP (already done): $seq_id  [$OUT_XTC]"
        ((skipped++))
        continue
    fi

    # ── Run gmx trjconv ───────────────────────────────────────────────────────
    echo "Processing: $seq_id  [$seq_type]"

    NDX_FLAG=""
    [[ -f "$NDX" ]] && NDX_FLAG="-n ${NDX}"

    printf '%s\n%s\n' "$CENTER_GROUP" "$OUTPUT_GROUP" | \
        gmx trjconv \
            -s  "$REF"    \
            -f  "$IN_XTC" \
            -o  "$OUT_XTC" \
            $NDX_FLAG      \
            -center        \
            -pbc mol       \
            -ur compact    \
            2>&1

    if [[ ${PIPESTATUS[1]} -eq 0 ]]; then
        echo "  → OK: $OUT_XTC"
        ((processed++))
    else
        echo "  → FAILED: $seq_id"
        rm -f "$OUT_XTC"   # remove partial output
        ((failed++))
    fi

done < "$SEQ_LIST"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================="
echo " Processed : $processed"
echo " Skipped   : $skipped"
echo " Failed    : $failed"
echo "============================="
