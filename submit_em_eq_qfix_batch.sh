#!/bin/bash
# submit_em_eq_qfix_batch.sh
# ─────────────────────────────────────────────────────────────────────────────
# For each sequence in a seq_ids_qfix_remaining_89.txt-style list (name,
# prefix, id, dir_type, tab-separated), submits em_PYR1_LCA_qfix.sh and
# chains equil_PYR1_LCA_qfix.sh (NVT+NPT) as a dependent job that only runs
# if EM finishes successfully (--dependency=afterok). Idempotent: skips a
# sequence entirely if EM_qfix+NVT_qfix+NPT_qfix are all already done, and
# skips straight to submitting just equil if EM_qfix is done but NVT/NPT
# aren't -- safe to rerun.
#
# Usage:
#   bash submit_em_eq_qfix_batch.sh [seq_ids_qfix_remaining_89.txt]
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EM_SCRIPT="${SCRIPT_DIR}/em_PYR1_LCA_qfix.sh"
EQ_SCRIPT="${SCRIPT_DIR}/equil_PYR1_LCA_qfix.sh"
LOG_DIR="${SCRIPT_DIR}/logs/em_eq_qfix"
mkdir -p "$LOG_DIR"

SEQ_LIST="${1:-${SCRIPT_DIR}/seq_ids_qfix_remaining_89.txt}"
if [[ ! -f "$SEQ_LIST" ]]; then
    echo "ERROR: seq list not found: $SEQ_LIST"
    exit 1
fi

BASE="/scratch/alpine/ivta1597/LCA_boltz_models"

n_skipped_done=0
n_eq_only=0
n_em_and_eq=0

while IFS=$'\t' read -r name prefix id dir_type; do
    [[ -z "$name" || "$name" == \#* ]] && continue

    case "$dir_type" in
        binders)       suffix="binder"    ;;
        nonbinders)    suffix="nb"        ;;
        neg_low_pkt)   suffix="low_pkt"   ;;
        neg_fail_gate) suffix="fail_gate" ;;
        *) echo "ERROR: unknown dir_type '$dir_type' for $name"; continue ;;
    esac

    seq_dir="${BASE}/${dir_type}/${prefix}_${id}_${suffix}"
    em_done="${seq_dir}/EM_qfix/em.gro"
    nvt_done="${seq_dir}/NVT_qfix/nvt.gro"
    npt_done="${seq_dir}/NPT_qfix/npt.gro"

    if [[ -f "$em_done" && -f "$nvt_done" && -f "$npt_done" ]]; then
        n_skipped_done=$((n_skipped_done + 1))
        continue
    fi

    if [[ -f "$em_done" ]]; then
        # EM already done, just chain equil on its own
        eq_out=$(sbatch --job-name="eq_qfix_${name}" \
                        --output="${LOG_DIR}/output_eq_${name}_%j.out" \
                        --error="${LOG_DIR}/error_eq_${name}_%j.err" \
                        "$EQ_SCRIPT" "$id" "$dir_type" "$prefix" 2>&1)
        eq_id=$(echo "$eq_out" | grep -oE '[0-9]+$')
        echo "EM already done, submitted EQ only: $name  (eq job $eq_id)"
        n_eq_only=$((n_eq_only + 1))
        continue
    fi

    em_out=$(sbatch --parsable --job-name="em_qfix_${name}" \
                    --output="${LOG_DIR}/output_em_${name}_%j.out" \
                    --error="${LOG_DIR}/error_em_${name}_%j.err" \
                    "$EM_SCRIPT" "$id" "$dir_type" "$prefix" 2>&1)
    if ! [[ "$em_out" =~ ^[0-9]+$ ]]; then
        echo "EM SUBMIT FAILED: $name  -- $em_out"
        continue
    fi
    em_id="$em_out"

    eq_out=$(sbatch --dependency=afterok:"$em_id" \
                    --job-name="eq_qfix_${name}" \
                    --output="${LOG_DIR}/output_eq_${name}_%j.out" \
                    --error="${LOG_DIR}/error_eq_${name}_%j.err" \
                    "$EQ_SCRIPT" "$id" "$dir_type" "$prefix" 2>&1)
    eq_id=$(echo "$eq_out" | grep -oE '[0-9]+$')
    echo "Submitted: $name  (em job $em_id -> eq job $eq_id, afterok)"
    n_em_and_eq=$((n_em_and_eq + 1))
done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Already fully done (skipped) : $n_skipped_done"
echo "  EQ-only submitted (EM was already done) : $n_eq_only"
echo "  EM+EQ chained and submitted  : $n_em_and_eq"
