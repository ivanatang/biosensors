#!/bin/bash
# submit_prod_qfix_batch.sh
# ─────────────────────────────────────────────────────────────────────────────
# Submits the initial prod_md_PYR1_LCA_qfix.sh chunk (up to 24h) for every
# sequence in a seq_ids_qfix_remaining_89.txt-style list (name, prefix, id,
# dir_type, tab-separated -- see audit_qfix_remaining_89.sh). Each job only
# runs after em_PYR1_LCA_qfix.sh / equil_PYR1_LCA_qfix.sh have already
# completed for that sequence (NPT_qfix/npt.gro must exist).
#
# A single 24h chunk will not reach the full 500ns -- once these jobs finish,
# check how far each got and chain xtnd_prod_PYR1_LCA_qfix.sh <ID> <SEQ_TYPE>
# <PREFIX> for whichever sequences need more time, same as was done for the
# original 95-sequence production and the 6-sequence qfix pilot.
#
# Usage:
#   bash submit_prod_qfix_batch.sh [seq_ids_qfix_remaining_89.txt]
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROD_SCRIPT="${SCRIPT_DIR}/prod_md_PYR1_LCA_qfix.sh"
LOG_DIR="${SCRIPT_DIR}/logs/prod_qfix"
mkdir -p "$LOG_DIR"

SEQ_LIST="${1:-${SCRIPT_DIR}/seq_ids_qfix_remaining_89.txt}"
if [[ ! -f "$SEQ_LIST" ]]; then
    echo "ERROR: seq list not found: $SEQ_LIST"
    exit 1
fi

FAILED_LIST="${SCRIPT_DIR}/submit_prod_qfix_batch_failed.txt"
> "$FAILED_LIST"

submitted=0
failed=0

while IFS=$'\t' read -r name prefix id dir_type; do
    [[ -z "$name" || "$name" == \#* ]] && continue

    sbatch_out=$(sbatch --job-name="prod_qfix_${name}" \
                        --output="${LOG_DIR}/output_${name}_%j.out" \
                        --error="${LOG_DIR}/error_${name}_%j.err" \
                        "$PROD_SCRIPT" "$id" "$dir_type" "$prefix" 2>&1)
    if [[ $? -ne 0 ]]; then
        echo "SUBMIT FAILED: $name  -- $sbatch_out"
        printf '%s\t%s\t%s\t%s\n' "$name" "$prefix" "$id" "$dir_type" >> "$FAILED_LIST"
        failed=$((failed + 1))
        continue
    fi

    echo "Submitting: $name  ($sbatch_out)"
    submitted=$((submitted + 1))
done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Failed    : $failed sequences (sbatch itself rejected the submission)"
if (( failed > 0 )); then
    echo "  $failed failed submissions written to: $FAILED_LIST"
fi
