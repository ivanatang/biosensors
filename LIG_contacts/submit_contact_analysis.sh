#!/bin/bash
# submit_contact_analysis.sh
# ---------------------------
# Submits one SLURM job per sequence in seq_list
#
# Usage:
#   bash submit_contact_analysis.sh                                   # seq_ids_orig.txt, full 40–500 ns, whole ligand
#   bash submit_contact_analysis.sh seq_ids_missing_water.txt         # custom list, full 40–500 ns
#   bash submit_contact_analysis.sh seq_ids_orig.txt 40 250           # 250 ns window
#   bash submit_contact_analysis.sh seq_ids_orig.txt 40 300           # 300 ns window
#   bash submit_contact_analysis.sh seq_ids_orig.txt 40 500 core      # steroid core only
#   bash submit_contact_analysis.sh seq_ids_orig.txt 40 500 tail      # carboxylate tail only
#
# Arguments:
#   $1  seq_list       - path to seq_ids.txt-style sequence list
#                         (default: /projects/ivta1597/biosensors/seq_ids_orig.txt)
#   $2  start_ns        (default: 40)
#   $3  end_ns          (default: 500)
#   $4  ligand_region   - whole | core | tail (default: whole)

BASE="/projects/ivta1597/biosensors/LIG_contacts"
WORKER="${BASE}/contact_type_worker.sh"

SEQ_IDS_FILE="${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}"
START_NS="${2:-40}"
END_NS="${3:-500}"
LIGAND_REGION="${4:-whole}"

FAILED_LIST="${BASE}/submit_contact_analysis_failed_${START_NS}_${END_NS}ns_${LIGAND_REGION}.txt"
> "$FAILED_LIST"

echo "============================================================"
echo "  Contact type analysis submission"
echo "  Seq list : $SEQ_IDS_FILE"
echo "  Window   : ${START_NS}–${END_NS} ns"
echo "  Region   : ${LIGAND_REGION}"
echo "  Output   : contact_type_results_${START_NS}_${END_NS}ns$([ "$LIGAND_REGION" != "whole" ] && echo "_${LIGAND_REGION}")/"
echo "============================================================"

submitted=0
failed=0

while read -r SEQ_ID _rest || [[ -n "${SEQ_ID}" ]]; do
    [[ -z "${SEQ_ID}" || "${SEQ_ID}" == \#* ]] && continue
    sbatch_out=$(sbatch --job-name="ct_${SEQ_ID}" \
                        --output="${BASE}/output_ct_${SEQ_ID}_%j.out" \
                        --error="${BASE}/error_ct_${SEQ_ID}_%j.err" \
                        --export=SEQ_ID="${SEQ_ID}",START_NS="${START_NS}",END_NS="${END_NS}",LIGAND_REGION="${LIGAND_REGION}" \
                        "${WORKER}" 2>&1)
    if [[ $? -ne 0 ]]; then
        echo "SUBMIT FAILED: ${SEQ_ID}  -- $sbatch_out"
        printf '%s\n' "${SEQ_ID}" >> "$FAILED_LIST"
        ((failed++))
        continue
    fi
    echo "Submitting: ${SEQ_ID}  window=${START_NS}–${END_NS}ns  region=${LIGAND_REGION}  ($sbatch_out)"
    ((submitted++))
done < "${SEQ_IDS_FILE}"

echo ""
echo "=== Done: submitted ${submitted} jobs, ${failed} failed ==="
if (( failed > 0 )); then
    echo "  $failed failed submissions written to: $FAILED_LIST"
    echo "  Re-submit just those: bash ${BASE}/submit_contact_analysis.sh $FAILED_LIST $START_NS $END_NS $LIGAND_REGION"
fi
