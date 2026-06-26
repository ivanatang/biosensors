#!/bin/bash
# submit_contact_analysis.sh
# ---------------------------
# Submits one SLURM job per sequence in seq_ids.txt
#
# Usage:
#   bash submit_contact_analysis.sh                  # full 40–500 ns (default)
#   bash submit_contact_analysis.sh 40 250           # 250 ns window
#   bash submit_contact_analysis.sh 40 300           # 300 ns window
#
# Arguments:
#   $1  start_ns  (default: 40)
#   $2  end_ns    (default: 500)

BASE="/scratch/alpine/ivta1597/LCA_boltz_models/LIG_contacts_flex"
SEQ_IDS_FILE="${BASE}/seq_ids.txt"
WORKER="${BASE}/contact_type_worker.sh"

START_NS="${1:-40}"
END_NS="${2:-500}"

echo "============================================================"
echo "  Contact type analysis submission"
echo "  Seq list : $SEQ_IDS_FILE"
echo "  Window   : ${START_NS}–${END_NS} ns"
echo "  Output   : contact_type_results_${START_NS}_${END_NS}ns/"
echo "============================================================"

submitted=0

while read -r SEQ_ID _rest || [[ -n "${SEQ_ID}" ]]; do
    [[ -z "${SEQ_ID}" || "${SEQ_ID}" == \#* ]] && continue
    echo "Submitting: ${SEQ_ID}  window=${START_NS}–${END_NS}ns"
    sbatch --export=SEQ_ID="${SEQ_ID}",START_NS="${START_NS}",END_NS="${END_NS}" \
           "${WORKER}"
    ((submitted++))
done < "${SEQ_IDS_FILE}"

echo ""
echo "=== Done: submitted ${submitted} jobs ==="
