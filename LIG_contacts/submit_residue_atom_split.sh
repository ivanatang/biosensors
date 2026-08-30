#!/bin/bash
# submit_residue_atom_split.sh
# ------------------------------
# Submits one SLURM job per sequence in seq_list, splitting a single
# residue's ligand contact into backbone vs side chain.
#
# Usage:
#   bash submit_residue_atom_split.sh                                    # seq_ids_orig.txt, res 116, core, 40-500ns
#   bash submit_residue_atom_split.sh seq_ids_orig.txt 116 core           # explicit residue + region
#   bash submit_residue_atom_split.sh seq_ids_orig.txt 116 whole 40 500   # whole-ligand region
#
# Arguments:
#   $1  seq_list       - path to seq_ids.txt-style sequence list
#                         (default: /projects/ivta1597/biosensors/seq_ids_orig.txt)
#   $2  resseq          (default: 116)
#   $3  ligand_region   - whole | core | tail (default: core)
#   $4  start_ns         (default: 40)
#   $5  end_ns           (default: 500)

BASE="/projects/ivta1597/biosensors/LIG_contacts"
WORKER="${BASE}/residue_atom_split_worker.sh"

# SLURM opens --output/--error at launch time using the worker's own
# hardcoded absolute path, and does not create missing parent directories --
# without this, every submission fails immediately with nowhere to write,
# and no log file anywhere to explain why.
mkdir -p "${BASE}/logs"

SEQ_IDS_FILE="${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}"
RESSEQ="${2:-116}"
LIGAND_REGION="${3:-core}"
START_NS="${4:-40}"
END_NS="${5:-500}"

echo "============================================================"
echo "  Residue backbone/side-chain contact split submission"
echo "  Seq list : $SEQ_IDS_FILE"
echo "  Residue  : ${RESSEQ}"
echo "  Region   : ${LIGAND_REGION}"
echo "  Window   : ${START_NS}-${END_NS} ns"
echo "  Output   : residue_atomsplit_results_${START_NS}_${END_NS}ns$([ "$LIGAND_REGION" != "whole" ] && echo "_${LIGAND_REGION}")/"
echo "============================================================"

submitted=0

while read -r SEQ_ID _rest || [[ -n "${SEQ_ID}" ]]; do
    [[ -z "${SEQ_ID}" || "${SEQ_ID}" == \#* ]] && continue
    echo "Submitting: ${SEQ_ID}  res=${RESSEQ}  region=${LIGAND_REGION}  window=${START_NS}-${END_NS}ns"
    sbatch --export=SEQ_ID="${SEQ_ID}",RESSEQ="${RESSEQ}",LIGAND_REGION="${LIGAND_REGION}",START_NS="${START_NS}",END_NS="${END_NS}" \
           "${WORKER}"
    ((submitted++))
done < "${SEQ_IDS_FILE}"

echo ""
echo "=== Done: submitted ${submitted} jobs ==="
