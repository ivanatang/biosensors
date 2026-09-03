#!/bin/bash
# submit_xtnd_prod_qfix_batch.sh
# ─────────────────────────────────────────────────────────────────────────────
# For each sequence in a seq_ids_qfix_remaining_89.txt-style list (name,
# prefix, id, dir_type, tab-separated), looks up that sequence's currently
# running/queued initial production job and submits xtnd_prod_PYR1_LCA_qfix.sh
# with --dependency=afterany so it only starts once the initial chunk ends --
# afterany, not afterok, since a 24h chunk is expected to be cut off by the
# walltime (TIMEOUT), not exit cleanly, and the extension needs to run either
# way to pick up from the checkpoint.
#
# submit_prod_qfix_batch.sh submits the job as "prod_qfix_<name>", but
# prod_md_PYR1_LCA_qfix.sh renames itself mid-run via scontrol (its #SBATCH
# job-name is static text set before ID/SEQ_TYPE/PREFIX are known, so it
# relabels once they are) to "<prefix>_<id>_prod_qfix" -- sacct records the
# renamed name, not the one it was submitted under, so that's what we search
# for below.
#
# Usage:
#   bash submit_xtnd_prod_qfix_batch.sh [seq_ids_qfix_remaining_89.txt]
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
XTND_SCRIPT="${SCRIPT_DIR}/xtnd_prod_PYR1_LCA_qfix.sh"
LOG_DIR="${SCRIPT_DIR}/logs/xtnd_prod_qfix"
mkdir -p "$LOG_DIR"

SEQ_LIST="${1:-${SCRIPT_DIR}/seq_ids_qfix_remaining_89.txt}"
if [[ ! -f "$SEQ_LIST" ]]; then
    echo "ERROR: seq list not found: $SEQ_LIST"
    exit 1
fi

NOT_FOUND_LIST="${SCRIPT_DIR}/submit_xtnd_prod_qfix_batch_no_job_found.txt"
> "$NOT_FOUND_LIST"

submitted=0
not_found=0

while IFS=$'\t' read -r name prefix id dir_type; do
    [[ -z "$name" || "$name" == \#* ]] && continue

    # sacct sees both still-running and already-finished (e.g. TIMEOUT) jobs
    # by name, unlike squeue which drops a job once it leaves the queue --
    # -X restricts to the main job record, excluding .batch/.extern substeps.
    # sort -n + tail -1 picks the most recent submission if there were retries.
    #
    # -S is required: sacct's default start time is midnight of *today*, not
    # "all history" -- without it, a job that finished on an earlier day is
    # silently excluded regardless of name, which is what made every lookup
    # come back empty even though the runs clearly completed on disk.
    #
    # Search by "<prefix>_<id>_prod_qfix", not "prod_qfix_<name>" -- see the
    # header comment: prod_md_PYR1_LCA_qfix.sh renames the job mid-run, and
    # sacct only ever has the renamed JobName on record.
    jobid=$(sacct -u "$(whoami)" --name="${prefix}_${id}_prod_qfix" \
                  -S 2025-01-01 \
                  --format=JobID,State --noheader --parsable2 -X 2>/dev/null \
            | awk -F'|' '{print $1}' | grep -E '^[0-9]+$' | sort -n | tail -1)

    if [[ -z "$jobid" ]]; then
        echo "NO PROD JOB FOUND: $name -- skipping"
        printf '%s\t%s\t%s\t%s\n' "$name" "$prefix" "$id" "$dir_type" >> "$NOT_FOUND_LIST"
        not_found=$((not_found + 1))
        continue
    fi

    xtnd_out=$(sbatch --dependency=afterany:"$jobid" \
                      --job-name="xtnd_qfix_${name}" \
                      --output="${LOG_DIR}/output_${name}_%j.out" \
                      --error="${LOG_DIR}/error_${name}_%j.err" \
                      "$XTND_SCRIPT" "$id" "$dir_type" "$prefix" 2>&1)

    echo "Submitted: $name  (after prod job $jobid)  -- $xtnd_out"
    submitted=$((submitted + 1))
done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted (dependent on prod job) : $submitted"
echo "  No prod job found (skipped)       : $not_found"
if (( not_found > 0 )); then
    echo "  Written to: $NOT_FOUND_LIST"
    echo "  (a sequence's prod job may have already finished within 24h and"
    echo "   already produced 500ns -- check before resubmitting these manually)"
fi
