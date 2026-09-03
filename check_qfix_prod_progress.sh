#!/bin/bash
# check_qfix_prod_progress.sh
# ─────────────────────────────────────────────────────────────────────────────
# For each sequence in a seq_ids_qfix_remaining_89.txt-style list (name,
# prefix, id, dir_type, tab-separated), reads the last "Step Time" entry
# logged in prod_md_500ns.log to see how far the initial 24h
# prod_md_PYR1_LCA_qfix.sh chunk got before its walltime cut it off.
#
# Writes any sequence that has a log but hasn't reached 500ns to
# <SEQ_LIST_basename>_incomplete.txt in the same 4-column format, ready to
# hand to submit_xtnd_prod_qfix_batch.sh.
#
# Run this on Alpine (a login node is fine) -- it reads scratch, not OneDrive.
#
# Usage:
#   bash check_qfix_prod_progress.sh [seq_ids_qfix_remaining_89.txt]
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEQ_LIST="${1:-${SCRIPT_DIR}/seq_ids_qfix_remaining_89.txt}"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
PROD_DIR="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd_qfix"
TARGET_PS=500000  # 500ns, in ps (log times are ps)

if [[ ! -f "$SEQ_LIST" ]]; then
    echo "ERROR: seq list not found: $SEQ_LIST"
    exit 1
fi

INCOMPLETE_LIST="${SCRIPT_DIR}/$(basename "${SEQ_LIST%.txt}")_incomplete.txt"
> "$INCOMPLETE_LIST"

n_total=0
n_complete=0
n_incomplete=0
n_no_log=0

printf "%-22s %-10s %-14s\n" "NAME" "STATUS" "LAST_TIME_NS"
echo "------------------------------------------------------------"

while IFS=$'\t' read -r name prefix id dir_type; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    n_total=$((n_total + 1))

    case "$dir_type" in
        binders)       suffix="binder"    ;;
        nonbinders)    suffix="nb"        ;;
        neg_low_pkt)   suffix="low_pkt"   ;;
        neg_fail_gate) suffix="fail_gate" ;;
        *) echo "ERROR: unknown dir_type '$dir_type' for $name"; continue ;;
    esac

    seq_dir="${BASE}/${dir_type}/${prefix}_${id}_${suffix}"
    prod_log="${seq_dir}/${PROD_DIR}/prod_md_500ns.log"

    if [[ ! -f "$prod_log" ]]; then
        printf "%-22s %-10s %-14s\n" "$name" "NO_LOG" "--"
        n_no_log=$((n_no_log + 1))
        continue
    fi

    # GROMACS periodically prints a "Step  Time" header followed by a line
    # of values -- same pattern used in check_xtc_files.sh -- tail -1 picks
    # the most recent one.
    last_time_ps=$(grep -a -A1 "^ *Step *Time$" "$prod_log" \
                    | grep -v "Step\|--" | awk '{print $2}' | tail -1)

    if [[ -z "$last_time_ps" ]]; then
        printf "%-22s %-10s %-14s\n" "$name" "UNREADABLE" "--"
        n_no_log=$((n_no_log + 1))
        continue
    fi

    last_time_ns=$(awk -v t="$last_time_ps" 'BEGIN{printf "%.1f", t/1000}')

    if awk -v t="$last_time_ps" -v target="$TARGET_PS" 'BEGIN{exit !(t >= target)}'; then
        printf "%-22s %-10s %-14s\n" "$name" "COMPLETE" "$last_time_ns"
        n_complete=$((n_complete + 1))
    else
        printf "%-22s %-10s %-14s\n" "$name" "INCOMPLETE" "$last_time_ns"
        printf '%s\t%s\t%s\t%s\n' "$name" "$prefix" "$id" "$dir_type" >> "$INCOMPLETE_LIST"
        n_incomplete=$((n_incomplete + 1))
    fi
done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Total       : $n_total"
echo "  Complete    : $n_complete (already at/past 500ns, no extension needed)"
echo "  Incomplete  : $n_incomplete (timed out short of 500ns -- written to $INCOMPLETE_LIST)"
echo "  No/bad log  : $n_no_log (never started or crashed before first log entry -- check these manually)"
