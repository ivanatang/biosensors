#!/bin/bash
# seq_list_utils.sh -- shared seq-list line parser for the *_qfix_batch.sh
# submit drivers. Source it, then call parse_seq_line on each line.
#
# Accepts both list formats in use:
#   4-column  name prefix id dir_type      (seq_ids_qfix_remaining_89.txt)
#   2-column  seq_id label                 (seq_ids_ba.txt, standard format)
# For the 2-column form, prefix/id/dir_type are derived from the seq_id
# (<prefix>_<id>_<suffix>), so the BA set needs only one list for both MD
# submission and featurization.

# Parses one tab-separated seq-list line into globals name, prefix, id,
# dir_type. Returns 1 (caller should skip) for blank/comment lines or an
# unrecognized suffix.
parse_seq_line() {
    local line="$1" f2 f3 f4 suffix stem
    IFS=$'\t' read -r name f2 f3 f4 <<< "$line"
    [[ -z "$name" || "$name" == \#* ]] && return 1

    if [[ -n "$f4" ]]; then
        prefix="$f2"; id="$f3"; dir_type="$f4"
        return 0
    fi

    case "$name" in
        *_fail_gate) suffix="fail_gate"; dir_type="neg_fail_gate" ;;
        *_low_pkt)   suffix="low_pkt";   dir_type="neg_low_pkt"   ;;
        *_binder)    suffix="binder";    dir_type="binders"       ;;
        *_nb)        suffix="nb";        dir_type="nonbinders"    ;;
        *) echo "ERROR: unrecognized seq_id suffix: $name" >&2; return 1 ;;
    esac
    stem="${name%_${suffix}}"
    prefix="${stem%%_*}"
    id="${stem#*_}"
    return 0
}
