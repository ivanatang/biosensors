#!/bin/bash
# this needs debugging - the dependency job doesn't work as is
# submit_prod_jobs.sh — submits initial prod run + chained continuation
# Usage: submit_prod_jobs.sh <seq_type> <ID>
#   seq_type: binder | nonbinder | low_pkt | fail_gate

SEQ_TYPE=$1
ID=$2

if [[ -z "$SEQ_TYPE" || -z "$ID" ]]; then
    echo "Usage: $0 <seq_type> <ID>"
    echo "  seq_type: binder | nonbinder | low_pkt | fail_gate"
    exit 1
fi

case "$SEQ_TYPE" in
    binder)
        PROD_SCRIPT=prod_md_binders_PYR1_LCA.sh
        XTND_SCRIPT=xtnd_b_prod_PYR1_LCA.sh
        ;;
    nonbinder)
        PROD_SCRIPT=prod_md_nonbinders_PYR1_LCA.sh
        XTND_SCRIPT=xtnd_nb_prod_PYR1_LCA.sh
        ;;
    low_pkt)
        PROD_SCRIPT=prod_md_low_pkt_PYR1_LCA.sh
        XTND_SCRIPT=xtnd_low_pkt_prod_PYR1_LCA.sh
        ;;
    fail_gate)
        PROD_SCRIPT=prod_md_fail_gate_PYR1_LCA.sh
        XTND_SCRIPT=xtnd_fail_gate_prod_PYR1_LCA.sh
        ;;
    *)
        echo "Unknown seq_type: $SEQ_TYPE"
        echo "  must be one of: binder, nonbinder, low_pkt, fail_gate"
        exit 1
        ;;
esac

jobid1=$(sbatch --parsable $PROD_SCRIPT $ID)
jobid2=$(sbatch --parsable --dependency=afterany:$jobid1 $XTND_SCRIPT $ID)

echo "seq_type: $SEQ_TYPE   ID: $ID"
echo "initial: $jobid1   continuation: $jobid2"
