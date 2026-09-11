#!/bin/bash
set -e
module load anaconda
source /projects/ivta1597/software/shap_venv/bin/activate
export MPLBACKEND=Agg
cd /projects/ivta1597/biosensors
python qfix_shap_summary.py --qfix_table qfix_500ns_feat_table.csv \
    --target_seq_list seq_ids_qfix_all95.txt
echo "DONE_MARKER_OK"
