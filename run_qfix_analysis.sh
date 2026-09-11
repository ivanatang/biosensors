#!/bin/bash
set -e
module load anaconda
conda activate biosensors
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"
cd /projects/ivta1597/biosensors

echo "=== compare_qfix_vs_standard.py ==="
python compare_qfix_vs_standard.py --qfix_table qfix_500ns_feat_table.csv \
    --target_seq_list seq_ids_qfix_all95.txt \
    --out_long qfix_vs_standard_deltas_long_all95.csv \
    --out_summary qfix_vs_standard_summary_all95.csv

echo ""
echo "=== model_swap_eval.py ==="
python model_swap_eval.py --qfix_table qfix_500ns_feat_table.csv \
    --target_seq_list seq_ids_qfix_all95.txt \
    --out model_swap_eval_all95_results.csv
