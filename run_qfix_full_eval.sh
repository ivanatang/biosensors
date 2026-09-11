#!/bin/bash
set -e
module load anaconda
conda activate biosensors
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"
cd /projects/ivta1597/biosensors
python qfix_full_performance_eval.py --qfix_table qfix_500ns_feat_table.csv \
    --target_seq_list seq_ids_qfix_all95.txt
