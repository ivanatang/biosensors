#!/bin/bash
# submit_ba_md.sh -- runs the *_qfix_batch.sh drivers on the bile acid
# (LCA / LCA-3-S) arm set: points them at the BA scratch base, uses the
# two-stage restrained EM (em_PYR1_BA.sh), and defaults to seq_ids_ba.txt.
#
# Usage (on Alpine, from /projects/ivta1597/biosensors):
#   bash submit_ba_md.sh em_eq    [seq_list]   # restrained+unrestrained EM -> NVT/NPT
#   bash submit_ba_md.sh prod     [seq_list]   # first 24h production chunk
#   bash submit_ba_md.sh xtnd     [seq_list]   # extension chained after prod
#   bash submit_ba_md.sh progress [seq_list]   # ns reached per arm
# Pass a subset list (same 2-column format) to submit in waves.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BASE=/scratch/alpine/ivta1597/BA_boltz_models
export EM_SCRIPT="${SCRIPT_DIR}/em_PYR1_BA.sh"

STAGE="$1"
SEQ_LIST="${2:-${SCRIPT_DIR}/seq_ids_ba.txt}"

case "$STAGE" in
    em_eq)    exec bash "${SCRIPT_DIR}/submit_em_eq_qfix_batch.sh"     "$SEQ_LIST" ;;
    prod)     exec bash "${SCRIPT_DIR}/submit_prod_qfix_batch.sh"      "$SEQ_LIST" ;;
    xtnd)     exec bash "${SCRIPT_DIR}/submit_xtnd_prod_qfix_batch.sh" "$SEQ_LIST" ;;
    progress) exec bash "${SCRIPT_DIR}/check_qfix_prod_progress.sh"    "$SEQ_LIST" ;;
    *) echo "Usage: bash submit_ba_md.sh {em_eq|prod|xtnd|progress} [seq_list]" >&2; exit 1 ;;
esac
