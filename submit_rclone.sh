#!/bin/bash
# submit_rclone.sh
# ─────────────────────────────────────────────────────────────────────────────
# Submits one rclone SLURM job per sequence in a seq_ids_qfix_remaining_89.txt-
# style list (name, prefix, id, dir_type, tab-separated), copying each
# sequence's full scratch run directory to the OneDrive backup
# (Shirts Lab/LCA_boltz_models), per CLAUDE.md's storage architecture.
#
# This is the OneDrive backup only -- it does NOT populate the PetaLibrary
# archive (/pl/active/shirts_archive/...) that contact_type_analysis.py,
# salt_bridge_analysis.py (via config_qfix.yaml), and
# extract_gate_latch_rmsd_feats.py read from. That copy is handled
# separately via Globus.
#
# Usage:
#   bash submit_rclone.sh [seq_ids_qfix_remaining_89.txt]
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEQ_LIST="${1:-${SCRIPT_DIR}/seq_ids_qfix_remaining_89.txt}"
BASE_SCRATCH=/scratch/alpine/ivta1597/LCA_boltz_models
BASE_ONEDRIVE="onedrive_ivana_cu:Shirts Lab/LCA_boltz_models"
LOG_DIR="/projects/ivta1597/biosensors/rclone_logs"

if [[ ! -f "$SEQ_LIST" ]]; then
    echo "ERROR: seq list not found: $SEQ_LIST"
    exit 1
fi

mkdir -p "$LOG_DIR"

submitted=0
skipped=()

while IFS=$'\t' read -r name prefix id dir_type; do
    [[ -z "$name" || "$name" == \#* ]] && continue

    src="${BASE_SCRATCH}/${dir_type}/${name}/"
    if [[ ! -d "$src" ]]; then
        skipped+=("$name")
        continue
    fi

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=rclone_${name}
#SBATCH --output=${LOG_DIR}/output_${name}_%j.out
#SBATCH --error=${LOG_DIR}/error_${name}_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=FAIL

module load slurm/alpine
module load rclone/1.58.0

rclone copy \
    "${src}" \
    "${BASE_ONEDRIVE}/${dir_type}/${name}/" \
    --transfers 4 \
    --checkers 8 \
    --log-file ${LOG_DIR}/rclone_${name}_\${SLURM_JOB_ID}.log \
    --log-level INFO
EOF

    echo "Submitted job for $name"
    submitted=$((submitted + 1))
done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
if [[ ${#skipped[@]} -gt 0 ]]; then
    echo "  Skipped (directory not found): ${#skipped[@]}"
    for s in "${skipped[@]}"; do
        echo "    $s"
    done
fi
