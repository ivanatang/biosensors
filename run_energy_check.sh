#!/bin/bash
# run_energy_check.sh
# ─────────────────────────────────────────────────────────────────────────────
# Loops through seq_ids.txt and, for each sequence, runs (with SUFFIX=""
# by default, or e.g. "_qfix" for the bond-order/charge-fix systems):
#
#   echo "10 0" | gmx energy -f EM${SUFFIX}/em.edr   -o EM${SUFFIX}/em_potential${SUFFIX}.xvg   # Potential
#   echo "16 0" | gmx energy -f NVT${SUFFIX}/nvt.edr -o NVT${SUFFIX}/nvt_temp${SUFFIX}.xvg      # Temperature
#   echo "24 0" | gmx energy -f NPT${SUFFIX}/npt.edr -o NPT${SUFFIX}/npt_density${SUFFIX}.xvg   # Density
#   echo "23 0" | gmx energy -f NPT${SUFFIX}/npt.edr -o NPT${SUFFIX}/npt_volume${SUFFIX}.xvg    # Volume
#
# Usage:
#   bash run_energy_check.sh                          # seq_ids.txt, standard EM/NVT/NPT
#   bash run_energy_check.sh my_seq_list.txt           # different seq list
#   bash run_energy_check.sh seq_ids.txt _qfix         # EM_qfix/NVT_qfix/NPT_qfix
#   bash run_energy_check.sh <(grep -E "^cdca_001|^glca_001" seq_ids.txt)
#                                                       # scope to specific sequences without
#                                                       # a separate tracked seq-list file
#
# seq_ids.txt format (tab-separated):
#   seq_id       seq_type (display)      optional_custom_path
#   bind_043     Binder
#   nonb_046     False Positive
#   seq14_binder Binder                  /scratch/.../seq14_binder   ← custom path used as-is
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

SEQ_LIST=${1:-seq_ids.txt}
SUFFIX=${2:-}

# ── Path configuration ─────────────────────────────────────────────────────
BASE_DIR="/scratch/alpine/ivta1597/LCA_boltz_models"
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
# If EM/NVT/NPT live one level deeper than ${BASE_DIR}/${dir_type}/${seq_id}*/
# (e.g. under an HMR/dodecahedron subfolder), set SUBDIR accordingly.
SUBDIR=""
# ─────────────────────────────────────────────────────────────────────────────

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

if [ ! -x "$GMX" ]; then
    echo "WARNING: $GMX not found/executable — falling back to 'gmx' from PATH"
    GMX="gmx"
fi

# ── Map display seq_type → directory name used in the file system ──────────
get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"      ;;
        "False Positive") echo "nonbinders"   ;;
        "Low Confidence") echo "neg_low_pkt"  ;;
        "Fail Geometry")  echo "neg_fail_gate";;
        *)                echo "$1"           ;;   # fallback: use as-is
    esac
}

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    # Skip empty lines and comments
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    if [[ -n "$custom_path" ]]; then
        WORKDIR="$custom_path"
    else
        dir_type=$(get_dir_type "$seq_type")
        parent="${BASE_DIR}/${dir_type}"

        # Directories carry a type suffix appended to seq_id, e.g.
        # bind_043 -> bind_043_binder, nonb_046 -> nonb_046_nb — glob for it
        # rather than hardcoding every suffix variant.
        shopt -s nullglob
        candidates=("${parent}/${seq_id}"*/)
        shopt -u nullglob

        if [[ ${#candidates[@]} -eq 0 ]]; then
            echo "SKIP (no directory matching ${seq_id}* under $parent): $seq_id"
            continue
        elif [[ ${#candidates[@]} -gt 1 ]]; then
            echo "WARNING: multiple dirs match ${seq_id}* under $parent — using first: ${candidates[0]}"
        fi

        WORKDIR="${candidates[0]%/}"
    fi
    [[ -n "$SUBDIR" ]] && WORKDIR="${WORKDIR}/${SUBDIR}"

    echo "Running energy extraction: $seq_id  [$seq_type -> $WORKDIR]  (suffix='${SUFFIX}')"

    # Skip per-output-file if it already exists -- gmx energy doesn't
    # overwrite in place, it renames the existing file to #name.N# and
    # piles those up on every re-run.
    extract() {
        local edr="$1" terms="$2" out="$3"
        if [[ -f "${WORKDIR}/${out}" ]]; then
            echo "  SKIP (exists): $out"
        else
            echo "$terms" | "$GMX" energy -f "$edr" -o "$out"
        fi
    }

    (
        cd "$WORKDIR" || exit 1
        extract "EM${SUFFIX}/em.edr"   "10 0" "EM${SUFFIX}/em_potential${SUFFIX}.xvg"
        extract "NVT${SUFFIX}/nvt.edr" "16 0" "NVT${SUFFIX}/nvt_temp${SUFFIX}.xvg"
        extract "NPT${SUFFIX}/npt.edr" "24 0" "NPT${SUFFIX}/npt_density${SUFFIX}.xvg"
        extract "NPT${SUFFIX}/npt.edr" "23 0" "NPT${SUFFIX}/npt_volume${SUFFIX}.xvg"
    )

done < "$SEQ_LIST"
