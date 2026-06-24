#!/bin/bash
# run_energy_check.sh
# ─────────────────────────────────────────────────────────────────────────────
# Loops through seq_ids.txt and, for each sequence, runs:
#
#   echo "10 0" | gmx energy -f EM/em.edr  -o EM/em_potential.xvg   # Potential
#   echo "16 0" | gmx energy -f NVT/nvt.edr -o NVT/nvt_temp.xvg     # Temperature
#   echo "24 0" | gmx energy -f NPT/npt.edr -o NPT/npt_density.xvg  # Density
#   echo "23 0" | gmx energy -f NPT/npt.edr -o NPT/npt_volume.xvg   # Volume
#
# Usage:
#   bash run_energy_check.sh                  # uses seq_ids.txt in cwd by default
#   bash run_energy_check.sh my_seq_list.txt  # pass a different file
#
# seq_ids.txt format (tab-separated):
#   seq_id       seq_type (display)      optional_custom_path
#   bind_043     Binder
#   nonb_046     False Positive
#   seq14_binder Binder                  /scratch/.../seq14_binder   ← custom path used as-is
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

SEQ_LIST=${1:-seq_ids.txt}

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

    echo "Running energy extraction: $seq_id  [$seq_type -> $WORKDIR]"

    (
        cd "$WORKDIR" || exit 1
        echo "10 0" | "$GMX" energy -f EM/em.edr   -o EM/em_potential.xvg
        echo "16 0" | "$GMX" energy -f NVT/nvt.edr -o NVT/nvt_temp.xvg
        echo "24 0" | "$GMX" energy -f NPT/npt.edr -o NPT/npt_density.xvg
        echo "23 0" | "$GMX" energy -f NPT/npt.edr -o NPT/npt_volume.xvg
    )

done < "$SEQ_LIST"
