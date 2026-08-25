#!/bin/bash
# run_energy_check_qfix.sh
# ─────────────────────────────────────────────────────────────────────────────
# Same extraction as run_energy_check.sh, adapted for the bond-order/charge-fix
# ("_qfix") pilot systems: reads EM_qfix/NVT_qfix/NPT_qfix instead of
# EM/NVT/NPT, and is scoped to a fixed pilot list rather than looping over
# seq_ids.txt.
#
#   echo "10 0" | gmx energy -f EM_qfix/em.edr   -o EM_qfix/em_potential_qfix.xvg   # Potential
#   echo "16 0" | gmx energy -f NVT_qfix/nvt.edr -o NVT_qfix/nvt_temp_qfix.xvg      # Temperature
#   echo "24 0" | gmx energy -f NPT_qfix/npt.edr -o NPT_qfix/npt_density_qfix.xvg   # Density
#   echo "23 0" | gmx energy -f NPT_qfix/npt.edr -o NPT_qfix/npt_volume_qfix.xvg    # Volume
#
# Energy-term indices (10/16/24/23) are unchanged from run_energy_check.sh --
# they come from the .mdp integrator/thermostat/barostat settings, which are
# byte-identical between the original and _qfix systems (same nvt.mdp/npt.mdp),
# not from the topology, so they apply the same way here.
#
# Usage:
#   bash run_energy_check_qfix.sh
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

BASE_DIR="/scratch/alpine/ivta1597/LCA_boltz_models"
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
if [ ! -x "$GMX" ]; then
    echo "WARNING: $GMX not found/executable — falling back to 'gmx' from PATH"
    GMX="gmx"
fi

# type_dir seq_id (pilot: 3 binders, 3 false positives)
PILOT=(
    "binders bind_022_binder"
    "binders bind_019_binder"
    "binders bind_020_binder"
    "nonbinders nonb_006_nb"
    "nonbinders nonb_008_nb"
    "nonbinders nonb_009_nb"
)

for pair in "${PILOT[@]}"; do
    set -- $pair
    type_dir="$1"; seq_id="$2"
    WORKDIR="${BASE_DIR}/${type_dir}/${seq_id}"

    if [[ ! -d "$WORKDIR" ]]; then
        echo "SKIP (not found): $WORKDIR"
        continue
    fi
    if [[ ! -f "$WORKDIR/NPT_qfix/npt.edr" ]]; then
        echo "SKIP (NPT_qfix/npt.edr not found -- equil not finished?): $seq_id"
        continue
    fi

    echo "Running energy extraction: $seq_id  [$WORKDIR]"
    (
        cd "$WORKDIR" || exit 1
        echo "10 0" | "$GMX" energy -f EM_qfix/em.edr   -o EM_qfix/em_potential_qfix.xvg
        echo "16 0" | "$GMX" energy -f NVT_qfix/nvt.edr -o NVT_qfix/nvt_temp_qfix.xvg
        echo "24 0" | "$GMX" energy -f NPT_qfix/npt.edr -o NPT_qfix/npt_density_qfix.xvg
        echo "23 0" | "$GMX" energy -f NPT_qfix/npt.edr -o NPT_qfix/npt_volume_qfix.xvg
    )
done
