#!/bin/bash
# audit_qfix_remaining_89.sh
# ─────────────────────────────────────────────────────────────────────────────
# Checks, for each of the 89 ngs_observed sequences not in the qfix pilot,
# whether the _qfix topology, EM_qfix, NVT_qfix, and NPT_qfix stages are
# already done -- before launching production for all 89.
#
# Usage: bash audit_qfix_remaining_89.sh [seq_ids_qfix_remaining_89.txt]
# ─────────────────────────────────────────────────────────────────────────────

SEQ_LIST="${1:-seq_ids_qfix_remaining_89.txt}"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"

if [[ ! -f "$SEQ_LIST" ]]; then
    echo "ERROR: seq list not found: $SEQ_LIST"
    exit 1
fi

n_total=0
n_ready=0
missing_topology=()
missing_em=()
missing_nvt=()
missing_npt=()

while IFS=$'\t' read -r name prefix id dir_type; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    n_total=$((n_total + 1))

    case "$dir_type" in
        binders)      suffix="binder"   ;;
        nonbinders)   suffix="nb"       ;;
        neg_low_pkt)  suffix="low_pkt"  ;;
        neg_fail_gate) suffix="fail_gate" ;;
        *) echo "ERROR: unknown dir_type '$dir_type' for $name"; continue ;;
    esac

    seq_dir="${BASE}/${dir_type}/${prefix}_${id}_${suffix}"
    top="${seq_dir}/${prefix}_${id}_${suffix}_dodecahedron_HMR_qfix.top"
    gro="${seq_dir}/${prefix}_${id}_${suffix}_dodecahedron_HMR_qfix.gro"
    em="${seq_dir}/EM_qfix/em.gro"
    nvt="${seq_dir}/NVT_qfix/nvt.gro"
    npt="${seq_dir}/NPT_qfix/npt.gro"

    ok=true
    [[ -f "$top" && -f "$gro" ]] || { missing_topology+=("$name"); ok=false; }
    [[ -f "$em" ]]  || { missing_em+=("$name");  ok=false; }
    [[ -f "$nvt" ]] || { missing_nvt+=("$name"); ok=false; }
    [[ -f "$npt" ]] || { missing_npt+=("$name"); ok=false; }

    $ok && n_ready=$((n_ready + 1))
done < "$SEQ_LIST"

echo "============================================================"
echo "  Qfix pipeline readiness: $n_ready / $n_total ready for production"
echo "============================================================"
echo "Missing topology (${#missing_topology[@]}): ${missing_topology[*]}"
echo "Missing EM_qfix  (${#missing_em[@]}): ${missing_em[*]}"
echo "Missing NVT_qfix (${#missing_nvt[@]}): ${missing_nvt[*]}"
echo "Missing NPT_qfix (${#missing_npt[@]}): ${missing_npt[*]}"
