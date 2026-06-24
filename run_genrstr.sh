#!/bin/bash
# run_genrstr.sh
# ─────────────────────────────────────────────────────────────────────────────
# Loops through seq_ids.txt and, for each sequence:
#   1. runs:
#        echo 4 | gmx genrestr -f EM/em.gro -o posre_protein.itp -fc 1000 1000 1000
#      (selection "4" = Backbone in the default GROMACS index groups)
#   2. inserts the POSRES include block right after the "#include ... protein"
#      line in {seq_id}_dodecahedron_HMR.top:
#        #ifdef POSRES
#        #include "posre_protein.itp"
#        #endif
#
# Usage:
#   bash run_genrestr.sh                  # uses seq_ids.txt in cwd by default
#   bash run_genrestr.sh my_seq_list.txt  # pass a different file
#
# seq_ids.txt format (tab-separated):
#   seq_id       seq_type (display)      optional_custom_path
#   bind_043     Binder
#   nonb_046     False Positive
#   seq14_binder Binder                  /scratch/.../seq14_binder   ← custom path used as-is
#
# .top naming convention (lives alongside posre_protein.itp in WORKDIR):
#   {seq_type}_{ID}_dodecahedron_HMR.top  e.g. bind_043_dodecahedron_HMR.top
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

SEQ_LIST=${1:-seq_ids.txt}

# ── Path configuration ─────────────────────────────────────────────────────
# Adjust these if your directory layout differs.
BASE_DIR="/scratch/alpine/ivta1597/LCA_boltz_models"
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
# If EM/ lives one level deeper than ${BASE_DIR}/${dir_type}/${seq_id}/ (e.g.
# under an HMR/dodecahedron subfolder), set SUBDIR accordingly. Leave empty
# if EM/ sits directly under the seq_id directory.
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

# ── Insert the POSRES include block right after the protein #include line ──
# Assumes the .top file is in the same directory as posre_protein.itp
# (i.e. WORKDIR), since GROMACS resolves #include paths relative to the
# directory the .top file is run from.
add_posres_include() {
    local top_file="$1"

    if [[ ! -f "$top_file" ]]; then
        echo "  SKIP top edit (.top not found): $top_file"
        return 1
    fi

    # Idempotency: don't double-insert if already patched
    if grep -q "POSRES" "$top_file"; then
        echo "  .top already has POSRES block, skipping: $top_file"
        return 0
    fi

    # Find the first #include line that references the protein topology
    # (case-insensitive match on "protein" in the line)
    local line_num
    line_num=$(grep -in '^#include.*protein' "$top_file" | head -n1 | cut -d: -f1)

    if [[ -z "$line_num" ]]; then
        echo "  WARNING: no '#include ... protein' line found in: $top_file"
        return 1
    fi

    awk -v n="$line_num" '
        { print }
        NR == n {
            print "#ifdef POSRES"
            print "#include \"posre_protein.itp\""
            print "#endif"
        }
    ' "$top_file" > "${top_file}.tmp" && mv "${top_file}.tmp" "$top_file"

    echo "  Inserted POSRES include after line $line_num: $top_file"
    return 0
}

completed=0
skipped=0
failed=0
top_patched=0
top_failed=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    # Skip empty lines and comments
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    if [[ -n "$custom_path" ]]; then
        WORKDIR="$custom_path"
    else
        dir_type=$(get_dir_type "$seq_type")
        parent="${BASE_DIR}/${dir_type}"

        # Directories have a type suffix appended to seq_id, e.g.
        # bind_043 -> bind_043_binder, nonb_046 -> nonb_046_nb, and we don't
        # want to hardcode every suffix variant, so glob for it instead.
        shopt -s nullglob
        candidates=("${parent}/${seq_id}"*/)
        shopt -u nullglob

        if [[ ${#candidates[@]} -eq 0 ]]; then
            echo "SKIP (no directory matching ${seq_id}* under $parent): $seq_id"
            ((skipped++))
            continue
        elif [[ ${#candidates[@]} -gt 1 ]]; then
            echo "WARNING: multiple dirs match ${seq_id}* under $parent — using first: ${candidates[0]}"
        fi

        WORKDIR="${candidates[0]%/}"
    fi
    [[ -n "$SUBDIR" ]] && WORKDIR="${WORKDIR}/${SUBDIR}"

    if [[ ! -f "${WORKDIR}/EM/em.gro" ]]; then
        echo "SKIP (missing EM/em.gro): $seq_id   [$WORKDIR]"
        ((skipped++))
        continue
    fi

    echo "Running genrestr: $seq_id  [$seq_type -> $WORKDIR]"

    (
        cd "$WORKDIR" || exit 1
        echo 4 | "$GMX" genrestr -f EM/em.gro -o posre_protein.itp -fc 1000 1000 1000
    )

    if [[ $? -eq 0 ]]; then
        ((completed++))

        # .top naming pattern: {seq_type}_{ID}_dodecahedron_HMR.top
        # seq_id is already "{seq_type}_{ID}" (e.g. bind_043, nonb_046), so it
        # maps directly onto the filename.
        TOP_FILE="${WORKDIR}/${seq_id}_dodecahedron_HMR.top"

        if add_posres_include "$TOP_FILE"; then
            ((top_patched++))
        else
            ((top_failed++))
        fi
    else
        echo "FAILED: $seq_id"
        ((failed++))
    fi

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  genrestr completed : $completed"
echo "  genrestr skipped    : $skipped (missing EM/em.gro)"
echo "  genrestr failed     : $failed"
echo "  .top patched         : $top_patched"
echo "  .top patch failed    : $top_failed"
