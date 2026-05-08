#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# STEP_LD_00_make_sparse_beagle.sh  (v4)
#
# Thin a per-chromosome BEAGLE to ~TARGET_N evenly-spaced markers.
# Default target raised to 2000 for better whole-chromosome LD density.
# =============================================================================

IN_BEAGLE="${1:?ERROR: provide input beagle.gz}"
CHR_LEN="${2:?ERROR: provide chromosome length in bp}"
OUTPREFIX="${3:?ERROR: provide output prefix}"
TARGET_N="${4:-2000}"

[[ -s "${IN_BEAGLE}" ]] || { echo "[ERROR] Missing: ${IN_BEAGLE}" >&2; exit 1; }

SPACING=$(( CHR_LEN / TARGET_N ))

# Round to nice value
if   (( SPACING < 25000 )); then
  SPACING=$(( (SPACING / 5000 + 1) * 5000 ))
elif (( SPACING < 100000 )); then
  SPACING=$(( (SPACING / 10000 + 1) * 10000 ))
else
  SPACING=$(( (SPACING / 25000 + 1) * 25000 ))
fi

echo "[INFO] Chr length: ${CHR_LEN} bp, Target: ${TARGET_N}, Spacing: ${SPACING} bp"

zcat "${IN_BEAGLE}" | awk -v SP="${SPACING}" '
BEGIN{OFS="\t"}
NR==1 { print; next }
{
  id=$1; gsub(/:/, "_", id)
  n=split(id, parts, "_"); pos=parts[n]+0
  bin=int(pos/SP)
  if (!(bin in seen)) { seen[bin]=1; kept++; print }
}
END { print kept " markers retained" > "/dev/stderr" }
' | gzip > "${OUTPREFIX}.sparse.beagle.gz"

gzip -t "${OUTPREFIX}.sparse.beagle.gz"

zcat "${OUTPREFIX}.sparse.beagle.gz" | awk '
BEGIN{OFS="\t"} NR==1{next}
{ id=$1; gsub(/:/, "_", id); n=split(id,p,"_")
  chr=""; for(i=1;i<n;i++) chr=(chr==""?p[i]:chr"_"p[i])
  print chr, p[n]+0
}' > "${OUTPREFIX}.sparse.pos"

N_RETAINED=$(zcat "${OUTPREFIX}.sparse.beagle.gz" | awk 'END{print NR-1}')
N_ORIG=$(zcat "${IN_BEAGLE}" | awk 'END{print NR-1}')

cat > "${OUTPREFIX}.sparse.stats" <<EOF
input_beagle=${IN_BEAGLE}
chr_len_bp=${CHR_LEN}
target_n=${TARGET_N}
spacing_bp=${SPACING}
markers_original=${N_ORIG}
markers_retained=${N_RETAINED}
thinning_ratio=$(awk "BEGIN{printf \"%.4f\", ${N_RETAINED}/${N_ORIG}}")
mean_spacing_bp=$(awk "BEGIN{printf \"%.0f\", ${CHR_LEN}/${N_RETAINED}}")
EOF

echo "[DONE] ${OUTPREFIX}.sparse.beagle.gz (${N_RETAINED}/${N_ORIG} markers, spacing ~${SPACING} bp)"
