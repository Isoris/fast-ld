#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J ngsLD_zoom
#SBATCH -o logs/ngsLD_zoom.%A_%a.out
#SBATCH -e logs/ngsLD_zoom.%A_%a.err

set -euo pipefail

# =============================================================================
# STEP_LD_00c_run_ngsLD_invzoom.slurm
#
# For each candidate inversion, extract a REGIONAL thinned BEAGLE covering
# the candidate + flanks, then run ngsLD with max_kb_dist = 0 (unlimited)
# on that regional subset. This gives FULL long-range LD coverage within
# the zoom region, which is what's needed for the zoom triangle/square plots.
#
# The key insight: the old localdense product used max_kb_dist=500 which
# only gives ~500 kb near-diagonal LD. A 12 Mb inversion needs 12+ Mb
# distance coverage. We fix this by:
#   1. Extract regional markers from the full BEAGLE
#   2. Thin to ~2000 markers within the zoom region
#   3. Run ngsLD with unlimited distance on this small set
#   4. ~2000 markers → ~2M pairs → very manageable
#
# Usage:
#   sbatch --array=0-64 STEP_LD_00c_run_ngsLD_invzoom.slurm candidate_ids.list
# =============================================================================

source ~/.bashrc
mamba activate assembly

#module load gcc-native/14.2
#module load PrgEnv-gnu/8.6.0
#module load cpeGNU/25.03
module load PrgEnv-gnu
module load gcc-native/14.2
module load GSL/2.8-cpeGNU-25.03


CANDIDATE_LIST="${1:?ERROR: provide candidate_ids.list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))
CID="$(sed -n "${LINE}p" "${CANDIDATE_LIST}")"
[[ -n "${CID}" ]] || { echo "[ERROR] No candidate at line ${LINE}" >&2; exit 1; }

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
NGSLD="${BASE}/ngsTools/ngsLD/ngsLD"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5C_Inversion_LD"
FAI="${BASE}/00-samples/fClaHyb_Gar_LG.fa.fai"
BEAGLEDIR="${BASE}/inversion_localpca_v7/pruned81_beagle_by_chr"
CAND_TABLE="${BASE}/inversion_localpca_v7/supplementary_tables/supp_candidate_master_table.tsv"
OUTDIR="${BASE}/inversion_localpca_v7/08_ld_invzoom"
THREADS="${SLURM_CPUS_PER_TASK}"
PAD_FRAC=0.05
TARGET_MARKERS=2000

mkdir -p "${OUTDIR}" logs

# ── Get candidate info ──────────────────────────────────────────────────
read CHROM INV_START INV_END < <(
  awk -F'\t' -v cid="${CID}" '
  NR==1 { for(i=1;i<=NF;i++) { if($i=="candidate_id") ci=i; if($i=="chrom") cc=i; if($i=="start_bp") cs=i; if($i=="end_bp") ce=i } next }
  $ci==cid { print $cc, $cs, $ce; exit }
  ' "${CAND_TABLE}"
)

[[ -n "${CHROM}" ]] || { echo "[ERROR] Candidate ${CID} not found in table" >&2; exit 1; }

CHR_LEN=$(awk -v c="${CHROM}" '$1==c {print $2}' "${FAI}")
[[ -n "${CHR_LEN}" ]] || { echo "[ERROR] ${CHROM} not in FAI" >&2; exit 1; }

# Compute zoom region
INV_LEN=$((INV_END - INV_START))
PAD=$(awk "BEGIN{v=int(${CHR_LEN}*${PAD_FRAC}); if(v<${INV_LEN}/2) v=int(${INV_LEN}/2); print v}")
ZOOM_START=$((INV_START - PAD))
(( ZOOM_START < 1 )) && ZOOM_START=1
ZOOM_END=$((INV_END + PAD))
(( ZOOM_END > CHR_LEN )) && ZOOM_END=${CHR_LEN}
ZOOM_SPAN=$((ZOOM_END - ZOOM_START))

echo "[INFO] Candidate ${CID}: ${CHROM}:${INV_START}-${INV_END} (${INV_LEN} bp)"
echo "[INFO] Zoom: ${ZOOM_START}-${ZOOM_END} (${ZOOM_SPAN} bp)"

PREFIX="${OUTDIR}/${CHROM}.cand${CID}.invzoom"

# Skip if already done
if [[ -s "${PREFIX}.pairs_r2.tsv" ]]; then
  echo "[INFO] Output exists, skipping: ${PREFIX}.pairs_r2.tsv"
  exit 0
fi

# ── Step 1: Extract + thin regional BEAGLE ─────────────────────────────
IN_BEAGLE="${BEAGLEDIR}/main_pruned81.${CHROM}.beagle.gz"
[[ -s "${IN_BEAGLE}" ]] || { echo "[ERROR] Missing BEAGLE: ${IN_BEAGLE}" >&2; exit 1; }

# Compute thinning spacing
SPACING=$(( ZOOM_SPAN / TARGET_MARKERS ))
# Round to nice value
if   (( SPACING < 5000 )); then SPACING=5000
elif (( SPACING < 25000 )); then SPACING=$(( (SPACING / 5000) * 5000 )); (( SPACING < 5000 )) && SPACING=5000
elif (( SPACING < 100000 )); then SPACING=$(( (SPACING / 10000) * 10000 ))
else SPACING=$(( (SPACING / 25000) * 25000 ))
fi

echo "[INFO] Thinning spacing: ${SPACING} bp (target ~${TARGET_MARKERS} markers)"

# Extract markers in zoom region, keep one per spacing bin
zcat "${IN_BEAGLE}" | awk -v ZS="${ZOOM_START}" -v ZE="${ZOOM_END}" -v SP="${SPACING}" '
BEGIN{OFS="\t"}
NR==1 { print; next }
{
  id=$1
  gsub(/:/, "_", id)
  n=split(id, parts, "_")
  pos=parts[n]+0
  if (pos < ZS || pos > ZE) next
  bin=int(pos/SP)
  if (!(bin in seen)) {
    seen[bin]=1
    print
  }
}' | gzip > "${PREFIX}.regional.beagle.gz"

# Build pos file
zcat "${PREFIX}.regional.beagle.gz" | awk '
BEGIN{OFS="\t"}
NR==1 { next }
{
  id=$1
  gsub(/:/, "_", id)
  n=split(id, parts, "_")
  chr=""
  for(i=1;i<n;i++) chr = (chr=="" ? parts[i] : chr"_"parts[i])
  pos=parts[n]+0
  print chr, pos
}' > "${PREFIX}.regional.pos"

N_SITES=$(zcat "${PREFIX}.regional.beagle.gz" | awk 'END{print NR-1}')
N_IND=$(zcat "${PREFIX}.regional.beagle.gz" | awk 'NR==1{n=(NF-3)/3} END{print n}')

echo "[INFO] Regional markers: ${N_SITES}, Individuals: ${N_IND}"

if (( N_SITES < 10 )); then
  echo "[WARN] Too few markers (${N_SITES}), skipping candidate ${CID}"
  exit 0
fi

# ── Step 2: Run ngsLD with UNLIMITED distance ──────────────────────────
"${NGSLD}" \
  --geno "${PREFIX}.regional.beagle.gz" \
  --probs \
  --n_ind "${N_IND}" \
  --n_sites "${N_SITES}" \
  --pos "${PREFIX}.regional.pos" \
  --max_kb_dist 0 \
  --max_snp_dist 0 \
  --rnd_sample 1.0 \
  --n_threads "${THREADS}" \
  --out "${PREFIX}.ld"

[[ -s "${PREFIX}.ld" ]] || { echo "[ERROR] Missing LD output" >&2; exit 1; }

# ── Step 3: Extract pairs_r2 ──────────────────────────────────────────
tail -n +2 "${PREFIX}.ld" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' \
  > "${PREFIX}.pairs_r2.tsv"

N_PAIRS=$(wc -l < "${PREFIX}.pairs_r2.tsv")

# Write stats
cat > "${PREFIX}.stats" <<EOF
candidate_id=${CID}
chrom=${CHROM}
inv_start=${INV_START}
inv_end=${INV_END}
zoom_start=${ZOOM_START}
zoom_end=${ZOOM_END}
spacing_bp=${SPACING}
markers_retained=${N_SITES}
individuals=${N_IND}
pairs=${N_PAIRS}
EOF

echo "[DONE] Candidate ${CID}: ${N_SITES} markers, ${N_PAIRS} pairs"
echo "  ${PREFIX}.pairs_r2.tsv"
