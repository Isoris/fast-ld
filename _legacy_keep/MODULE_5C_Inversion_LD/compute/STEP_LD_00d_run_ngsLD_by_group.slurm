#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J ngsLD_grp
#SBATCH -o logs/ngsLD_grp.%A_%a.out
#SBATCH -e logs/ngsLD_grp.%A_%a.err

set -euo pipefail

# =============================================================================
# STEP_LD_00d_run_ngsLD_by_group.slurm
#
# For each candidate, use STEP17c-exported sample lists to subset the BEAGLE
# into per-group BEAGLEs (FULL_A, FULL_B, HALF), thin markers, and run ngsLD
# with unlimited distance for each group separately.
#
# This enables group-specific LD heatmaps and comparative LD plots.
#
# Usage:
#   sbatch --array=0-64 STEP_LD_00d_run_ngsLD_by_group.slurm candidate_ids.list
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
FAI="${BASE}/00-samples/fClaHyb_Gar_LG.fa.fai"
BEAGLEDIR="${BASE}/inversion_localpca_v7/pruned81_beagle_by_chr"
SAMPLE_ORDER="${BEAGLEDIR}/bamlist226.sample_order.txt"
CAND_TABLE="${BASE}/inversion_localpca_v7/supplementary_tables/supp_candidate_master_table.tsv"
GROUP_DIR="${BASE}/inversion_localpca_v7/13_scoring_table/contrast_groups"
OUTDIR="${BASE}/inversion_localpca_v7/08_ld_by_group"
THREADS="${SLURM_CPUS_PER_TASK}"
PAD_FRAC=0.05
TARGET_MARKERS=2000
MIN_SAMPLES=5

mkdir -p "${OUTDIR}" logs

# ── Get candidate info ──────────────────────────────────────────────────
read CHROM INV_START INV_END < <(
  awk -F'\t' -v cid="${CID}" '
  NR==1 { for(i=1;i<=NF;i++) { if($i=="candidate_id") ci=i; if($i=="chrom") cc=i; if($i=="start_bp") cs=i; if($i=="end_bp") ce=i } next }
  $ci==cid { print $cc, $cs, $ce; exit }
  ' "${CAND_TABLE}"
)
[[ -n "${CHROM}" ]] || { echo "[ERROR] Candidate ${CID} not found" >&2; exit 1; }

CHR_LEN=$(awk -v c="${CHROM}" '$1==c {print $2}' "${FAI}")
INV_LEN=$((INV_END - INV_START))
PAD=$(awk "BEGIN{v=int(${CHR_LEN}*${PAD_FRAC}); if(v<${INV_LEN}/2) v=int(${INV_LEN}/2); print v}")
ZOOM_START=$((INV_START - PAD)); (( ZOOM_START < 1 )) && ZOOM_START=1
ZOOM_END=$((INV_END + PAD)); (( ZOOM_END > CHR_LEN )) && ZOOM_END=${CHR_LEN}
ZOOM_SPAN=$((ZOOM_END - ZOOM_START))

SPACING=$(( ZOOM_SPAN / TARGET_MARKERS ))
if   (( SPACING < 5000 )); then SPACING=5000
elif (( SPACING < 25000 )); then SPACING=$(( (SPACING / 5000) * 5000 )); (( SPACING < 5000 )) && SPACING=5000
elif (( SPACING < 100000 )); then SPACING=$(( (SPACING / 10000) * 10000 ))
else SPACING=$(( (SPACING / 25000) * 25000 )); fi

IN_BEAGLE="${BEAGLEDIR}/main_pruned81.${CHROM}.beagle.gz"
[[ -s "${IN_BEAGLE}" ]] || { echo "[ERROR] Missing BEAGLE: ${IN_BEAGLE}" >&2; exit 1; }

echo "[INFO] Candidate ${CID}: ${CHROM}:${INV_START}-${INV_END}, zoom ${ZOOM_START}-${ZOOM_END}, spacing ${SPACING}"

# ── Read pruned81 sample order from BEAGLE header ──────────────────────
BEAGLE_SAMPLES="${OUTDIR}/.tmp_beagle_samples_${CHROM}.txt"
zcat "${IN_BEAGLE}" | head -1 | tr '\t' '\n' | tail -n +4 | awk 'NR%3==1' > "${BEAGLE_SAMPLES}"
N_BEAGLE_SAMPLES=$(wc -l < "${BEAGLE_SAMPLES}")

# ── Process each group ─────────────────────────────────────────────────
for GROUP in FULL_A FULL_B HALF; do
  SAMPLE_LIST="${GROUP_DIR}/candidate_${CID}.${GROUP}.samples.txt"
  [[ -s "${SAMPLE_LIST}" ]] || { echo "[INFO] No samples for ${GROUP}, skipping"; continue; }

  N_GROUP=$(wc -l < "${SAMPLE_LIST}")
  if (( N_GROUP < MIN_SAMPLES )); then
    echo "[INFO] ${GROUP}: only ${N_GROUP} samples (min=${MIN_SAMPLES}), skipping"
    continue
  fi

  PREFIX="${OUTDIR}/${CHROM}.cand${CID}.${GROUP}"

  # Skip if already done
  [[ -s "${PREFIX}.pairs_r2.tsv" ]] && { echo "[INFO] ${PREFIX}.pairs_r2.tsv exists, skipping"; continue; }

  # Build column indices for this group's samples in the BEAGLE
  # BEAGLE columns: marker allele1 allele2 [sample1 sample1 sample1] [sample2 ...]
  # For sample at 0-based index i in the BEAGLE sample order:
  #   columns = 3*i+1, 3*i+2, 3*i+3 (1-based)
  # Plus columns 1,2,3 for marker/alleles

  CUTCOLS="1,2,3"
  N_KEPT=0
  while read -r sid; do
    # Find this sample's 1-based line number in the beagle sample order
    IDX_LINE=$(grep -n -x "${sid}" "${BEAGLE_SAMPLES}" | head -1 | cut -d: -f1)
    if [[ -z "${IDX_LINE}" ]]; then continue; fi
    I0=$((IDX_LINE - 1))  # 0-based
    C1=$((3 * I0 + 4)); C2=$((C1 + 1)); C3=$((C2 + 1))
    CUTCOLS="${CUTCOLS},${C1},${C2},${C3}"
    N_KEPT=$((N_KEPT + 1))
  done < "${SAMPLE_LIST}"

  if (( N_KEPT < MIN_SAMPLES )); then
    echo "[INFO] ${GROUP}: only ${N_KEPT} matched in BEAGLE (min=${MIN_SAMPLES}), skipping"
    continue
  fi

  echo "[INFO] ${GROUP}: ${N_KEPT} samples, subsetting + thinning..."

  # Subset BEAGLE to group samples + thin to zoom region with spacing
  zcat "${IN_BEAGLE}" | cut -f${CUTCOLS} | awk -v ZS="${ZOOM_START}" -v ZE="${ZOOM_END}" -v SP="${SPACING}" '
  BEGIN{OFS="\t"}
  NR==1 { print; next }
  {
    id=$1; gsub(/:/, "_", id)
    n=split(id, parts, "_"); pos=parts[n]+0
    if (pos < ZS || pos > ZE) next
    bin=int(pos/SP)
    if (!(bin in seen)) { seen[bin]=1; print }
  }' | gzip > "${PREFIX}.regional.beagle.gz"

  # POS file
  zcat "${PREFIX}.regional.beagle.gz" | awk '
  BEGIN{OFS="\t"} NR==1{next}
  { id=$1; gsub(/:/, "_", id); n=split(id,p,"_")
    chr=""; for(i=1;i<n;i++) chr=(chr==""?p[i]:chr"_"p[i])
    print chr, p[n]+0
  }' > "${PREFIX}.regional.pos"

  N_SITES=$(zcat "${PREFIX}.regional.beagle.gz" | awk 'END{print NR-1}')

  if (( N_SITES < 10 )); then
    echo "[WARN] ${GROUP}: too few markers (${N_SITES}), skipping"
    continue
  fi

  echo "[INFO] ${GROUP}: ${N_SITES} markers, ${N_KEPT} individuals, running ngsLD..."

  "${NGSLD}" \
    --geno "${PREFIX}.regional.beagle.gz" \
    --probs \
    --n_ind "${N_KEPT}" \
    --n_sites "${N_SITES}" \
    --pos "${PREFIX}.regional.pos" \
    --max_kb_dist 0 \
    --max_snp_dist 0 \
    --rnd_sample 1.0 \
    --n_threads "${THREADS}" \
    --out "${PREFIX}.ld"

  [[ -s "${PREFIX}.ld" ]] || { echo "[ERROR] ${GROUP}: ngsLD failed" >&2; continue; }

  tail -n +2 "${PREFIX}.ld" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' \
    > "${PREFIX}.pairs_r2.tsv"

  N_PAIRS=$(wc -l < "${PREFIX}.pairs_r2.tsv")
  echo "[DONE] ${GROUP}: ${N_SITES} markers, ${N_PAIRS} pairs → ${PREFIX}.pairs_r2.tsv"
done

rm -f "${BEAGLE_SAMPLES}"
echo "[DONE] Candidate ${CID} group LD complete"
