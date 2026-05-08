#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J ngsLD_sparse
#SBATCH -o logs/ngsLD_sparse.%A_%a.out
#SBATCH -e logs/ngsLD_sparse.%A_%a.err

set -euo pipefail

# =============================================================================
# STEP_LD_00b_run_ngsLD_chrsparse.slurm
#
# Run ngsLD on sparse-thinned BEAGLEs (~1000 markers per chromosome)
# with UNLIMITED max_kb_dist to capture chromosome-wide long-range LD.
#
# This produces the data needed for whole-chromosome LD heatmaps that
# show long-range haploblock structure (like Todesco et al. sunflower Fig. 4e).
#
# Requires: STEP_LD_00_make_sparse_beagle.sh to have been run first.
#
# Usage:
#   sbatch --array=0-27%7 STEP_LD_00b_run_ngsLD_chrsparse.slurm chrom.list
# =============================================================================

source ~/.bashrc
mamba activate assembly

#module load gcc-native/14.2
#module load PrgEnv-gnu/8.6.0
#module load cpeGNU/25.03
module load PrgEnv-gnu
module load gcc-native/14.2
module load GSL/2.8-cpeGNU-25.03


CHROM_LIST="${1:?ERROR: provide chromosome list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))
CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"
[[ -n "${CHROM}" ]] || { echo "[ERROR] No chromosome at line ${LINE}" >&2; exit 1; }

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
NGSLD="${BASE}/ngsTools/ngsLD/ngsLD"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5C_Inversion_LD"
FAI="${BASE}/00-samples/fClaHyb_Gar_LG.fa.fai"
BEAGLEDIR="${BASE}/inversion_localpca_v7/pruned81_beagle_by_chr"
OUTDIR="${BASE}/inversion_localpca_v7/08_ld_chrsparse"
THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "${OUTDIR}" logs

# ── Get chromosome length ──────────────────────────────────────────────
CHR_LEN=$(awk -v c="${CHROM}" '$1==c {print $2}' "${FAI}")
[[ -n "${CHR_LEN}" ]] || { echo "[ERROR] Chromosome ${CHROM} not in FAI" >&2; exit 1; }

# ── Step 1: Make sparse BEAGLE if not already done ─────────────────────
SPARSE_PREFIX="${OUTDIR}/${CHROM}.chrsparse"
if [[ ! -s "${SPARSE_PREFIX}.sparse.beagle.gz" ]]; then
  echo "[INFO] Building sparse BEAGLE for ${CHROM}..."
  bash "${SCRIPTDIR}/STEP_LD_00_make_sparse_beagle.sh" \
    "${BEAGLEDIR}/main_pruned81.${CHROM}.beagle.gz" \
    "${CHR_LEN}" \
    "${SPARSE_PREFIX}" \
    2000
fi

[[ -s "${SPARSE_PREFIX}.sparse.beagle.gz" ]] || { echo "[ERROR] Missing sparse BEAGLE" >&2; exit 1; }
[[ -s "${SPARSE_PREFIX}.sparse.pos" ]] || { echo "[ERROR] Missing sparse POS" >&2; exit 1; }

# ── Step 2: Count markers/individuals ──────────────────────────────────
N_SITES=$(zcat "${SPARSE_PREFIX}.sparse.beagle.gz" | awk 'END{print NR-1}')
N_IND=$(zcat "${SPARSE_PREFIX}.sparse.beagle.gz" | awk 'NR==1{n=(NF-3)/3} END{print n}')

echo "[INFO] CHROM=${CHROM}"
echo "[INFO] Sparse markers: ${N_SITES}"
echo "[INFO] Individuals: ${N_IND}"

# ── Step 3: Run ngsLD with UNLIMITED distance ──────────────────────────
# max_kb_dist=0 means unlimited in ngsLD
# max_snp_dist=0 means unlimited too
# With ~1000 markers this produces ~500K pairs — very manageable
"${NGSLD}" \
  --geno "${SPARSE_PREFIX}.sparse.beagle.gz" \
  --probs \
  --n_ind "${N_IND}" \
  --n_sites "${N_SITES}" \
  --pos "${SPARSE_PREFIX}.sparse.pos" \
  --max_kb_dist 0 \
  --max_snp_dist 0 \
  --rnd_sample 1.0 \
  --n_threads "${THREADS}" \
  --out "${SPARSE_PREFIX}.ld"

[[ -s "${SPARSE_PREFIX}.ld" ]] || { echo "[ERROR] Missing LD output" >&2; exit 1; }

# ── Step 4: Extract pairs_r2 table ────────────────────────────────────
tail -n +2 "${SPARSE_PREFIX}.ld" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' \
  > "${SPARSE_PREFIX}.pairs_r2.tsv"

N_PAIRS=$(wc -l < "${SPARSE_PREFIX}.pairs_r2.tsv")

echo "[DONE] ${CHROM}"
echo "  Sparse BEAGLE: ${SPARSE_PREFIX}.sparse.beagle.gz"
echo "  LD output:     ${SPARSE_PREFIX}.ld"
echo "  Pairs r2:      ${SPARSE_PREFIX}.pairs_r2.tsv (${N_PAIRS} pairs)"
