#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_ld_v4
#SBATCH -o logs/inv_ld_v4.%j.out
#SBATCH -e logs/inv_ld_v4.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5C_Inversion_LD"

CHR_LEN="${INVDIR}/08_ld_pruned81/chr_lengths.tsv"
if [[ ! -s "${CHR_LEN}" ]]; then
  printf "chrom\tchr_len\n" > "${CHR_LEN}"
  awk -F'\t' '{print $1"\t"$2}' "${BASE}/00-samples/fClaHyb_Gar_LG.fa.fai" >> "${CHR_LEN}"
fi

mkdir -p "${INVDIR}/ld_cache_v4" "${INVDIR}/ld_plots_v4" logs

# ─── STEP 1: Run STEP_LD_00b for chromosome-sparse (if not done) ────────
# sbatch --array=0-27%7 ${SCRIPTDIR}/STEP_LD_00b_run_ngsLD_chrsparse.slurm ${INVDIR}/chrom.list

# ─── STEP 2: Run STEP_LD_00c for per-candidate invzoom (if not done) ────
# seq 1 65 > ${INVDIR}/candidate_ids.list
# sbatch --array=0-64 ${SCRIPTDIR}/STEP_LD_00c_run_ngsLD_invzoom.slurm ${INVDIR}/candidate_ids.list

# ─── STEP 3: Build caches + plot ────────────────────────────────────────
python3 "${SCRIPTDIR}/run_ld_candidate_suite_v4.py" \
  --candidate-tsv "${INVDIR}/supplementary_tables/supp_candidate_master_table.tsv" \
  --chrsparse-dir "${INVDIR}/08_ld_chrsparse" \
  --invzoom-dir "${INVDIR}/08_ld_invzoom" \
  --chr-len-tsv "${CHR_LEN}" \
  --cache-dir "${INVDIR}/ld_cache_v4" \
  --outdir "${INVDIR}/ld_plots_v4" \
  --target-bins 250 \
  --pad-frac 0.05 \
  --transform power --power 0.5 \
  --cmap apex \
  --dpi 400 \
  --script-dir "${SCRIPTDIR}"

echo "[DONE] LD suite v4"
