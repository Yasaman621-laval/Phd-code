#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-sduchesn
#SBATCH --job-name=mixer_snps_EAS
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-22
# 22 chromosomes


module load StdEnv/2023 python/3.11
source ~/scratch/ENmixer311/bin/activate

set -euo pipefail

# get chromosome from SLURM
CHR=${SLURM_ARRAY_TASK_ID}

# --- Paths ---
MIXER_ROOT="/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer"
MIXER_PY="${MIXER_ROOT}/precimed/mixer.py"
LIBBG="${MIXER_ROOT}/src/build/lib/libbgmg.so"

SC_DIR="/lustre06/project/6005709/yatah3/simulation/SimuGenotype/sim_hsq0.5_rho0.4_train10000-10000-80000"
REFDIR="${SC_DIR}/Chr${CHR}/PlinkFormat"
LD_DIR="${SC_DIR}/ld_mixer"

echo "[$(date)] Generating EAS .snps for chr${CHR}"

srun python3 "${MIXER_PY}" snps \
  --lib "${LIBBG}" \
  --bim-file "${REFDIR}/EASchr${CHR}bedformat.bim" \
  --ld-file  "${LD_DIR}/chr${CHR}.EAS.ld" \
  --chr2use  "${CHR}" \
  --out      "${LD_DIR}/EAS_chr${CHR}.prune_maf0p05_rand2M_r2p8.snps" \
  --maf 0.05 \
  --r2 0.8 
 

echo "[$(date)] Done EAS .snps for chr${CHR}"


#If two SNPs have: r² > 0.8 ? They give almost the same information.? Keep only one. r² <= 0.8 ? Keep both; they are not too redundant.