#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-thchlava
#SBATCH --job-name=mixer_ld_AFR
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

# Load Python and activate environment
module load StdEnv/2023 python/3.11
source ~/scratch/ENmixer311/bin/activate


# Get chromosome index from SLURM array
CHR=${SLURM_ARRAY_TASK_ID}

# Set directories
SC_DIR="/lustre06/project/6005709/yatah3/simulation/SimuGenotype/sim_hsq0.5_rho0.4_train10000-10000-80000/"
REFDIR="${SC_DIR}/Chr${CHR}/PlinkFormat/"
LD_DIR="${SC_DIR}/ld_mixer"

# Mixer paths
MIXER_PY="/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer/precimed/mixer.py"
LIBBG="/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer/src/build/lib/libbgmg.so"

# Create output and log directories
mkdir -p "${LD_DIR}" logs

# Run Mixer LD command
python3 "$MIXER_PY" ld \
  --bfile "${REFDIR}/AFRchr${CHR}bedformat" \
  --out   "${LD_DIR}/chr${CHR}.AFR.ld" \
  --lib   "$LIBBG" \
  --r2min 0.05 \
  --ldscore-r2min 0.0001 \
  --ld-window-kb 10000

