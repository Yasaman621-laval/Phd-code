#!/bin/bash
#SBATCH --account=def-thchlava
#def-thchlava
#def-sduchesn
#rrg-thchlava
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-66
#66



module load StdEnv/2020  r/4.3.1



Rscript create_ld_block_prccsx.r $SLURM_ARRAY_TASK_ID 


echo "Job finished with exit code $? at: `date`"




