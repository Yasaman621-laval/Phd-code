#!/bin/bash
#SBATCH --account=def-sduchesn
#def-thchlava
#def-sduchesn
#rrg-thchlava
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20g
#SBATCH --array=2
#5



module load StdEnv/2020  plink/1.9b_6.21-x86_64  r/4.3.1 


Rscript --vanilla step2_non_meta.r $SLURM_ARRAY_TASK_ID 



echo "Job finished with exit code $? at: `date`"
