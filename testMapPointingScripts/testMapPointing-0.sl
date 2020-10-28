#!/usr/bin/zsh
#
#SBATCH --job-name=testMapPointing
#SBATCH --output=log/testMapPointing-part_%a.log
#SBATCH --partition=big
#SBATCH --nodelist=compute-0
#SBATCH --array=49-63
##SBATCH --cpus-per-task=4


##export OMP_NUM_THREADS=1

srun ./build/testMapPointing $SLURM_ARRAY_TASK_ID
#srun root -l -b -q testMapPointing.C\($SLURM_ARRAY_TASK_ID\)