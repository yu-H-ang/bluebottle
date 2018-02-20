#!/bin/sh
#
#SBATCH --partition=gpu
#SBATCH --time=24:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --job-name=test
#SBATCH --output=bluebottle.out
#SBATCH --open-mode=truncate
sleep 5
./build-restart.sh $SLURM_JOB_NAME $SLURM_JOB_ID
sleep 5
sbatch restart.sh
sleep 5
srun ./bluebottle -r
