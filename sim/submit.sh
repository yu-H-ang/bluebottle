#!/bin/sh
#
#SBATCH --partition=gpu
#SBATCH --time=24:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --job-name=test
#SBATCH --open-mode=truncate

srun ./bluebottle
