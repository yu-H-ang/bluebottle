#!/bin/sh
# submit.sh
#
#SBATCH --partition=devel
#SBATCH --gres=gpu:1
#SBATCH --job-name=bb-Euler

srun ./bluebottle
