#!/bin/bash

#SBATCH --job-name=strong_landau
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=70G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=krasny0
#SBATCH --partition=gpu
#SBATCH  --gpus=1

module load python ffmpeg/6.0.0
python FS_int.py -rd ./ --gpu -r
