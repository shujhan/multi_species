#!/bin/bash

#SBATCH --job-name=strong_landau_amr_movie
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --account=krasny0
#SBATCH --partition=standard

module load python ffmpeg/6.0.0
python FS_int.py -rd ./ -ps
