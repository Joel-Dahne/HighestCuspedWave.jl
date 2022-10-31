#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
# Run with sbatch --partition=shared PDC/slurm-test.sh
#SBATCH --job-name tune
#SBATCH --account snic2022-22-821
#SBATCH --time 00:60:00

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=4

#SBATCH -o PDC/logs/tune.o
#SBATCH -e PDC/logs/tune.e

julia --project=. PDC/scripts/tune.jl "$@"
