#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
# Run with sbatch --partition=shared PDC/slurm-test.sh
#SBATCH --job-name slurm-test
#SBATCH --account snic2022-22-821
#SBATCH --time 00:05:00

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2

#SBATCH -e PDC/logs/slurm_test.e
#SBATCH -o PDC/logs/slurm_test.o

julia --project=. PDC/scripts/slurm-test.jl
