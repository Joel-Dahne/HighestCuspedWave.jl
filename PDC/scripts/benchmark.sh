#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
# Run with sbatch --partition=shared PDC/slurm-test.sh
#SBATCH --job-name benchmark
#SBATCH --account snic2022-22-821
#SBATCH --time 00:30:00

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=4

#SBATCH -o PDC/logs/benchmark.o
#SBATCH -e PDC/logs/benchmark.e

# Use precompiled sysimage if it exists
if test -f $HOME/.julia/sysimages/HighestCuspedWave.so ; then
    time julia --sysimage=$HOME/.julia/sysimages/HighestCuspedWave.so --project=. PDC/scripts/benchmark.jl "$@"
else
    time julia --project=. PDC/scripts/benchmark.jl "$@"
fi
