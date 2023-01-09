#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
# Run with sbatch --partition=shared PDC/slurm-test.sh
#SBATCH --job-name run_proof
#SBATCH --account snic2022-22-821
#SBATCH --mail-type=ALL
#SBATCH --time 04:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=4

#SBATCH -o PDC/logs/run_proof.o
#SBATCH -e PDC/logs/run_proof.e

# Use precompiled sysimage if it exists
if test -f $HOME/.julia/sysimages/HighestCuspedWave.so ; then
    time julia --sysimage=$HOME/.julia/sysimages/HighestCuspedWave.so --project=. PDC/scripts/run_proof.jl "$@"
else
    time julia --project=. PDC/scripts/run_proof.jl "$@"
fi
