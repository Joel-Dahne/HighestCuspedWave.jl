#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
#SBATCH --job-name run_proof
#SBATCH --account naiss2023-22-649
#SBATCH --mail-type=ALL
#SBATCH --time 04:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=4

#SBATCH -o Dardel/logs/run_proof.o
#SBATCH -e Dardel/logs/run_proof.e

time julia --project=. Dardel/scripts/run_proof.jl "$@"
