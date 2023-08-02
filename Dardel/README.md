# Computations on Dardel
Some of the computations for the computer assisted part of the proof
were done on the
[Dardel<](https://www.pdc.kth.se/hpc-services/computing-systems/about-the-dardel-hpc-system-1.1053338)
computer cluster, This directory contains the scripts used for those
computations.

Access to Dardel was provided by the National Academic Infrastructure
for Supercomputing in Sweden (NAISS) and the Swedish National
Infrastructure for Computing (SNIC) at PDC partially funded by the
Swedish Research Council through grant agreements no. 2022-06725 and
no. 2018-05973.

The Dardel cluster makes use of Slurm for job scheduling. The scripts
for running the computations are written in terms of Slurm scripts.
The scripts assumes that Julia is available and that the directory
where Julia stores configurations is the default `$HOME/.julia/`. From
a fresh clone of this repository an initial setup is done with the
following commands, executed from the root of this repository.

``` shell
# Install all required packages for Julia - this takes some time
julia --project=. --eval 'using Pkg; Pkg.instantiate()'
```

The computations are run using the script `Dardel/scripts/run_proof.sh
start stop m`. Here `start`, `stop` and `m` are argument that indicate
which part of the interval to handle. The interval `[-0.9999,
-0.0012]` is split into 32 subintervals according to
`HighestCupsedWave.proof_interval_subdivisions`. The argument `start`
determines which of these 32 subintervals to start at and `stop`
determines which to stop at, if `stop` is not given then it computes
up to the last one. The optional argument `m` is given to
`HighestCuspedWave.proof_interval_subdivisions_mince` and allows for
doing the computations on a part of the subintervals, for the full
proof it should be omitted to cover the whole subintervals.

The number of nodes, tasks and CPUs per task to use is set in
`Dardel/scripts/run_proof.sh`. It is tuned for Dardel, which has 256
threads per node, and uses 64 tasks and 4 CPUs per task in most cases.
By default Julia makes use of all available tasks and CPUs, it can be
manually set using the environmental variables `HCW_WORKERS` and
`HCW_THREADS` respectively, see the function `create_workers` in
`Dardel/scripts/helper.jl`. Near -1 the computations use a lot of
memory and to reduce the peak memory usage some parameters are
modified. In this case it uses a `HCW_THREADS` value that is half of
the number of CPUs per task, this only has a minor impact on
performance due to simultaneous multithreading. Very close to -1 it
also uses 32 tasks with 8 CPUs per task instead.

We split the full computation into 7 jobs, each handling a subset of
the 32 subintervals. These jobs can be started with the following
commands, executed from the root of this repository. The data is
written to a timestamped subdirectory of `Dardel/data/proof/`. Note
that the data that the proof in the paper is based on is available at
`proofs/data/Dardel/`.

``` shell
# Intervals 1:1
# Takes around 7 hours in total.
HCW_THREADS=4 sbatch --ntasks=32 --cpus-per-task=8 --time=9:00:00 -p main --job-name=run_proof_1 -o Dardel/logs/run_proof_1.o -e Dardel/logs/run_proof_1.e Dardel/scripts/run_proof.sh 1 1

# Intervals 2:2
# Takes around 5 hours in total.
HCW_THREADS=4 sbatch --ntasks=32 --cpus-per-task=8 --time=7:00:00 -p main --job-name=run_proof_2 -o Dardel/logs/run_proof_2.o -e Dardel/logs/run_proof_2.e Dardel/scripts/run_proof.sh 2 2

# Intervals 3:3
# Takes around 8 hours in total.
HCW_THREADS=4 sbatch --ntasks=32 --cpus-per-task=8 --time=10:00:00 -p main --job-name=run_proof_3 -o Dardel/logs/run_proof_3.o -e Dardel/logs/run_proof_3.e Dardel/scripts/run_proof.sh 3 3

# Intervals 4:4
# Takes around 6 hours in total.
HCW_THREADS=2 sbatch --time=8:00:00 -p main --job-name=run_proof_4 -o Dardel/logs/run_proof_4.o -e Dardel/logs/run_proof_4.e Dardel/scripts/run_proof.sh 4 4

# Intervals 5:11
# Takes around 9 hours in total.
sbatch --time=11:00:00 -p main --job-name=run_proof_5 -o Dardel/logs/run_proof_5.o -e Dardel/logs/run_proof_5.e Dardel/scripts/run_proof.sh 5 11

# Intervals 12:15
# Takes around 8 hours in total.
sbatch --time=10:00:00 -p main --job-name=run_proof_6 -o Dardel/logs/run_proof_6.o -e Dardel/logs/run_proof_6.e Dardel/scripts/run_proof.sh 12 15

# Intervals 16:end
# Takes around 7 hours in total.
sbatch --time=9:00:00 -p main --job-name=run_proof_7 -o Dardel/logs/run_proof_7.o -e Dardel/logs/run_proof_7.e Dardel/scripts/run_proof.sh 16
```

It is also possible to run the computations on a non-SLURM system
(running Linux). In this case the number of workers and threads needs
to be set explicitly using `HCW_WORKERS` and `HCW_THREADS`
respectively. For example the following code could be used on a system
with 256 threads.

``` shell
# 1:1
HCW_WORKERS=32 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 1 1

# 2:2
HCW_WORKERS=32 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 2 2

# 3:3
HCW_WORKERS=32 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 3 3

# 4:4
HCW_WORKERS=64 HCW_THREADS=2 sh Dardel/scripts/run_proof.sh 4 4

# 5:11
HCW_WORKERS=64 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 5 11

# 12:15
HCW_WORKERS=64 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 12 15

# 16:end
HCW_WORKERS=64 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 16
```
