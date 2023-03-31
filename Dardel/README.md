# Computations on Dardel
Some of the computations for the computer assisted part of the proof
were done on the computer cluster [Dardel
cluser](https://www.pdc.kth.se/hpc-services/computing-systems/about-the-dardel-hpc-system-1.1053338),
This directory contains the scripts used for those computations.

Access to Dardel was provided by the Swedish National Infrastructure
for Computing (SNIC) at the PDC Center for High Performance Computing,
KTH Royal Institute of Technology, partially funded by the Swedish
Research Council through grant agreement no. 2018-05973.

The Dardel cluster makes use of Slurm for job scheduling. The scripts
for running the computations are written in terms of Slurm scripts.
The scripts assumes that Julia is available and that the directory
where Julia stores configurations is the default `$HOME/.julia/`. From
a fresh clone of this repository an initial setup is done with the
following commands, executed from the root of this repository.

``` shell
# Install all required packages for Julia
julia --project=. --eval 'using Pkg; Pkg.instantiate()'

# Optionally precompile the code using PackageCompiler.jl. This
# reduces the total runtime since things only have to be compiled
# once. Note that if you later change the code you then need to
# compile it again.
julia --eval 'using Pkg; Pkg.add("PackageCompiler")' # Install PackageCompiler.jl
bash Dardel/scripts/compile.sh
```

The computations are run using the script `Dardel/scripts/run_proof.sh
start stop m`. Here `start`, `stop` and `m` are argument that indicate
which part of the interval to handle. The interval `[-0.9999,
-0.0012]` is split into 31 subintervals according to
`HighestCupsedWave.proof_interval_subdivisions`. The argument `start`
determines which of these 31 subintervals to start at and `stop`
determines which to stop at, if `stop` is not given then it computes
up to the last one. The optional argument `m` is given to
`HighestCuspedWave.proof_interval_subdivisions_mince` and allows for
doing the computations on a part of the subintervals, for the full
proof it should be omitted to cover the whole subintervals.

The number of nodes, tasks and CPUs per task to use is set in
`Dardel/scripts/run_proof.sh`. It is tuned for Dardel, which has 256
threads per node, and uses 64 tasks and 4 CPUs per task. By default
Julia makes use of all available tasks and CPUs, it can be manually
set using the environmental variables `HCW_WORKERS` and `HCW_THREADS`
respectively, see the function `create_workers` in
`Dardel/scripts/helper.jl`. In some cases it is beneficial to use a
`HCW_THREADS` value that is half of the number of CPUs per task, this
reduces the peak memory usage but only has a minor impact on
performance due to simultaneous multithreading.

We split the full computation into 6 jobs, each handling a subset of
the 31 subintervals. These jobs can be started with the following
commands, executed from the root of this repository. The data is
written to a timestamped subdirectory of `Dardel/data/proof/`. Note
that the data that the proof in the paper is based on is available at
`proofs/data/`.

``` shell
# The first three intervals have a higher memory usage and for that
# reason we use fewer threads to avoid using too much memory. This is
# slower than if we would have used all threads. However it is not
# that much slower due to SMT, only about 20-30% slower it seems.

# Intervals 1:1
# Takes around 12 hours in total.
HCW_THREADS=2 sbatch --time=14:00:00 -p main -o Dardel/logs/run_proof_1.o -e Dardel/logs/run_proof_1.e Dardel/scripts/run_proof.sh 1 1

# Intervals 2:2
# Takes around 9 hours in total.
HCW_THREADS=2 sbatch --time=11:00:00 -p main -o Dardel/logs/run_proof_2.o -e Dardel/logs/run_proof_2.e Dardel/scripts/run_proof.sh 2 2

# Intervals 3:3
# Takes around 6 hours in total.
HCW_THREADS=2 sbatch --time=8:00:00 -p main -o Dardel/logs/run_proof_3.o -e Dardel/logs/run_proof_3.e Dardel/scripts/run_proof.sh 3 3

# Intervals 4:10
# Takes around 10 hours in total.
sbatch --time=12:00:00 -p main -o Dardel/logs/run_proof_4.o -e Dardel/logs/run_proof_4.e Dardel/scripts/run_proof.sh 4 10

# Intervals 11:14
# Takes around 8 hours in total.
sbatch --time=10:00:00 -p main -o Dardel/logs/run_proof_5.o -e Dardel/logs/run_proof_5.e Dardel/scripts/run_proof.sh 11 14

# Intervals 15:end
# Takes around 7 hours in total.
sbatch --time=9:00:00 -p main -o Dardel/logs/run_proof_6.o -e Dardel/logs/run_proof_6.e Dardel/scripts/run_proof.sh 15
```

It is also possible to run the computations on a non-SLURM system
(running Linux). In this case the number of workers and threads needs
to be set explicitly using `HCW_WORKERS` and `HCW_THREADS`
respectively. For example the following code could be used on a system
with 256 threads.

``` shell
# 1:1
HCW_WORKERS=64 HCW_THREADS=2 sh Dardel/scripts/run_proof.sh 1 1

# 2:2
HCW_WORKERS=64 HCW_THREADS=2 sh Dardel/scripts/run_proof.sh 2 2

# 3:3
HCW_WORKERS=64 HCW_THREADS=2 sh Dardel/scripts/run_proof.sh 3 3

# 4:10
HCW_WORKERS=64 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 4 10

# 11:14
HCW_WORKERS=64 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 11 14

# 15:end
HCW_WORKERS=64 HCW_THREADS=4 sh Dardel/scripts/run_proof.sh 15
```
