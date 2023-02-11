# Highest Cusped Waves for the fractional KdV equations

This repository contains the code for the computer assisted parts of
the proofs for the paper [Highest Cusped Waves for the fractional KdV
equations](). It also contains code related to the paper [Highest
Cusped Waves for the Burgers-Hilbert
equation](https://arxiv.org/abs/2205.00802), though the code is
updated compared to that paper.

The proof handles the half-open interval [-1, 0) and is split into
four parts
1. -1
2. (-1, -0.9997)
3. [-0.9997, -0.0012]
4. (-0.0012, 0)

If you are interested in seeing the results of the computations
without running anything yourself you can look at the html-file in the
`proofs` directory. You will likely have to download the files and
open them locally in your browser. It contains the following files
1. `bh.jl.html` - contains the proof for -1.
2. `bhkdv.jl.html` - contains the proof for (-1, -0.9997).
3. `data-analysis.jl.html` - contains an analysis of the data
   generated for the proof for [-0.9997, -0.0012]. See further down
   for how the data is generated.
4. `kdvzero.jl.html` - contains the proof for (-0.0012, 0).

## Reproducing the proof

The html-files mentioned above are generated from
[Pluto](https://github.com/fonsp/Pluto.jl) notebooks. You can
reproduce the proof by running the notebooks yourself as described
below. The exception is the interval [-0.9997, -0.0012], for which the
data needs to be computed separately, see the next section for more
details.

The proofs were run with Julia version 1.8.5 but should likely work
with later versions as well. This repository contains the same
`Manifest.toml` file as was used when running the proofs, this allows
us to install exactly the same versions of the packages. To do this
start by downloading this repository, enter the directory and start
Julia. Now run the code

``` julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will likely take some time the first time you run it since it has
to download all the packages. You can see if it seems to work by
running `Pkg.test()`. This should give you some output related to all
the installed packages and then ending in

```

```

To run the notebooks run

``` julia
using Pluto
Pluto.run()
```

which should open a Pluto window in your browser. Now you can open the
notebooks inside the `proofs` directory through this and it should run
the proof.

### Reproducing the proof for [-0.9997, -0.0012]
The proof for [-0.9997, -0.0012] takes significantly longer time than
the other parts. For that reason the necessary data is computed
separately and then analysed in the notebook `proof/data-analysis.jl`.

The data used for the proof can be found in the directory
`proof/data/`. The computations were done on the [Dardel]() for which
access was granted through [SNIC](). The full computations took around
13 000 core hours.

Dardel uses SLURM and the full computations can be run with the script
`PDC/scripts/run_proof_all.sh`. This makes use of the scripts
`PDC/scripts/run_proof.sh` and `PDC/scripts/run_proof.jl`. The number
of tasks and threads to use is set in `PDC/scripts/run_proof.sh`, it
is tuned for Dardel (which has 256 threads per node) and uses 64 tasks
and 4 cpus per task. When running on a different configuration these
numbers would need to be updated.

To run the computations on a non-SLURM system (running Linux) it
should be possible to run the part argument to all the cases in
`PDC/scripts/run_proof_all` by themselves, the number of workers and
threads to use can be set with the environmental variables
`HCW_WORKERS` and `HCW_THREADS` respectively. For example the
following code could be used on a system with 256 threads, in which
case the computations should take roughly 50 hours. **TODO:** Consider
having a script handling this.

``` sh
# The first three intervals have a higher memory usage and for that
# reason we use fewer threads to avoid using too much memory. This is
# slower than if we would have used all threads. However it is not
# that much slower due to SMT, only about 20-30% slower it seems.

# 1:1
HCW_WORKERS=64 HCW_THREADS=2 sh PDC/scripts/run_proof.sh 1 1

# 2:2
HCW_WORKERS=64 HCW_THREADS=2 sh PDC/scripts/run_proof.sh 2 2

# 3:3
HCW_WORKERS=64 HCW_THREADS=2 sh PDC/scripts/run_proof.sh 3 3

# 4:10
HCW_WORKERS=64 HCW_THREADS=4 sh PDC/scripts/run_proof.sh 4 10

# 11:14
HCW_WORKERS=64 HCW_THREADS=4 sh PDC/scripts/run_proof.sh 11 14

# 14:end
HCW_WORKERS=64 HCW_THREADS=4 sh PDC/scripts/run_proof.sh 15
```

The data is written to a timestamped subdirectory of `PDC/data/proof/`.

# Notes about the implementation
The code uses [Arblib.jl](https://github.com/kalmarek/Arblib.jl),
which is an interface to [Arb](https://www.arblib.org/), for rigorous
numerics. Some of the methods used for computing bounds are from
[ArbExtras.jl](https://github.com/Joel-Dahne/ArbExtras.jl), most
notably `ArbExtras.enclose_maximum` which is used for computing most
of the bounds.

The code includes implementations of some special functions, these are
found in `src/special-functions/`. Most notably are the
implementations of the Clausen functions, these are found in
`src/special-functions/clausenc.jl` and
`src/special-functions/clausens.jl`. The other functions are mostly
wrappers for Arb implementations, sometimes with slight modifications.

Depending on the interval the approximation u₀ is represented by
different type.
1. For -1 it is represented by the type `BHAnsatz`, for which the
   corresponding methods are found in `src/BurgersHilbert/`.
2. For (-1, -0.9997) it is represented by the type `BHKdVAnsatz`, for
   which the corresponding methods are found in `src/BHKdV/`.
3. For [-0.9997, -0.0012] it is either represented by the type
   `FractionalKdVAnsatz` or the type `KdVZero`. The methods
   corresponding to `FractionalKdV` are found in `src/FractionalKdV/`
   and the ones corresponding to `KdVZero` in `src/KdVZero/`
4. For (-0.0012, 0) it is represented by the type `KdVZero`, for which
   the corresponding methods are found in `src/KdVZero/`.

Some other notable files are
- `src/TaylorModel.jl` - contains code for computing with Taylor
  models
- `src/proof.jl` - contains code for bisecting the interval [-0.9997,
  -0.0012] into smaller subintervals and the starting point for
  proving the inequality on for each subinterval.
- `src/data-handling.jl` - contains code for handling the data
  produced for the interval [-0.9997, -0.0012].
