# Highest Cusped Waves for the fractional KdV equations

This repository contains the code for the computer assisted parts of
the proofs for the paper [Highest Cusped Waves for the fractional KdV
equations](). It also contains code related to the paper [Highest
Cusped Waves for the Burgers-Hilbert
equation](https://arxiv.org/abs/2205.00802), though the code is
updated compared to that paper (for the original version see
[BurgersHilbertWave.jl](https://github.com/Joel-Dahne/BurgersHilbertWave.jl)).

The proof handles the half-open interval [-1, 0) and is split into
four parts
1. -1
2. (-1, -0.9999)
3. [-0.9999, -0.0012]
4. (-0.0012, 0)

If you are interested in seeing the results of the computations
without running anything yourself you can look at the html-file in the
`proofs` directory. You will likely have to download the files and
open them locally in your browser. It contains the following files
1. `bh.jl.html` - contains the proof for -1.
2. `bhkdv.jl.html` - contains the proof for (-1, -0.9999).
3. `kdv.jl.html` - illustrates the idea of the proof for [-0.9999,
   -0.0012].
4. `kdv-data-analysis.jl.html` - contains an analysis of the data
   generated for the proof for [-0.9999, -0.0012]. See further down
   for how the data is generated.
5. `kdvzero.jl.html` - contains the proof for (-0.0012, 0).

## Reproducing the proof

The html-files mentioned above are generated from
[Pluto](https://github.com/fonsp/Pluto.jl) notebooks. You can
reproduce the proof by running the notebooks yourself as described
below. The exception is the interval [-0.9999, -0.0012], for which the
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

### Reproducing the proof for [-0.9999, -0.0012]
The proof for [-0.9999, -0.0012] takes significantly longer time than
the other parts. For that reason the necessary data is computed
separately and then analysed in the notebook
`proof/kdv-data-analysis.jl`.

The computations for generating the data used for the proof were done
on the computer cluster
[Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/about-the-dardel-hpc-system-1.1053338).
See `Dardel/README.md` for more details about how the computations
were run.

The data used for the proof can be found in the directory
`proof/data/`. The full computations took around 13 000 core hours.

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
2. For (-1, -0.9999) it is represented by the type `BHKdVAnsatz`, for
   which the corresponding methods are found in `src/BHKdV/`.
3. For [-0.9999, -0.0012] it is either represented by the type
   `FractionalKdVAnsatz` or the type `KdVZero`. The methods
   corresponding to `FractionalKdV` are found in `src/FractionalKdV/`
   and the ones corresponding to `KdVZero` in `src/KdVZero/`
4. For (-0.0012, 0) it is represented by the type `KdVZero`, for which
   the corresponding methods are found in `src/KdVZero/`.

Some other notable files are
- `src/TaylorModel.jl` - contains code for computing with Taylor
  models
- `src/proof.jl` - contains code for bisecting the interval [-0.9999,
  -0.0012] into smaller subintervals and the starting point for
  proving the inequality on for each subinterval.
- `src/data-handling.jl` - contains code for handling the data
  produced for the interval [-0.9999, -0.0012].

## Naming conventions
Many of the variables in the paper have an index $\alpha$, in the code
this is most of the time replaced with a `0`. We have the following
correspondence between the name in the paper and in the code.

| Paper                          | Code             |
|--------------------------------|------------------|
| $n_\alpha$                     | `n0` or `n₀`     |
| $\delta_\alpha$                | `delta0` or `δ₀` |
| $D_\alpha$                     | `D0` or `D₀`     |
| $u_\alpha$                     | `u0`             |
| $F_\alpha$                     | `F0`             |
| $\mathcal{T}_\alpha$           | `T0`             |
| $\mathcal{H}^\alpha[u_\alpha]$ | `H(u0)`          |

## Asymptotic and non-asymptotic evaluation
Most functions can be evaluated in two ways. One using direct ball
arithmetic and one optimized for `x` close to zero that uses
expansions at `x = 0`. In general the method of evaluation is set by
giving either `Ball()` or `Asymptotic()` as an argument to the
function, with `Ball()` being the default.

```
julia> using HighestCuspedWave, Arblib

julia> setprecision(Arb, 100)
100

julia> u0 = FractionalKdVAnsatz(Arb(-0.6))
FractionalKdVAnsatz{Arb} N₀ = 8, N₁ = 16
α = [-0.599999999999999977795539507497 +/- 1.31e-31], p = [0.799999999999999988897769753748 +/- 4.35e-31]

julia> x = Arb(0.1)
[0.100000000000000005551115123126 +/- 2.18e-31]

julia> u0(x, Ball()) # Direct evaluation with ball arithmetic
[0.318537784850320154969 +/- 1.08e-22]

julia> u0(x) # The defaults is ball arithmetic
[0.318537784850320154969 +/- 1.08e-22]

julia> u0(x, Asymptotic()) # Evaluation using the asymptotic expansion
[0.31853778485 +/- 6.48e-12]

julia> F0(u0)(x)
[-1.99183775964e-8 +/- 5.40e-20]

julia> F0(u0, Asymptotic())(x)
[-1.99e-8 +/- 4.54e-11]

julia> T0(u0)(x)
[0.87 +/- 6.64e-3]

julia> T0(u0, Asymptotic())(x)
[0.8675738 +/- 5.30e-8]

```
