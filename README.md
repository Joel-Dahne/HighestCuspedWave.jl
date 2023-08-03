# Highest Cusped Waves for the fractional KdV equations

This repository contains the code for the computer assisted parts of
the proofs for the paper [Highest Cusped Waves for the Fractional KdV
Equations](). It also contains code related to the paper [Highest
Cusped Waves for the Burgers-Hilbert
Equation](https://doi.org/10.1007/s00205-023-01904-6), though the code
is updated compared to that paper (for the original version see
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

The proofs were run with Julia version 1.9.2. This repository contains
the same `Manifest.toml` file as was used when running the proofs,
this allows us to install exactly the same versions of the packages.
The computation of the approximation is very sensitive to rounding for
some values of $\alpha$, for this reason using a different Julia
version or other versions for packages might produce slightly
different bounds or fail to compute bounds.

To reproduce the results downloading this repository, enter the
directory and start Julia. Now run the code

``` julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will likely take some time the first time you run it since it has
to download and compile all the packages. You can see if it seems to
work by running `Pkg.test()`. This should give you some output related
to all the installed packages and then ending in

```
     Testing Running tests...
Test Summary:                |   Pass   Total     Time
HighestCuspedWave            | 311905  311905  5m43.5s
  fx_div_x                   | 188757  188757     3.0s
  TaylorModel                |  17291   17291    13.7s
  Special functions          |  97617   97617    32.8s
    _reduce_argument_clausen |   4235    4235     0.3s
    clausenc                 |  41701   41701     7.7s
    clausenc_expansion       |   3132    3132     4.2s
    clausencmzeta            |   8149    8149     1.2s
    clausens                 |  35701   35701     6.1s
    clausens_expansion       |   3000    3000     3.8s
    special-functions        |   1699    1699     0.4s
  FractionalKdV              |   1408    1408    32.4s
    evaluation               |   1248    1248    14.0s
    T0                       |    160     160    18.1s
  BurgersHilbert             |    730     730  1m20.7s
    evaluation               |    699     699    54.5s
    T0                       |     31      31    24.2s
  KdVZero                    |   5986    5986    56.9s
    expansion                |    796     796     1.6s
    evaluation               |   4888    4888    22.2s
    T0                       |    298     298    20.8s
    proof                    |      4       4    11.7s
  BHKdV                      |    116     116  2m03.4s
    evaluation               |     85      85  1m02.8s
    T0                       |     31      31    59.5s
     Testing HighestCuspedWave tests passed
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
See [`Dardel/README.md`](Dardel/README.md) for more details about how
the computations were run.

The data used for the proof can be found in the directory
`proof/data/`. The full computation took around 12 000 core hours.

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
different types.
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

```julia
julia> using HighestCuspedWave, Arblib

julia> setprecision(Arb, 100)
100

julia> u0 = FractionalKdVAnsatz(Arb(-0.6))
FractionalKdVAnsatz{Arb} N₀ = 6, N₁ = 16
α = [-0.599999999999999977795539507497 +/- 1.31e-31], p = [0.799999999999999988897769753748 +/- 4.35e-31]

julia> x = Arb(0.1)
[0.100000000000000005551115123126 +/- 2.18e-31]

julia> u0(x, Ball()) # Direct evaluation with ball arithmetic
[0.318537048012475859459769 +/- 6.40e-25]

julia> u0(x) # The defaults is ball arithmetic
[0.318537048012475859459769 +/- 6.40e-25]

julia> u0(x, Asymptotic()) # Evaluation using the asymptotic expansion
[0.318537048 +/- 2.25e-10]

julia> F0(u0)(x)
[-2.58138884562876e-7 +/- 3.59e-22]

julia> F0(u0, Asymptotic())(x)
[-2.6e-7 +/- 3.33e-9]

julia> T0(u0)(x)
[0.87 +/- 6.63e-3]

julia> T0(u0, Asymptotic())(x)
[0.8675758 +/- 4.67e-8]

```
