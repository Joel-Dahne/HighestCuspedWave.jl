module HighestCuspedWave

using ArbExtras
using Arblib
using NLsolve
using OffsetArrays
using OrderedCollections
using SpecialFunctions
using DataFrames
using Tables
using CSV
import Distributed
import ProgressLogging

include("fmpz.jl")
include("arb.jl")
include("series.jl")
include("TaylorModel.jl")
include("special-functions/special-functions.jl")
include("special-functions/clausenc.jl")
include("special-functions/clausens.jl")
include("types.jl")
include("evaluation.jl")

include("bounds.jl")
include("estimates.jl")
include("T0.jl")
include("proof.jl")
include("data-handling.jl")

include("FractionalKdV/FractionalKdVAnsatz.jl")
include("FractionalKdV/determination.jl")
include("FractionalKdV/evaluation.jl")
include("FractionalKdV/n0.jl")
include("FractionalKdV/delta0.jl")
include("FractionalKdV/T0.jl")
include("FractionalKdV/T01.jl")
include("FractionalKdV/T02.jl")
include("FractionalKdV/T0_asymptotic_bhkdv.jl")
include("FractionalKdV/D0.jl")
include("FractionalKdV/proof.jl")

include("BurgersHilbert/BHAnsatz.jl")
include("BurgersHilbert/determination.jl")
include("BurgersHilbert/evaluation.jl")
include("BurgersHilbert/n0.jl")
include("BurgersHilbert/delta0.jl")
include("BurgersHilbert/T0.jl")
include("BurgersHilbert/T01.jl")
include("BurgersHilbert/T02.jl")
include("BurgersHilbert/D0.jl")

include("BHKdV/BHKdVAnsatz.jl")
include("BHKdV/evaluation.jl")
include("BHKdV/n0.jl")
include("BHKdV/delta0.jl")
include("BHKdV/T0.jl")
include("BHKdV/T0_asymptotic.jl")
include("BHKdV/T01.jl")
include("BHKdV/T02.jl")
include("BHKdV/D0.jl")

include("KdVZero/KdVZeroAnsatz.jl")
include("KdVZero/expansion.jl")
include("KdVZero/evaluation.jl")
include("KdVZero/n0.jl")
include("KdVZero/delta0.jl")
include("KdVZero/T0.jl")
include("KdVZero/D0.jl")
include("KdVZero/proof.jl")

include("equations.jl")
include("lemmas.jl")
include("lemmas_bh.jl") # Lemmas from the Burgers-Hilbert paper

using PrecompileTools

@compile_workload begin
    setprecision(Arb, 100) do
        αs = [
            HighestCuspedWave.proof_interval_subdivisions_mince(15, thin = true)[1]
            HighestCuspedWave.proof_interval_subdivisions_mince(25, thin = true)[1]
        ]

        HighestCuspedWave.prove(αs, only_estimate_D0 = true)

        # Close to α = -1 it takes to long to run the full computations
        u0 = FractionalKdVAnsatz(
            HighestCuspedWave.proof_interval_subdivisions_mince(10)[end],
        )
        F0(u0)(Arb(1))
        F0(u0)(ArbSeries((1, 1)))
        F0(u0, Asymptotic())(Arb(0.1))
        T0(u0)(Arb(1))
        T0(u0, Asymptotic())(Arb(0.1))
    end
end
end # module
