module HighestCuspedWave

using ArbExtras
using Arblib
using NLsolve
using OffsetArrays
using OrderedCollections
using SpecialFunctions

import Nemo

include("arb.jl")
include("special-functions/polylog.jl")
include("special-functions/clausenc.jl")
include("special-functions/clausens.jl")
include("types.jl")
include("evaluation.jl")

include("alpha.jl")
include("delta.jl")
include("T0.jl")
include("CB.jl")

include("FractionalKdV/FractionalKdVAnsatz.jl")
include("FractionalKdV/determination.jl")
include("FractionalKdV/coefficients.jl")
include("FractionalKdV/errorterms.jl")
include("FractionalKdV/evaluation.jl")
include("FractionalKdV/alpha.jl")
include("FractionalKdV/delta.jl")
include("FractionalKdV/T0.jl")
include("FractionalKdV/T01.jl")
include("FractionalKdV/T02.jl")
include("FractionalKdV/CB.jl")
include("FractionalKdV/proof.jl")

include("BurgersHilbert/BHAnsatz.jl")
include("BurgersHilbert/determination.jl")
include("BurgersHilbert/evaluation.jl")
include("BurgersHilbert/alpha.jl")
include("BurgersHilbert/delta.jl")
include("BurgersHilbert/T0.jl")
include("BurgersHilbert/T01.jl")
include("BurgersHilbert/T02.jl")
include("BurgersHilbert/CB.jl")

include("BHKdV/BHKdVAnsatz.jl")
include("BHKdV/evaluation.jl")
include("BHKdV/alpha.jl")
include("BHKdV/delta.jl")
include("BHKdV/T0.jl")
include("BHKdV/T0_asymptotic.jl")
include("BHKdV/T01.jl")
include("BHKdV/T02.jl")
include("BHKdV/CB.jl")
include("BHKdV/lemmas.jl")

include("KdVZero/KdVZeroAnsatz.jl")
include("KdVZero/expansion_with_remainder.jl")
include("KdVZero/expansion.jl")
include("KdVZero/evaluation.jl")
include("KdVZero/alpha0.jl")
include("KdVZero/delta.jl")
include("KdVZero/T0.jl")
include("KdVZero/CB.jl")

end # module
