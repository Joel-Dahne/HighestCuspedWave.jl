module HighestCuspedWave

using ArbTools
using NLsolve
using Nemo
using OffsetArrays
using Printf
using OrderedCollections
using Arblib
using ArbExtras

import SpecialFunctions

import ArbTools: getinterval
import Nemo: midpoint, radius

include("arb.jl")
include("clausian.jl")
include("special-functions/polylog.jl")
include("special-functions/clausenc.jl")
include("special-functions/clausens.jl")
include("types.jl")
include("bounded_by.jl")
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

end # module
