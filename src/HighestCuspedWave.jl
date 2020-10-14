module HighestCuspedWave

using ArbTools
using NLsolve
using Nemo
using OffsetArrays
using Printf

import SpecialFunctions

include("arb.jl")
include("types.jl")
include("bounded_by.jl")

include("FractionalKdV/FractionalKdVAnsatz.jl")
include("FractionalKdV/determination.jl")
include("FractionalKdV/evaluation.jl")
include("FractionalKdV/alpha.jl")
include("FractionalKdV/delta.jl")
include("FractionalKdV/T0.jl")
include("FractionalKdV/CB.jl")
include("FractionalKdV/proof.jl")

end # module
