using Test
using HighestCuspedWave
using Nemo
using ArbTools

@testset "HighestCuspedWave" begin
    include("arb-test.jl")
    include("clausian-test.jl")
    include("evaluation-test.jl")

    include("T0-test.jl")

    include("alpha-test.jl")
    include("delta-test.jl")
    #include("CB-test.jl")
end
