using Test
using HighestCuspedWave
using Nemo
using ArbTools
using Arblib

setprecision(Arb, 64) do
    @testset "HighestCuspedWave" begin
        include("arb-test.jl")
        include("clausian-test.jl")
        include("evaluation-test.jl")

        #include("T0-test.jl") # Not finished and takes a long time

        include("alpha-test.jl")
        include("delta-test.jl")
        #include("CB-test.jl")
    end
end
