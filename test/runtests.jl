using Arblib, HighestCuspedWave, SpecialFunctions, Test

setprecision(Arb, 128) do
    @testset "HighestCuspedWave" begin
        @testset "Special functions" begin
            include("clausenc.jl")
            include("clausens.jl")
        end
    end
end
