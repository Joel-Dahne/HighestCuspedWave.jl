using Arblib, HighestCuspedWave, OffsetArrays, SpecialFunctions, Test

setprecision(Arb, 128) do
    @testset "HighestCuspedWave" begin
        @testset "Special functions" begin
            include("clausenc.jl")
            include("clausens.jl")
        end

        @testset "KdVZero" begin
            include("KdVZero/expansion.jl")
            include("KdVZero/evaluation.jl")
        end
    end
end
