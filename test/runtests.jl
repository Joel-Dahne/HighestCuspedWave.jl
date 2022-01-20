using Arblib, HighestCuspedWave, OffsetArrays, SpecialFunctions, Test

setprecision(Arb, 128) do
    @testset "HighestCuspedWave" begin
        @testset "Special functions" begin
            include("clausenc.jl")
            include("clausens.jl")
        end

        @testset "KdVZero" begin
            include("KdVZero/with_remainder.jl")
            include("KdVZero/expansion.jl")
            include("KdVZero/evaluation.jl")
            include("KdVZero/T0.jl")
        end
    end
end
