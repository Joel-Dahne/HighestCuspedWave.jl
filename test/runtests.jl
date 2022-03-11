using Arblib, HighestCuspedWave, OffsetArrays, SpecialFunctions, Test

setprecision(Arb, 128) do
    @testset "HighestCuspedWave" begin
        @testset "Special functions" begin
            include("arb.jl")
            include("clausenc.jl")
            include("clausens.jl")
        end

        @testset "BurgersHilbert" begin
            include("KdVZero/evaluation.jl")
            include("KdVZero/T0.jl")
        end

        @testset "KdVZero" begin
            include("KdVZero/expansion.jl")
            include("KdVZero/evaluation.jl")
            include("KdVZero/T0.jl")
        end
    end
end
