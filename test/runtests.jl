using Arblib, HighestCuspedWave, OffsetArrays, SpecialFunctions, Test

setprecision(Arb, 128) do
    @testset "HighestCuspedWave" begin
        @testset "Special functions" begin
            include("arb.jl")
            include("clausenc.jl")
            include("clausens.jl")
            include("special-functions.jl")
        end

        @testset "FractionalKdV" begin
            include("FractionalKdV/evaluation.jl")
            include("FractionalKdV/T0.jl")
        end

        @testset "BurgersHilbert" begin
            include("BurgersHilbert/evaluation.jl")
            include("BurgersHilbert/T0.jl")
        end

        @testset "KdVZero" begin
            include("KdVZero/expansion.jl")
            include("KdVZero/evaluation.jl")
            include("KdVZero/T0.jl")
        end
    end
end
