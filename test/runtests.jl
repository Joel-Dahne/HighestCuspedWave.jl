using Arblib, HighestCuspedWave, OffsetArrays, SpecialFunctions, Test

setprecision(Arb, 128) do
    @testset "HighestCuspedWave" verbose = true begin
        include("series.jl")
        include("TaylorModel.jl")

        @testset "Special functions" verbose = true begin
            include("clausenc.jl")
            include("clausens.jl")
            include("special-functions.jl")
        end

        @testset "FractionalKdV" verbose = true begin
            include("FractionalKdV/evaluation.jl")
            include("FractionalKdV/T0.jl")
        end

        @testset "BurgersHilbert" verbose = true begin
            include("BurgersHilbert/evaluation.jl")
            include("BurgersHilbert/T0.jl")
        end

        @testset "KdVZero" verbose = true begin
            include("KdVZero/expansion.jl")
            include("KdVZero/evaluation.jl")
            include("KdVZero/T0.jl")
        end
    end
end
