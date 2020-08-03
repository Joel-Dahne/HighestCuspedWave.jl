# Mainly test that it actually runs.
# TODO: So far this only works for a very limited range of parameters,
# this needs to be improved.
@testset "CB" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 6)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 0, 0)
        C_B = HighestCuspedWave.CB(u0)
        C_B_low = HighestCuspedWave.CB_estimate(u0)
        @test isfinite(C_B) && C_B > 0
        @test !(C_B_low > C_B)
    end
end
