# Mainly test that it actually runs.
# TODO: So far this only works for a very limited range of parameters,
# this needs to be improved.
@testset "alpha0" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 6)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 0, 0)
        α0 = HighestCuspedWave.alpha0(u0)
        α0_low = HighestCuspedWave.alpha0_estimate(u0)
        @test isfinite(α0) && α0 > 0
        @test !(α0_low > α0)
    end
end
