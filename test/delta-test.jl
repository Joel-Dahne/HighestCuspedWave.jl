# Mainly test that it actually runs.
@testset "delta0" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 4)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 0, 0)
        δ0 = HighestCuspedWave.delta0(u0)
        δ0_low = HighestCuspedWave.delta0_estimate(u0)
        @test isfinite(δ0)
        @test !(δ0_low > δ0)
    end
end
