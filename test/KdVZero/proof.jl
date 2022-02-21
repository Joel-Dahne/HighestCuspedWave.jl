@testset "proof" begin
    verbose = false
    @testset "α0 = 0" begin
        ϵ = Arb(-1.2e-3)
        u0 = KdVZeroAnsatz(Arb((ϵ, 0)))

        α0 = alpha0(u0; verbose)

        @test 1.56 < α0 < 1.58

        A = delta0(u0, maxevals = 10000, atol = Arb(0.0025); verbose)[2]

        @test -9e-3 < abs_ubound(A) < 9e-3

        B = CB(u0; verbose)[1]

        @test 0.23 < B < 0.25

        @test A < B^2 / 4α0
    end

    @testset "α0 = -1 / 6" begin
        interval = Arb((-1 // 6, -0.1))
        n = 1000
        α = HighestCuspedWave.mince(interval, n)[1]
        u0 = KdVZeroAnsatz(α, midpoint(Arb, α))

        proof_data = HighestCuspedWave.prove(u0; verbose)

        @test proof_data.proved

        @test 1.348 < proof_data.α₀ < 1.350

        @test 5e-5 < proof_data.δ₀ < 7e-5

        @test 0.97 < proof_data.C_B_estimate < 0.99
    end
end
