@testset "proof" begin
    verbose = true
    @testset "α0 = 0" begin
        setprecision(Arb, 100) do
            ϵ = Arb(-1.2e-3)
            u0 = KdVZeroAnsatz(Arb((ϵ, 0)))

            n₀ = n0_bound(u0; verbose)

            @test 1.56 < n₀ < 1.58

            A = delta0_bound(u0, maxevals = 10000, atol = Arb(0.0025); verbose).p[2]

            @test -9e-3 < abs_ubound(A) < 9e-3

            B = D0_bound(u0; verbose).p[1]

            @test 0.23 < B < 0.25

            @test A < B^2 / 4n₀
        end
    end

    @testset "α0 = -1 / 6" begin
        setprecision(Arb, 100) do
            interval = Arb((-1 // 6, -0.1))
            n = 1000
            α = HighestCuspedWave.mince(interval, n)[1]
            u0 = KdVZeroAnsatz(α, midpoint(Arb, α))

            proof_data = HighestCuspedWave.prove(u0; verbose)

            @test proof_data.proved

            @test 1.348 < proof_data.n₀ < 1.350

            @test 5e-5 < proof_data.δ₀ < 7e-5

            @test 0.97 < proof_data.D₀_estimate < 0.99
        end
    end
end
