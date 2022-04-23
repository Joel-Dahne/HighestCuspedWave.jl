@testset "proof" begin
    verbose = true

    @time u0 = BHAnsatz{Arb}()

    @time n₀ = n0_bound(u0; verbose)

    @test 0.53681 < n₀ < 0.53682

    @time δ₀ = delta0_bound(u0; verbose)

    @test 0.0008496 < δ₀ < 0.0008498

    @time D₀ = D0_bound(u0; verbose)

    @test 0.944 < D₀ < 0.946

    @test δ₀ < (1 - D₀)^2 / 4n₀
end
