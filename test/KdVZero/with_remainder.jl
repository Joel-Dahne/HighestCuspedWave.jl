@testset "with_remainder" begin
    # Expansion of gamma(2α) / gamma(α) at a point close to zero.
    # Computed both by direct evaluation and by expanding at zero.

    α = ArbSeries((-1e-6, 1), degree = 5)

    # Direct evaluation
    res1 = gamma(2α) / gamma(α)

    # Expansion at 0 with remainder term
    w = ArbSeries(α, degree = Arblib.degree(α) + 1)
    w[0] = 0
    interval = union(α[0], Arb(0))

    # Expansion of rgamma(α)
    rgamma_expansion_1 =
        HighestCuspedWave.compose_with_remainder(HighestCuspedWave.rgamma, w, interval)

    # Expansion of rgamma(2α)
    rgamma_expansion_2 = HighestCuspedWave.compose_with_remainder(
        α -> HighestCuspedWave.rgamma(2α),
        w,
        interval,
    )

    # Expansion of rgamma(α) / rgamma(2α) = gamma(2α) / gamma(α)
    gammadivgamma_expansion = HighestCuspedWave.div_with_remainder(
        rgamma_expansion_1 << 1,
        rgamma_expansion_2 << 1,
        interval,
    )

    # Evaluation from expansion
    res2 = HighestCuspedWave.collapse_from_remainder(gammadivgamma_expansion, interval)

    # Check α * gamma(α)
    r1 = HighestCuspedWave.collapse_from_remainder(
        HighestCuspedWave.compose_with_remainder(inv, rgamma_expansion_1 << 1, interval),
        interval,
    )
    r2 = α * gamma(α)
    @test_broken Arblib.overlaps(r1, r2)

    r1 = HighestCuspedWave.collapse_from_remainder(
        HighestCuspedWave.compose_with_remainder(inv, rgamma_expansion_2 << 1, interval),
        interval,
    )
    r2 = α * gamma(2α)
    @test_broken Arblib.overlaps(r1, r2)

    @test_broken Arblib.overlaps(res1, res2)
end
