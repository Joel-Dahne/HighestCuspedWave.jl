@testset "clausenc" begin
    # Check that the evaluation with polylog and zeta agree on (0, 2π)
    for s in [range(Arb(-4), Arb(4), length = 10); Arb.(-3:3)]
        for x in range(Arb(0), 2Arb(π), length = 100)[2:end-1]
            res1 = HighestCuspedWave._clausenc_polylog(x, s)
            res2 = HighestCuspedWave._clausenc_zeta(x, s)
            @test isfinite(res1)
            @test isfinite(res2)
            @test Arblib.overlaps(res1, res2)
        end
    end

    # Check that the evaluation with polylog and zeta agree for
    # complex values
    for s in range(Arb(-4), Arb(4), length = 10)[2:end-1]
        for x in range(Arb(0), 2Arb(π), length = 10)[2:end-1]
            for y in range(Arb(-10), Arb(10), length = 10)
                z = Acb(x, y)
                res1 = (polylog(Acb(s), exp(im * z)) + polylog(Acb(s), exp(-im * z))) / 2
                res2 = HighestCuspedWave._clausenc_zeta(z, s)
                @test isfinite(res1)
                @test isfinite(res2)
                @test Arblib.overlaps(res1, res2)
            end
        end
    end

    # Check that _clausen_zeta throws an error outside of the domain
    @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(-1), Arb(2.5))
    @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(7), Arb(2.5))

    # Check evaluation with integer s
    for s = -4:4
        @test isfinite(clausenc(one(Arb), s))
        @test isfinite(clausenc(one(Arb), Arb(s)))
    end

    # Check evaluation on wide intervals
    for s in (-2.5, 2.5)
        for lower in range(Arb(-10), 8, length = 10)
            for upper in range(lower + 1, 10, length = 10)
                interval = Arb((lower, upper))
                y = clausenc(interval, s)
                for x in range(lower, upper, length = 10)
                    @test Arblib.overlaps(clausenc(x, s), y)
                end
            end
        end
    end

    # Check evaluation with other types
    @test clausenc(1.5, 2) == Float64(clausenc(Arb(1.5), 2))
    @test clausenc(2, 2) == Float64(clausenc(Arb(2), 2))

    # TODO: Add tests for wide values of s

    # TODO: Add tests for x::ArbSeries

    # s::ArbSeries
    # Check that the evaluation with polylog and zeta agree on (0, 2π)
    for s in [range(Arb(-4), Arb(4), length = 6); Arb.(-3:3)]
        s_series = ArbSeries((s, 1), degree = 2)
        for x in range(Arb(0), 2Arb(π), length = 20)[2:end-1]
            res1 = HighestCuspedWave._clausenc_polylog(x, s_series)
            res2 = HighestCuspedWave._clausenc_zeta(x, s_series)
            @test isfinite(res1)
            @test isfinite(res2)
            @test Arblib.overlaps(res1, res2)
        end
    end

    # TODO: Add tests for clausenc(x, s, β)
end

@testset "clausenc_expansion" begin
    s = Arb(0.5)
    for M in [3, 6]
        for x in range(Arb(0), 2Arb(π), length = 100)
            C, e, P, E = HighestCuspedWave.clausenc_expansion(x, s, M)
            @test Arblib.overlaps(C * abs(x)^e + P(x) + E * x^2M, clausenc(x, s))
        end
    end
end

@testset "clausencmzeta" begin
    # In this case we only really care about s > 1

    # Check that the evaluation with zeta and naive implementation
    # agree on (0, 2π)
    for s in range(Arb(1), Arb(4), length = 5)
        for x in range(Arb(0), 2Arb(π), length = 100)[2:end-1]
            res1 = HighestCuspedWave._clausencmzeta_zeta(x, s)
            res2 = clausenc(x, s) - zeta(s)
            @test Arblib.overlaps(res1, res2)
        end
    end

    # Check that _clausen_zeta throws an error outside of the domain
    @test_throws DomainError HighestCuspedWave._clausencmzeta_zeta(Arb(-1), Arb(2.5))
    @test_throws DomainError HighestCuspedWave._clausencmzeta_zeta(Arb(7), Arb(2.5))

    # We only really care about s > 1 in this case
    for s = 2:5
        @test isfinite(clausencmzeta(one(Arb), s))
        @test isfinite(clausencmzeta(one(Arb), Arb(s)))
    end

    # Check evaluation on wide intervals
    for s in (1.5, 2.5)
        for lower in range(Arb(-10), 8, length = 10)
            for upper in range(lower + 1, 10, length = 10)
                interval = Arb((lower, upper))
                y = clausencmzeta(interval, s)
                for x in range(lower, upper, length = 10)
                    @test Arblib.overlaps(clausencmzeta(x, s), y)
                end
            end
        end
    end

    # Check that promotion seems to work properaly
    @test isequal(clausencmzeta(Arb(1.5), 2), clausencmzeta(Arb(1.5), Arb(2)))

    # Test with other types
    @test clausencmzeta(1.5, 2) ≈ Float64(clausencmzeta(Arb(1.5), 2))
    @test clausencmzeta(2, 2) ≈ Float64(clausencmzeta(Arb(2), 2))

    # TODO: Add tests for wide values of s

    # TODO: Add tests for x::ArbSeries

    # TODO: Add tests for s::ArbSeries

    # TODO: Add tests for clausenc(x, s, β)
    for s = 2:3
        for x in range(Arb(0), 2Arb(π), length = 100)[2:end-1]
            res1 = clausencmzeta(x, s, 1)
            res2 = clausenc(x, s, 1) - HighestCuspedWave.dzeta(Arb(s))
            @test Arblib.overlaps(res1, res2)
        end
    end
end
