@testset "clausens" begin
    # Check that the evaluation with polylog and zeta agree on (0, 2π)
    for s in range(Arb(-4), Arb(4), length = 10)
        for x in range(Arb(0), 2Arb(π), length = 100)[2:end-1]
            res1 = HighestCuspedWave._clausens_polylog(x, s)
            res2 = HighestCuspedWave._clausens_zeta(x, s)
            @test Arblib.overlaps(res1, res2)
        end
    end

    # Check that _clausen_zeta throws an error outside of the domain
    @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(-1), Arb(2.5))
    @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(7), Arb(2.5))

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
                y = clausens(interval, s)
                for x in range(lower, upper, length = 10)
                    @test Arblib.overlaps(clausens(x, s), y)
                end
            end
        end
    end

    # Check evaluation with other types
    @test clausens(1.5, 2) == Float64(clausens(Arb(1.5), 2))
    @test clausens(2, 2) == Float64(clausens(Arb(2), 2))

    # TODO: Add tests for wide values of s

    # TODO: Add tests for x::ArbSeries

    # TODO: Add tests for s::ArbSeries

    # TODO: Add tests for clausenc(x, s, β)
end

@testset "clausens_expansion" begin
    s = Arb(0.5)
    for M in [3, 6]
        for x in range(Arb(0), 2Arb(π), length = 100)
            C, e, P, E = HighestCuspedWave.clausens_expansion(x, s, M)
            @test Arblib.overlaps(C * abs(x)^e + P(x) + E * x^(2M + 1), clausens(x, s))
        end
    end
end
