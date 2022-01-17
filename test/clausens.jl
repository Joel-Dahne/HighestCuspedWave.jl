@testset "clausens" begin
    @testset "clausens(x, s)" begin
        # Check that the evaluation with polylog and zeta agree on (0, 2π)
        for s in [range(Arb(-4), Arb(4), length = 10); Arb.(-3:3)]
            for x in range(Arb(0), 2Arb(π), length = 100)[2:end-1]
                res1 = HighestCuspedWave._clausens_polylog(x, s)
                res2 = HighestCuspedWave._clausens_zeta(x, s)
                @test isfinite(res1)
                @test isfinite(res2)
                @test Arblib.overlaps(res1, res2)
            end
        end

        # Check that _clausen_zeta throws an error outside of the domain
        @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(-1), Arb(2.5))
        @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(7), Arb(2.5))

        # Check evaluation with integer s
        for s in Arb.(-4:4)
            @test Arblib.overlaps(
                HighestCuspedWave._clausens_polylog(one(s), s),
                HighestCuspedWave._clausens_zeta(one(s), s),
            )
            @test isfinite(clausenc(zero(s), s)) == (s > 1)
            @test isfinite(clausenc(one(s), s))
        end

        # Check evaluation on wide intervals
        for s in Arb[-2.5, 2.5]
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

        # Test s overlapping integers
        for x in range(Arb(0), 2Arb(π), length = 10)[2:end-1]
            for s in Arb.(-1:4)
                s_interval = Arblib.add_error!(copy(s), Arb(0.0001))
                y1 = clausens(x, s_interval)
                @test isfinite(y1)
                for ss in [s; range(getinterval(Arb, s_interval)..., length = 10)]
                    y2 = clausens(x, ss)
                    @test isfinite(y2)
                    @test Arblib.overlaps(y1, y2)
                end
            end
        end

        # x::ArbSeries
        # Very simple tests - just check that it runs and is finite
        @test isfinite(clausens(ArbSeries((2, 1, 0)), Arb(2)))
        @test isfinite(clausens(ArbSeries((2, 1, 0)), Arb(2.5)))

        # s::ArbSeries
        # Check that the evaluation with polylog and zeta agree on (0, 2π)
        for s in [range(Arb(-4), Arb(4), length = 6); Arb.(-3:3)]
            s_series = ArbSeries((s, 1), degree = 2)
            for x in range(Arb(0), 2Arb(π), length = 20)[2:end-1]
                res1 = HighestCuspedWave._clausens_polylog(x, s_series)
                res2 = HighestCuspedWave._clausens_zeta(x, s_series)
                @test isfinite(res1)
                @test isfinite(res2)
                @test Arblib.overlaps(res1, res2)
            end
        end
    end

    @testset "clausens(x, s, β)" begin
        # Check that the evaluation with polylog and zeta agree on (0, 2π)
        for β in [1, 2]
            for s in [range(Arb(-4), Arb(4), length = 10); Arb.(-3:3)]
                for x in range(Arb(0), 2Arb(π), length = 10)[2:end-1]
                    res1 = HighestCuspedWave._clausens_polylog(x, s, β)
                    res2 = HighestCuspedWave._clausens_zeta(x, s, β)
                    @test isfinite(res1)
                    @test isfinite(res2)
                    @test Arblib.overlaps(res1, res2)
                end
            end
        end

        # Check evaluation on wide intervals
        for β in [1, 2]
            for s in Arb[-2.5, 2.5]
                for lower in range(Arb(-10), 8, length = 10)
                    for upper in range(lower + 0.01, 10, length = 10)
                        interval = Arb((lower, upper))
                        y = clausens(interval, s, β)
                        @test isfinite(y) || !(s > 1)
                        for x in range(lower, upper, length = 10)
                            @test Arblib.overlaps(clausens(x, s, β), y)
                        end
                    end
                end
            end
        end

        # Check that _clausen_zeta throws an error outside of the domain
        @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(-1), Arb(2.5), 1)
        @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(7), Arb(2.5), 1)

        # Check evaluation with integer s
        for β in [1, 2]
            for s in Arb.(-4:4)
                @test isfinite(clausens(zero(s), s, β)) == (s > 1)
                @test isfinite(clausens(one(s), s, β))
            end
        end

        # Check evaluation with other types
        @test clausens(1.5, 2) == Float64(clausens(Arb(1.5), 2))
        @test clausens(2, 2) == Float64(clausens(Arb(2), 2))

        # TODO: Add tests for wide values of s

        # TODO: Add tests for x::ArbSeries
    end
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
