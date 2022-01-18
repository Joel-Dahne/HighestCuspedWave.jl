@testset "_reduce_argument_clausen" begin
    pi = Arb(π)
    twopi = let tmp = Arb(π)
        Arblib.mul_2exp!(tmp, tmp, 1)
    end


    xs1 = range(Arb(-20), 20, length = 100)
    xs2 = Arblib.add_error!.(deepcopy(xs1), Arb(0.1))
    xs3 = Arb(π) * collect(-100:100)
    xs4 = [Arb((i * Arb(π), (i + 1) * Arb(π))) for i = -100:100]
    xs5 = [
        Arblib.nonnegative_part!(Arb(), Arb((0, π))),
        Arblib.nonnegative_part!(Arb(), Arb((0, 2Arb(π)))),
        Arb("[-2.9387358858138341719519480363289240901e-38 +/- 8.65e-77]"), # Failed before
    ]

    for x in [xs1; xs2; xs3; xs4; xs5]
        y, haszero, haspi, has2pi = HighestCuspedWave._reduce_argument_clausen(x)

        @test Arblib.contains_int((x - y) / twopi)

        @test haszero == Arblib.contains_zero(y)
        @test haspi == (Arblib.overlaps(y, pi) || Arblib.overlaps(y, -pi))
        @test has2pi == Arblib.overlaps(y, twopi)
        @test !Arblib.overlaps(y, -twopi)

        @test !has2pi || (has2pi && haszero)

        # Check that radius is not too much higher, the value used
        # here is somewhat random
        @test radius(Arb, y) <= 2radius(Arb, x) + sqrt(eps(Arb))
    end

end

@testset "clausenc" begin
    @testset "clausenc(x, s)" begin
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
                    res1 =
                        (polylog(Acb(s), exp(im * z)) + polylog(Acb(s), exp(-im * z))) / 2
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
        for s in Arb.(-4:4)
            @test Arblib.overlaps(
                HighestCuspedWave._clausenc_polylog(one(s), s),
                HighestCuspedWave._clausenc_zeta(one(s), s),
            )
            @test isfinite(clausenc(zero(s), s)) == (s > 1)
            @test isfinite(clausenc(one(s), s))
        end

        # Check evaluation on wide intervals
        for s in Arb[-2.5, 2.5]
            for lower in range(Arb(-10), 8, length = 10)
                for upper in range(lower + 0.01, 10, length = 10)
                    interval = Arb((lower, upper))
                    y = clausenc(interval, s)
                    @test isfinite(y) || !(s > 1)
                    for x in range(lower, upper, length = 10)
                        @test Arblib.overlaps(clausenc(x, s), y)
                    end
                end
            end
        end

        # Check evaluation with other types
        @test clausenc(1.5, 2) == Float64(clausenc(Arb(1.5), 2))
        @test clausenc(2, 2) == Float64(clausenc(Arb(2), 2))

        # Test s overlapping integers
        for x in range(Arb(0), 2Arb(π), length = 10)[2:end-1]
            for s in Arb.(-1:4)
                s_interval = Arblib.add_error!(copy(s), Arb(0.0001))
                y1 = clausenc(x, s_interval)
                @test isfinite(y1)
                for ss in [s; range(getinterval(Arb, s_interval)..., length = 10)]
                    y2 = clausenc(x, ss)
                    @test isfinite(y2)
                    @test Arblib.overlaps(y1, y2)
                end
            end
        end

        # x::ArbSeries
        # Very simple tests - just check that it runs and is finite
        @test isfinite(clausenc(ArbSeries((2, 1, 0)), Arb(2)))
        @test isfinite(clausenc(ArbSeries((2, 1, 0)), Arb(2.5)))

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

        # Wide s around integer
        for x in range(Arb(0), 2Arb(π), length = 10)[2:end-1]
            for s0 in Arb.(-1:4)
                s0_interval = Arblib.add_error!(copy(s0), Arb(1e-10))
                s = ArbSeries((s0_interval, 1), degree = 2)
                y1 = clausenc(x, s)
                @test isfinite(y1)
                for ss in [s0; range(getinterval(Arb, s0_interval)..., length = 10)]
                    y2 = clausenc(x, ArbSeries((ss, 1), degree = 2))
                    @test isfinite(y2)
                    @test Arblib.overlaps(y1, y2)
                end
            end
        end
    end

    @testset "clausenc(x, s, β)" begin
        # Check that the evaluation with polylog and zeta agree on (0, 2π)
        for β in [1, 2]
            for s in [range(Arb(-4), Arb(4), length = 10); Arb.(-3:3)]
                for x in range(Arb(0), 2Arb(π), length = 10)[2:end-1]
                    res1 = HighestCuspedWave._clausenc_polylog(x, s, β)
                    res2 = HighestCuspedWave._clausenc_zeta(x, s, β)
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
                        y = clausenc(interval, s, β)
                        @test isfinite(y) || !(s > 1)
                        for x in range(lower, upper, length = 10)
                            @test Arblib.overlaps(clausenc(x, s, β), y)
                        end
                    end
                end
            end
        end

        # Check that _clausen_zeta throws an error outside of the domain
        @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(-1), Arb(2.5), 1)
        @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(7), Arb(2.5), 1)

        # Check evaluation with integer s
        for β in [1, 2]
            for s in Arb.(-4:4)
                @test isfinite(clausenc(zero(s), s, β)) == (s > 1)
                @test isfinite(clausenc(one(s), s, β))
            end
        end

        # Check evaluation with other types
        @test clausenc(1.5, 2) == Float64(clausenc(Arb(1.5), 2))
        @test clausenc(2, 2) == Float64(clausenc(Arb(2), 2))

        # TODO: Add tests for wide values of s

        # TODO: Add tests for x::ArbSeries
    end
end

@testset "clausenc_expansion" begin
    s = Arb(0.5)
    for M in [3, 6]
        for x in range(Arb(0), 2Arb(π), length = 100)[1:end-1]
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
            res1 = clausencmzeta(x, s)
            res2 = clausenc(x, s) - zeta(s)
            @test Arblib.overlaps(res1, res2)
        end
    end

    # We only really care about s > 1 in this case
    for s = 2:5
        @test isfinite(clausencmzeta(one(Arb), s))
        @test isfinite(clausencmzeta(one(Arb), Arb(s)))
    end

    # Check evaluation on wide intervals
    for s in (1.5, 2.5)
        for lower in range(Arb(-10), 8, length = 10)
            for upper in range(lower + 0.01, 10, length = 10)
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
