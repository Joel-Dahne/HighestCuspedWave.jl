@testset "clausens" begin
    xs = range(Arb(0), 2Arb(π), 12)
    # Combination of integers and non-integers
    ss = [range(Arb(-4), Arb(4), 8); Arb.(-3:3)]
    @testset "clausens(x, s)" begin
        @testset "x::Arb, s::Arb" begin
            # Check that different evaluations agree
            for x in xs[2:end-1]
                for s in ss
                    res1 = HighestCuspedWave._clausens_polylog(x, s)
                    res2 = HighestCuspedWave._clausens_zeta(x, s)
                    res3 = clausens(x, s)
                    @test isfinite(res1)
                    @test isfinite(res2)
                    @test isfinite(res3)
                    @test Arblib.overlaps(res1, res2)
                    @test Arblib.overlaps(res1, res3)
                    @test Arblib.overlaps(res2, res3)
                end
            end

            # Wide x
            for lower in range(Arb(-10), 8, 10)
                for upper in range(lower + 0.01, 10, 10)
                    x_interval = Arb((lower, upper))
                    for s in ss
                        y = clausens(x_interval, s)
                        @test isfinite(y) || !(s > 1)
                        for x in range(lower, upper, length = 5)
                            @test Arblib.overlaps(clausens(x, s), y)
                        end
                    end
                end
            end

            # Wide s
            for x in xs[2:end-1]
                for s in ss
                    s_interval = add_error(s, Arb(0.0001))
                    y1 = clausens(x, s_interval)
                    @test isfinite(y1)
                    for t in [s; range(getinterval(Arb, s_interval)..., 4)]
                        y2 = clausens(x, t)
                        @test isfinite(y2)
                        @test Arblib.overlaps(y1, y2)
                    end
                end
            end

            # Check that _clausen_zeta throws an error outside of the domain
            @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(-1), Arb(2.5))
            @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(7), Arb(2.5))
        end

        @testset "x::ArbSeries, s::Arb" begin
            for s in ss
                # Check that it agrees with clausens
                @test Arblib.overlaps(
                    clausens(ArbSeries((2, 1, 0)), s),
                    Arblib.derivative(-clausenc(ArbSeries((2, 1, 0, 0)), s + 1)),
                )

                # Check that it works as Taylor expansion for clausens(sin(x), s)
                for degree = 1:4
                    x0 = Arb(1)
                    x = Arb((0.9, 1.1))
                    p = clausens(sin(ArbSeries((x0, 1); degree)), s)
                    R =
                        (x - x0)^(degree + 1) *
                        clausens(sin(ArbSeries((x, 1), degree = degree + 1)), s)[degree+1]
                    @test Arblib.overlaps(p(0.9 - x0) + R, clausens(sin(Arb(0.9)), s))
                    @test Arblib.overlaps(p(1.1 - x0) + R, clausens(sin(Arb(1.1)), s))
                end
            end
        end

        @testset "x::Arb, s::ArbSeries" begin
            # Check that different evaluations agree
            for x in xs[2:end-1]
                for s in ss
                    s_series = ArbSeries((s, 1), degree = 2)
                    res1 = HighestCuspedWave._clausens_polylog(x, s_series)
                    res2 = HighestCuspedWave._clausens_zeta(x, s_series)
                    res3 = clausens(x, s_series)
                    @test isfinite(res1)
                    @test isfinite(res2)
                    @test isfinite(res3)
                    @test Arblib.overlaps(res1, res2)
                    @test Arblib.overlaps(res1, res3)
                    @test Arblib.overlaps(res2, res3)
                end
            end

            # Wide s
            for x in xs[2:end-1]
                for s0 in ss
                    s0_interval = add_error(s0, Arb(0.0001))
                    s = ArbSeries((s0_interval, 1), degree = 2)
                    y1 = clausens(x, s)
                    @test isfinite(y1)
                    for t in [s0; range(getinterval(Arb, s0_interval)..., 4)]
                        y2 = clausens(x, ArbSeries((t, 1), degree = 2))
                        @test isfinite(y2)
                        @test Arblib.overlaps(y1, y2)
                    end
                end
            end
        end

        @testset "x::Any, s::Any" begin
            @test clausens(1.5, 2.4) ≈ Float64(clausens(Arb(1.5), 2.4))
            @test clausens(1.5, 2) == Float64(clausens(Arb(1.5), 2))
            @test clausens(2, 2.4) ≈ Float64(clausens(Arb(2), 2.4))
            @test clausens(2, 2) == Float64(clausens(Arb(2), 2))
        end
    end

    @testset "clausens(x, s, β)" begin
        @testset "x::Arb, s::Arb" begin
            # Check that different evaluations agree
            for β in [1, 2]
                for x in xs[2:end-1]
                    for s in ss
                        res1 = HighestCuspedWave._clausens_polylog(x, s, β)
                        res2 = HighestCuspedWave._clausens_zeta(x, s, β)
                        res3 = clausens(x, s, β)
                        @test isfinite(res1)
                        @test isfinite(res2)
                        @test isfinite(res3)
                        @test Arblib.overlaps(res1, res2)
                        @test Arblib.overlaps(res1, res3)
                        @test Arblib.overlaps(res2, res3)
                    end
                end
            end

            # Wide x
            for β in [1, 2]
                for lower in range(Arb(-10), 8, 10)
                    for upper in range(lower + 0.01, 10, 10)
                        x_interval = Arb((lower, upper))
                        for s in ss
                            y = clausens(x_interval, s, β)
                            @test isfinite(y) || !(s > 1)
                            for x in range(lower, upper, length = 5)
                                @test Arblib.overlaps(clausens(x, s, β), y)
                            end
                        end
                    end
                end
            end

            # Wide s
            for x in xs[2:end-1]
                for s in ss
                    s_interval = add_error(s, Arb(0.0001))
                    y1 = clausens(x, s_interval)
                    @test isfinite(y1)
                    for t in [s; range(getinterval(Arb, s_interval)..., 4)]
                        y2 = clausens(x, t)
                        @test isfinite(y2)
                        @test Arblib.overlaps(y1, y2)
                    end
                end
            end

            # Check that _clausen_zeta throws an error outside of the domain
            @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(-1), Arb(2.5))
            @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(7), Arb(2.5))
        end

        @testset "x::Any, s::Any" begin
            for β in [1, 2]
                @test clausens(1.5, 2.4, β) ≈ Float64(clausens(Arb(1.5), 2.4, β))
                @test clausens(1.5, 2, β) == Float64(clausens(Arb(1.5), 2, β))
                @test clausens(2, 2.4, β) ≈ Float64(clausens(Arb(2), 2.4, β))
                @test clausens(2, 2, β) == Float64(clausens(Arb(2), 2, β))
            end
        end
    end
end

@testset "clausens_expansion" begin
    s = Arb(0.5)
    for M in [3, 6]
        for x in range(Arb(0), 2Arb(π), length = 100)[1:end-1]
            C, e, P, E = HighestCuspedWave.clausens_expansion(x, s, M)
            @test Arblib.overlaps(C * abs(x)^e + P(x) + E * x^(2M + 1), clausens(x, s))
        end
    end
end
