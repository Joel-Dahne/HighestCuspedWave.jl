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
    xs = range(Arb(0), 2Arb(π), 12)
    # Combination of integers and non-integers
    ss = [range(Arb(-4), Arb(4), 8); Arb.(-3:3)]

    @testset "clausenc(x, s)" begin
        @testset "x::Arb, s::Arb" begin
            # Check that different evaluations agree
            for x in xs[2:end-1]
                for s in ss
                    res1 = HighestCuspedWave._clausenc_polylog(x, s)
                    res2 = HighestCuspedWave._clausenc_zeta(x, s)
                    res3 = clausenc(x, s)
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
                        y = clausenc(x_interval, s)
                        @test isfinite(y) || !(s > 1)
                        for x in range(lower, upper, length = 5)
                            @test Arblib.overlaps(clausenc(x, s), y)
                        end
                    end
                end
            end

            # Wide s
            for x in xs[2:end-1]
                for s in ss
                    s_interval = add_error(s, Arb(0.0001))
                    y1 = clausenc(x, s_interval)
                    @test isfinite(y1)
                    for t in [s; range(getinterval(Arb, s_interval)..., 4)]
                        y2 = clausenc(x, t)
                        @test isfinite(y2)
                        @test Arblib.overlaps(y1, y2)
                    end
                end
            end

            # Check that _clausen_zeta throws an error outside of the domain
            @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(-1), Arb(2.5))
            @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(7), Arb(2.5))
        end

        @testset "x::Acb, s::Arb" begin
            # Check that different evaluations agree
            for x in range(Arb(0), 2Arb(π), length = 10)[2:end-1]
                for y in range(Arb(-10), Arb(10), length = 8)
                    z = Acb(x, y)
                    for s in ss
                        res1 =
                            (
                                HighestCuspedWave.polylog(Acb(s), exp(im * z)) +
                                HighestCuspedWave.polylog(Acb(s), exp(-im * z))
                            ) / 2
                        res2 = HighestCuspedWave._clausenc_zeta(z, s)
                        res3 = clausenc(z, s)
                        @test isfinite(res1)
                        @test isfinite(res2) || Arblib.contains_int(s)
                        @test isfinite(res3)
                        @test Arblib.overlaps(res1, res2)
                        @test Arblib.overlaps(res1, res3)
                        @test Arblib.overlaps(res2, res3)
                    end
                end
            end
        end

        @testset "x::ArbSeries, s::Arb" begin
            for s in ss
                # Check that it agrees with clausens
                @test Arblib.overlaps(
                    clausenc(ArbSeries((2, 1, 0)), s),
                    Arblib.derivative(clausens(ArbSeries((2, 1, 0, 0)), s + 1)),
                )

                # Check that it works as Taylor expansion for clausenc(sin(x), s)
                for degree = 1:4
                    x0 = Arb(1)
                    x = Arb((0.9, 1.1))
                    p = clausenc(sin(ArbSeries((x0, 1); degree)), s)
                    R =
                        (x - x0)^(degree + 1) *
                        clausenc(sin(ArbSeries((x, 1), degree = degree + 1)), s)[degree+1]
                    @test Arblib.overlaps(p(0.9 - x0) + R, clausenc(sin(Arb(0.9)), s))
                    @test Arblib.overlaps(p(1.1 - x0) + R, clausenc(sin(Arb(1.1)), s))
                end
            end
        end

        @testset "x::Arb, s::ArbSeries" begin
            # Check that different evaluations agree
            for x in xs[2:end-1]
                for s in ss
                    s_series = ArbSeries((s, 1), degree = 2)
                    res1 = HighestCuspedWave._clausenc_polylog(x, s_series)
                    res2 = HighestCuspedWave._clausenc_zeta(x, s_series)
                    res3 = clausenc(x, s_series)
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
                    y1 = clausenc(x, s)
                    @test isfinite(y1)
                    for t in [s0; range(getinterval(Arb, s0_interval)..., 4)]
                        y2 = clausenc(x, ArbSeries((t, 1), degree = 2))
                        @test isfinite(y2)
                        @test Arblib.overlaps(y1, y2)
                    end
                end
            end
        end

        @testset "x::Any, s::Any" begin
            @test clausenc(1.5, 2.4) ≈ Float64(clausenc(Arb(1.5), 2.4))
            @test clausenc(1.5, 2) == Float64(clausenc(Arb(1.5), 2))
            @test clausenc(2, 2.4) ≈ Float64(clausenc(Arb(2), 2.4))
            @test clausenc(2, 2) == Float64(clausenc(Arb(2), 2))
        end
    end

    @testset "clausenc(x, s, β)" begin
        @testset "x::Arb, s::Arb" begin
            # Check that different evaluations agree
            for β in [1, 2]
                for x in xs[2:end-1]
                    for s in ss
                        res1 = HighestCuspedWave._clausenc_polylog(x, s, β)
                        res2 = HighestCuspedWave._clausenc_zeta(x, s, β)
                        res3 = clausenc(x, s, β)
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
                            y = clausenc(x_interval, s, β)
                            @test isfinite(y) || !(s > 1)
                            for x in range(lower, upper, length = 5)
                                @test Arblib.overlaps(clausenc(x, s, β), y)
                            end
                        end
                    end
                end
            end

            # Wide s
            for x in xs[2:end-1]
                for s in ss
                    s_interval = add_error(s, Arb(0.0001))
                    y1 = clausenc(x, s_interval)
                    @test isfinite(y1)
                    for t in [s; range(getinterval(Arb, s_interval)..., 4)]
                        y2 = clausenc(x, t)
                        @test isfinite(y2)
                        @test Arblib.overlaps(y1, y2)
                    end
                end
            end

            # Check that _clausen_zeta throws an error outside of the domain
            @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(-1), Arb(2.5))
            @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(7), Arb(2.5))
        end

        @testset "x::ArbSeries, s::Arb" begin
            for β in [1, 2]
                for s in ss
                    # Check that it works as Taylor expansion for clausenc(sin(x), s, β)
                    for degree = 1:4
                        x0 = Arb(1)
                        x = Arb((0.9, 1.1))
                        p = clausenc(sin(ArbSeries((x0, 1); degree)), s, β)
                        R =
                            (x - x0)^(degree + 1) *
                            clausenc(sin(ArbSeries((x, 1), degree = degree + 1)), s, β)[degree+1]
                        @test Arblib.overlaps(
                            p(0.9 - x0) + R,
                            clausenc(sin(Arb(0.9)), s, β),
                        )
                        @test Arblib.overlaps(
                            p(1.1 - x0) + R,
                            clausenc(sin(Arb(1.1)), s, β),
                        )
                    end
                end
            end
        end

        @testset "x::Any, s::Any" begin
            for β in [1, 2]
                @test clausenc(1.5, 2.4, β) ≈ Float64(clausenc(Arb(1.5), 2.4, β))
                @test clausenc(1.5, 2, β) == Float64(clausenc(Arb(1.5), 2, β))
                @test clausenc(2, 2.4, β) ≈ Float64(clausenc(Arb(2), 2.4, β))
                @test clausenc(2, 2, β) == Float64(clausenc(Arb(2), 2, β))
            end
        end
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

    # TODO: Add tests for wide s

    # TODO: Add tests for β
end

@testset "clausenc_expansion_odd_s_singular_K1_K2" begin
    for s in Arb[1.5, 2.9, 3.1]
        for m in [1, 2, 3]
            K1, K2 = HighestCuspedWave.clausenc_expansion_odd_s_singular_K1_K2(s, m)
            K3 = HighestCuspedWave.clausenc_expansion_odd_s_singular_K3(m)
            for x in (Arb(-0.5), Arb(0.1), ArbSeries((-0.5, 2, 3)), ArbSeries((0.1, 3, 2)))
                r1 =
                    K1 * abs(x)^(s - 1) +
                    K2 * x^2m +
                    K3 *
                    HighestCuspedWave.x_pow_s_x_pow_t_m1_div_t(x, Arb(2m), s - (2m + 1))
                r2 =
                    gamma(1 - s) * sinpi(s / 2) * abs(x)^(s - 1) +
                    (-1)^m * zeta(s - 2m) * x^2m / factorial(2m)
                @test Arblib.overlaps(r1, r2)
            end
        end
    end
end

@testset "clausencmzeta" begin
    xs = range(Arb(0), 2Arb(π), 20)
    # Combination of integers and non-integers In this case we only
    # really care about s > 1
    ss = [range(Arb(1), Arb(4), 8)[2:end]; Arb.(2:3)]

    @testset "x::Arb, s::Arb" begin
        # Check that different evaluations agree
        for x in xs[2:end-1]
            for s in ss
                res1 = clausencmzeta(x, s)
                res2 = clausenc(x, s) - zeta(s)
                @test isfinite(res1)
                @test isfinite(res2)
                @test Arblib.overlaps(res1, res2)
            end
        end

        # Wide x
        for lower in range(Arb(-10), 8, 10)
            for upper in range(lower + 0.01, 10, 10)
                x_interval = Arb((lower, upper))
                for s in ss
                    y = clausencmzeta(x_interval, s)
                    @test isfinite(y) || !(s > 1)
                    for x in range(lower, upper, length = 5)
                        @test Arblib.overlaps(clausencmzeta(x, s), y)
                    end
                end
            end
        end

        # Wide s
        for x in xs[2:end-1]
            for s in ss
                s_interval = add_error(s, Arb(0.0001))
                y1 = clausencmzeta(x, s_interval)
                @test isfinite(y1)
                for t in [s; range(getinterval(Arb, s_interval)..., 4)]
                    y2 = clausencmzeta(x, t)
                    @test isfinite(y2)
                    @test Arblib.overlaps(y1, y2)
                end
            end
        end
    end

    @testset "x::ArbSeries, s::Arb" begin
        for s in ss
            # Check that it agrees with clausens
            @test Arblib.overlaps(
                Arblib.derivative(clausencmzeta(ArbSeries((2, 1, 0, 0)), s)),
                Arblib.derivative(clausenc(ArbSeries((2, 1, 0, 0)), s)),
            )

            # Check that it works as Taylor expansion for clausenc(sin(x), s)
            for degree = 1:4
                x0 = Arb(1)
                x = Arb((0.9, 1.1))
                p = clausencmzeta(sin(ArbSeries((x0, 1); degree)), s)
                R =
                    (x - x0)^(degree + 1) *
                    clausencmzeta(sin(ArbSeries((x, 1), degree = degree + 1)), s)[degree+1]
                @test Arblib.overlaps(p(0.9 - x0) + R, clausencmzeta(sin(Arb(0.9)), s))
                @test Arblib.overlaps(p(1.1 - x0) + R, clausencmzeta(sin(Arb(1.1)), s))
            end
        end
    end

    @testset "x::Arb, s::ArbSeries" begin
        for s0 in ss
            # Check that it works as Taylor expansion for clausencmzeta(x, sin(s))
            s = add_error(s0, Mag(0.1))
            x = Arb(1)
            for degree = 1:4
                p = clausencmzeta(x, sin(ArbSeries((s0, 1); degree)))
                R =
                    (s - s0)^(degree + 1) *
                    clausencmzeta(x, sin(ArbSeries((s, 1), degree = degree + 1)))[degree+1]
                @test Arblib.overlaps(
                    p(lbound(s) - s0) + R,
                    clausencmzeta(x, sin(lbound(Arb, s))),
                )
                @test Arblib.overlaps(
                    p(ubound(s) - s0) + R,
                    clausencmzeta(x, sin(ubound(Arb, s))),
                )
            end
        end
    end

    @testset "x::Arb, s::Arb, β::Int" begin
        for β in [1, 2]
            for x in xs[2:end-1]
                for s in ss
                    res1 = clausencmzeta(x, s, β)
                    res2 =
                        clausenc(x, s, β) -
                        zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
                    @test Arblib.overlaps(res1, res2)
                end
            end
        end
    end

    @testset "x::Any, s::Any" begin
        @test clausencmzeta(1.5, 2.4) ≈ Float64(clausencmzeta(Arb(1.5), 2.4))
        @test clausencmzeta(1.5, 2) ≈ Float64(clausencmzeta(Arb(1.5), 2))
        @test clausencmzeta(2, 2.4) ≈ Float64(clausencmzeta(Arb(2), 2.4))
        @test clausencmzeta(2, 2) ≈ Float64(clausencmzeta(Arb(2), 2))
    end
end
