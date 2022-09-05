@testset "TaylorModel" begin
    fs = [sin, cos, exp, atan]

    # Test if M looks like it is a Taylor model of f. We check so that
    # M overlaps with f at the midpoint, the endpoints and two points
    # close to the midpoint.
    test_f = (M, f) -> begin
        Arblib.overlaps(M.p[0], f(M.x0)) || return false
        if !Arblib.isexact(M.I)
            x_l, x_u = getinterval(Arb, M.I)
            Arblib.overlaps(M(x_l), f(x_l)) || return false
            Arblib.overlaps(M(x_u), f(x_u)) || return false
            # Points close to the midpoint
            x1 = (4x_l + 6x_u) / 10
            x2 = (6x_l + 4x_u) / 10
            Arblib.overlaps(M(x1), f(x1)) || return false
            Arblib.overlaps(M(x2), f(x2)) || return false
        end
        return true
    end

    @testset "Construction" begin
        for f in fs
            for x0 in Arb[-2, 0, 1]
                for r in Mag[0, 1e-10, 1e-5, 1e0]
                    I = add_error(x0, r)
                    for n = 0:3
                        M1 = TaylorModel(f, I, x0, degree = n, enclosure_degree = -1)
                        M2 = TaylorModel(f, I, x0, degree = n, enclosure_degree = 0)
                        M3 = TaylorModel(f, I, x0, degree = n, enclosure_degree = 1)
                        @test Arblib.overlaps(M1, M2)
                        @test Arblib.overlaps(M1, M3)
                        @test Arblib.overlaps(M2, M3)
                        @test test_f(M1, f)
                        @test test_f(M2, f)
                        @test test_f(M3, f)
                    end
                end
            end
        end

        M = zero(TaylorModel(ArbSeries((1, 2)), Arb((-1, 1)), Arb(0.5)))
        @test izero(M)
        @test !isone(M)
        @test isequal(M.I, Arb((-1, 1)))
        @test isequal(M.x0, 0.5)

        M = one(TaylorModel(ArbSeries((1, 2)), Arb((-1, 1)), Arb(0.5)))
        @test one(M)
        @test !iszero(M)
        @test isequal(M.I, Arb((-1, 1)))
        @test isequal(M.x0, 0.5)
    end

    @testset "checkcompatible" begin
        @test HighestCuspedWave.checkcompatible(
            Bool,
            TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
            TaylorModel(cos, Arb((-1, 1)), Arb(0), degree = 3),
        )

        @test !HighestCuspedWave.checkcompatible(
            Bool,
            TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
            TaylorModel(cos, Arb((-1, 2)), Arb(0), degree = 3),
        )

        @test !HighestCuspedWave.checkcompatible(
            Bool,
            TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
            TaylorModel(cos, Arb((-1, 1)), Arb(1), degree = 3),
        )

        @test !HighestCuspedWave.checkcompatible(
            Bool,
            TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
            TaylorModel(cos, Arb((-1, 1)), Arb(0), degree = 4),
        )

        @test_throws ErrorException HighestCuspedWave.checkcompatible(
            TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
            TaylorModel(cos, Arb((-1, 2)), Arb(0), degree = 3),
        )

        @test_throws ErrorException HighestCuspedWave.checkcompatible(
            TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
            TaylorModel(cos, Arb((-1, 1)), Arb(1), degree = 3),
        )

        @test_throws ErrorException HighestCuspedWave.checkcompatible(
            TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
            TaylorModel(cos, Arb((-1, 1)), Arb(0), degree = 4),
        )
    end

    @testset "truncate" begin
        for f in fs
            for x0 in Arb[-2, 0, 1]
                for r in Mag[0, 1e-10, 1e-5, 1e0]
                    I = add_error(x0, r)
                    for n = 0:3
                        M = TaylorModel(f, I, x0, degree = n)
                        for m = 0:n
                            M_truncated = HighestCuspedWave.truncate(M, degree = m)
                            @test test_f(M_truncated, f)
                        end
                    end
                end
            end
        end
    end

    @testset "compose" begin
        for f in fs
            for g in fs
                for x0 in Arb[-2, 0, 1]
                    for r in Mag[0, 1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for n = 0:3
                            Mf = TaylorModel(f, I, x0, degree = n)
                            Mgf = TaylorModel(g ∘ f, I, x0, degree = n)

                            @test Arblib.overlaps(Mgf, HighestCuspedWave.compose(g, Mf))
                            @test test_f(HighestCuspedWave.compose(g, Mf), g ∘ f)
                        end
                    end
                end
            end
        end
    end

    @testset "Arithmetic" begin
        for f in fs
            for g in fs
                for x0 in Arb[-2, 0, 1]
                    for r in Mag[0, 1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for n = 0:3
                            Mf = TaylorModel(f, I, x0, degree = n)
                            Mg = TaylorModel(g, I, x0, degree = n)

                            Mf_plus_g = TaylorModel(x -> f(x) + g(x), I, x0, degree = n)
                            Mf_minus_g = TaylorModel(x -> f(x) - g(x), I, x0, degree = n)
                            Mf_mul_g = TaylorModel(x -> f(x) * g(x), I, x0, degree = n)
                            Mneg_f = TaylorModel(x -> -f(x), I, x0, degree = n)

                            @test Arblib.overlaps(Mf + Mg, Mf_plus_g)
                            @test Arblib.overlaps(Mf - Mg, Mf_minus_g)
                            @test Arblib.overlaps(Mf * Mg, Mf_mul_g)
                            @test Arblib.overlaps(-Mf, Mneg_f)

                            @test test_f(Mf + Mg, x -> f(x) + g(x))
                            @test test_f(Mf - Mg, x -> f(x) - g(x))
                            @test test_f(Mf * Mg, x -> f(x) * g(x))
                            @test test_f(-Mf, x -> -f(x))

                            if !Arblib.contains_zero(Mg.p[0])
                                Mf_div_g = TaylorModel(x -> f(x) / g(x), I, x0, degree = n)

                                @test Arblib.overlaps(Mf / Mg, Mf_div_g)
                                @test test_f(Mf / Mg, x -> f(x) / g(x))
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Scalar arithmetic" begin
        for f in fs
            for x0 in Arb[-2, 0, 1]
                for r in Mag[0, 1e-10, 1e-5, 1e0]
                    I = add_error(x0, r)
                    for n = 0:3
                        Mf = TaylorModel(f, I, x0, degree = n)

                        c = Arb(2)
                        Mf_add_c = TaylorModel(x -> f(x) + c, I, x0, degree = n)
                        Mf_sub_c = TaylorModel(x -> f(x) - c, I, x0, degree = n)
                        Mf_mul_c = TaylorModel(x -> f(x) * c, I, x0, degree = n)
                        Mf_div_c = TaylorModel(x -> f(x) / c, I, x0, degree = n)
                        c_add_Mf = TaylorModel(x -> c + f(x), I, x0, degree = n)
                        c_sub_Mf = TaylorModel(x -> c - f(x), I, x0, degree = n)
                        c_mul_Mf = TaylorModel(x -> c * f(x), I, x0, degree = n)
                        c_div_Mf = TaylorModel(x -> c / f(x), I, x0, degree = n)

                        @test Arblib.overlaps(Mf + c, Mf_add_c)
                        @test Arblib.overlaps(Mf - c, Mf_sub_c)
                        @test Arblib.overlaps(Mf * c, Mf_mul_c)
                        @test Arblib.overlaps(Mf / c, Mf_div_c)
                        @test Arblib.overlaps(c + Mf, Mf_add_c)
                        @test Arblib.overlaps(c - Mf, c_sub_Mf)
                        @test Arblib.overlaps(c * Mf, Mf_mul_c)
                        @test Arblib.overlaps(c / Mf, c_div_Mf)
                    end
                end
            end
        end
    end

    @testset "Shift" begin
        for f in fs
            for x0 in Arb[-2, 0, 1]
                for r in Mag[0, 1e-10, 1e-5, 1e0]
                    I = add_error(x0, r)
                    for n = 0:3
                        Mf = TaylorModel(f, I, x0, degree = n)
                        Mf_mul_x = TaylorModel(x -> (x - x0) * f(x), I, x0, degree = n + 1)
                        Mf_mul_x2 =
                            TaylorModel(x -> (x - x0)^2 * f(x), I, x0, degree = n + 2)

                        @test Arblib.overlaps(Mf, Mf_mul_x << 1)
                        @test Arblib.overlaps(Mf, Mf_mul_x2 << 2)
                        @test Arblib.overlaps(Mf_mul_x, Mf_mul_x2 << 1)
                        @test Arblib.overlaps(Mf >> 1, Mf_mul_x)
                        @test Arblib.overlaps(Mf >> 2, Mf_mul_x2)
                        @test Arblib.overlaps(Mf_mul_x >> 1, Mf_mul_x2)
                    end
                end
            end
        end
    end

    @testset "div_removable" begin
        for f in [sin, x -> exp(x) - 1, atan]
            for g in [sin, x -> exp(x) - 1, atan]
                for x0 in Arb[-2, 0, 1]
                    for r in Mag[0, 1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for order = 1:3
                            for n = 0:3
                                Mf = TaylorModel(
                                    x -> f(x - x0)^order,
                                    I,
                                    x0,
                                    degree = n + order,
                                )
                                Mg = TaylorModel(
                                    x -> g(x - x0)^order,
                                    I,
                                    x0,
                                    degree = n + order,
                                )

                                Mf_div_g = TaylorModel(I, x0, degree = n) do x
                                    # fx_div_x only uses ArbSeries in
                                    # cases when it overlaps zero.
                                    if x isa ArbSeries || Arblib.contains_zero(x - x0)
                                        HighestCuspedWave.fx_div_x(f, x - x0)^order /
                                        HighestCuspedWave.fx_div_x(g, x - x0)^order
                                    else
                                        f(x - x0)^order / g(x - x0)^order
                                    end
                                end

                                @test Arblib.overlaps(
                                    HighestCuspedWave.div_removable(Mf, Mg, order),
                                    Mf_div_g,
                                )
                                @test test_f(
                                    Mf_div_g,
                                    x -> if Arblib.contains_zero(x - x0)
                                        HighestCuspedWave.fx_div_x(f, x - x0)^order /
                                        HighestCuspedWave.fx_div_x(g, x - x0)^order
                                    else
                                        f(x - x0)^order / g(x - x0)^order
                                    end,
                                )
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "clausenc" begin
        for x in Arb[0.1, 1, 3]
            for f in [identity, sin, exp]
                for s0 in Arb[0.5, 1.6, 2.7]
                    for r in Mag[0, 1e-10, 1e-5]
                        I = add_error(s0, r)
                        for n = 0:3
                            Ms = TaylorModel(f, I, s0, degree = n)
                            M = clausencmzeta(x, Ms)
                            @test test_f(M, s -> clausencmzeta(x, f(s)))
                        end
                    end
                end
            end
        end
    end

    @testset "abspow" begin
        for x in Arb[0, 0.5, 1.5]
            for f in [identity, sin, exp, y -> -y]
                for y0 in Arb[0.1, 1.6, 2.7]
                    for r in Mag[0, 1e-10, 1e-5]
                        I = add_error(y0, r)
                        for n = 0:3
                            My = TaylorModel(f, I, y0, degree = n)
                            M = HighestCuspedWave.abspow(x, My)
                            @test test_f(M, y -> HighestCuspedWave.abspow(x, f(y)))
                        end
                    end
                end
            end
        end
    end
end
