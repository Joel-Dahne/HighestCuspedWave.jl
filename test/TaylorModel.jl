@testset "TaylorModel" begin
    fs = [sin, cos, exp, atan]

    @testset "Construction" begin
        for f in fs
            for x0 in Arb[-2, 0, 1]
                for r in Mag[1e-10, 1e-5, 1e0]
                    I = add_error(x0, r)
                    for n = 0:5
                        M1 = TaylorModel(f, I, x0, degree = n, enclosure_degree = -1)
                        M2 = TaylorModel(f, I, x0, degree = n, enclosure_degree = 0)
                        M3 = TaylorModel(f, I, x0, degree = n, enclosure_degree = 1)
                        @test Arblib.overlaps(M1, M2)
                        @test Arblib.overlaps(M1, M3)
                        @test Arblib.overlaps(M2, M3)
                        for x in range(x0 - r, x0 + r, 100)
                            fx = f(x)
                            @test Arblib.overlaps(fx, M1(x))
                            @test Arblib.overlaps(fx, M2(x))
                            @test Arblib.overlaps(fx, M3(x))
                        end
                    end
                end
            end
        end
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
                for r in Mag[1e-10, 1e-5, 1e0]
                    I = add_error(x0, r)
                    for n = 0:5
                        M = TaylorModel(f, I, x0, degree = n)
                        for m = 0:n
                            M_truncated = HighestCuspedWave.truncate(M, degree = m)
                            for x in range(x0 - r, x0 + r, 10)
                                @test Arblib.overlaps(M_truncated(x), M(x))
                            end
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
                    for r in Mag[1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for n = 0:5
                            Mf = TaylorModel(f, I, x0, degree = n)
                            Mgf = TaylorModel(g ∘ f, I, x0, degree = n)

                            @test Arblib.overlaps(Mgf, HighestCuspedWave.compose(g, Mf))
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
                    for r in Mag[1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for n = 0:5
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

                            if !Arblib.contains_zero(Mg.p[0])
                                Mf_div_g = TaylorModel(x -> f(x) / g(x), I, x0, degree = n)

                                @test Arblib.overlaps(Mf / Mg, Mf_div_g)
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
                for r in Mag[1e-10, 1e-5, 1e0]
                    I = add_error(x0, r)
                    for n = 0:5
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
                for r in Mag[1e-10, 1e-5, 1e0]
                    I = add_error(x0, r)
                    for n = 0:5
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
                    for r in Mag[1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for order = 1:3
                            for n = 0:5
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
                                        f(y - x0)^order / g(y - x0)^order
                                    end
                                end

                                @test Arblib.overlaps(
                                    HighestCuspedWave.div_removable(Mf, Mg, order),
                                    Mf_div_g,
                                )
                            end
                        end
                    end
                end
            end
        end
    end
end
