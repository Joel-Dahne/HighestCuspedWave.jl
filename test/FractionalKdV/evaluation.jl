@testset "evaluation" begin
    u0s = [
        FractionalKdVAnsatz(Arb(-0.9)),
        FractionalKdVAnsatz(add_error(Arb(-0.9), Arb(1e-6))),
        FractionalKdVAnsatz(Arb(-0.6)),
        FractionalKdVAnsatz(add_error(Arb(-0.6), Arb(1e-6))),
        FractionalKdVAnsatz(Arb(-0.3)),
        FractionalKdVAnsatz(add_error(Arb(-0.3), Arb(1e-6))),
        FractionalKdVAnsatz(Arb(-0.9), use_bhkdv = true),
        FractionalKdVAnsatz(add_error(Arb(-0.9), Arb(1e-6)), use_bhkdv = true),
        FractionalKdVAnsatz(Arb(-0.995), use_bhkdv = true),
        FractionalKdVAnsatz(add_error(Arb(-0.995), Arb(1e-8)), use_bhkdv = true),
        # These have s values that overlap integers
        FractionalKdVAnsatz(Arb("[-0.60572740288 +/- 4.35e-11]")),
        FractionalKdVAnsatz(Arb("[-0.6268228870 +/- 8.17e-11]")),
    ]

    xs = exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), 5))
    # In theory the asymptotic an non-asymptotic results for ArbSeries
    # might not overlap. The reason being that the remainder term is
    # not computed using ArbSeries. In practice they do however
    # overlap.
    xs_series = [ArbSeries((x, 1, 0)) for x in xs]

    @testset "u0" begin
        for u0 in u0s
            for x in xs
                y1 = u0(x)
                y2 = u0(x, Asymptotic())
                @test isfinite(y1)
                @test isfinite(y2)
                @test Arblib.overlaps(y1, y2)
            end

            @test Arblib.overlaps(
                u0(Arb((0, xs[1])), Asymptotic()),
                u0(Arb(xs[1]), Asymptotic()),
            )
        end
    end

    @testset "u0_derivative" begin
        for u0 in u0s
            for x in xs
                y1 = HighestCuspedWave.u0_derivative(u0, x)
                y2 = ArbExtras.derivative_function(x -> u0(x))(x)
                @test isfinite(y1)
                @test isfinite(y2)
                @test Arblib.overlaps(y1, y2)
            end
        end
    end

    @testset "H(u0)" begin
        for u0 in u0s
            f = H(u0)
            g = H(u0, Asymptotic())
            for x in xs
                y1 = f(x)
                y2 = g(x)
                @test isfinite(y1)
                @test isfinite(y2)
                @test Arblib.overlaps(y1, y2)
            end

            @test Arblib.overlaps(g(Arb((0, xs[1]))), g(xs[1]))
        end
    end

    @testset "F0(u0)" begin
        for u0 in u0s
            f = F0(u0)
            g = F0(u0, Asymptotic(), ϵ = Arb("6e-1"))
            for x in ifelse(u0.use_bhkdv, xs, [xs; xs_series])
                y1 = f(x)
                y2 = g(x)
                @test isfinite(y1)
                @test isfinite(y2)
                @test Arblib.overlaps(y1, y2)
            end

            @test Arblib.overlaps(g(Arb((0, xs[1]))), g(xs[1]))
        end
    end

    @testset "inv_u0_normalised(u0)" begin
        for u0 in u0s
            f = x -> x^(-u0.α) / u0(x)
            g = HighestCuspedWave.inv_u0_normalised(u0)
            for x in [xs; xs_series]
                y1 = f(x)
                y2 = g(x)
                @test isfinite(y1)
                @test isfinite(y2)
                @test Arblib.overlaps(y1, y2)
            end

            @test Arblib.overlaps(g(Arb((0, xs[1]))), g(xs[1]))
        end
    end
end
