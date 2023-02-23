@testset "evaluation" begin
    u0 = BHKdVAnsatz(Arb(1e-5 * (1 - 0.9997)), BHAnsatz{Arb}())

    xs = exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), 4))
    xs_series = [ArbSeries((x, 2, 3)) for x in xs]

    @testset "u0" begin
        for x in xs
            y1 = u0(x)
            y2 = u0(x, Asymptotic())
            y3 = u0.v0(x)

            @test isfinite(y1)
            @test isfinite(y2)
            @test isfinite(y3)
            @test Arblib.overlaps(y1, y2)
            @test Arblib.overlaps(y1, y3)
            @test Arblib.overlaps(y2, y3)
        end

        for x in xs_series
            y1 = u0(x)
            y2 = u0.v0(x)

            @test isfinite(y1)
            @test isfinite(y2)
            @test Arblib.overlaps(y1, y2)
        end
    end

    @testset "H(u0)" begin
        for x in [xs; xs_series]
            y1 = H(u0)(x)
            y2 = H(u0.v0)(x)

            @test isfinite(y1)
            @test isfinite(y2)
            @test Arblib.overlaps(y1, y2)
        end
    end

    @testset "F0(u0)" begin
        f = F0(u0)
        g = HighestCuspedWave.F0_bound(u0)
        h = F0(u0, Asymptotic(), ϵ = Arb("6e-1"))
        for x in xs
            y1 = f(x)
            y2 = g(x)
            y3 = h(x)
            @test isfinite(y1)
            @test isfinite(y2)
            @test isfinite(y3)
            # These are not strictly guaranteed to overlap since g and
            # h only compute upper bounds. In practice the upper bound
            # is however good enough for them to overlap.
            @test Arblib.overlaps(y1, y2)
            @test Arblib.overlaps(abs(y1), y3) # h computes absolute value
            @test Arblib.overlaps(abs(y2), y3) # h computes absolute value
        end

        @test Arblib.overlaps(g(Arb((0, xs[1]))), g(xs[1]))
    end
end
