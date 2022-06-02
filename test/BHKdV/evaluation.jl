@testset "evaluation" begin
    # IMPROVE: This takes a long time to compute, consider computing
    # it elsewhere.
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
        for x in xs
            y1 = H(u0)(x)
            #y2 = H(u0, Asymptotic())(x)
            y2 = H(u0.v0)(x)

            @test isfinite(y1)
            @test isfinite(y2)
            @test Arblib.overlaps(y1, y2)
        end

        for x in xs_series
            y1 = H(u0)(x)
            y2 = H(u0.v0)(x)

            @test isfinite(y1)
            @test isfinite(y2)
            @test Arblib.overlaps(y1, y2)
        end
    end
end
