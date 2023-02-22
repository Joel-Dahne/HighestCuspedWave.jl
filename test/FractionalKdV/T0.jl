@testset "T0" begin
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
    ]

    xs = exp.(range(log(Arb("1e-3")), log(Arb("5e-1")), 5))

    for u0 in u0s
        f = T0(u0)
        g = T0(u0, Asymptotic())
        for x in xs
            y1 = f(x)
            y2 = g(x)
            @test isfinite(y1)
            @test isfinite(y2)
            # The asymptotic version only computes an upper bound and
            # not an enclosure
            @test y1 < y2 || Arblib.overlaps(y1, y2)
        end

        @test Arblib.overlaps(g(Arb((0, xs[1]))), g(xs[1]))
    end
end
