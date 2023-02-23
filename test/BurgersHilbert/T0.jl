@testset "T0" begin
    u0 = BHAnsatz{Arb}()

    f = T0(u0)
    g = T0(u0, Asymptotic())
    for x in exp.(range(log(Arb("1e-4")), log(Arb("1e-1")), length = 10))
        y1 = f(x)
        y2 = g(x)
        @test isfinite(y1)
        @test isfinite(y2)
        # The asymptotic version only computes an upper bound and
        # not an enclosure
        @test y1 < y2 || Arblib.overlaps(y1, y2)
    end

    @test Arblib.overlaps(g(Arb((0, 1e-4))), g(Arb(1e-4)))
end
