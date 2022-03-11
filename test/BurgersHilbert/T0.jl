@testset "T0" begin
    u0 = BHAnsatz{Arb}()

    f = T0(u0)
    g = T0(u0, Asymptotic())
    for x in exp.(range(log(Arb("1e-4")), log(Arb("1e-1")), length = 50))
        y1 = f(x)
        y2 = g(x)
        @test y1 < y2 || Arblib.overlaps(y1, y2)
    end
end
