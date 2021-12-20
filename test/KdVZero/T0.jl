@testset "T0" begin
    # Set up interval of α to test the expansions on and take a number
    # of points on this interval
    α_lower = Arb(-0.001)
    α_upper = Arb(0)
    α_interval = Arb((α_lower, α_upper))
    αs = range(α_lower, α_upper, length = 100)[1:end-1]

    # Construct KdVZeroansatz for the interval
    u0 = KdVZeroAnsatz(α_interval)

    # Construct FractionalKdvansatz for each α
    u0s = map(αs) do α
        a = OffsetVector([HighestCuspedWave.finda0(α), HighestCuspedWave._finda1a2(α)...], 0:2)
        FractionalKdVAnsatz{Arb}(
            α,
            HighestCuspedWave.findp0(α),
            a,
            Arb[],
            one(Arb),
            Set{NTuple{3,Int}}([(2, 0, 0)]),
        )
    end

    @testset "Ball" begin
        x = Arb(1.5)
        p = T0(u0)(x)
        for (i, α) in enumerate(αs)
            @test Arblib.overlaps(p(α), T0(u0s[i])(x))
            #@show p(α) - T0(u0s[i])(x)
        end
    end

    @testset "Asymptotic" begin
        x = Arb(0.01)
        p = T0(u0, Asymptotic())(x)
        for (i, α) in enumerate(αs)
            @test Arblib.overlaps(p(α), T0(u0s[i], Asymptotic())(x))
            #@show p(α) - T0(u0s[i], Asymptotic())(x)
        end
    end

    @testset "Ball vs Asymptotic" begin
        f = T0(u0)
        g = T0(u0, Asymptotic())
        for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 50))
            p1 = f(x)
            p2 = g(x)
            @test all(Arblib.overlaps.(Arblib.coeffs(p1), Arblib.coeffs(p2)))
        end
    end
end
