@testset "evaluation" begin
    # Set up interval of α to test the expansions on and take a number
    # of points on this interval
    α_lower = Arb(-0.001)
    α_upper = Arb(0)
    α_interval = Arb((α_lower, α_upper))
    αs = range(α_lower, α_upper, length = 100)[1:end-1]

    # Construct KdVZeroansatz for the interval
    u0 = KdVZeroAnsatz(α_interval)

    # Ansatz on a much smaller interval
    u0_tight = KdVZeroAnsatz(Arb((-1e-10, 0)))

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

    @testset "u0" begin
        @testset "Ball" begin
            x = Arb(1.5)
            p = u0(x)
            for (i, α) in enumerate(αs)
                @test Arblib.overlaps(p(α), u0s[i](x))
                #@show p(α) - u0s[i](x)
            end
        end

        @testset "Asymptotic" begin
            x = Arb(1e-10)
            p = eval_expansion(u0, u0(x, AsymptoticExpansion()), x)
            for (i, α) in enumerate(αs)
                @test Arblib.overlaps(p(α), u0s[i](x, Asymptotic()))
                #@show p(α) - u0s[i](x, Asymptotic())
            end
        end

        @testset "Ball vs Asymptotic" begin
            expansion = u0_tight(Arb(1), AsymptoticExpansion())
            for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                p1 = u0_tight(x)
                p2 = eval_expansion(u0_tight, expansion, x)
                @test all(Arblib.overlaps.(Arblib.coeffs(p1), Arblib.coeffs(p2)))
                #@show p1 - p2
            end
        end
    end

    @testset "H(u0)" begin
        @testset "Ball" begin
            x = Arb(1.5)
            p = H(u0)(x)
            for (i, α) in enumerate(αs)
                @test Arblib.overlaps(p(α), H(u0s[i])(x))
                #@show p(α) - H(u0s[i])(x)
            end
        end

        @testset "Asymptotic" begin
            x = Arb(1e-10)
            p = eval_expansion(u0, H(u0, AsymptoticExpansion())(x), x)
            for (i, α) in enumerate(αs)
                @test Arblib.overlaps(p(α), H(u0s[i], Asymptotic())(x))
                #@show p(α) - H(u0s[i], Asymptotic())(x)
            end
        end

        @testset "Ball vs Asymptotic" begin
            f = H(u0_tight)
            expansion = H(u0_tight, AsymptoticExpansion())(Arb(1))
            g = x -> eval_expansion(u0_tight, expansion, x)
            for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                p1 = f(x)
                p2 = g(x)
                @test all(Arblib.overlaps.(Arblib.coeffs(p1), Arblib.coeffs(p2)))
                #@show p1 - p2
            end
        end
    end

    @testset "D" begin
        @testset "Ball" begin
            x = Arb(1.5)
            p = D(u0)(x)
            for (i, α) in enumerate(αs)
                @test Arblib.overlaps(p(α), D(u0s[i])(x))
                #@show p(α) - D(u0s[i])(x)
            end
        end

        @testset "Asymptotic" begin
            x = Arb(1e-10)
            p = eval_expansion(u0, D(u0, AsymptoticExpansion())(x), x)
            for (i, α) in enumerate(αs)
                @test Arblib.overlaps(p(α), D(u0s[i], Asymptotic())(x))
                #@show p(α)  D(u0s[i], Asymptotic())(x)
            end
        end


        @testset "Ball vs Asymptotic" begin
            f = D(u0_tight)
            expansion = D(u0_tight, AsymptoticExpansion())(Arb(1))
            g = x -> eval_expansion(u0_tight, expansion, x)
            for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                p1 = f(x)
                p2 = g(x)
                @test all(Arblib.overlaps.(Arblib.coeffs(p1), Arblib.coeffs(p2)))
                #@show p1 - p2
            end
        end
    end

    @testset "F0" begin
        @testset "Ball" begin
            x = Arb(1.5)
            p = F0(u0)(x)
            for (i, α) in enumerate(αs)
                @test Arblib.overlaps(p(α), F0(u0s[i])(x))
                #@show p(α) - F0(u0s[i])(x)
            end
        end

        @testset "Asymptotic" begin
            x = Arb(1e-10)
            p = F0(u0, Asymptotic(); ϵ = x)(x)
            for (i, α) in enumerate(αs)
                @test Arblib.overlaps(p(α), F0(u0s[i], Asymptotic())(x))
                #@show p(α) - F0(u0s[i], Asymptotic())(x)
            end
        end

        @testset "Ball vs Asymptotic" begin
            f = F0(u0_tight)
            g = F0(u0_tight, Asymptotic(), ϵ = Arb(1))
            for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                p1 = f(x)
                p2 = g(x)
                @test all(Arblib.overlaps.(Arblib.coeffs(p1), Arblib.coeffs(p2)))
                #@show p1 - p2
            end
        end
    end
end
