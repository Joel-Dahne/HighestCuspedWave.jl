@testset "evaluation" begin
    α0s = [Arb(0), Arb((-1e-3))]

    intervals = [
        -Arblib.nonnegative_part!(zero(Arb), Arb((0, 1e-3))),
        Arblib.add_error!(Arb(-1e-3), Arb(1e-6)),
    ]

    @testset "α0 = $α0" for (α0, interval) in zip(α0s, intervals)
        αs = range(getinterval(Arb, interval)..., 100)
        if iszero(α0)
            αs = αs[1:end-1]
        end

        # Ansatz on the full interval
        u0 = KdVZeroAnsatz(interval, α0)
        # Ansatz on a smaller interval
        if iszero(α0)
            interval_narrow = Arb((-1e-8, 0))
        else
            interval_narrow = Arblib.add_error!(Arb(α0), Arb(1e-10))
        end
        u0_narrow = KdVZeroAnsatz(interval_narrow, α0)
        # Ansatz on just the point
        u0_tight = KdVZeroAnsatz(α0, α0)

        # Construct FractionalKdvansatz for each α
        u0s = map(αs) do α
            a = OffsetVector(
                [HighestCuspedWave.finda0(α), HighestCuspedWave._finda1a2(α)...],
                0:2,
            )
            FractionalKdVAnsatz{Arb}(
                α,
                HighestCuspedWave.findp0(α),
                a,
                Arb[],
                one(Arb),
                Set{NTuple{3,Int}}([(2, 0, 0)]),
            )
        end

        # Take x exactly equal to 1.5 and x equal to 1.5 but with a radius
        # large enough to make iswide(x) true
        xs = [Arb(1.5), Arb(1.5)]
        Arblib.add_error!(xs[2], ubound(sqrt(eps(xs[2]))))
        @test HighestCuspedWave.iswide(xs[2])

        # Take x exactly equal to 1e-10 and x equal to 1e-10 but with a
        # radius large enough to make iswide(x) true
        xs_asym = [Arb(1e-10), Arb(1e-10)]
        Arblib.add_error!(xs_asym[2], ubound(sqrt(eps(xs_asym[2]))))
        @test HighestCuspedWave.iswide(xs_asym[2])

        @testset "u0" begin
            @testset "Ball" begin
                for x in xs
                    p = u0(x)
                    for (i, α) in enumerate(αs)
                        @test Arblib.overlaps(p(α), u0s[i](x))
                    end
                end
            end

            @testset "Asymptotic" begin
                for x in xs_asym
                    p = eval_expansion(u0, u0(x, AsymptoticExpansion()), x)
                    for (i, α) in enumerate(αs)
                        @test Arblib.overlaps(p(α), u0s[i](x, Asymptotic()))
                    end
                end
            end

            @testset "Ball vs Asymptotic - $type" for (u0, type) in (
                (u0_narrow, "narrow"),
                (u0_tight, "tight"),
            )
                expansion = u0(Arb(1), AsymptoticExpansion())
                for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                    p1 = u0_tight(x)
                    p2 = eval_expansion(u0, expansion, x)
                    @test Arblib.overlaps(p1, p2, require_compatible = false)
                end
            end
        end

        @testset "H(u0)" begin
            @testset "Ball" begin
                for x in xs
                    p = H(u0)(x)
                    for (i, α) in enumerate(αs)
                        @test Arblib.overlaps(p(α), H(u0s[i])(x))
                    end
                end
            end

            @testset "Asymptotic" begin
                for x in xs_asym
                    p = eval_expansion(u0, H(u0, AsymptoticExpansion())(x), x)
                    for (i, α) in enumerate(αs)
                        @test Arblib.overlaps(p(α), H(u0s[i], Asymptotic())(x))
                    end
                end
            end

            @testset "Ball vs Asymptotic - $type" for (u0, type) in (
                (u0_narrow, "narrow"),
                (u0_tight, "tight"),
            )
                f = H(u0)
                expansion = H(u0, AsymptoticExpansion())(Arb(1))
                g = x -> eval_expansion(u0, expansion, x)
                for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                    p1 = f(x)
                    p2 = g(x)
                    @test Arblib.overlaps(p1, p2, require_compatible = false)
                end
            end
        end

        @testset "D" begin
            @testset "Ball" begin
                for x in xs
                    p = D(u0)(x)
                    for (i, α) in enumerate(αs)
                        @test Arblib.overlaps(p(α), D(u0s[i])(x))
                    end
                end
            end

            @testset "Asymptotic" begin
                for x in xs_asym
                    p = eval_expansion(u0, D(u0, AsymptoticExpansion())(x), x)
                    for (i, α) in enumerate(αs)
                        @test Arblib.overlaps(p(α), D(u0s[i], Asymptotic())(x))
                    end
                end
            end


            @testset "Ball vs Asymptotic - $type" for (u0, type) in (
                (u0_narrow, "narrow"),
                (u0_tight, "tight"),
            )
                f = D(u0)
                expansion = D(u0, AsymptoticExpansion())(Arb(1))
                g = x -> eval_expansion(u0, expansion, x)
                for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                    p1 = f(x)
                    p2 = g(x)
                    @test Arblib.overlaps(p1, p2, require_compatible = false)
                end
            end
        end

        @testset "F0" begin
            @testset "Ball" begin
                for x in xs
                    p = F0(u0)(x)
                    for (i, α) in enumerate(αs)
                        @test Arblib.overlaps(p(α), F0(u0s[i])(x))
                    end
                end
            end

            @testset "Asymptotic" begin
                for x in xs_asym
                    p = F0(u0, Asymptotic())(x)
                    for (i, α) in enumerate(αs)
                        @test Arblib.overlaps(p(α), F0(u0s[i], Asymptotic())(x))
                    end
                end
            end

            @testset "Ball vs Asymptotic - $type" for (u0, type) in (
                (u0_narrow, "narrow"),
                (u0_tight, "tight"),
            )
                f = F0(u0)
                g = F0(u0, Asymptotic())
                for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                    p1 = f(x)
                    p2 = g(x)
                    @test Arblib.overlaps(p1, p2, require_compatible = false)
                end
            end

            if !iszero(u0.α0)
                @testset "F0_direct" begin
                    f = F0(u0, Asymptotic())
                    g = HighestCuspedWave.F0_direct(u0, Asymptotic())
                    for x in exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), length = 100))
                        @test Arblib.overlaps(f(x)(u0.α), g(x))
                    end
                end
            end
        end
    end
end
