@testset "evaluation" begin
    u0 = BHAnsatz{Arb}()

    xs = exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), 10))

    @testset "u0" begin
        for x in xs
            y1 = u0(x)
            y2 = u0(x, Asymptotic())
            @test isfinite(y1)
            @test isfinite(y2)
            @test Arblib.overlaps(y1, y2)
        end
    end

    @testset "H(u0)" begin
        f = H(u0)
        g = H(u0, Asymptotic())
        for x in xs
            y1 = f(x)
            y2 = g(x)
            @test isfinite(y1)
            @test isfinite(y2)
            @test Arblib.overlaps(y1, y2)
        end
    end

    @testset "F0(u0)" begin
        @testset "Ball vs asymptotic" begin
            f = F0(u0)
            g = F0(u0, Asymptotic(), ϵ = Arb("6e-1"))
            for x in xs
                y1 = f(x)
                y2 = g(x)
                @test isfinite(y1)
                @test isfinite(y2)
                @test Arblib.overlaps(y1, y2)
            end
        end

        @testset "Exponent limit" begin
            f1 = F0(u0, Asymptotic())
            f2 = F0(u0, Asymptotic(), exponent_limit = Arb(1 // 4))
            f3 = F0(u0, Asymptotic(), exponent_limit = Arb("1e-2"))
            for x in exp.(range(log(Arb("1e-1000")), log(Arb("1e-10")), 100))
                y1 = f1(x)
                y2 = f2(x)
                y3 = f3(x)
                @test isfinite(y1)
                @test isfinite(y2)
                @test isfinite(y3)
                @test Arblib.overlaps(y1, y2)
                @test Arblib.overlaps(y1, y3)
                @test Arblib.overlaps(y2, y3)
            end
        end
    end

    @testset "inv_u0_normalised(u0)" begin
        f = x -> -x * log(x) / u0(x)
        g = HighestCuspedWave.inv_u0_normalised(u0, ϵ = Arb("6e-1"))
        for x in xs
            y1 = f(x)
            y2 = g(x)
            @test isfinite(y1)
            @test isfinite(y2)
            @test Arblib.overlaps(y1, y2)
        end
    end
end
