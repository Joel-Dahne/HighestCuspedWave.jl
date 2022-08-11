@testset "evaluation" begin
    u0s = [
        FractionalKdVAnsatz(Arb(-0.9)),
        FractionalKdVAnsatz(add_error(Arb(-0.9), Arb(1e-6))),
        FractionalKdVAnsatz(Arb(-0.6)),
        FractionalKdVAnsatz(add_error(Arb(-0.6), Arb(1e-6))),
        FractionalKdVAnsatz(Arb(-0.3)),
        FractionalKdVAnsatz(add_error(Arb(-0.3), Arb(1e-6))),
    ]

    xs = exp.(range(log(Arb("1e-5")), log(Arb("5e-1")), 10))

    @testset "u0" begin
        for u0 in u0s
            for x in xs
                y1 = u0(x)
                y2 = u0(x, Asymptotic())
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
        end
    end

    @testset "F0(u0)" begin
        for u0 in u0s
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
    end

    @testset "inv_u0_normalised(u0)" begin
        for u0 in u0s
            f = x -> x^(-u0.α) / u0(x)
            g = HighestCuspedWave.inv_u0_normalised(u0)
            for x in xs
                y1 = f(x)
                y2 = g(x)
                @test isfinite(y1)
                @test isfinite(y2)
                @test Arblib.overlaps(y1, y2)
            end
        end
    end
end
