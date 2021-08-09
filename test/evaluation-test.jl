@testset "FractionalKdV asymptotic u0" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 10)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 3, 3)
        for x in range(-π, stop = π, length = 100)
            x = RR(x)
            y1 = u0(x)
            y2 = u0(x, HighestCuspedWave.Asymptotic())
            @test overlaps(y1, y2)
        end
    end
end

@testset "FractionalKdV asymptotic H(u0)" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 10)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 3, 3)
        for x in range(-π, stop = π, length = 100)
            x = RR(x)
            y1 = HighestCuspedWave.H(u0)(x)
            y2 = HighestCuspedWave.H(u0, HighestCuspedWave.Asymptotic())(x)
            @test overlaps(y1, y2)
        end
    end
end

@testset "FractionalKdV asymptotic D(u0)" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 10)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 3, 3)
        for x in range(-π, stop = π, length = 100)
            x = RR(x)
            y1 = HighestCuspedWave.D(u0)(x)
            y2 = HighestCuspedWave.D(u0, HighestCuspedWave.Asymptotic())(x)
            @test overlaps(y1, y2)
        end
    end
end

@testset "FractionalKdV asymptotic F(u0)" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 10)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 3, 3)
        for x in range(-0.1, stop = 0.1, length = 100)
            x = RR(x)
            y1 = HighestCuspedWave.F0(u0)(x)
            y2 = HighestCuspedWave.F0(u0, HighestCuspedWave.Asymptotic())(x)
            @test overlaps(y1, y2)
        end
    end
end

@testset "FractionalKdV asymptotic expansion u0" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 6)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 3, 3)
        for x in range(-π, stop = π, length = 50)
            x = RR(x)
            y1 = u0(x, HighestCuspedWave.Asymptotic())
            y2 = HighestCuspedWave.eval_expansion(
                u0,
                u0(x, HighestCuspedWave.AsymptoticExpansion()),
                x,
            )
            @test overlaps(y1, y2)
        end
    end
end

@testset "FractionalKdV asymptotic expansion H(u0)" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 6)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 3, 3)
        for x in range(-π, stop = π, length = 50)
            x = RR(x)
            y1 = HighestCuspedWave.H(u0, HighestCuspedWave.Asymptotic())(x)
            y2 = HighestCuspedWave.eval_expansion(
                u0,
                HighestCuspedWave.H(u0, HighestCuspedWave.AsymptoticExpansion())(x),
                x,
            )
            @test overlaps(y1, y2)
        end
    end
end

@testset "FractionalKdV asymptotic expansion D(u0)" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 6)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 3, 3)
        for x in range(-π, stop = π, length = 50)
            x = RR(x)
            y1 = HighestCuspedWave.D(u0, HighestCuspedWave.Asymptotic())(x)
            y2 = HighestCuspedWave.eval_expansion(
                u0,
                HighestCuspedWave.D(u0, HighestCuspedWave.AsymptoticExpansion())(x),
                x,
            )
            @test overlaps(y1, y2)
        end
    end
end

@testset "FractionalKdV bound hat(u0)" begin
    RR = RealField(64)
    for α in range(-1, stop = 0, length = 10)[2:end-1]
        u0 = FractionalKdVAnsatz(RR(α), 2, 3)
        for ϵ in RR.([1e-2, 1e-1])
            c = HighestCuspedWave.c(u0, ϵ)
            for x in range(0, stop = Float64(ϵ), length = 50)
                x = RR(x)
                y1 = HighestCuspedWave.hat(u0)(x)
                y2 = c * x^u0.p0
                # Due to y1 not being evaluated exactly we might not
                # have y1 <= y2. Instead we make sure that y1 is never
                # strictly greater than y2
                @test !(y1 > y2)
            end
        end
    end
end
