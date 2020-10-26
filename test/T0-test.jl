@testset "T0" begin
    @testset "T01 asymptotic gives bound" begin
        RR = RealField(64)
        for α in range(-0.8, stop = -0.2, length = 6)
            u0 = FractionalKdVAnsatz(RR(α))
            f = HighestCuspedWave.T01(u0, Ball())
            g = HighestCuspedWave.T01(u0, Asymptotic(), nonasymptotic_u0 = true)
            h = HighestCuspedWave.T01(u0, Asymptotic(), nonasymptotic_u0 = false)
            for x in [1e-6; 1e-4; collect(range(0, stop = π, length = 10)[2:end])]
                x = RR(x)
                y1 = f(x)
                y2 = g(x)
                y3 = h(x)
                @test !(y1 > y2)
                @test !(y1 > y3)
                @test overlaps(y2, y3)
            end
        end
    end

    @testset "T02 asymptotic gives bound" begin
        RR = RealField(64)
        for α in range(-0.8, stop = -0.2, length = 6)
            u0 = FractionalKdVAnsatz(RR(α))
            f = HighestCuspedWave.T02(u0, Ball(), δ2 = RR(1e-6))
            g = HighestCuspedWave.T02(u0, Asymptotic(), nonasymptotic_u0 = true)
            h = HighestCuspedWave.T02(u0, Asymptotic(), nonasymptotic_u0 = false)
            for x in [1e-6; 1e-4; collect(range(0, stop = π, length = 10)[2:end])]
                x = RR(x)
                y1 = f(x)
                y2 = g(x)
                y3 = h(x)
                @test !(y1 > y2)
                @test !(y1 > y3)
                @test overlaps(y2, y3)
            end
        end
    end
end
