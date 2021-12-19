@testset "expansion" begin
    # Set up interval of α to test the expansions on and take a number
    # of points on this interval
    α_lower = Arb(-0.01)
    α_upper = Arb(0)
    α_interval = Arb((α_lower, α_upper))
    αs = range(α_lower, α_upper, length = 100)[1:end-1]

    # Construct KdVZeroansatz for the interval
    u0 = KdVZeroAnsatz(α_interval)

    @testset "p0" begin
        for α in αs
            @test Arblib.overlaps(u0.p0(α), HighestCuspedWave.findp0(α))
        end
    end

    # Construct FractionalKdvansatz for each αs
    #u0s = FractionalKdVAnsatz.(αs, 2, 0)

    @testset "as" begin
        for α in αs
            @test Arblib.overlaps(u0.a[0](α), HighestCuspedWave.finda0(α))

            a1, a2 = HighestCuspedWave._finda1a2(α)
            @test Arblib.overlaps(u0.a[1](α), a1)
            @test Arblib.overlaps(u0.a[2](α), a2)
        end
    end
end
