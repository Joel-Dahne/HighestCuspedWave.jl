@testset "expansion" begin
    α0s = [Arb(0), Arb((-1e-3))]

    intervals = [
        -Arblib.nonnegative_part!(zero(Arb), Arb((0, 1e-3))),
        Arblib.add_error!(Arb(-1e-3), Arb(1e-6)),
    ]

    for (α0, interval) in zip(α0s, intervals)
        αs = range(getinterval(Arb, interval)..., 100)
        if iszero(α0)
            αs = αs[1:end-1]
        end

        u0 = KdVZeroAnsatz(interval, α0)

        @testset "p0" begin
            for α in αs
                @test Arblib.overlaps(u0.p0(α), HighestCuspedWave.findp0(α))
            end
        end

        @testset "as" begin
            for α in αs
                @test Arblib.overlaps(u0.a[0](α), HighestCuspedWave.finda0(α))

                a1, a2 = HighestCuspedWave._finda1a2(α)
                @test Arblib.overlaps(u0.a[1](α), a1)
                @test Arblib.overlaps(u0.a[2](α), a2)
            end
        end
    end
end
