@testset "Clausen functions" begin
    @testset "polylog" begin
        # TODO: Add tests for wide values of s
    end

    @testset "clausenc" begin
        # Check that the evaluation with polylog and zeta agree on (0, π)
        for s in range(Arb(1), 5, length = 10)[2:end-1]
            for x in range(Arb(0), π, length = 100)[2:end-1]
                res1 = HighestCuspedWave._clausenc_polylog(x, s)
                res2 = HighestCuspedWave._clausenc_zeta(x, s)
                @test Arblib.overlaps(res1, res2)
            end
        end

        # Check that _clausen_zeta throws an error outside of the domain
        @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(-1), Arb(2.5))
        @test_throws DomainError HighestCuspedWave._clausenc_zeta(Arb(4), Arb(2.5))

        # Check evaluation with integer s
        @test isfinite(clausenc(Arb(1), 2))
        @test isfinite(clausenc(Arb(1), 3))
        @test isfinite(clausenc(Arb(1), Arb(2)))
        @test isfinite(clausenc(Arb(1), Arb(3)))

        # Check evaluation on wide intervals
        let interval = Arb((1, 2)), s = Arb(2.5)
            y = clausenc(interval, s)
            for x in range(getinterval(Arb, interval)..., length = 100)
                @test Arblib.overlaps(clausenc(x, s), y)
            end
        end

        # Check evaluation with other types
        @test clausenc(1.5, 2) == Float64(clausenc(Arb(1.5), 2))
        @test clausenc(2, 2) == Float64(clausenc(Arb(2), 2))

        # TODO: Add tests for wide values of s
    end

    @testset "clausens" begin
        # Check that the evaluation with polylog and zeta agree on (0, π)
        for s in range(Arb(1), 5, length = 10)[2:end-1]
            for x in range(Arb(0), π, length = 100)[2:end-1]
                res1 = HighestCuspedWave._clausens_polylog(x, s)
                res2 = HighestCuspedWave._clausens_zeta(x, s)
                @test Arblib.overlaps(res1, res2)
            end
        end

        # Check that _clausen_zeta throws an error outside of the domain
        @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(-1), Arb(2.5))
        @test_throws DomainError HighestCuspedWave._clausens_zeta(Arb(4), Arb(2.5))

        # Check evaluation with integer s
        @test isfinite(clausens(Arb(1), 2))
        @test isfinite(clausens(Arb(1), 3))
        @test isfinite(clausens(Arb(1), Arb(2)))
        @test isfinite(clausens(Arb(1), Arb(3)))

        # Check evaluation on wide intervals
        let interval = Arb((1, 2)), s = Arb(2.5)
            y = clausens(interval, s)
            for x in range(getinterval(Arb, interval)..., length = 100)
                @test Arblib.overlaps(clausens(x, s), y)
            end
        end

        # Check evaluation with other types
        @test clausens(1.5, 2) == Float64(clausens(Arb(1.5), 2))
        @test clausens(2, 2) == Float64(clausens(Arb(2), 2))

        # TODO: Add tests for wide values of s
    end
end
