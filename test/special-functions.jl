@testset "special-functions" begin
    @testset "_sinc" begin
        # Test around removable singularity
        for ϵ in Arb(2) .^ (-50:5:0)
            x0 = Arb((-ϵ, ϵ))
            x = ArbSeries((x0, 1, 1))
            y1 = HighestCuspedWave._sinc(x)
            for xx0 in range(-ϵ, ϵ, 10)
                xx = ArbSeries((xx0, 1, 1))
                y2 = sinc(xx)
                @test Arblib.overlaps(y1, y2)
            end
        end
    end

    @testset "zeta_deflated" begin
        for a in range(Arb(0), 2, 10)[2:end]
            for s0 in [1; range(Arb(0), 2, 50)]
                s = ArbSeries((s0, 1, 1))
                z1 = Arblib.zeta_series!(zero(s), s, a, 1, length(s))
                z2 = HighestCuspedWave._zeta_deflated(s, a)
                @test Arblib.overlaps(z1, z2)
            end
        end

        # Test removable singularity
        for a in range(Arb(0), 2, 10)[2:end]
            for ϵ in Arb(2) .^ (-50:5:0)
                s0 = Arb((1 - ϵ, 1 + ϵ))
                s = ArbSeries((s0, 1, 1))
                z1 = HighestCuspedWave._zeta_deflated(s, a)
                for ss0 in range(1 - ϵ, 1 + ϵ, 10)
                    ss = ArbSeries((ss0, 1, 1))
                    z2 = Arblib.zeta_series!(zero(ss), ss, a, 1, length(ss))
                    @test Arblib.overlaps(z1, z2)
                end
            end
        end

        # Test indeterminate for a <= 0
        s = ArbSeries((0.5, 1))
        @test all(isnan, Arblib.coeffs(HighestCuspedWave._zeta_deflated(s, Arb(0))))
        @test all(isnan, Arblib.coeffs(HighestCuspedWave._zeta_deflated(s, Arb(-1))))
        @test all(isnan, Arblib.coeffs(HighestCuspedWave._zeta_deflated(s, Arb((-1, 1)))))

        # Check that the non-series version also works
        s = ArbSeries(0.5)
        a = Arb(0.5)
        @test isequal(
            HighestCuspedWave._zeta_deflated(s[0], a),
            HighestCuspedWave._zeta_deflated(s, a)[0],
        )

    end
end
