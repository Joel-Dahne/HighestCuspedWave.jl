@testset "fx_div_x" begin
    @testset "x::Arb" begin
        for f in [sin, HighestCuspedWave.rgamma, tan, x -> x^2]
            for radius in Arb(2) .^ (-10:1)
                x = Arblib.add_error!(zero(Arb), radius)
                if isfinite(f(ArbSeries((x, 1)))) # No point to test otherwise
                    for extra_degree = 0:2
                        enclosure = HighestCuspedWave.fx_div_x(f, x; extra_degree)
                        @test isfinite(enclosure)
                        for y in range(getinterval(Arb, x)..., 100)
                            @test Arblib.overlaps(enclosure, f(y) / y)
                        end
                    end
                end
            end
        end
    end

    @testset "x::ArbSeries" begin
        # Pure derivative
        for f in [sin, HighestCuspedWave.rgamma, tan, x -> x^2]
            for radius in Arb(2) .^ (-10:1)
                for degree = 0:4
                    x = ArbSeries((Arblib.add_error!(zero(Arb), radius), 1); degree)
                    if isfinite(f(ArbSeries((x[0], 1)))) # No point to test otherwise
                        for extra_degree = 0:2
                            enclosure = HighestCuspedWave.fx_div_x(f, x; extra_degree)
                            @test isfinite(enclosure)
                            for y in range(getinterval(Arb, x[0])..., 100)
                                y_s = ArbSeries(x)
                                y_s[0] = y
                                @test Arblib.overlaps(enclosure, f(y_s) / y_s)
                            end
                        end
                    end
                end
            end
        end

        # General series
        for f in [sin, HighestCuspedWave.rgamma, tan, x -> x^2]
            for radius in Arb(2) .^ (-10:1)
                for degree = 0:4
                    x = ArbSeries((Arblib.add_error!(zero(Arb), radius), 2, 3); degree)
                    if isfinite(f(ArbSeries((x[0], 1)))) # No point to test otherwise
                        for extra_degree = 0:2
                            enclosure = HighestCuspedWave.fx_div_x(f, x; extra_degree)
                            @test isfinite(enclosure)
                            for y in range(getinterval(Arb, x[0])..., 100)
                                y_s = ArbSeries(x)
                                y_s[0] = y
                                @test Arblib.overlaps(enclosure, f(y_s) / y_s)
                            end
                        end
                    end
                end
            end
        end
    end
end
