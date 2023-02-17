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

    @testset "x::Acb" begin
        for f in [sin, tan, x -> x^2]
            for radius in Arb(2) .^ (-10:1)
                x = Arblib.add_error!(zero(Acb), radius)
                if isfinite(f(AcbSeries((x, 1)))) # No point to test otherwise
                    for extra_degree = 0:2
                        enclosure = HighestCuspedWave.fx_div_x(f, x; extra_degree)
                        @test isfinite(enclosure)
                        y_reals = range(getinterval(Arb, real(x))..., 100)
                        y_imags = range(getinterval(Arb, imag(x))..., 100)
                        for y_real in y_reals
                            y1 = Acb(y_real, y_imags[1])
                            y2 = Acb(y_real, y_imags[end])
                            @test Arblib.overlaps(enclosure, f(y1) / y1)
                            @test Arblib.overlaps(enclosure, f(y2) / y2)
                        end
                        for y_imag in y_imags
                            y1 = Acb(y_reals[1], y_imag)
                            y2 = Acb(y_reals[end], y_imag)
                            @test Arblib.overlaps(enclosure, f(y1) / y1)
                            @test Arblib.overlaps(enclosure, f(y2) / y2)
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
