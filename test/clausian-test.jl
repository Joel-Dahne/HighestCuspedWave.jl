@testset "Ci" begin
    RR = RealField(64)
    PP = ArbPolyRing(RR, :x)
    N = 8
    s = RR(0.5)

    for f in [x -> Ci(x, s), x -> Ci(sin(x), s)]
        for a in range(0, stop = π, length = 5)[2:end]
            a = RR(a)
            poly_a = arb_series(PP([a, one(a)]), N)
            poly1 = f(poly_a)
            for x in range(0, stop = π, length = 25)[2:end]
                x = RR(x)

                # Real value
                y = f(x)

                # Approximation
                y_approx = evaluate(poly1.poly, x - a)

                # Error term
                poly_ball = arb_series(PP([setunion(x, a), one(x)]), N + 1)
                poly2 = f(poly_ball)
                restterm = abs(x - a)^N*coeff(poly2.poly, N)

                @test overlaps(y_approx + restterm, y)
            end
        end
    end
end

@testset "Ci_expansion" begin
    RR = RealField(64)
    s = RR(0.5)
    for M in [3, 7]
        for x in range(-π, stop = π, length = 100)
            x = RR(x)
            (C, e, P, E) = HighestCuspedWave.Ci_expansion(x, s, M)
            @test overlaps(C*abs(x)^e + evaluate(P.poly, x) + E*x^(2M), Ci(x, s))
        end
    end
end

@testset "Si_expansion" begin
    RR = RealField(64)
    s = RR(0.5)
    for M in [3, 7]
        for x in range(-π, stop = π, length = 100)
            x = RR(x)
            (C, e, P, E) = HighestCuspedWave.Si_expansion(x, s, M)
            @test overlaps(C*ifelse(x < 0, -1, 1)*abs(x)^e + evaluate(P.poly, x) + E*x^(2M + 1), Si(x, s))
        end
    end
end
