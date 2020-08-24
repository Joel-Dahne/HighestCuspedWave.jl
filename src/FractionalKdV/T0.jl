export T0

"""
    T0(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that T0(u0)(x) is the function whose supremum
on [0, π] gives C_B. The strategy for evaluation depends on type of
evaltype.
"""
T0(u0::FractionalKdVAnsatz; kwargs...) = T0(u0, Ball(); kwargs...)

function T0(u0::FractionalKdVAnsatz{arb},
            evaltype::Ball;
            rtol = -1.0,
            show_trace = false,
            )
    δ0 = parent(u0.α)("1e-3")
    δ1 = parent(u0.α)("1e-3")
    δ2 = parent(u0.α)("1e-6")

    return x -> begin
        ## Integral on [0, x] - Change to t = y/x

        T01 = (
            T011(u0, evaltype, δ0 = δ0)(x)
            + T012(u0, evaltype, δ0 = δ0, δ1 = δ1, rtol = rtol, show_trace = show_trace)(x)
            + T013(u0, evaltype, δ1 = δ1)(x)
        )

        if isnan(T01)
            # Short circuit on NaN
            return T01
        end

        ## Integral on [x, π]

        # Integral on [x, x + δ2]
        # TODO: Implement this
        T021 = zero(u0.α)/(parent(u0.α)(π)*u0.w(x)*u0(x))

        T02 = (
            T021
            + T022(u0, δ2 = δ2, rtol = rtol, show_trace = show_trace)(x)
        )

        return T01 + T02
    end
end

function T0(u0::FractionalKdVAnsatz{T}, evaltype::Asymptotic) where {T}
    error("not yet implemented")

    return x -> begin
        return zero(u0.α)
    end
end

"""
    T011(u0::FractionalKdVAnstaz{arb}, evaltype; δ0)
Returns a function such that T011(u0, δ0)(x) computes the integral
T_{0,1,1} from the paper.

Uses an expansion with `N` terms (of degree `N - 1`).
"""
T011(u0::FractionalKdVAnsatz{arb}; δ0::arb = parent(u0.α)("1e-6")) = T011(u0, Ball(), δ0 = δ0)

function T011(u0::FractionalKdVAnsatz{arb},
              evaltype::Ball;
              δ0::arb = parent(u0.α)("1e-6"),
              N::Integer = 9,
              )
    Γ = Nemo.gamma
    α = u0.α

    PP = ArbPolyRing(parent(α), :x)

    return x -> begin
        # Prove: that the expression inside the absolute value of the
        # integrand is negative

        # The first two terms are well behaved at t = 0, we compute their
        # Taylor expansions at the point
        t_series = arb_series(PP([0, 1]), N)
        part1_series = Ci(x*(1 - t_series), -α) + Ci(x*(1 + t_series), -α)

        t_series_restterm = arb_series(PP([setinterval(zero(α), δ0), one(α)]), N + 1)
        part1_series_restterm = Ci(x*(1 - t_series_restterm), -α) + Ci(x*(1 + t_series_restterm), -α)
        part1_restterm = ball(zero(α), δ0^N*part1_series_restterm[N])

        # The third term has a singularity at t = 0 so we expand and
        # handle the singular term by itself
        singular_exponent = -α - 1
        singular_coefficient = gamma(1 + α)*sinpi(-α/2)*abspow(x, singular_exponent)

        part2_series = arb_series(PP(), N)
        M = div(N, 2) + 1
        part2_series[0] = zeta(-α)
        for m = 1:M-1
            part2_series[2m] = (-1)^m*zeta(-α - 2m)/factorial(fmpz(2m))*abspow(x, 2m)
        end
        part2_restterm = ball(zero(α),
                              2(2parent(α)(π))^(1 - α - 2M)
                              *zeta(2M + 1 + α)*δ0^(2M)/(4parent(α)(π)^2 - abspow(x*δ0, 2))
                              *abspow(x*δ0, 2M),
                              )

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        res -= 2singular_coefficient*δ0^(singular_exponent + u0.p + 1)/(singular_exponent + u0.p + 1)

        # Integrate the terms in the series termwise on the interval [0, δ0]
        full_series = part1_series - 2part2_series
        for i = 0:N-1
            res += full_series[i]*δ0^(i + u0.p + 1)/(i + u0.p + 1)
        end

        # Add the error given by the restterm times the length of the interval
        res += δ0*(part1_restterm - 2part2_restterm)

        return -res*x/(parent(α)(π)*u0(x))
    end
end

"""
    T012(u0::FractionalKdVAnsatz{arb}, δ0, δ1, evaltype)
Returns a function such that T012(u0, δ0, δ1)(x) computes the integral
T_{0,1,2} from the paper.
"""
T012(u0::FractionalKdVAnsatz;
     δ0::arb = parent(u0.α)("1e-6"),
     δ1::arb = parent(u0.α)("1e-6"),
     rtol = -1.0,
     atol = -1.0,
     show_trace = false,
     ) = T012(u0, Ball(), δ0 = δ0, δ1 = δ1, rtol = rtol, atol = atol, show_trace = show_trace)

function T012(u0::FractionalKdVAnsatz{arb},
              evaltype::Ball;
              δ0::arb = parent(u0.α)("1e-6"),
              δ1::arb = parent(u0.α)("1e-6"),
              rtol = -1.0,
              atol = -1.0,
              show_trace = false,
              )
    return x -> begin
        CC = ComplexField(prec(parent(u0.α)))
        mα = CC(-u0.α)
        a = CC(δ0)
        b = CC(1 - δ1)

        # PROVE: That there are no branch cuts that interact with the
        # integral
        F(t, analytic = false) = ArbTools.real_abs(
            Ci(x*(1 - t), mα) + Ci(x*(1 + t), mα) - 2Ci(x*t, mα),
            analytic = analytic,
        )*t^u0.p
        res = real(ArbTools.integrate(CC, F, a, b,
                                      rel_tol = rtol,
                                      abs_tol = atol,
                                      eval_limit = 3000,
                                      verbose = Int(show_trace),
                                      ))

        return res*x/(parent(u0.α)(π)*u0(x))
    end
end

"""
    T013(u0::FractionalKdVAnstaz{arb}, evaltype; δ1)
Returns a function such that T013(u0, δ0)(x) computes the integral
T_{0,1,3} from the paper.

Uses an expansion with `N` terms (of degree `N - 1`).
"""
T013(u0::FractionalKdVAnsatz{arb}; δ1::arb = parent(u0.α)("1e-6")) = T013(u0, Ball(), δ1 = δ1)

function T013(u0::FractionalKdVAnsatz{arb},
              evaltype::Ball;
              δ1::arb = parent(u0.α)("1e-6"),
              N::Integer = 9,
              )
    Γ = Nemo.gamma
    α = u0.α

    PP = ArbPolyRing(parent(α), :x)

    return x -> begin
        # Prove: that the expression inside the absolute value of the
        # integrand is positive

        # The first two terms are well behaved at t = 1, we compute their
        # Taylor expansions at the point
        t_series = arb_series(PP([1, 1]), N)
        part1_series = Ci(x*(1 + t_series), -α) - 2Ci(x*t_series, -α)

        t_series_restterm = arb_series(PP([setinterval(1 - δ1, one(α)), one(α)]), N + 1)
        part1_series_restterm = Ci(x*(1 + t_series_restterm), -α) - 2Ci(x*t_series_restterm, -α)
        part1_restterm = ball(zero(α), δ1^N*part1_series_restterm[N])

        # The first term has a singularity at t = 1 so we expand and
        # handle the singular term by itself
        singular_exponent = -α - 1
        singular_coefficient = gamma(1 + α)*sinpi(-α/2)*abspow(x, singular_exponent)

        part2_series = arb_series(PP(), N)
        M = div(N, 2) + 1
        part2_series[0] = zeta(-α)
        for m = 1:M-1
            part2_series[2m] = (-1)^m*zeta(-α - 2m)/factorial(fmpz(2m))*abspow(x, 2m)
        end
        part2_restterm = ball(zero(α),
                              2(2parent(α)(π))^(1 - α - 2M)
                              *zeta(2M + 1 + α)*(x*(1 - δ1))^(2M)/(4parent(α)(π)^2 - (x*(1 - δ1))^2),
                              )

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        # We are integrating ∫_(1 - δ1)^1 |t - 1|^s*t^(p) dt
        # = ∫_(1 - δ1)^1 (1 - t)^s*t^(p) dt which is given by
        # Γ(1 + s)*Γ(1 + p)/Γ(2 + s + p) - B(1 + p, 1 + s; 1 - δ1)
        # Where B(a, b; z) is the incomplete Beta-function
        res += singular_coefficient*(
            Γ(1 + singular_exponent)*Γ(1 + u0.p)/Γ(2 + singular_exponent + u0.p)
            - beta_inc(1 + u0.p, 1 + singular_exponent, 1 - δ1)
        )

        # Integrate the terms in the series termwise on the interval [1 - δ1, 1]
        # We are integrating ∫_(1 - δ1)^1 (t - 1)^i*t^(p) dt which is given by
        # (-1)^i (Γ(1 + i)*Γ(1 + p)/Γ(2 + i + p) - B(1 + p, 1 + i; 1 - δ1))
        # Where B(a, b; z) is the incomplete Beta-function
        full_series = part1_series + part2_series
        for i = 0:N-1
            res += full_series[i] * (-1)^i * (
                Γ(parent(α)(1 + i))*Γ(1 + u0.p)/Γ(2 + i + u0.p)
                - beta_inc(1 + u0.p, parent(α)(1 + i), 1 - δ1)
            )
        end

        # Add the error given by the restterm times the length of the interval
        res += δ1*(part1_restterm + part2_restterm)


        return res*x/(parent(u0.α)(π)*u0(x))
    end
end

"""
    T022(u0::FractionalKdVAnsatz;
Returns a function such that T022(u0, δ0 = δ0)(x) computes the (not yet
existing) integral T_{0,2,2} from the paper. That is, the integral
from x + δ0 to π. If x + δ0 >= π the function returns zero, note that
the integral T02 might be far from zero in this case.
)
"""
T022(u0::FractionalKdVAnsatz;
     δ2::arb = parent(u0.α)("1e-6"),
     rtol = -1.0,
     atol = -1.0,
     show_trace = false,
     ) = T022(u0, Ball(), δ2 = δ2, rtol = rtol, atol = atol, show_trace = show_trace)

function T022(u0::FractionalKdVAnsatz{arb},
              evaltype::Ball;
              δ2::arb = parent(u0.α)("1e-6"),
              rtol = -1.0,
              atol = -1.0,
              show_trace = false,
              )
    return x -> begin
        if x + δ2 >= parent(u0.α)(π)
            return zero(u0.α)
        end

        CC = ComplexField(prec(parent(u0.α)))
        mα = CC(-u0.α)
        a = CC(x + δ2)
        b = CC(parent(u0.α)(π))

        # PROVE: That there are no branch cuts that interact with the
        # integral
        F(y, analytic = false) = ArbTools.real_abs(
            Ci(x - y, mα) + Ci(x + y, mα) - 2Ci(y, mα),
            analytic = analytic,
        )*y^u0.p
        res = real(ArbTools.integrate(CC, F, a, b,
                                      rel_tol = rtol,
                                      abs_tol = atol,
                                      eval_limit = 2000,
                                      verbose = Int(show_trace),
                                      ))

        return res/(parent(u0.α)(π)*u0.w(x)*u0(x))
    end
end
