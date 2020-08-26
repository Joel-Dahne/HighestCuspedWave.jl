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
            atol = -1.0,
            show_trace = false,
            )
    δ0 = parent(u0.α)(1e-4)
    δ1 = parent(u0.α)(1e-4)
    δ2 = parent(u0.α)(1e-4)

    return x -> begin
        ## Integral on [0, x] - Change to t = y/x
        part1 = (
            T011(u0, evaltype, δ0 = δ0)(x)
            + T012(u0, evaltype, δ0 = δ0, δ1 = δ1,
                   rtol = rtol, atol = atol, show_trace = show_trace)(x)
            + T013(u0, evaltype, δ1 = δ1)(x)
        )

        if isnan(part1)
            # Short circuit on NaN
            return part1
        end

        ## Integral on [x, π]
        part2 = T02(u0, δ2 = δ2, rtol = rtol, atol = atol, show_trace = show_trace)(x)

        return part1 + part2
    end
end

function T0(u0::FractionalKdVAnsatz{T}, ::Asymptotic) where {T}
    @warn "T0(u0) is not yet implemented"

    return x -> begin
        return zero(u0.α)
    end
end

"""
    T011(u0::FractionalKdVAnstaz{arb}; δ0)
Returns a function such that T011(u0, δ0 = δ0)(x) computes the
integral T_{0,1,1} from the paper.

The strategy for evaluation is to compute an expansion for the
integrand which is then integrated termwise on the interval. The first
two terms inside the absolute value are analytic and their Taylor
expansions of degree `N - 1` are computed at `t = 0` and the error
term is enclose. The last term inside the absolute value is not
analytic and instead we use the expansion for Clausians from the paper
and also here enclose the error term.

The value inside the absolute value has constant sign so we can remove
it. Switching integration and summation gives us terms of the form
`∫_0^δ0 t^s*t^p dt` where `s` depends on the term and `p = u0.p`. This
is easily calculated to be `δ0^(s + p + 1*/(s + p + 1)`. The error
handled as constant values which are just multiplied by the length of
the interval.
"""
T011(u0::FractionalKdVAnsatz{arb}; kwargs...) = T011(u0, Ball(); kwargs...)

function T011(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ0::arb = parent(u0.α)(1e-4),
              N::Integer = 3,
              )
    Γ = Nemo.gamma
    α = u0.α

    PP = ArbPolyRing(parent(α), :x)

    return x -> begin
        # Analytic terms
        t_series = arb_series(PP([0, 1]), N)
        part1_series = Ci(x*(1 - t_series), -α) + Ci(x*(1 + t_series), -α)

        t_series_restterm = arb_series(PP([setinterval(zero(α), δ0), one(α)]), N + 1)
        part1_series_restterm = Ci(x*(1 - t_series_restterm), -α) + Ci(x*(1 + t_series_restterm), -α)
        part1_restterm = ball(zero(α), δ0^N*part1_series_restterm[N])

        # Singular term
        singular_exponent = -α - 1
        singular_coefficient = gamma(1 + α)*sinpi(-α/2)*abspow(x, singular_exponent)

        part2_series = arb_series(PP(), N)
        M = div(N, 2) + 1
        part2_series[0] = zeta(-α)
        for m = 1:M-1
            part2_series[2m] = (-1)^m*zeta(-α - 2m)/factorial(fmpz(2m))*x^(2m)
        end
        part2_restterm = ball(zero(α),
                              2(2parent(α)(π))^(1 - α - 2M)
                              *zeta(2M + 1 + α)*δ0^(2M)/(4parent(α)(π)^2 - (x*δ0)^2)
                              *(x*δ0)^(2M),
                              )

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        res -= 2singular_coefficient*δ0^(singular_exponent + u0.p + 1)/(singular_exponent + u0.p + 1)

        # Integrate the analytic terms
        full_series = part1_series - 2part2_series
        for i = 0:N-1
            res += full_series[i]*δ0^(i + u0.p + 1)/(i + u0.p + 1)
        end

        # Add the error term
        res += δ0*(part1_restterm - 2part2_restterm)

        # Prove: that the expression inside the absolute value of the
        # integrand is negative
        return -res*x/(parent(α)(π)*u0(x))
    end
end

"""
    T012(u0::FractionalKdVAnsatz{arb}; δ0, δ1)
Returns a function such that T012(u0, δ0 = δ0, δ1 = δ1)(x) computes
the integral T_{0,1,2} from the paper.
"""
T012(u0::FractionalKdVAnsatz{arb}; kwargs...) = T012(u0, Ball(); kwargs...)

function T012(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ0::arb = parent(u0.α)(1e-4),
              δ1::arb = parent(u0.α)(1e-4),
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
    T013(u0::FractionalKdVAnstaz{arb}; δ1)
Returns a function such that T013(u0, δ0 = δ0)(x) computes the integral
T_{0,1,3} from the paper.

The strategy for evaluation is the same as for T011 except that the
first term is singular and the last two are analytic and their Taylor
expansion is computed at `t = 1`.

The integral that needs to be computed in this case is `∫_(1 - δ1)^1
(1 - t)^s*t^p dt` which is given by `Γ(1 + s)*Γ(1 + p)/Γ(2 + s + p) -
B(1 + p, 1 + s; 1 - δ1)` where B(a, b; z) is the incomplete
Beta-function.
"""
T013(u0::FractionalKdVAnsatz{arb}; kwargs...) = T013(u0, Ball(); kwargs...)

function T013(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ1::arb = parent(u0.α)(1e-4),
              N::Integer = 3,
              )
    Γ = Nemo.gamma
    α = u0.α

    PP = ArbPolyRing(parent(α), :x)

    return x -> begin
        # Analytic terms
        t_series = arb_series(PP([1, 1]), N)
        part1_series = Ci(x*(1 + t_series), -α) - 2Ci(x*t_series, -α)

        t_series_restterm = arb_series(PP([setinterval(1 - δ1, one(α)), one(α)]), N + 1)
        part1_series_restterm = Ci(x*(1 + t_series_restterm), -α) - 2Ci(x*t_series_restterm, -α)
        part1_restterm = ball(zero(α), δ1^N*part1_series_restterm[N])

        # Singular terms
        singular_exponent = -α - 1
        singular_coefficient = gamma(1 + α)*sinpi(-α/2)*abspow(x, singular_exponent)

        part2_series = arb_series(PP(), N)
        M = div(N, 2) + 1
        part2_series[0] = zeta(-α)
        for m = 1:M-1
            part2_series[2m] = (-1)^m*zeta(-α - 2m)/factorial(fmpz(2m))*x^(2m)
        end
        part2_restterm = ball(zero(α),
                              2(2parent(α)(π))^(1 - α - 2M)
                              *zeta(2M + 1 + α)*(x*(1 - δ1))^(2M)/(4parent(α)(π)^2 - (x*(1 - δ1))^2),
                              )

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        # Using ∫_(1 - δ1)^1 |t - 1|^s*t^(p) dt = ∫_(1 - δ1)^1 (1 - t)^s*t^(p) dt
        res += singular_coefficient*(
            Γ(1 + singular_exponent)*Γ(1 + u0.p)/Γ(2 + singular_exponent + u0.p)
            - beta_inc(1 + u0.p, 1 + singular_exponent, 1 - δ1)
        )

        # Integrate the analytic part
        # Using ∫_(1-δ1)^1 (t-1)^i*t^(p) dt = (-1)^i ∫_(1-δ1)^1 (t-1)^i*t^(p) dt
        full_series = part1_series + part2_series
        for i = 0:N-1
            res += full_series[i] * (-1)^i * (
                Γ(parent(α)(1 + i))*Γ(1 + u0.p)/Γ(2 + i + u0.p)
                - beta_inc(1 + u0.p, parent(α)(1 + i), 1 - δ1)
            )
        end

        # Add the error term
        res += δ1*(part1_restterm + part2_restterm)

        # Prove: that the expression inside the absolute value of the
        # integrand is positive
        return res*x/(parent(u0.α)(π)*u0(x))
    end
end

"""
    T02(u0::FractionalKdVAnsatz; δ2)
Returns a function such that T02(u0, δ2 = δ2)(x) computes the integral
T_{0,2} from the paper.

If `x + δ2 > π` the function shows a warning and returns zero.

The split between T021 and T022 is not done exactly at x + δ2 but at
an exact floating point value slightly above this.
"""
T02(u0; kwargs...) = T02(u0, Ball(); kwargs...)

function T02(u0::FractionalKdVAnsatz{arb},
             ::Ball;
             δ2::arb = parent(u0.α)(1e-4),
             rtol = -1.0,
             atol = -1.0,
             show_trace = false,
             )
    return x -> begin
        a = ArbTools.ubound(x + δ2)

        if !(a < parent(u0.α)(π))
            @warn "Evaluating T02 to close to x = π - NOT rigorous"
            return zero(u0.α)
        end

        part1 = T021(u0, Ball(), a, x)

        part2 = T022(u0, Ball(), a, x, rtol = rtol, atol = atol, show_trace = show_trace)

        return part1 + part2
    end
end

"""
    T021(u0::FractionalKdVAnstaz{arb}, a::arb, x::arb)
Computes the (not yet existing) integral T_{0,2,1} from the paper.

The strategy for evaluation is the same as for T011 except that the
first term is singular and the last two are analytic and their Taylor
expansion is computed at `t = x`.

The integral that needs to be computed in this case is `∫_x^a (y -
x)^s*y^p dy` which is given by `x^(s + p + 1)*(Γ(1 + s)*Γ(-1 - s -
p)/Γ(-p) - B(-1 - s - p, 1 + s; x/a)` where B(a, b; z) is the
incomplete Beta-function. For `p = 1` we instead use the expression
`(a - x)^(1+s)*((a - x)/(2 + s) + x/(1 + s))`.
"""
T021(u0::FractionalKdVAnsatz{arb}, a, x) = T021(u0, Ball(), a, x)

function T021(u0::FractionalKdVAnsatz{arb},
              ::Ball,
              a::arb,
              x::arb;
              N::Integer = 3,
              )
    Γ = Nemo.gamma
    α = u0.α
    δ2 = a - x

    PP = ArbPolyRing(parent(α), :x)

    # Analytic terms
    y_series = arb_series(PP([x, one(α)]), N)
    part1_series = Ci(x + y_series, -α) - 2Ci(y_series, -α)

    y_series_restterm = arb_series(PP([setinterval(x, a), one(α)]), N + 1)
    part1_series_restterm = Ci(x + y_series_restterm, -α) - 2Ci(y_series_restterm, -α)
    part1_restterm = ball(zero(α), δ2^N*part1_series_restterm[N])

    # Singular term
    singular_exponent = -α - 1
    singular_coefficient = gamma(1 + α)*sinpi(-α/2)

    part2_series = arb_series(PP(), N)
    M = div(N, 2) + 1
    part2_series[0] = zeta(-α)
    for m = 1:M-1
        part2_series[2m] = (-1)^m*zeta(-α - 2m)/factorial(fmpz(2m))
    end
    part2_restterm = ball(zero(α),
                          2(2parent(α)(π))^(1 - α - 2M)
                          *zeta(2M + 1 + α)*δ2^(2M)/(4parent(α)(π)^2 - δ2^2),
                          )

    # Compute the integral
    res = zero(α)
    # Integrate the singular term
    # Using ∫_x^a |x - y|^s*y^p dy = ∫_x^a (y - x)^s*y^p dy
    if u0.p == 1
        res += singular_coefficient*δ2^(1+singular_exponent)*(
            δ2/(2+singular_exponent)+x/(1+singular_exponent)
        )
    else
        res += singular_coefficient*x^(1 + singular_exponent + u0.p)*(
            Γ(1 + singular_exponent)*Γ(-1 - singular_exponent - u0.p)/Γ(-u0.p)
            - beta_inc(-1 -singular_exponent - u0.p, 1 + singular_exponent, x/a)
        )
    end

    # Integrate the analytic terms
    full_series = part1_series + part2_series
    for i = 0:N-1
        if u0.p == 1
            res += full_series[i]*δ2^(1 + i)*(
                δ2/(2 + i) + x/(1 + i)
            )
        else
            res += full_series[i] * x^(1 + i + u0.p) * (
                Γ(parent(α)(1 + i))*Γ(-1 -i - u0.p)/Γ(-u0.p)
                - beta_inc(-1 - i - u0.p, parent(α)(1 + i), x/a)
            )
        end
    end

    # Add the error term
    res += δ2*(part1_restterm + part2_restterm)

    # Prove: that the expression inside the absolute value of the
    # integrand is positive
    return res/(parent(u0.α)(π)*u0.w(x)*u0(x))
end

"""
    T022(u0::FractionalKdVAnsatz{arb}, a::arb, x::arb)
Computes the (not yet existing) integral T_{0,2,2} from the paper.
)
"""
T022(u0::FractionalKdVAnsatz{arb}, a, x; kwargs...) = T022(u0, Ball(), a, x; kwargs...)

function T022(u0::FractionalKdVAnsatz{arb},
              ::Ball,
              a::arb,
              x::arb;
              rtol = -1.0,
              atol = -1.0,
              show_trace = false,
              )
    CC = ComplexField(prec(parent(u0.α)))
    mα = CC(-u0.α)
    a = CC(a)
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
