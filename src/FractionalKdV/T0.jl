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
        part1 = T01(u0, evaltype, δ0 = δ0, δ1 = δ1,
                    rtol = rtol, atol = atol, show_trace = show_trace)(x)

        if isnan(part1)
            # Short circuit on NaN
            return part1
        end

        ## Integral on [x, π]
        part2 = T02(u0, evaltype, δ2 = δ2, rtol = rtol, atol = atol, show_trace = show_trace)(x)

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
    T01(u0::FractionalKdVAnsatz; δ1, δ2)
Returns a function such that T01(u0, δ1 = δ1, δ2 = δ2)(x) computes the integral
T_{0,1} from the paper.
"""
T01(u0; kwargs...) = T01(u0, Ball(); kwargs...)

function T01(u0::FractionalKdVAnsatz{arb},
             evaltype::EvalType;
             δ0::arb = parent(u0.α)(1e-4),
             δ1::arb = parent(u0.α)(1e-4),
             rtol = -1.0,
             atol = -1.0,
             show_trace = false,
             )
    return x -> begin
        return (
            T011(u0, evaltype, δ0 = δ0)(x)
            + T012(
                u0,
                evaltype,
                δ0 = δ0,
                δ1 = δ1,
                rtol = rtol,
                atol = atol,
                show_trace = show_trace
            )(x)
            + T013(u0, evaltype, δ1 = δ1)(x)
        )
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
        (P, P_E) = taylor_with_error(zero(α), setunion(zero(α), δ0), N) do t
            Ci(x*(1 - t), -α) + Ci(x*(1 + t), -α)
        end
        P_restterm = ball(zero(α), P_E*δ0^N)

        # Singular term
        M = div(N, 2) + 1
        (C, e, P2, P2_E) = Ci_expansion(x*δ0, -α, M)
        C *= x^e
        for m = 1:M-1
            P2[2m] *= x^(2m)
        end
        P2_restterm = P2_E*(x*δ0)^(2M)

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        res -= 2C*δ0^(e + u0.p + 1)/(e + u0.p + 1)

        # Integrate the analytic terms
        full_series = P - 2P2
        for i = 0:N-1
            res += full_series[i]*δ0^(i + u0.p + 1)/(i + u0.p + 1)
        end

        # Add the error term
        res += δ0*(P_restterm - P2_restterm)

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

If `x` is equal or very close to π (determined by `ϵ`) then the Taylor
expansion gives a very poor approximation for `Ci(x*(t + 1), -α)`. In
this case we make use of the fact that it's 2π periodic and even, so
that `Ci(x*(t + 1), -α) = Ci(x*(t + 1) - 2π, -α) = Ci(2π - x*(t + 1),
-α)`, to be able to use the asymptotic expansion instead. That gives
us the integral ```∫_(1 - δ1)^1 (2π - x*(t + 1))^s*t^p dt` which is
given by `(2π - x)^(1 + p + s)*x^(-1 - p)*(B(1 + p, 1 + s, x/(2π - x))
- B(1 + p, 1 + s, (x - δ1*x)/(2π - x)))```. The value of `x/(2π - x)`
will always be less than or equal to 1 for `x` less than or equal to
π, however due to overestimation the enclosing ball might contain
values greater than one, we therefore have to use `beta_inc_zeroone`
to be able to get finite results in that case.
"""
T013(u0::FractionalKdVAnsatz{arb}; kwargs...) = T013(u0, Ball(); kwargs...)

function T013(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ1::arb = parent(u0.α)(1e-4),
              ϵ::arb = parent(u0.α)(1e-2),
              N::Integer = 3,
              )
    Γ = Nemo.gamma
    α = u0.α
    π = parent(α)(pi)

    PP = ArbPolyRing(parent(α), :x)

    return x -> begin
        # Determine if the asymptotic expansion or the Taylor
        # expansion should be used for the second term
        use_asymptotic = π - x < ϵ

        # Analytic terms
        (P, E) = taylor_with_error(one(α), setunion(1 - δ1, one(α)), N) do t
            if !use_asymptotic
                return Ci(x*(1 + t), -α) - 2Ci(x*t, -α)
            else
                return -2Ci(x*t, -α)
            end
        end
        P_restterm = ball(zero(α), E*δ1^N)

        # Singular term
        M = div(N, 2) + 1
        (C, e, P2, P2_E) = Ci_expansion(x*δ1, -α, M)
        C *= x^e
        for m = 1:M-1
            P2[2m] *= x^(2m)
        end
        P2_restterm = P2_E*(x*δ1)^(2M)

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        # Using ∫_(1 - δ1)^1 |t - 1|^s*t^(p) dt = ∫_(1 - δ1)^1 (1 - t)^s*t^(p) dt
        res += C*(Γ(1 + e)*Γ(1 + u0.p)/Γ(2 + e + u0.p) - beta_inc(1 + u0.p, 1 + e, 1 - δ1))

        # Integrate the analytic part
        # Using ∫_(1-δ1)^1 (t-1)^i*t^(p) dt = (-1)^i ∫_(1-δ1)^1 (t-1)^i*t^(p) dt
        full_series = P + P2
        for i = 0:N-1
            res += full_series[i] * (-1)^i * (
                Γ(parent(α)(1 + i))*Γ(1 + u0.p)/Γ(2 + i + u0.p)
                - beta_inc(1 + u0.p, parent(α)(1 + i), 1 - δ1)
            )
        end

        # Add the error term
        res += δ1*(P_restterm + P2_restterm)

        if use_asymptotic
            # Handle asymptotic expansion of Ci(x*(t + 1), -α)
            (C, e, P3, P3_E) = Ci_expansion(2π - x*(2 - δ1), -α, M)
            P3_restterm = P2_E*(2π - x*(2 - δ1))^(2M)

            # Add the singular part to the integral
            res += C*(2π - x)^(1 + u0.p + e)*x^(-1 - u0.p)*(
                beta_inc_zeroone(1 + u0.p, 1 + e, x/(2π - x))
                - beta_inc_zeroone(1 + u0.p, 1 + e, (x - δ1*x)/(2π - x))
            )

            for i = 0:2:N-1
                # Only even terms
                res += P3[i]*(2π - x)^(1 + u0.p + i)*x^(-1 - u0.p)*(
                beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), x/(2π - x))
                - beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), (x - δ1*x)/(2π - x))
                )
            end

            # Add error term
            res += δ1*P3_restterm
        end

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

If `x` is close to π (`π - x < ϵ`) then use only the asymptotic
expansion for the full integral.
"""
T02(u0; kwargs...) = T02(u0, Ball(); kwargs...)

function T02(u0::FractionalKdVAnsatz{arb},
             ::Ball;
             δ2::arb = parent(u0.α)(1e-4),
             rtol = -1.0,
             atol = -1.0,
             show_trace = false,
             ϵ::arb = parent(u0.α)(1e-1),
             )
    return x -> begin
        a = ArbTools.ubound(x + δ2)

        if parent(u0.α)(π) - x < ϵ
            # Use only the asymptotic expansion on the whole interval
            return T021(u0, Ball(), parent(u0.α)(π), x, ϵ = ϵ)
        end

        part1 = T021(u0, Ball(), a, x, ϵ = ϵ)

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

If `x` is equal or very close to π (determined by `ϵ`) then the Taylor
expansion gives a very poor approximation for `Ci(x + y, -α)`. In this
case we make use of the fact that it's 2π periodic and even, so that
`Ci(x + y, -α) = Ci(x + y - 2π, -α) = Ci(2π - (x + y), -α)`, to be
able to use the asymptotic expansion instead. That gives us the
integral ```∫_x^a (2π - (x + y))^s*t^p dy` which is given by `(2π -
x)^(1 + p + s)*(B(1 + p, 1 + s, a/(2π - x)) - B(1 + p, 1 + s, x/(2π -
x)))```. The value of `x/(2π - x)` will always be less than or equal
to 1 for `x` less than or equal to π, however due to overestimation
the enclosing ball might contain values greater than one, we therefore
have to use `beta_inc_zeroone` to be able to get finite results in
that case.
"""
T021(u0::FractionalKdVAnsatz{arb}, a, x; kwargs...) = T021(u0, Ball(), a, x; kwargs...)

function T021(u0::FractionalKdVAnsatz{arb},
              ::Ball,
              a::arb,
              x::arb;
              ϵ::arb = parent(u0.α)(1e-1),
              N::Integer = 3,
              )
    Γ = Nemo.gamma
    α = u0.α
    δ2 = a - x
    π = parent(α)(pi)

    PP = ArbPolyRing(parent(α), :x)

    # Determine if the asymptotic expansion or the Taylor
    # expansion should be used for the second term
    use_asymptotic = π - x < ϵ

    # Analytic terms
    (P, E) = taylor_with_error(x, setunion(x, a), N) do y
        if !use_asymptotic
            return Ci(x + y, -α) - 2Ci(y, -α)
        else
            return -2Ci(y, -α)
        end
    end
    P_restterm = ball(zero(α), E*δ2^N)

    # Singular term
    M = div(N, 2) + 1
    (C, e, P2, P2_E) = Ci_expansion(δ2, -α, M)
    P2_restterm = P2_E*(δ2)^(2M)

    # Compute the integral
    res = zero(α)
    # Integrate the singular term
    # Using ∫_x^a |x - y|^s*y^p dy = ∫_x^a (y - x)^s*y^p dy
    if u0.p == 1
        res += C*abspow(δ2, 1 + e)*(δ2/(2 + e) + x/(1 + e))
    else
        res += C*x^(1 + e + u0.p)*(
            Γ(1 + e)*Γ(-1 - e - u0.p)/Γ(-u0.p)
            - beta_inc_zeroone(-1 -e - u0.p, 1 + e, x/a)
        )
    end

    # Integrate the analytic terms
    full_series = P + P2
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
    res += δ2*(P_restterm + P2_restterm)

    if use_asymptotic
        # Handle asymptotic expansion of Ci(x + y, -α)
        # The furthest away from 2π we are is at y = x
        (C, e, P3, P3_E) = Ci_expansion(2π - 2x, -α, M)
        P3_restterm = P2_E*(2π - 2x)^(2M)

        # Add the singular part to the integral
        res += C*(2π - x)^(1 + e + u0.p)*(
            beta_inc_zeroone(1 + u0.p, 1 + e, a/(2π - x))
            - beta_inc_zeroone(1 + u0.p, 1 + e, x/(2π - x))
        )

        for i = 0:2:N-1
            # Only even terms
            res += P3[i]*(2π - x)^(1 + i + u0.p)*(
            beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), a/(2π - x))
            - beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), x/(2π - x))
        )
        end

        # Add error term
        res += δ2*P3_restterm
    end

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
