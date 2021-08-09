"""
    T01(u0::FractionalKdVAnsatz, ::Ball; δ1, δ2)
Returns a function such that `T01(u0, Ball(); δ1, δ2)(x)` computes the
integral T_{0,1} from the paper.

The integral is split into three parts depending on `δ1` and `δ2`.
These are given by `T011`, `T012` and `T013`.
"""
function T01(
    u0::FractionalKdVAnsatz{arb},
    evaltype::Ball;
    δ0::arb = parent(u0.α)(1e-2),
    δ1::arb = parent(u0.α)(1e-2),
    rtol = -1.0,
    atol = -1.0,
    show_trace = false,
)
    f = T011(u0, evaltype; δ0)
    g = T012(u0, evaltype; δ0, δ1, rtol, atol, show_trace)
    h = T013(u0, evaltype; δ1)
    return x -> begin
        return f(x) + g(x) + h(x)
    end
end

function T01(
    u0::FractionalKdVAnsatz{Arb},
    evaltype::Ball;
    δ0::Arf = Arf(1e-2),
    δ1::Arf = Arf(1e-2),
    rtol = 0,
    atol = 0,
    show_trace = false,
)
    f = x -> zero(x) # T011(u0, evaltype; δ0)
    g = T012(u0, evaltype; δ0, δ1, rtol, atol, show_trace)
    h = x -> zero(x) # T013(u0, evaltype; δ1)
    return x -> begin
        return f(x) + g(x) + h(x)
    end
end

"""
    T01(u0::FractionalKdVAnsatz, ::Asymptotic)
Returns a function such that `T01(u0, Asymptotic())(x)` computes an
**upper bound** of the integral T_{0,1} from the paper using an
evaluation strategy that works asymptotically as `x` goes to 0.

It splits the Clausians into the main singular part and the analytic
expansion. The integral of the singular part is computed by finding
where the integrand is positive respectively negative and then
integrating explicitly. The expansion is integrated term wise and the
resulting sum is bounded.

FIXME: There is something wrong somewhere. The result doesn't match up
completely with the non-asymptotic version for small values of `x`. I
don't know which one is the correct one...
```
u0 = FractionalKdVAnsatz(RR(-0.6), pp = RR(1))
x = RR(1e-6)
a = HighestCuspedWave.T01(u0, Ball())(x)
b = HighestCuspedWave.T01(u0, Asymptotic())(x)
overlaps(a, b)
```
"""
function T01(
    u0::FractionalKdVAnsatz{arb},
    ::Asymptotic;
    nonasymptotic_u0 = false, # Mainly for testing
)
    @warn "T01(u0, Asymptotic()) is not yet complete - the tail of the sum is not bounded"
    Γ = Nemo.gamma
    α = u0.α
    p = u0.p
    RR = parent(α)
    CC = ComplexField(precision(RR))
    π = RR(pi)

    # Compute
    # c_α = |Γ(1 + α)*sin(πα/2)|∫0^1 |(1 - t)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)|tᵖ dt
    # Find the unique zero of the integrand on [0, 1]
    # PROVE: That there is at most one zero on [0, 1]
    f = t -> (1 - t)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)
    roots, flags =
        isolateroots(f, parent(α)(0.1), parent(α)(0.9), refine = true, evaltype = :taylor)
    @assert only(flags)
    s = setunion(only(roots)...)

    c_α =
        abs(Γ(1 + α) * sinpi(α / 2)) * (
            (
                hypgeom_2f1(1 + p, 1 + α, 2 + p, -one(α)) -
                2s^(1 + p) * hypgeom_2f1(1 + p, 1 + α, 2 + p, -s)
            ) / (1 + p) - 2beta_inc(1 + p, -α, s) +
            (2 - 4s^(p - α)) / (α - p) +
            Γ(-α) * Γ(1 + p) / Γ(1 - α + p)
        )

    return x -> begin
        if p == 1
            # TODO: Bound the tail - the terms go to zero extremely
            # fast so it should be negligible
            c_ϵ = zero(α)
            for m = 1:10
                m = fmpz(m)
                c_ϵ +=
                    (-one(α))^m / (factorial(2m)) * zeta(-α - 2m) * m * (RR(4)^m - 1) /
                    ((m + 1) * (2m + 1)) * abspow(x, RR(2m - 2))
            end
            c_ϵ *= 2
        else
            # TODO: Bound the tail - the terms go to zero extremely
            # fast so it should be negligible
            c_ϵ = zero(α)
            for m = 1:10
                m = fmpz(m)
                c_ϵ +=
                    (-one(α))^m / (factorial(2m)) *
                    zeta(-α - 2m) *
                    (sum(binom(RR(2m), unsigned(2k)) / (2k + p + 1) for k = 0:Int(m)-1)) *
                    abspow(x, RR(2m - 2))
            end
            c_ϵ *= 2
        end

        if nonasymptotic_u0
            # Version without asymptotically expanding u0(x)
            res = c_α * abspow(x, -α) + c_ϵ * abspow(x, 3)
            return res / (π * u0(x))
        end

        res = c_α + c_ϵ * abspow(x, 3 + α)
        # Ball containing 1 + hat(u0)(x)
        L = ball(parent(α)(1), c(u0, ArbTools.abs_ubound(x)) * abspow(x, u0.p0))

        return L / (π * a0(u0, 0)) * res
    end
end

"""
    T011(u0::FractionalKdVAnstaz{arb}, evaltype = Ball(); δ0)
Returns a function such that T011(u0; δ0)(x) computes the integral
T_{0,1,1} from the paper.

The strategy for evaluation is to compute an expansion for the
integrand which is then integrated termwise on the interval. The first
two terms inside the absolute value are analytic and their Taylor
expansions of degree `N - 1` are computed at `t = 0` and the error
term is enclosed. The last term inside the absolute value is not
analytic and instead we use the expansion for Clausians from the paper
and also here enclose the error term.

The value inside the absolute value has constant sign so we can remove
it. Switching integration and summation gives us terms of the form
`∫_0^δ0 t^s*t^p dt` where `s` depends on the term and `p = u0.p`. This
is easily calculated to be `δ0^(s + p + 1*/(s + p + 1)`. The errors
are handled as constant values which are just multiplied by the length
of the interval.
"""
function T011(
    u0::FractionalKdVAnsatz{arb},
    ::Ball = Ball();
    δ0::arb = parent(u0.α)(1e-2),
    N::Integer = 3,
)
    Γ = Nemo.gamma
    α = u0.α

    PP = ArbPolyRing(parent(α), :x)

    M = div(N, 2) + 1

    return x -> begin
        # Analytic terms
        (P, P_E) = taylor_with_error(zero(α), setunion(zero(α), δ0), N) do t
            Ci(x * (1 - t), -α) + Ci(x * (1 + t), -α)
        end
        P_restterm = ball(zero(α), P_E * δ0^N)

        # Singular term
        (C, e, P2, P2_E) = Ci_expansion(x * δ0, -α, M)
        C *= x^e
        for m = 1:M-1
            P2[2m] *= x^(2m)
        end
        P2_restterm = P2_E * (x * δ0)^(2M)

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        res -= 2C * δ0^(e + u0.p + 1) / (e + u0.p + 1)

        # Integrate the analytic terms
        full_series = P - 2P2
        for i = 0:N-1
            res += full_series[i] * δ0^(i + u0.p + 1) / (i + u0.p + 1)
        end

        # Add the error term
        res += δ0 * (P_restterm - P2_restterm)

        # Prove: that the expression inside the absolute value of the
        # integrand is negative
        return -res * x / (parent(α)(π) * u0(x))
    end
end

"""
    T012(u0::FractionalKdVAnsatz{arb}, evaltype = Ball(); δ0, δ1)
Returns a function such that T012(u0; δ0, δ1)(x) computes the integral
T_{0,1,2} from the paper.
"""
function T012(
    u0::FractionalKdVAnsatz{arb},
    ::Ball = Ball();
    δ0::arb = parent(u0.α)(1e-2),
    δ1::arb = parent(u0.α)(1e-2),
    rtol = -1.0,
    atol = -1.0,
    show_trace = false,
)
    α = u0.α
    CC = ComplexField(precision(parent(α)))
    mα = CC(-α)
    a = CC(δ0)
    b = CC(1 - δ1)

    return x -> begin
        # PROVE: That there are no branch cuts that interact with the
        # integral
        F(t, analytic = false) = begin
            if isreal(t)
                res = CC(
                    Ci(x * (1 - real(t)), -α) + Ci(x * (1 + real(t)), -α) -
                    2Ci(x * real(t), -α),
                )
            else
                res = Ci(x * (1 - t), mα) + Ci(x * (1 + t), mα) - 2Ci(x * t, mα)
            end

            return ArbTools.real_abs(res, analytic = analytic) * t^u0.p
        end
        res = real(
            ArbTools.integrate(
                CC,
                F,
                a,
                b,
                rel_tol = rtol,
                abs_tol = atol,
                eval_limit = 3000,
                verbose = Int(show_trace),
            ),
        )

        return res * x / (parent(u0.α)(π) * u0(x))
    end
end

function T012(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball = Ball();
    δ0::Arf = Arf(1e-2),
    δ1::Arf = Arf(1e-2),
    rtol = 0,
    atol = 0,
    show_trace = false,
)
    α = u0.α
    mα = Acb(-α)
    a = Acb(δ0)
    b = 1 - Acb(δ1)

    return x -> begin
        # PROVE: That there are no branch cuts that interact with the
        # integral
        f(t; analytic = false) = begin
            t = Acb(t)
            if isreal(t)
                res = Acb(
                    Ci(x * (1 - real(t)), -α) + Ci(x * (1 + real(t)), -α) -
                    2Ci(x * real(t), -α),
                )
            else
                res = Ci(x * (1 - t), mα) + Ci(x * (1 + t), mα) - 2Ci(x * t, mα)
            end

            return Arblib.real_abs!(res, res, analytic) * t^u0.p
        end

        # TODO: Increase maximum number of evaluations
        res = real(Arblib.integrate(f, a, b, check_analytic = true; rtol, atol))

        return res * x / (Arb(π) * u0(x))
    end
end

"""
    T013(u0::FractionalKdVAnstaz{arb}, evaltype = Ball(); δ1)
Returns a function such that T013(u0; δ1)(x) computes the integral
T_{0,1,3} from the paper.

The strategy for evaluation is the same as for [`T011`](@ref) except
that the first term is singular and the last two are analytic and
their Taylor expansion is computed at `t = 1`.

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

TODO: We could precompute some of the values, in particular the
beta_inc functions can be precomputed.
"""
function T013(
    u0::FractionalKdVAnsatz{arb},
    ::Ball = Ball();
    δ1::arb = parent(u0.α)(1e-2),
    ϵ::arb = parent(u0.α)(1e-2),
    N::Integer = 3,
)
    Γ = Nemo.gamma
    α = u0.α
    π = parent(α)(pi)

    PP = ArbPolyRing(parent(α), :x)

    M = div(N, 2) + 1

    return x -> begin
        # Determine if the asymptotic expansion or the Taylor
        # expansion should be used for the second term
        use_asymptotic = π - x < ϵ

        # Analytic terms
        (P, E) = taylor_with_error(one(α), setunion(1 - δ1, one(α)), N) do t
            if !use_asymptotic
                return Ci(x * (1 + t), -α) - 2Ci(x * t, -α)
            else
                return -2Ci(x * t, -α)
            end
        end
        P_restterm = ball(zero(α), E * δ1^N)

        # Singular term
        (C, e, P2, P2_E) = Ci_expansion(x * δ1, -α, M)
        C *= x^e
        for m = 1:M-1
            P2[2m] *= x^(2m)
        end
        P2_restterm = P2_E * (x * δ1)^(2M)

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        # Using ∫_(1 - δ1)^1 |t - 1|^s*t^(p) dt = ∫_(1 - δ1)^1 (1 - t)^s*t^(p) dt
        res +=
            C * (
                Γ(1 + e) * Γ(1 + u0.p) / Γ(2 + e + u0.p) -
                beta_inc(1 + u0.p, 1 + e, 1 - δ1)
            )

        # Integrate the analytic part
        # Using ∫_(1-δ1)^1 (t-1)^i*t^(p) dt = (-1)^i ∫_(1-δ1)^1 (t-1)^i*t^(p) dt
        full_series = P + P2
        for i = 0:N-1
            res +=
                full_series[i] *
                (-1)^i *
                (
                    Γ(parent(α)(1 + i)) * Γ(1 + u0.p) / Γ(2 + i + u0.p) -
                    beta_inc(1 + u0.p, parent(α)(1 + i), 1 - δ1)
                )
        end

        # Add the error term
        res += δ1 * (P_restterm + P2_restterm)

        if use_asymptotic
            # Handle asymptotic expansion of Ci(x*(t + 1), -α)
            (C, e, P3, P3_E) = Ci_expansion(2π - x * (2 - δ1), -α, M)
            P3_restterm = P2_E * (2π - x * (2 - δ1))^(2M)

            # Add the singular part to the integral
            res +=
                C *
                (2π - x)^(1 + u0.p + e) *
                x^(-1 - u0.p) *
                (
                    beta_inc_zeroone(1 + u0.p, 1 + e, x / (2π - x)) -
                    beta_inc_zeroone(1 + u0.p, 1 + e, (x - δ1 * x) / (2π - x))
                )

            for i = 0:2:N-1
                # Only even terms
                res +=
                    P3[i] *
                    (2π - x)^(1 + u0.p + i) *
                    x^(-1 - u0.p) *
                    (
                        beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), x / (2π - x)) -
                        beta_inc_zeroone(
                            1 + u0.p,
                            parent(α)(1 + i),
                            (x - δ1 * x) / (2π - x),
                        )
                    )
            end

            # Add error term
            res += δ1 * P3_restterm
        end

        # Prove: that the expression inside the absolute value of the
        # integrand is positive
        return res * x / (parent(u0.α)(π) * u0(x))
    end
end
