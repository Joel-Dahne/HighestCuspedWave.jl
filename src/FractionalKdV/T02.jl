"""
    T02(u0::FractionalKdVAnsatz; δ2)
Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral T_{0,2} from the paper.

If `u0.p == 1` it uses a closed form expression for the integral.
Otherwise it computes the integral directly for most of the interval
and uses and asymptotic expansion for `y < x + δ2`.

If `x` is close to π (`π - x < ϵ`) then use only the asymptotic
expansion for the full integral.

The choice of both `δ2` and `ϵ` can be tuned a lot. For `δ2` it
depends on both `x` and `α`, whereas for `ϵ` it only depends on `α`.
For `δ2` it should be larger when both `x` and `α` are larger. In
theory we could compute with several values and take the best result,
but that would likely be to costly. For `ϵ` it's mainly a question of
cost, we don't want to compute the expensive integral when we believe
the asymptotic expansion will be the best anyway, larger values of `α`
should give a lower value of `ϵ`. The choice of `ϵ` is partially based
on /figures/optimal-epsilon-choice.png, but can likely be tuned
further.

TODO: Look closer at computing with the asymptotic expansion and using
the best result. Consider rewriting `T021` and `T022` to be more like
the other methods here.
"""
function T02(
    u0::FractionalKdVAnsatz{arb},
    ::Ball;
    δ2::arb = parent(u0.α)(1e-2),
    ϵ::arb = 1 + u0.α,
    rtol = -1.0,
    atol = -1.0,
    show_trace = false,
)
    π = parent(u0.α)(pi)

    if u0.p == 1
        # Use the closed form expression
        α = u0.α
        p = u0.p
        return x -> begin
            # TODO: Handle the case when x contains π. Then Si(2x, 1 -
            # α) evaluates to NaN. Si has issues as soon as the
            # argument is a ball containing a multiple of 2π.
            res =
                Ci(x + π, 2 - α) - Ci(π, 2 - α) + Ci(x, 2 - α) -
                (Ci(2x, 2 - α) + zeta(2 - α)) / 2 + x * Si(x, 1 - α)
            if π - x < 1e-4
                # When 2x is close to 2π direct evaluation
                # of Si fails. Use that Si(2x, 1 - α) = Si(2x - 2π, 1
                # - α) and the asymptotic expansion.
                y = 2x - 2π
                M = 3
                C, e, P, E = Si_expansion(y, 1 - α, M)
                res -=
                    x / 2 * (-C * abspow(y, e) + evaluate(P.poly, y) + E * abs(y)^(2M + 1))
            else
                res -= x / 2 * Si(2x, 1 - α)
            end
            return 2 / (π * u0.w(x) * u0(x)) * res
        end
    else
        return x -> begin
            a = ArbTools.ubound(x + δ2)

            # Compute with the asymptotic expansion on the whole interval
            res_asymptotic = T021(u0, Ball(), π, x, ϵ = π)

            if π < a || π - x < ϵ
                return res_asymptotic
            end

            part1 = T021(u0, Ball(), a, x; ϵ)

            part2 = T022(u0, Ball(), a, x; rtol, atol, show_trace)

            res = part1 + part2
            if radius(res) < radius(res_asymptotic)
                return res
            else
                return res_asymptotic
            end
        end
    end
end

function T02(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball;
    δ2::Arf = Arf(1e-2),
    ϵ::Arb = 1 + u0.α,
    skip_div_u0 = false,
)
    if u0.p == 1
        # Use the closed form expression
        α = u0.α
        p = u0.p
        return x -> begin
            # TODO: Handle the case when x contains π. Then Si(2x, 1 -
            # α) evaluates to NaN. Si has issues as soon as the
            # argument is a ball containing a multiple of 2π.
            res =
                Ci(x + π, 2 - α) - Ci(Arb(π), 2 - α) + Ci(x, 2 - α) -
                (Ci(2x, 2 - α) + zeta(2 - α)) / 2 + x * Si(x, 1 - α)
            if π - x < 1e-4
                # When 2x is close to 2π direct evaluation
                # of Si fails. Use that Si(2x, 1 - α) = Si(2x - 2π, 1
                # - α) and the asymptotic expansion.
                y = 2x - 2Arb(π)
                M = 3
                C, e, P, E = Si_expansion(y, 1 - α, M)
                res -= x / 2 * (-C * abspow(y, e) + P(y) + E * abs(y)^(2M + 1))
            else
                res -= x / 2 * Si(2x, 1 - α)
            end
            if skip_div_u0
                return 2 / (π * u0.w(x)) * res
            else
                return 2 / (π * u0.w(x) * u0(x)) * res
            end
        end
    else
        return x -> begin
            a = Arblib.ubound(Arb, x + δ2)

            # Compute with the asymptotic expansion on the whole interval
            res_asymptotic = T021(u0, Ball(), Arb(π), x, ϵ = Arb(π), skip_div_u0 = true)

            if π < a || π - x < ϵ
                return res_asymptotic
            end

            part1 = T021(u0, Ball(), a, x, skip_div_u0 = true; ϵ)

            part2 = T022(u0, Ball(), a, x, skip_div_u0 = true)

            res = intersect(part1 + part2, res_asymptotic)
            if skip_div_u0
                return res
            else
                return res / u0(x)
            end
        end
    end
end

"""
    T02(u0::FractionalKdVAnsatz, ::Asymptotic)
Returns a function such that `T02(u0, Asymptotic())(x)` computes an
**upper bound** of the integral T_{0,2} from the paper using an
evaluation strategy that works asymptotically as `x` goes to 0.
"""
function T02(
    u0::FractionalKdVAnsatz{arb},
    ::Asymptotic;
    N::Integer = 100,
    nonasymptotic_u0 = false, # Mainly for testing
)
    α = u0.α
    p = u0.p
    RR = parent(α)
    π = RR(pi)

    if p == 1
        return x -> begin
            (A, _, P1, E1) = Ci_expansion(x, 2 - α, 2)
            (_, _, _, E1p) = Ci_expansion(2x, 2 - α, 2)
            (B, _, P2, E2) = Si_expansion(x, 1 - α, 1)
            (_, _, _, E2p) = Si_expansion(2x, 1 - α, 1)
            # PROVE: That P3[1] == P3[3] == 0
            (P3, E3) = taylor_with_error(π, setunion(π, π + x), 4) do y
                return Ci(y, 2 - α)
            end

            c_α = A * (1 - parent(α)(2)^(-α)) + B * (1 - parent(α)(2)^(-α - 1))
            K = 2 / (π * a0(u0, 0))
            # Ball containing 1 + hat(u0)(x)
            L = ball(parent(α)(1), c(u0, ArbTools.abs_ubound(x)) * abspow(x, u0.p0))

            res = (
                K * c_α * L +
                K / 2 * L * (zeta(-α) - Ci(π, -α)) * abspow(x, 1 + α) +
                K * L * abspow(x, α + 3) * (E3 + E1 - 8 * E1p + E2 - 8 * E2p)
            )

            return res
        end
    end

    @warn "T02(u0, Asymptotic()) is not yet rigorous when u0.p != 1 - the tail of the sum is not bounded"
    Γ = Nemo.gamma
    CC = ComplexField(precision(RR))

    # Compute
    # c_α = |Γ(1 + α)*sin(πα/2)|∫₁^∞ |(t - 1)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)|tᵖ dt
    # PROVE: That this is real
    c_α =
        abs(Γ(1 + α) * sinpi(α / 2)) * (
            real(
                CC(-1)^CC(-p) * (CC(-1)^CC(α) + CC(-1)^CC(p)) * Γ(α - p) * Γ(1 + p) /
                Γ(1 + α) + CC(-1)^CC(α) * Γ(-α) * Γ(1 + p) / (Γ(1 - α + p)),
            ) - 2 / (α - p) - hypgeom_2f1(1 + α, 1 + p, 2 + p, RR(-1)) / (1 + p)
        )

    # We evaluate the sum in the paper at x = 0, the expression
    # coming from the interior integral then simplifies a lot.
    # TODO: Bound the tail - the terms go to zero extremely
    # fast so it should be negligible
    # PROVE: That the value at x = 0 is indeed a bound of the
    # sum. Could it be that this is not the case?
    c_ϵ = zero(α)
    for m = 1:30
        integral = 2 * π^(2m - 1 + p) * binom(RR(2m), unsigned(2m - 2)) / (2m - 1 + p)
        c_ϵ += (-one(α))^m / factorial(2fmpz(m)) * zeta(-α - 2m) * integral
    end

    return x -> begin
        if nonasymptotic_u0
            # Version without asymptotically expanding u0(x)
            res = c_α * abspow(x, -α) + c_ϵ * abspow(x, 2 - p)
            return res / (π * u0(x))
        end

        res = c_α + c_ϵ * abspow(x, 2 - p + α)
        # Ball containing 1 + hat(u0)(x)
        L = ball(parent(α)(1), c(u0, ArbTools.abs_ubound(x)) * abspow(x, u0.p0))
        return L / (π * a0(u0, 0)) * res
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
T021(u0::FractionalKdVAnsatz, a, x; kwargs...) = T021(u0, Ball(), a, x; kwargs...)

function T021(
    u0::FractionalKdVAnsatz{arb},
    ::Ball;
    δ2::arb = parent(u0.α)(1e-4),
    ϵ::arb = parent(u0.α)(1e-1),
    N::Integer = 3,
)
    return x -> begin
        T021(u0, Ball(), ArbTools.ubound(x + δ2), x, ϵ = ϵ, N = N)
    end
end

function T021(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball;
    δ2::Arf = Arf(1e-4),
    ϵ::Arb = Arb(1e-1),
    N::Integer = 3,
)
    return x -> begin
        T021(u0, Ball(), Arblib.ubound(Arb, x + δ2), x; ϵ, N)
    end
end

function T021(
    u0::FractionalKdVAnsatz{arb},
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
    P_restterm = ball(zero(α), E * δ2^N)

    # Singular term
    M = div(N, 2) + 1
    (C, e, P2, P2_E) = Ci_expansion(δ2, -α, M)
    P2_restterm = P2_E * (δ2)^(2M)

    # Compute the integral
    res = zero(α)
    # Integrate the singular term
    # Using ∫_x^a |x - y|^s*y^p dy = ∫_x^a (y - x)^s*y^p dy
    if u0.p == 1
        res += C * abspow(δ2, 1 + e) * (δ2 / (2 + e) + x / (1 + e))
    else
        res +=
            C *
            x^(1 + e + u0.p) *
            (
                Γ(1 + e) * Γ(-1 - e - u0.p) / Γ(-u0.p) -
                beta_inc_zeroone(-1 - e - u0.p, 1 + e, x / a)
            )
    end

    # Integrate the analytic terms
    full_series = P + P2
    for i = 0:N-1
        if u0.p == 1
            res += full_series[i] * δ2^(1 + i) * (δ2 / (2 + i) + x / (1 + i))
        else
            res +=
                full_series[i] *
                x^(1 + i + u0.p) *
                (
                    Γ(parent(α)(1 + i)) * Γ(-1 - i - u0.p) / Γ(-u0.p) -
                    beta_inc(-1 - i - u0.p, parent(α)(1 + i), x / a)
                )
        end
    end

    # Add the error term
    res += δ2 * (P_restterm + P2_restterm)

    if use_asymptotic
        # Handle asymptotic expansion of Ci(x + y, -α)
        # The furthest away from 2π we are is at y = x
        (C, e, P3, P3_E) = Ci_expansion(2π - 2x, -α, M)
        P3_restterm = P2_E * (2π - 2x)^(2M)

        # Add the singular part to the integral
        res +=
            C *
            (2π - x)^(1 + e + u0.p) *
            (
                beta_inc_zeroone(1 + u0.p, 1 + e, a / (2π - x)) -
                beta_inc_zeroone(1 + u0.p, 1 + e, x / (2π - x))
            )

        for i = 0:2:N-1
            # Only even terms
            res +=
                P3[i] *
                (2π - x)^(1 + i + u0.p) *
                (
                    beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), a / (2π - x)) -
                    beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), x / (2π - x))
                )
        end

        # Add error term
        res += δ2 * P3_restterm
    end

    # Prove: that the expression inside the absolute value of the
    # integrand is positive
    return res / (parent(u0.α)(π) * u0.w(x) * u0(x))
end

function T021(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball,
    a::Arb,
    x::Arb;
    ϵ::Arb = Arb(1e-1),
    N::Integer = 3,
    skip_div_u0 = false,
)
    Γ = SpecialFunctions.gamma
    α = u0.α
    δ2 = a - x

    # Determine if the asymptotic expansion or the Taylor
    # expansion should be used for the second term
    use_asymptotic = π - x < ϵ

    # Analytic terms
    (P, E) = taylor_with_error(x, union(x, a), N) do y
        if !use_asymptotic
            return Ci(x + y, -α) - 2Ci(y, -α)
        else
            return -2Ci(y, -α)
        end
    end
    P_restterm = Arblib.add_error!(zero(α), E * δ2^N)

    # Singular term
    M = N ÷ 2 + 1
    (C, e, P2, P2_E) = Ci_expansion(δ2, -α, M)
    P2_restterm = P2_E * δ2^(2M)

    # Compute the integral
    res = zero(α)
    # Integrate the singular term
    # Using ∫_x^a |x - y|^s*y^p dy = ∫_x^a (y - x)^s*y^p dy
    if u0.p == 1
        res += C * abspow(δ2, 1 + e) * (δ2 / (2 + e) + x / (1 + e))
    else
        res +=
            C *
            x^(1 + e + u0.p) *
            (
                Γ(1 + e) * Γ(-1 - e - u0.p) / Γ(-u0.p) -
                beta_inc_zeroone(-1 - e - u0.p, 1 + e, x / a)
            )
    end

    # Integrate the analytic terms
    full_series = P + P2
    for i = 0:N-1
        if u0.p == 1
            res += full_series[i] * δ2^(1 + i) * (δ2 / (2 + i) + x / (1 + i))
        else
            res +=
                full_series[i] *
                x^(1 + i + u0.p) *
                (
                    Γ(Arb(1 + i)) * Γ(-1 - i - u0.p) / Γ(-u0.p) -
                    beta_inc(-1 - i - u0.p, Arb(1 + i), x / a)[1]
                )
        end
    end

    # Add the error term
    res += δ2 * (P_restterm + P2_restterm)

    if use_asymptotic
        # Handle asymptotic expansion of Ci(x + y, -α)
        # The furthest away from 2π we are is at y = x
        (C, e, P3, P3_E) = Ci_expansion(2Arb(π) - 2x, -α, M)
        P3_restterm = P2_E * (2Arb(π) - 2x)^(2M)

        # Add the singular part to the integral
        res +=
            C *
            (2Arb(π) - x)^(1 + e + u0.p) *
            (
                beta_inc_zeroone(1 + u0.p, 1 + e, a / (2Arb(π) - x)) -
                beta_inc_zeroone(1 + u0.p, 1 + e, x / (2Arb(π) - x))
            )

        for i = 0:2:N-1
            # Only even terms
            res +=
                P3[i] *
                (2Arb(π) - x)^(1 + i + u0.p) *
                (
                    beta_inc_zeroone(1 + u0.p, Arb(1 + i), a / (2Arb(π) - x)) -
                    beta_inc_zeroone(1 + u0.p, Arb(1 + i), x / (2Arb(π) - x))
                )
        end

        # Add error term
        res += δ2 * P3_restterm
    end

    # Prove: that the expression inside the absolute value of the
    # integrand is positive
    if skip_div_u0
        return res / (Arb(π) * u0.w(x))
    else
        return res / (Arb(π) * u0.w(x) * u0(x))
    end
end

"""
    T022(u0::FractionalKdVAnsatz{arb}, a::arb, x::arb)
Computes the (not yet existing) integral T_{0,2,2} from the paper.

TODO: Write about how this is done.
"""
T022(u0::FractionalKdVAnsatz, a, x; kwargs...) = T022(u0, Ball(), a, x; kwargs...)

function T022(
    u0::FractionalKdVAnsatz{arb},
    ::Ball;
    δ2::arb = parent(u0.α)(1e-4),
    rtol = -1.0,
    atol = -1.0,
    show_trace = false,
)
    return x -> begin
        T022(
            u0,
            Ball(),
            ArbTools.ubound(x + δ2),
            x,
            rtol = rtol,
            atol = atol,
            show_trace = show_trace,
        )
    end
end

function T022(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball;
    δ2::Arf = Arf(1e-4),
    skip_div_u0 = false,
)
    return x -> begin
        T022(u0, Ball(), Arblib.ubound(Arb, x + δ2), x; skip_div_u0)
    end
end

function T022(
    u0::FractionalKdVAnsatz{arb},
    ::Ball,
    a::arb,
    x::arb;
    rtol = -1.0,
    atol = -1.0,
    show_trace = false,
)
    CC = ComplexField(precision(parent(u0.α)))
    mα = CC(-u0.α)
    a = CC(a)
    b = CC(parent(u0.α)(π))

    # PROVE: That there are no branch cuts that interact with the
    # integral
    F(y; analytic = false) = begin
        if isreal(y)
            res = CC(
                Ci(x - real(y), real(mα)) + Ci(x + real(y), real(mα)) -
                2Ci(real(y), real(mα)),
            )
        else
            res = Ci(x - y, mα) + Ci(x + y, mα) - 2Ci(y, mα)
        end

        return ArbTools.real_abs(res, analytic = analytic) * y^u0.p
    end
    res = real(
        ArbTools.integrate(
            CC,
            F,
            a,
            b,
            rel_tol = rtol,
            abs_tol = atol,
            checkanalytic = true,
            verbose = Int(show_trace),
        ),
    )

    return res / (parent(u0.α)(π) * u0.w(x) * u0(x))

end

function T022(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball,
    a::Arb,
    x::Arb;
    skip_div_u0 = false,
)
    mα = Acb(-u0.α)
    a = Acb(a)
    b = Acb(π)

    # PROVE: That there are no branch cuts that interact with the
    # integral
    f(y; analytic = false) = begin
        y = Acb(y)
        if isreal(y)
            res = Acb(Ci(x - real(y), -u0.α) + Ci(x + real(y), -u0.α) - 2Ci(real(y), -u0.α))
        else
            res = Ci(x - y, mα) + Ci(x + y, mα) - 2Ci(y, mα)
        end

        return Arblib.real_abs!(res, res, analytic) * y^u0.p
    end

    # TODO: Increase maximum number of evaluations
    res = real(
        Arblib.integrate(
            f,
            a,
            b,
            check_analytic = true,
            rtol = 1e-5,
            atol = 1e-5,
            warn_on_no_convergence = false,
        ),
    )

    if skip_div_u0
        return res / (Arb(π) * u0.w(x))
    else
        return res / (Arb(π) * u0.w(x) * u0(x))
    end
end
