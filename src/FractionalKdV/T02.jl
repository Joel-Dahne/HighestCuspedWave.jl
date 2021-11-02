"""
    T02(u0::FractionalKdVAnsatz{Arb}; δ2)

Returns a function such that `T02(u0; δ2, ϵ)(x)` computes the integral
\$T_{0,2}\$ from the paper.

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
on [/figures/optimal-epsilon-choice.png], but can likely be tuned
further.

TODO: Look closer at computing with the asymptotic expansion and using
the best result. Consider rewriting `T021` and `T022` to be more like
the other methods here.
"""
function T02(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball;
    δ2::Arf = Arf(1e-2),
    ϵ::Arb = 1 + u0.α,
    skip_div_u0 = false,
)
    if isone(u0.p)
        # Use the closed form expression
        α = u0.α
        p = u0.p
        return x -> begin
            # TODO: Handle the case when x contains π. Then
            # clausens(2x, 1 - α) evaluates to NaN. clausens has
            # issues as soon as the argument is a ball containing a
            # multiple of 2π.
            res =
                clausenc(x + π, 2 - α) - clausenc(Arb(π), 2 - α) + clausenc(x, 2 - α) -
                (clausenc(2x, 2 - α) + zeta(2 - α)) / 2 + x * clausens(x, 1 - α)
            if π - x < 1e-4
                # When 2x is close to 2π direct evaluation of clausens
                # fails. Use that clausens(2x, 1 - α) = clausens(2x -
                # 2π, 1 - α) and the asymptotic expansion.
                y = 2x - 2Arb(π)
                M = 3
                C, e, P, E = clausens_expansion(y, 1 - α, M)
                res -= x / 2 * (-C * abspow(y, e) + P(y) + E * abs(y)^(2M + 1))
            else
                res -= x / 2 * clausens(2x, 1 - α)
            end

            if skip_div_u0
                return 2res / (π * u0.w(x))
            else
                return 2res / (π * u0.w(x) * u0(x))
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
    T02(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic)

Returns a function such that `T02(u0, Asymptotic())(x)` computes an
**upper bound** of the integral \$T_{0,2}\$ from the paper using an
evaluation strategy that works asymptotically as `x` goes to 0.

It splits the Clausen functions into the main singular part and the
analytic expansion. The integral of the singular part is computed
directly. The expansion is integrated term wise and the resulting sum
is bounded.

For the singular part we get the integral
```
c_α = abs(gamma(1 + α) * sinpi(α / 2)) * ∫ abs((t - 1)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)) * t^p dt
```
from `1` to `∞`. The expression inside the absolute value is positive
and after removing it we can integrate explicitly. Skipping the factor
in front we have
```
∫ ((t - 1)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)) * t^p dt =
    (hypgeom_2f1(1 + α, α - p, 1 + α - p, Arb(-1)) - 2) / (α - p) +
    gamma(-α) * gamma(α - p) / gamma(-p)
```
This was computed using Mathematica with the following code
```
Integrate[((t - 1)^(-1 - a) + (t + 1)^(-1 - a) - 2 t^(-1 - a))*t^p, {t, 1, Infinity}, Assumptions -> -1 < a < 0]
```
- **FIXME:** Mathematica gives the condition `-1 < Re[p] < 0 && a >
  Re[p]` which is not satisfied in our case. However the expression
  seems to give the correct values so we use it for now.
- **PROVE:** That the expression is positive.

- **TODO:** Write down the expression for the tail and the formulas used
  for the integration. They are in the paper.
- **TODO:** Bound the error terms for the tails.
- **TODO:** Write down the last parts of the computations.
- **TODO:** Explain evaluation strategy for the case `p == 1`.
"""
function T02(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    nonasymptotic_u0 = false, # Mainly for testing
)
    α = u0.α
    p = u0.p
    π = Arb(Irrational{:π}())

    if isone(p)
        return x -> begin
            (A, _, P1, E1) = clausenc_expansion(x, 2 - α, 2)
            (_, _, _, E1p) = clausenc_expansion(2x, 2 - α, 2)
            (B, _, P2, E2) = clausens_expansion(x, 1 - α, 1)
            (_, _, _, E2p) = clausens_expansion(2x, 1 - α, 1)

            # PROVE: That P3[1] == P3[3] == 0
            (P3, E3) = taylor_with_error(π, setunion(π, π + x), 4) do y
                return clausenc(y, 2 - α)
            end

            c_α = A * (1 - 2^(-α)) + B * (1 - 2^(-α - 1))

            K = 2 / (π * a0(u0, 0))

            # Ball containing 1 + hat(u0)(x)
            L = Arblib.add_error!(
                one(α),
                c(u0, Arblib.abs_ubound(Arb, x)) * abspow(x, u0.p0),
            )

            res = (
                K * c_α * L +
                K / 2 * L * (zeta(-α) - clausenc(π, -α)) * abspow(x, 1 + α) +
                K * L * abspow(x, α + 3) * (E3 + E1 - 8 * E1p + E2 - 8 * E2p)
            )

            return res
        end
    end

    @warn "T02(u0, Asymptotic()) is not yet rigorous - tail of sum not bounded"

    # Compute c_α
    c_α = let
        integral =
            (hypgeom_2f1(1 + α, α - p, 1 + α - p, Arb(-1)) - 2) / (α - p) +
            gamma(-α) * gamma(α - p) / gamma(-p)

        abs(gamma(1 + α) * sinpi(α / 2)) * integral
    end

    # We evaluate the sum in the paper at x = 0, the expression
    # coming from the interior integral then simplifies a lot.
    # TODO: Bound the tail - the terms go to zero extremely
    # fast so it should be negligible
    # PROVE: That the value at x = 0 is indeed a bound of the
    # sum. Could it be that this is not the case?
    c_ϵ = zero(α)
    for m = 1:20
        integral = 2 * binomial(2m, 2m - 2) * π^(2m - 1 + p) / (2m - 1 + p)
        term = (-1)^m * zeta(-α - 2m) * integral / factorial(Arb(2m))

        c_ϵ += term
    end

    return x -> begin
        if nonasymptotic_u0
            # Version without asymptotically expanding u0(x)
            res = c_α * abspow(x, -α) + c_ϵ * abspow(x, 2 - p)
            return res / (π * u0(x))
        end

        res = c_α + c_ϵ * abspow(x, 2 - p + α)
        # Ball containing 1 + hat(u0)(x)

        L = Arblib.add_error!(one(α), c(u0, Arblib.abs_ubound(Arb, x)) * abspow(x, u0.p0))

        return L / (π * a0(u0, 0)) * res
    end
end

"""
    T021(u0::FractionalKdVAnstaz{Arb}, a::Arb, x::Arb)
Computes the (not yet existing) integral T_{0,2,1} from the paper.

That is the integral
```
∫abs(clausenc(y - x, -α) + clausenc(y + x, -α) - 2clausenc(y, -α)) * y^p dy
```
from `x` to `a`.

The strategy for evaluation is the same as for [`T011`](@ref) except
that the first term is singular and the last two are analytic and
their Taylor expansion is computed at `t = x`.

The integral that needs to be computed in this case is
```
∫_x^a (y - x)^s*y^p dy
```
which is given by
```
x^(s + p + 1) * (gamma(1 + s) * gamma(-1 - s - p) / gamma(-p) - B(-1 - s - p, 1 + s; x / a))
```
where B(a, b; z) is the incomplete Beta-function. In the end we are
dividing by `w(x) = x^p` and hence we can simplify it by cancelling
the `p` in the `x`-exponent directly.

For `p = 1` we instead use the expression
```
(a - x)^(1+s) * ((a - x) / (2 + s) + x / (1 + s))
```
and we divide by the weight directly.

If `x` is equal or very close to π (determined by `ϵ`) then the Taylor
expansion gives a very poor approximation for `clausenc(x + y, -α)`.
In this case we make use of the fact that it's 2π periodic and even,
so that `clausenc(x + y, -α) = clausenc(x + y - 2π, -α) = clausenc(2π
- (x + y), -α)`, to be able to use the asymptotic expansion instead.
That gives us the integral
```
∫_x^a (2π - (x + y))^s*t^p dy
```
which is given by
```
(2π - x)^(1 + p + s)*(B(1 + p, 1 + s, a/(2π - x)) - B(1 + p, 1 + s, x/(2π - x)))
```
The value of `x/(2π - x)` will always be less than or equal to 1 for
`x` less than or equal to π, however due to overestimation the
enclosing ball might contain values greater than one, we therefore
have to use `beta_inc_zeroone` to be able to get finite results in
that case.
"""
T021(u0::FractionalKdVAnsatz, a, x; kwargs...) = T021(u0, Ball(), a, x; kwargs...)

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
    u0::FractionalKdVAnsatz{Arb},
    ::Ball,
    a::Arb,
    x::Arb;
    ϵ::Arb = Arb(1e-1),
    N::Integer = 3,
    skip_div_u0 = false,
)
    Γ = gamma
    α = u0.α
    δ2 = a - x

    # Determine if the asymptotic expansion or the Taylor
    # expansion should be used for the second term
    use_asymptotic = π - x < ϵ

    # Compute expansion

    # Singular term
    M = N ÷ 2 + 1
    (C, e, P1, P1_E) = clausenc_expansion(δ2, -α, M)
    P1_restterm = P1_E * δ2^(2M)

    # Analytic terms
    (P2, P2_E) = taylor_with_error(x, union(x, a), N) do y
        if !use_asymptotic
            return clausenc(x + y, -α) - 2clausenc(y, -α)
        else
            return -2clausenc(y, -α)
        end
    end
    P2_restterm = Arblib.add_error!(zero(α), P2_E * δ2^N)

    # Compute the integrals

    # Integrate the singular term and divide by w(x) = x^p
    if isone(u0.p)
        singular_term = C * abspow(δ2, 1 + e) * (δ2 / (2 + e) + x / (1 + e)) / u0.w(x)
    else
        singular_term =
            C *
            x^(1 + e) *
            (
                Γ(1 + e) * Γ(-1 - e - u0.p) / Γ(-u0.p) -
                beta_inc_zeroone(-1 - e - u0.p, 1 + e, x / a)
            )
    end

    # Integrate the analytic terms and divide them by w(x) = x^p
    full_series = P1 + P2

    analytic_term = zero(α)
    for i = 0:N-1
        if isone(u0.p)
            analytic_term +=
                full_series[i] * δ2^(1 + i) * (δ2 / (2 + i) + x / (1 + i)) / u0.w(x)
        else
            analytic_term +=
                full_series[i] *
                x^(1 + i) *
                (
                    Γ(Arb(1 + i)) * Γ(-1 - i - u0.p) / Γ(-u0.p) -
                    beta_inc(-1 - i - u0.p, Arb(1 + i), x / a)[1]
                )
        end
    end

    res = singular_term + analytic_term

    # Add rest term
    res += δ2 * (P1_restterm + P2_restterm)

    if use_asymptotic
        # Handle asymptotic expansion of clausenc(x + y, -α)
        # The furthest away from 2π we are is at y = x
        (C, e, P3, P3_E) = clausenc_expansion(2Arb(π) - 2x, -α, M)
        P3_restterm = P3_E * (2Arb(π) - 2x)^(2M)

        # Add the singular part divided by u0.w(x)
        singular_term_2 =
            C *
            (2Arb(π) - x)^(1 + e + u0.p) *
            (
                beta_inc_zeroone(1 + u0.p, 1 + e, a / (2Arb(π) - x)) -
                beta_inc_zeroone(1 + u0.p, 1 + e, x / (2Arb(π) - x))
            ) / u0.w(x)

        for i = 0:2:N-1
            # Only even terms
            singular_term_2 +=
                P3[i] *
                (2Arb(π) - x)^(1 + i + u0.p) *
                (
                    beta_inc_zeroone(1 + u0.p, Arb(1 + i), a / (2Arb(π) - x)) -
                    beta_inc_zeroone(1 + u0.p, Arb(1 + i), x / (2Arb(π) - x))
                ) / u0.w(x)
        end

        # Add rest term
        singular_term_2 += δ2 * P3_restterm

        # Add to res
        res += singular_term_2
    end

    # Prove: that the expression inside the absolute value of the
    # integrand is positive

    if skip_div_u0
        return res / π
    else
        return res / (π * u0(x))
    end
end

"""
    T022(u0::FractionalKdVAnsatz{Arb}, a::Arb, x::Arb)

Compute the (not yet existing) integral \$T_{0,2,2}\$ from the paper.

This is the integral
```
∫abs(clausenc(y - x, -α) + clausenc(y + x, -α) - 2clausenc(y, -α)) * y^p dy
```
from `a` to `π`. We should have `x < a < π`.

On this interval
```
clausenc(y - x, -α) + clausenc(y + x, -α) - 2clausenc(y, -α)
```
is positive so we can skip the absolute value.
- **PROVE:** That `clausenc(y - x, -α) + clausenc(y + x, -α) -
  2clausenc(y, -α)` is positive on the interval `[x, π]`.
"""
T022(u0::FractionalKdVAnsatz, a, x; kwargs...) = T022(u0, Ball(), a, x; kwargs...)

function T022(u0::FractionalKdVAnsatz{Arb}, ::Ball, a::Arb, x::Arb; skip_div_u0 = false)
    mα = -u0.α
    cp = Acb(u0.p)

    integrand(y) = begin
        # The integrand is singular at y = x so check that the real
        # part of y is greater than x or return an indeterminate
        # result. This also ensures that y^u0.p is analytic on all of
        # y.
        x < Arblib.realref(y) || return Arblib.indeterminate!(zero(y))

        if isreal(y)
            y = real(y)
            return Acb(
                (clausenc(y - x, mα) + clausenc(x + y, mα) - 2clausenc(y, mα)) * y^u0.p,
            )
        else
            y = Acb(y)
            return (clausenc(x - y, mα) + clausenc(x + y, mα) - 2clausenc(y, mα)) * y^cp
        end
    end

    res = real(
        Arblib.integrate(
            integrand,
            a,
            π,
            check_analytic = false,
            rtol = 1e-5,
            atol = 1e-5,
            warn_on_no_convergence = false,
            opts = Arblib.calc_integrate_opt_struct(0, 5_000, 0, 0, 0),
        ),
    )

    if skip_div_u0
        return res / (π * u0.w(x))
    else
        return res / (π * u0.w(x) * u0(x))
    end
end
