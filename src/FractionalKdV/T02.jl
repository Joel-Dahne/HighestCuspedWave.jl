"""
    T02(u0::FractionalKdVAnsatz{Arb}; δ2)

Returns a function such that `T02(u0; δ2, ϵ)(x)` computes the integral
``T_{0,2}`` from the paper.

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

**IMPROVE:** Look closer at computing with the asymptotic expansion
and using the best result. Consider rewriting `T021` and `T022` to be
more like the other methods here.
"""
function T02(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball;
    δ2::Arf = Arf(1e-2),
    ϵ::Arb = 1 + u0.α,
    skip_div_u0 = false,
)
    return x -> begin
        a = ubound(Arb, x + δ2)

        # Compute with the asymptotic expansion on the whole interval
        res_asymptotic = T021(u0, Ball(), Arb(π), x, ϵ = Arb(π), skip_div_u0 = true)

        if !(π < a || π - x < ϵ)
            part1 = T021(u0, Ball(), a, x, skip_div_u0 = true; ϵ)

            part2 = T022(u0, Ball(), a, x, skip_div_u0 = true)

            res = intersect(part1 + part2, res_asymptotic)
        else
            res = res_asymptotic
        end

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T02(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M = 5, ϵ = one(Arb))

Returns a function for computing an **upper bound** of the integral of
`T02(u0)`, using an evaluation strategy that works asymptotically as
`x` goes to zero.

# Arguments
- `M::Integer` determines the number of terms in the expansions.
- `ϵ::Arb` determines the interval ``[0, ϵ]`` on which the expansion
  is valid.

# Implementation
It first splits the function as
```
T01(u0)(x) = inv(π) * inv(u0(x) / x^-α) * (U02(x) / x^(-α + p))
```
where `α = u0.α`, `p = u0.p` and
```
U02(x) = x^(1 + p) * ∫ (clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t^u0.p
```
integrated from `1` to `π / x`. The factor `inv(u0(x) / x^-α)` is
computed using [`inv_u0_normalised`](@ref) so the remaining work is in
bounding the `U02 / x^(-α + p)` factor.

From the lemma in the paper we have
```
U02 / x^(-α + p) <= c + d * x^(2 + α - p)
```
with
```
c = gamma(1 + α) * sinpi(-α / 2) * (
    gamma(-α) * gamma(α - p) / gamma(-p) +
    (hypgeom_2f1(1 + α, α - p, 1 + α - p, -1) - 2) / (α - p)
)
```
and
```
d = 2π^(p - 1) * do m
        (-1)^m *
        zeta(-α - 2m) *
        Arb(π)^2m /
        factorial(2m) *
        sum(binomial(2m, 2k) * (ϵ / π)^(2(m - 1 - k)) / (2k + 1) for k = 0:m-1)
    end +
    6π^(p - 1) * sum((-1)^m * zeta(-α - 2m) * (3π / 2)^2m / factorial(2m) for m = N:Inf)
```
for any `N >= 1`.

To compute an enclosure of the tail of `d` we note that it is the same
as the sum in [`clausenc_expansion_remainder`](@ref) with `x = 3π / 2`.

# Notes
The expression for `c` was was computed using Mathematica with the
following code
```
Integrate[((t - 1)^(-1 - a) + (t + 1)^(-1 - a) - 2 t^(-1 - a))*t^p, {t, 1, Infinity}, Assumptions -> -1 < a < 0]
```

**TODO:** We could get a much better enclosure by using that we
only integrate from `1` to `π / x` and not to infinity. In particular
for `α` close to `-1` this makes a big difference. In that case `c`
would depend on `x` and be given by
```
gamma(1 + α) * sinpi(-α / 2) *
(
    gamma(-α) * gamma(α - p) / gamma(-p) +
    2(1 - (π / x)^(-α + p)) / (-α + p) +
    (-1)^(-α + p) * (beta_inc(α - p, -α, -1) - beta_inc(α - p, -α, -x / π))
    - beta_inc(α - p, -α, x / π)
)
```
computed using Mathematica with the following code
```
Integrate[((t - 1)^(-a - 1) + (1 + t)^(-a - 1) - 2 t^(-a - 1))*t^p, {t, 1, Pi/x}, Assumptions -> a < 0 && 0 < x < Pi]
```
There are some problems with this tough
- Evaluating the term with `(-1)^(-α + p)` requires complex
  arithmetic.
- We can't compute it when `x` is very small.
- For `α` close to `-1` it only gives reasonable enclosures for very
  thin balls for `α`. We might have to improve the computed
  enclosures.
"""
function T02(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M::Integer = 5, ϵ::Arb = Arb(1))
    @assert ϵ <= Arb(π) / 2 # This is required for the bound of the tail of d

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    α = u0.α
    p = u0.p

    c =
        gamma(1 + α) *
        sinpi(-α / 2) *
        (
            gamma(-α) * gamma(α - p) / gamma(-p) +
            (hypgeom_2f1(1 + α, α - p, 1 + α - p, -one(α)) - 2) / (α - p)
        )

    # c has a removable singularity for p = 1 and for some values of α
    # for p < 1. In practice we don't encounter these values because
    # those are covered by T0_p_one.
    isfinite(c) || @error "non-finite enclosure for c in T02" α p

    # Take a lot of terms since the remainder is computed at a large
    # value
    N = 30
    # Sum first N - 1 terms
    d =
        2Arb(π)^(p - 1) * sum(1:N-1) do m
            (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) *
            sum(binomial(2m, 2k) * (ϵ / π)^(2(m - 1 - k)) / (2k + 1) for k = 0:m-1)
        end
    # Enclose remainder
    d +=
        6Arb(π)^(p - 1) *
        (3Arb(π) / 2)^(2N) *
        clausenc_expansion_remainder(3Arb(π) / 2, -α, N)

    return x::Arb -> begin
        @assert x <= ϵ

        # TODO: This is not covered in the paper yet. It doesn't work
        # very well for wide values of α close to -1, the enclosures
        # are bad. This will have to be improved if we don't come up
        # with a different solution.

        # Attempt to compute a better bound for c by only integrating
        # up to π / x. Since the integral is increasing in x it is
        # enough to evaluate at a lower bound of x to get an upper
        # bound.
        cx = let x = abs_lbound(Arb, x)
            cx_complex_part =
                Acb(-1)^(-α + p) * (
                    beta_inc(Acb(α - p), -Acb(α), Acb(-1)) -
                    beta_inc(Acb(α - p), -Acb(α), -Acb(x / π))
                )
            @assert Arblib.contains_zero(imag(cx_complex_part))

            cx =
                gamma(1 + α) *
                sinpi(-α / 2) *
                (
                    gamma(-α) * gamma(α - p) / gamma(-p) +
                    2(1 - (π / x)^(-α + p)) / (-α + p) +
                    real(cx_complex_part) - beta_inc(α - p, -α, x / π)
                )

            cx
        end

        if isfinite(cx) && cx < c
            U02 = (cx + d * abspow(x, 2 + α - p)) / π
        else
            U02 = (c + d * abspow(x, 2 + α - p)) / π
        end

        return inv_u0(x) * U02
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
        T021(u0, Ball(), ubound(Arb, x + δ2), x; ϵ, N)
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

Compute the (not yet existing) integral ``T_{0,2,2}`` from the paper.

This is the integral
```
∫abs(clausenc(y - x, -α) + clausenc(y + x, -α) - 2clausenc(y, -α)) * y^p dy
```
from `a` to `π`. We should have `x < a < π`.

Notice that due to lemma [`lemma_integrand_2`](@ref) the expression
inside the absolute value is always positive, so we can remove the
absolute value.
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
        x < Arblib.realref(y) || return indeterminate(y)

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
