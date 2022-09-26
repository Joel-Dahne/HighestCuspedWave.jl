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
    f = T021(u0, Ball(), skip_div_u0 = true)

    return x -> begin
        a = ubound(Arb, x + δ2)

        # Compute with T021 on the whole interval
        res = f(x, Arb(π))

        if !(π < a || π - x < ϵ)
            part1 = f(x, a)

            part2 = T022(u0, Ball(), a, x, skip_div_u0 = true)

            res = intersect(part1 + part2, res)
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
- `return_enclosure::Bool` if true it returns an enclosure instead of
  an upper bound by returning the interval between zero and the
  computer upper bound.

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
function T02(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
    return_enclosure::Bool = false,
)
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
        # with a different solution. It is not used at the moment. It
        # makes the computed enclosure non-monotone in x which leads
        # to issues in the way ϵ is currently chose in D0_bounded_by.

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

        if false #isfinite(cx) && cx < c
            U02 = (cx + d * abspow(x, 2 + α - p)) / π
        else
            U02 = (c + d * abspow(x, 2 + α - p)) / π
        end

        res = inv_u0(x) * U02

        if return_enclosure
            return union(zero(res), res)
        else
            return res
        end
    end
end

"""
    T021(u0::FractionalKdVAnstaz{Arb}, ::Ball; skip_div_u0 = false)

Return a function such that `T021(u0)(x, a)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * ∫ (clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * u0.w(y) dy
```
where the integration is taken from `x` to `a`.

Since `u0.w` is bounded on the interval of integration we can factor
it out as
```
inv(π * u0(x) * u0.w(x)) * u0.w(Arb((x, a))) * ∫ clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α) dy
```
and compute the integral explicitly.
```
∫ clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α) dy = (
    -clausens(x - a, 1 - α) +
    clausens(x + a, 1 - α) -
    2clausens(a, 1 - α) -
    clausens(2x, 1 - α) +
    2clausens(x, 1 - α)
)
```
Using that `Arb((x, a)) = x * Arb((1, a / x))` and `u0.w(x * y) =
u0.w(x) * u0.w(y)` we can simplify the result to
```
inv(π * u0(x)) * u0.w(Arb((1, a / x))) * (
    -clausens(x - a, 1 - α) +
    clausens(x + a, 1 - α) -
    2clausens(a, 1 - α) -
    clausens(2x, 1 - α) +
    2clausens(x, 1 - α)
)
```

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.

**IMPROVE:** An earlier version of this method worked by expanding the
integrand and integrating the expansion. This was more complicated and
in general it gave worse results. However, for `x` close to zero it
could give better results. If need be we could reintroduce this
method, though it is most likely a better idea to do that in the
`::Asymptotic` version instead.
"""
function T021(u0::FractionalKdVAnsatz{Arb}, ::Ball = Ball(); skip_div_u0 = false)
    return (x::Arb, a::Arb) -> begin
        weight_factor = u0.w(Arb((1, a / x)))

        # Compute a tighter enclosure by expanding in x. If x is
        # closer to zero it is beneficial to use a higher degree when
        # computing the enclosure.
        degree = ifelse(x < 0.5, 4, 1)
        s = 1 - u0.α
        clausen_factor = ArbExtras.enclosure_series(x; degree) do x
            -clausens(x - a, s) + clausens(x + a, s) - 2clausens(a, s) - clausens(2x, s) + 2clausens(x, s)
        end

        res = weight_factor * clausen_factor / π

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
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
