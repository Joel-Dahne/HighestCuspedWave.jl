"""
    T02(u0::FractionalKdVAnsatz{Arb}; δ2::Arb = Arb(1e-2), ϵ::Arb = 1 + u0.α, skip_div_u0 = false)

Return a function such that `T02(u0; δ2)(x)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * ∫ (clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * u0.w(y) dy
```
where the integration is taken from `x` to `π`.

The interval of integration is split into two parts, one from `x` to
`x + a` and one from `x + a` to `π`. The first part is handled using
[`T021`](@ref) and the second one using [`T022`](@ref). The
computation is done two times with different values of `a` and the
intersection of the results is taken.

The result is first computed using `a = π`, in which case only the
first interval is used. It is then done using `a = ubound(x + δ2)`. If
`x` is close to `π`, as determined by `π - x < ϵ`, or if `ubound(x +
δ2)` is larger than `π` then only a = π` is used.

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.

**IMPROVE:** The choice of both `δ2` and `ϵ` can be tuned. For `δ2` it
depends on both `x` and `α`, whereas for `ϵ` it only depends on `α`.
For `δ2` it should be larger when both `x` and `α` are larger. In
theory we could compute with several values and take the best result,
but that would likely be to costly. For `ϵ` it's mainly a question of
cost, we don't want to compute the expensive integral when we believe
the asymptotic expansion will be the best anyway, larger values of `α`
should give a lower value of `ϵ`. The choice of `ϵ` is partially based
on [/figures/optimal-epsilon-choice.png], but can likely be tuned
further.
"""
function T02(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball;
    δ2::Arb = Arb(1e-2),
    ϵ::Arb = Arb(0.1),
    skip_div_u0 = false,
)
    f = T021(u0, Ball(), skip_div_u0 = true)
    g = T022(u0, Ball(), skip_div_u0 = true)

    return x -> begin
        a = ubound(Arb, x + δ2)

        # Compute with T021 on the whole interval
        res = f(x, Arb(π))

        if !(π < a || π - x < ϵ)
            part1 = f(x, a)

            part2 = g(x, a)

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
d = -gamma(1 + α) * sinpi(-α / 2) * (1 + α) * (2 + α) / (2 + α - p) / π^(2 + α - p) +
    2π^(p - 1) * do m
        (-1)^m *
        zeta(-α - 2m) *
        Arb(π)^2m /
        factorial(2m) *
        sum(binomial(2m, 2k) * (ϵ / π)^(2(m - 1 - k)) / (2k + 1 + p) for k = 0:m-1)
    end +
    6π^(p - 1) * sum((-1)^m * zeta(-α - 2m) * (3π / 2)^2m / factorial(2m) for m = N:Inf)
```
for any `N >= 1`.

To compute an enclosure of the tail of `d` we note that it is the same
as the sum in [`clausenc_expansion_remainder`](@ref) with `x = 3π / 2`.

# Notes
The expression for the integral used when determining `c` and `d`
where computed with Mathematica using the following code
```
Integrate[((t - 1)^(-a - 1) + (1 + t)^(-a - 1) - 2 t^(-a - 1))*t^p, {t, 1, Pi/x}, Assumptions -> a < 0 && 0 < x < Pi]
```
and using that
```
beta_inc(α - p, -α, z) = z^(α - p) / (α - p) * hypgeom_2f1(1 + α, α - p, 1 + α - p, z)
```

**IMPROVE:** We could get a slightly better enclosure by keeping more
terms from the expansions of the 2F1 function, though in practice this
would likely not help much.
"""
function T02(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
    return_enclosure::Bool = false,
)
    # This is required for the bound of the tail of d
    ϵ <= Arb(π) / 2 || throw(ArgumentError("we require that ϵ <= π / 2"))

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    α = u0.α
    p = u0.p

    # This is a requirement in the lemma giving the bound. In practice
    # we would just get an indeterminate result if it doesn't hold.
    Arblib.overlaps(1 + α, p) && throw(ArgumentError("we require that 1 + α != p"))

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

    d =
        -gamma(1 + α) * sinpi(-α / 2) * (1 + α) * (2 + α) / (2 + α - p) /
        (Arb(π)^(2 + α - p))

    # Take a lot of terms since the remainder is computed at a large
    # value
    N = 30
    # Sum first N - 1 terms
    d +=
        2Arb(π)^(p - 1) * sum(1:N-1) do m
            (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) *
            sum(binomial(2m, 2k) * (ϵ / π)^(2(m - 1 - k)) / (2k + 1 + p) for k = 0:m-1)
        end
    # Enclose remainder
    d +=
        6Arb(π)^(p - 1) *
        (3Arb(π) / 2)^(2N) *
        clausenc_expansion_remainder(3Arb(π) / 2, -α, N)

    return x::Arb -> begin
        @assert x <= ϵ

        U02 = c + d * abspow(x, 2 + α - p)

        res = inv_u0(x) * U02 / π

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

**IMPROVE:** For `x` overlapping or very close to `π` the Clausen
functions are not differentiable and the computed enclosure is very
wide. This could possibly be improved by expanding them at zero.
"""
function T021(u0::FractionalKdVAnsatz{Arb}, ::Ball = Ball(); skip_div_u0 = false)
    return (x::Arb, a::Arb) -> begin
        weight_factor = u0.w(Arb((1, a / x)))

        # Compute a tighter enclosure by expanding in x. If x is
        # closer to zero or π it is beneficial to use a higher degree
        # when computing the enclosure.
        degree = ifelse(0.5 < x < 3, 1, 4)
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
    T022(u0::FractionalKdVAnsatz{Arb}, ::Ball; skip_div_u0 = true)

Returns a functions such that `T022(u0)(x, a)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * ∫ (clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * u0.w(y) dy
```
where the integration is taken from `a` to `π`.

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.
"""
function T022(u0::FractionalKdVAnsatz{Arb}, ::Ball = Ball(); skip_div_u0 = false)
    return (x::Arb, a::Arb) -> begin
        integrand(y) = begin
            ry = real(y)

            # The integrand is singular at y = x so check that the real
            # part of y is greater than x or return an indeterminate
            # result. This also ensures that y^u0.p is analytic on all of
            # y.
            x < ry || return indeterminate(y)

            if isreal(y)
                return (
                    clausenc(ry - x, -u0.α) + clausenc(x + ry, -u0.α) -
                    2clausenc(ry, -u0.α)
                ) * y^u0.p
            else
                return (
                    clausenc(x - y, -u0.α) + clausenc(x + y, -u0.α) -
                    2clausenc(Acb(y), -u0.α)
                ) * y^u0.p
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
            return res / (π * u0(x) * u0.w(x))
        end
    end
end
