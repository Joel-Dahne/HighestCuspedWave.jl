"""
    T0_bhkdv(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic)

Compute `T0(x)` when `u0.use_bhkdv` is true.

It uses that the weight is given by `x^u0.p * log(2ℯ + inv(x)).

# Implementation
It first splits the function as
```
T0(x) = inv(π) * inv(u0(x) / x^-α) * (log(inv(x)) * x^p / u0.w(x)) * (U0(x) / (log(inv(x)) * x^(-α + p)))
```
where
```
U0(x) = x * ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * (x * t)^p * log(2ℯ + inv(x * t)) dt
```
integrated from `0` to `π / x`. It computes `inv(u0(x) / x^-u0.α)`
using [`inv_u0_normalised`](@ref) and
```
log(inv(x)) * x^u0.p / u0.w(x) = log(inv(x)) / log(2ℯ + inv(x))
```
using that it is decreasing in `x` and `1` for `x = 0`. What remains
is to bound `U0 / (log(inv(x)) * x^(-α + p))`.

# Bounding `U0 / (log(inv(x)) * x^(-α + p))`
We denote `U0 / (log(inv(x)) * x^(-α + p))` by `I(x)`.

We start by using that
```
log(2ℯ + inv(x * t)) = log(1 + 2ℯ * x * t) + log(inv(x)) + log(inv(t))
```
to split `U0` into
```
U01(x) = x^(1 + p) * ∫ abs(_integrand_I_hat(x, t, α)) * t^p * log(1 + 2ℯ * x * t) dt
U02(x) = log(inv(x)) * x^(1 + p) *∫ abs(_integrand_I_hat(x, t, α)) * t^p dt
U03(x) = x^(1 + p) *∫ abs(_integrand_I_hat(x, t, α)) * t^p * log(inv(t)) dt
```
with `U0(x) = U01(x) + U02(x) + U03(x)`. This allows us
to split `I(x)` into
```
I1(x) = x^(1 + α) / log(inv(x)) * ∫ abs(_integrand_I_hat(x, t, α)) * t^p * log(1 + 2ℯ * x * t) dt
I2(x) = x^(1 + α) * ∫ abs(_integrand_I_hat(x, t, α)) * t^p dt
I3(x) = x^(1 + α) / log(inv(x)) *∫ abs(_integrand_I_hat(x, t, α)) * t^p * log(inv(t)) dt
```
Satisfying `I(x) = I1(x) + I2(x) + I3(x)`. Bounds for these are
implemented in [`_T0_bhkdv_I1`](@ref), [`_T0_bhkdv_I2`](@ref) and
[`_T0_bhkdv_I3`](@ref) respectively.
"""
function T0_bhkdv(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
)
    # Function for enclosing inv(u0(x) / x^-α)
    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    # Function for enclosing
    # log(inv(x)) * x^u0.p / u0.w(x) = log(inv(x)) / log(2ℯ + inv(x))
    f = x -> if iszero(x)
        one(x)
    elseif Arblib.contains_zero(x)
        lower = let xᵤ = ubound(Arb, x)
            log(inv(xᵤ)) / log(2Arb(ℯ) + inv(xᵤ))
        end
        upper = one(x)
        Arb((lower, upper))
    else
        log(inv(x)) / log(2Arb(ℯ) + inv(x))
    end

    I1 = _T0_bhkdv_I1(u0.α, u0.p, ϵ)
    I2 = _T0_bhkdv_I2(u0.α, u0.p, ϵ)
    I3 = _T0_bhkdv_I3(u0.α, u0.p, ϵ)

    return x::Arb -> begin
        return inv(π) * inv_u0(x) * f(x) * (I1(x) + I2(x) + I3(x))
    end
end

"""
    _T0_bhkdv_I1(α, p, ϵ)

Return a function for computing a bound of
```
x^(1 + α) / log(inv(x)) * ∫ abs(_integrand_I_hat(x, t, α)) * t^p * log(1 + 2ℯ * x * t) dt
```
integrated from `0` to `π / x` for `x < ϵ`.

# Implementation
The factor `log(1 + 2ℯ * x * t)` is bounded on the interval of
integration and the idea is to factor it out and compute the remaining
integral.

Focusing on the part
```
x^(1 + α) * ∫ abs(_integrand_I_hat(x, t, α)) * t^p * log(1 + 2ℯ * x * t) dt
```
We first split it as
```
x^(1 + α) * ∫_0^1 abs(_integrand_I_hat(x, t, α)) * t^p * log(1 + 2ℯ * x * t) dt
+ x^(1 + α) * ∫_1^(π / x)( abs(_integrand_I_hat(x, t, α)) * t^p * log(1 + 2ℯ * x * t) dt
```
Factoring out the logarithm we get
```
log(1 + 2ℯ * x * Arb((0, 1))) * x^(1 + α) * ∫_0^1 abs(_integrand_I_hat(x, t, α)) * t^p dt
+ log(1 + 2ℯ * Arb((x, π)) * x^(1 + α) * ∫_1^(π / x) abs(_integrand_I_hat(x, t, α)) * t^p dt
```
We can compute
```
x^(1 + α) * ∫_0^1 abs(_integrand_I_hat(x, t, α)) * t^p dt
```
and
```
x^(1 + α) * ∫_1^(π / x) abs(_integrand_I_hat(x, t, α)) * t^p dt
```
using [`_T0_bhkdv_I2`](@ref) with `return_parts` set to true.

For the integral from `0` to `1` this gives a good enough enclosure.
For the integration from `1` to `π / x` this only give a good enough
enclosure for small values of `x`, around `x < 1e-5`. For larger
values of `x` we have to work a bit harder to get a good enough
enclosure.

# Handling `x > 1e-5`
In this case we integrate numerically. From `1` to `a` for some `a`
close to `1` we factor out the log-factor and integrate explicitly.
From `a` to `π / x` we integrate numerically.

In this case `_integrand_I_hat` is positive and we can remove the
absolute value. For the numerical integration we thus want to compute
```
∫ _integrand_I_hat(x, t, α) * t^p * log(1 + 2ℯ * x * t) dt
```

For the integration from `1` to `a` we factor out the log-factor as
well as the factor `t^(1 - p)`, leaving us with
```
Arb((1, a))^(1 - p) * log(1 + 2ℯ * x * Arb((1, a))) * ∫ _integrand_I_hat(x, t, α) * t^p dt
```
The integral can be computed explicitly using that the primitive
function is
```
(clausenc(x * (1 - t), 2 - α) / x^2 - t * clausens(x * (1 - t), 1 - α) / x) +
    (clausenc(x * (1 + t), 2 - α) / x^2 + t * clausens(x * (1 + t), 1 - α) / x) -
    2(clausenc(x * t, 2 - α) / x^2 + t * clausens(x * t, 1 - α) / x)
```
Which can be simplified to
```
(
    (
        clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) -
        2clausenc(x * t, 2 - α)
    ) / x + t * (
        -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) -
        2clausens(x * t, 1 - α)
    )
) / x
```
"""
function _T0_bhkdv_I1(α, p, ϵ)
    I2 = _T0_bhkdv_I2(α, p, ϵ, return_parts = true)

    # Precompute expansion to use when computing integrand from a to π
    M = 10
    C, _, P, _ = clausenc_expansion(Arb(0), -α, M, skip_constant = true)

    return x -> begin
        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        I21, I22 = I2(x)

        # Enclosure of x^(1 + α) * ∫_0^1 abs(_integrand_I_hat(x, t, α)) * t^p * log(1 + 2ℯ * x * t) dt
        part1 = log(1 + 2Arb(ℯ) * x * Arb((0, 1))) * I21

        # Enclosure of x^(1 + α) * ∫_1^(π / x) _integrand_I_hat(x, t, α) * t^p * log(1 + 2ℯ * x * t) dt
        if x > 1e-5
            a = Arb(1.1)

            primitive(x, t) =
                (
                    (
                        clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) -
                        2clausenc(x * t, 2 - α)
                    ) / x +
                    t * (
                        -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) -
                        2clausens(x * t, 1 - α)
                    )
                ) / x


            # Integral from 1 to a
            part21 =
                Arb((1, a))^(1 - p) *
                log(1 + 2Arb(ℯ) * x * Arb((1, a))) *
                ArbExtras.enclosure_series(x, degree = 4) do x
                    primitive(x, a) - primitive(x, one(a))
                end

            # Integral from a to π / x
            integrand(t) = begin
                rt = real(t)

                # The integrand is singular at t = 1 so check that the
                # real part of t is greater than 1 or return an
                # indeterminate result. This also ensures that the
                # integrand is analytic everywhere
                rt > 1 || return indeterminate(t)

                if isreal(t)
                    if x * (1 + rt) < 2Arb(π)
                        E = clausenc_expansion_remainder(x * (1 + rt), -α, M)
                        _integrand_I_hat_series(x, rt, α, C, P, E) *
                        rt^p *
                        log(1 + 2Arb(ℯ) * x * rt)
                    else
                        # In some cases when x and t is wide we get
                        # that x * (1 + real(t)) < 2Arb(π) and the
                        # remainder term is unbounded. In this cases
                        # fall back to direct evaluation.
                        ArbExtras.enclosure_series(x, degree = 8) do x
                            _integrand_I_hat(x, rt, α) * rt^p * log(1 + 2Arb(ℯ) * x * rt)
                        end
                    end
                else
                    _integrand_I_hat(x, t, α) * t^p * log(1 + 2Arb(ℯ) * x * t)
                end
            end

            part22_lower = real(
                Arblib.integrate(
                    integrand,
                    a,
                    lbound(Arb, π / x),
                    rtol = 1e-5,
                    warn_on_no_convergence = false,
                    opts = Arblib.calc_integrate_opt_struct(0, 500, 0, 0, 0),
                ),
            )
            part22_upper = real(
                Arblib.integrate(
                    integrand,
                    a,
                    ubound(Arb, π / x),
                    rtol = 1e-5,
                    warn_on_no_convergence = false,
                    opts = Arblib.calc_integrate_opt_struct(0, 500, 0, 0, 0),
                ),
            )
            part22 = Arb((part22_lower, part22_upper))

            part2 = x^(1 + α) * (part21 + part22)
        else
            part2 = log(1 + 2Arb(ℯ) * Arb((x, π))) * I22
        end

        return invloginvx * (part1 + part2)
    end
end

"""
    _T0_bhkdv_I2(α, p, ϵ; return_parts = false)

Return a function for computing a bound of
```
x^(1 + α) * ∫ abs(_integrand_I_hat(x, t, α)) * t^p dt
```
integrated from `0` to `π / x` for `x < ϵ`.

If `return_parts` is true then return the integration from `0` to `1`
and the integration from `1` to `π / x` separately. This is used in
[`_T0_bhkdv_I1`](@ref).

# Implementation
This is the same integral as when the weight is `x^p` and we can
therefore use same method to compute it. This is a combination of the
asymptotic version of `T01` and `T02` for that weight.

# Improvements
To get better enclosures we do some changes to better handle `α` close
to `-1`. The main problem is the sum
```
c2 + d2 * abspow(x, 2 + α - p)
```
which has large cancellations. To handle this we split `d2` into two
parts
```
d2_part1 = -gamma(1 + α) * sinpi(-α / 2) * (1 + α) * (2 + α) / (2 + α - p) / (Arb(π)^(2 + α - p))
```
and `d2_part2` is the rest. We then factor out `gamma(1 + α) *
sinpi(-α / 2)` from both `c2` and `d2_part1`, giving us
```
c2_modified = (
    gamma(-α) * gamma(α - p + 1) / gamma(-p) +
    hypgeom_2f1(1 + α, α - p, 1 + α - p, -one(α)) -
    2
) / (α - p)
d2_part1_modified = -(1 + α) * (2 + α) / (2 + α - p) / (Arb(π)^(2 + α - p))
```
and
```
c2 + d2_part1 * abspow(x, 2 + α - p) =
    gamma(1 + α) * sinpi(-α / 2) * (c2_modified + d2_part1_modified * abspow(x, 2 + α - p))
```
"""
function _T0_bhkdv_I2(α, p, ϵ; return_parts = false)
    # This is required for the bound of the tail of d2
    ϵ <= Arb(π) / 2 || throw(ArgumentError("we require that ϵ <= π / 2"))

    # This is a requirement in the lemma giving the bound. In practice
    # we would just get an indeterminate result if it doesn't hold.
    Arblib.overlaps(1 + α, p) && throw(ArgumentError("we require that 1 + α != p"))

    r = _integrand_compute_root(FractionalKdVAnsatz, zero(α), α)

    # It is important to get good enclosures of both c1 and c2 for
    # wide values of α. For c1 this is done by rewriting it to work
    # better for α close to -1 and for c2 by bisection.
    c1 =
        gamma(1 + α) *
        sinpi(-α / 2) *
        (
            2 / (α - p) +
            gamma(-α) * gamma(1 + p) / gamma(1 - α + p) +
            hypgeom_2f1(1 + α, 1 + p, 2 + p, -one(α)) / (1 + p) -
            2r^(1 + p) * (
                2r^(-α - 1) / (α - p) +
                (
                    hypgeom_2f1(1 + α, 1 + p, 2 + p, -r) +
                    hypgeom_2f1(1 + α, 1 + p, 2 + p, r)
                ) / (1 + p)
            )
        )

    c2_modified =
        ArbExtras.extrema_enclosure(getinterval(α)..., degree = -1) do α
            (
                gamma(-α) * gamma(α - p + 1) / gamma(-p) +
                hypgeom_2f1(1 + α, α - p, 1 + α - p, -one(α)) - 2
            ) / (α - p)
        end |> Arb

    # c2 has a removable singularity for p = 1 and for some values of
    # α for p < 1. In practice we don't encounter these values because
    # those are covered by T0_p_one.
    isfinite(c2_modified) ||
        @error "non-finite enclosure for c2_modified in _T0_bhkdv_12" α p

    d1 = let N = 10, d1 = zero(α)
        # Sum first N - 1 terms
        d1 +=
            2sum(1:N-1) do m
                (-1)^m * zeta(-α - 2m) * ϵ^(2m - 2) / factorial(2m) *
                sum(binomial(2m, 2k) / (2k + 1 + p) for k = 0:m-1)
            end

        # Enclose remainder
        # Note that 4(2ϵ)^(2M - 2) = (2ϵ)^(2M) / ϵ^2
        d1 += 4(2ϵ)^(2N - 2) * clausenc_expansion_remainder(2ϵ, -α, N)

        d1
    end

    d2 = let N = 30, d2 = zero(α)
        d2 += ArbExtras.enclosure_series(α) do α
            -gamma(1 + α) * sinpi(-α / 2) * (1 + α) * (2 + α) / (2 + α - p) /
            (Arb(π)^(2 + α - p))
        end

        # Sum first N - 1 terms
        d2 +=
            2Arb(π)^(p - 1) * sum(1:N-1) do m
                (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) * sum(
                    binomial(2m, 2k) * (ϵ / π)^(2(m - 1 - k)) / (2k + 1 + p) for k = 0:m-1
                )
            end

        # Enclose remainder
        d2 +=
            6Arb(π)^(p - 1) *
            (3Arb(π) / 2)^(2N) *
            clausenc_expansion_remainder(3Arb(π) / 2, -α, N)

        d2
    end

    d2_part1_modified = ArbExtras.enclosure_series(α) do α
        -(1 + α) * (2 + α) / (2 + α - p) / (Arb(π)^(2 + α - p))
    end

    d2_part2 = let N = 30, d2 = zero(α)
        # Sum first N - 1 terms
        d2 +=
            2Arb(π)^(p - 1) * sum(1:N-1) do m
                (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) * sum(
                    binomial(2m, 2k) * (ϵ / π)^(2(m - 1 - k)) / (2k + 1 + p) for k = 0:m-1
                )
            end

        # Enclose remainder
        d2 +=
            6Arb(π)^(p - 1) *
            (3Arb(π) / 2)^(2N) *
            clausenc_expansion_remainder(3Arb(π) / 2, -α, N)

        d2
    end

    return x::Arb -> begin
        @assert x <= ϵ

        part1 = c1 + d1 * abspow(x, 3 + α)
        part2 =
            gamma(1 + α) *
            sinpi(-α / 2) *
            (c2_modified + d2_part1_modified * abspow(x, 2 + α - p)) +
            d2_part2 * abspow(x, 2 + α - p)

        if return_parts
            return part1, part2
        else
            return part1 + part2
        end
    end
end

"""
    _T0_bhkdv_I3(α, p, ϵ)

Return a function for computing a bound of
```
x^(1 + α) / log(inv(x)) * ∫ abs(_integrand_I_hat(x, t, α)) * t^p * log(inv(t)) dt
```
integrated from `0` to `π / x` for `x < ϵ`.

# Implementation
Expanding the Clausen functions in `_integrand_I_hat(x, t, α)` at zero
allows us to write it as
```
_integrand_I_hat(x, t, α) =
    -gamma(1 + α) * sinpi(-α / 2) * (abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * x^(-α - 1)
    + sum((-1)^m * zeta(-α - 2m) * ((1 - t)^2m + (1 + t)^2m - 2t^2m) * x^2m / factorial(2m) for m = 1:Inf)
```
Using that
```
(1 - t)^2m + (1 + t)^2m - 2t^2m = 2sum(binomial(2m, 2k) * t^2k for k = 0:m-1)
```
we can rewrite that sum as
```
sum((-1)^m * zeta(-α - 2m) * ((1 - t)^2m + (1 + t)^2m - 2t^2m) * x^2m / factorial(2m) for m = 1:Inf) =
    2sum((-1)^m * zeta(-α - 2m) * x^2m / factorial(2m) * sum(binomial(2m, 2k) * t^2k for k = 0:m-1) for m = 1:Inf)
```

From this we get
```
I3(x) <= gamma(1 + α) * sinpi(-α / 2) / log(inv(x)) *
    ∫ abs((abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))) * t^p * log(inv(t)) dt +
    2x^(1 + α) / log(inv(x)) * sum((-1)^m * zeta(-α - 2m) * x^2m / factorial(2m) * sum(binomial(2m, 2k) * ∫ t^(2k + p) * log(inv(t)) dt for k = 0:m-1) for m = 1:Inf)
```
We denote the first term by `I3_M(x)` and the second by `I3_R(x)`,
bounds for them are implemented in [`_T0_bhkdv_I3_M`](@ref) and
[`_T0_bhkdv_I3_R`](@ref) respectively.
"""
function _T0_bhkdv_I3(α, p, ϵ)
    I3_M = _T0_bhkdv_I3_M(α, p, ϵ)
    I3_R = _T0_bhkdv_I3_R(α, p, ϵ)

    return x -> I3_M(x) + I3_R(x)
end

"""
    _T0_bhkdv_I3_M(α, p, ϵ)

Return a function for computing a bound of
```
gamma(1 + α) * sinpi(-α / 2) / log(inv(x)) * ∫ abs((abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))) * t^p * log(inv(t)) dt
```
integrated from `0` to `π / x` for `x < ϵ`.

# Implementation
We split the integration into two parts, one from `0` to `1` and one
from `1` to `π / x`. They are implemented in [`_T0_bhkdv_I3_M1`](@ref) and
[`_T0_bhkdv_I3_M2`](@ref) respectively.
"""
function _T0_bhkdv_I3_M(α, p, ϵ)
    I3_M1 = _T0_bhkdv_I3_M1(α, p, ϵ)
    I3_M2 = _T0_bhkdv_I3_M2(α, p, ϵ)

    return x -> I3_M1(x) + I3_M2(x)
end

"""
    _T0_bhkdv_I3_M1(α, p, ϵ)

Return a function for computing a bound of
```
gamma(1 + α) * sinpi(-α / 2) / log(inv(x)) * ∫ abs((abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))) * t^p * log(inv(t)) dt
```
integrated from `0` to `1` for `x < ϵ`.

# Implementation
To begin with we note that `1 - t` is positive and we can remove the
absolute value for `abs(1 - t)`.

To compute the integral we integrate numerically from `0` to `a` for
some `a` close to `1`. From `a` to `1` we use that `log(inv(t))` is
bounded to factor it out and integrate explicitly.

The only complication for the numerical integration is enclosing the
integrand close to zero. For this we rewrite it as
```
abs(((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))) * t^p * log(inv(t))
= -abs((t^(1 // 2) * ((1 - t)^(-α - 1) + (1 + t)^(-α - 1)) - 2t^(-α - 1 // 2))) * t^(p - 1 // 2) * log(t)
```
The factor with the absolute value can then be enclosed directly. For
the remaining factors we note that they are equal to `-logabspow(t, 1,
p - 1 // 2)`, which can then be used to bound them.

For the integration from `a` to `1` we factor out an enclosure of
`log(inv(t))`, giving us
```
log(inv(Arb((a, 1)))) * ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))) * t^p dt
```
Since
```
abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))
```
is increasing in `t` it is enough to check that it is positive at `t =
a` to be able to remove the absolute value. Doing that we are left integrating
```
∫ ((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p dt
```
The primitive function is given by
```
(beta_inc(1 + p, -α, t) - inv(α - p) * inv(t)^(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, -inv(t)) - 2))
```
Integrating from `a` to `1`we thus get
```
(beta_inc(1 + p, -α, 1) - inv(α - p) * inv(1)^(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, -inv(1)) - 2)) -
(beta_inc(1 + p, -α, a) - inv(α - p) * inv(a)^(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, -inv(a)) - 2))
```
which simplifies to
```
(beta_inc(1 + p, -α, 1) - inv(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, -1) - 2)) -
(beta_inc(1 + p, -α, a) - inv(α - p) * inv(a)^(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, -inv(a)) - 2))
```
"""
function _T0_bhkdv_I3_M1(α, p, ϵ)
    C = gamma(1 + α) * sinpi(-α / 2)

    J = begin
        integrand(t; analytic = false) = begin
            if isreal(t)
                t = real(t)
                @assert !analytic
                if Arblib.contains_zero(t)
                    -abs(
                        abspow(t, Arb(1 // 2)) * ((1 - t)^(-α - 1) + (1 + t)^(-α - 1)) -
                        2abspow(t, -α - 1 // 2),
                    ) * logabspow(t, 1, p - 1 // 2)
                else
                    return -abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                           t^p *
                           log(t)
                end
            else
                res = (1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
                Arblib.real_abs!(res, res, analytic)

                return -res * t^p * log(t)
            end
        end

        a = Arb(0.999)

        # Integral from 0 to a
        part1 =
            real(Arblib.integrate(integrand, 0, a, rtol = 1e-5, check_analytic = true))

        # Check that the integrand is positive after removing the
        # absolute value
        @assert Arblib.ispositive((1 - a)^(-α - 1) + (1 + a)^(-α - 1) - 2a^(-α - 1))

        # Integral from a to 1
        part2 =
            log(inv(Arb((a, 1)))) * (
                (
                    beta_inc(1 + p, -α, one(α)) -
                    inv(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, -one(α)) - 2)
                ) - (
                    beta_inc(1 + p, -α, a) -
                    inv(α - p) *
                    inv(a)^(α - p) *
                    (hypgeom_2f1(1 + α, α - p, 1 + α - p, -inv(a)) - 2)
                )
            )

        part1 + part2
    end

    return x::Arb -> begin
        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        return C * invloginvx * J
    end
end

"""
    _T0_bhkdv_I3_M2(α, p, ϵ)

Return a function for computing a bound of
```
gamma(1 + α) * sinpi(-α / 2) / log(inv(x)) * ∫ abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p * log(inv(t)) dt
```
integrated from `1` to `π / x` for `x < ϵ`.

# Implementation
To begin with we note that `1 - t` is negative and we can remove the
absolute value for `abs(1 - t)` and put a minus sign, furthermore the
outer absolute value can be removed. This gives us the integral
```
J = ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p * log(inv(t)) dt
```

Since `log(inv(t))` is negative for `t > 1` this is decreasing in `x`.
To get an upper bound it is thus enough to work with an upper bound of
`x`.

We split the integral into four parts. One from `1` to `a`, one from
`b` to `c` and one from `c` to `π / x`. Here `a = 1.001`, `b = min(π /
x, 1e10)` and `c = min(π / x, 1e50)`. We denote the four subintegrals
by `J1`, `J2`, `J3` and `J4`. On the different parts we use different
methods for computing the integral. Depending on the value of `x` it
might be that some of the intervals are empty.

# Primitive function
In several of the cases we factor out the logarithm and integrate
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p dt
```
explicitly. For this we use that the primitive function is given by
```
-inv(α - p) * inv(t)^(α - p) * (
    hypgeom_2f1(1 + α, α - p, 1 + α - p, inv(t)) +
    hypgeom_2f1(1 + α, α - p, 1 + α - p, -inv(t)) -
    2
)
```
We define
```
primitive_inv(a, b) = -inv(α - p) * (
    b^(α - p) * (
        hypgeom_2f1(1 + α, α - p, 1 + α - p, b) +
        hypgeom_2f1(1 + α, α - p, 1 + α - p, -b) -
        2
    )
) * (
    a^(α - p) * (
        hypgeom_2f1(1 + α, α - p, 1 + α - p, a) +
        hypgeom_2f1(1 + α, α - p, 1 + α - p, -a) -
        2
    )
)
```
such that `primitive_inv(inv(a), inv(b))` gives the integral from `a`
to `b`.

## Handling small `t`
For small values of `t` the function `primitive_inv(a, b)` has large
cancellations that needs to be accounted for. To handle this we
rewrite it using the series expansion.

From the series expansion of `hypgeom_2f1` we get
```
(hypgeom_2f1(1 + α, α - p, 1 + α - p, y) + hypgeom_2f1(1 + α, α - p, 1 + α - p, -y) - 2)
= 2sum(
    rising(1 + α, 2k) * rising(α - p, 2k) / rising(1 + α - p, 2k) / factorial(2k) * y^2k
    for k in 1:Inf
)
```
We can simplify the rising factorials to
```
rising(1 + α, 2k) * rising(α - p, 2k) / rising(1 + α - p, 2k) = rising(1 + α, 2k) * (α - p) / (2k + α - p)
```
Note that this is negative for `k >= 1` since `α - p` is negative and
the other two factors are positive. To enclose the sum we sum the
first `N-1` terms explicitly and then enclose the tail
```
2sum(
    rising(1 + α, 2k) * rising(α - p, 2k) / rising(1 + α - p, 2k) / factorial(2k) * y^2k
    for k in N:Inf
)
```
Since all the terms are negative a trivial upper bound is zero. For a
lower bound we add more terms to make it look like the tail for the
`hypgeom_2f1` function. Namely we use that for `y > 0`
```
sum(
    rising(1 + α, 2k) * rising(α - p, 2k) / rising(1 + α - p, 2k) / factorial(2k) * y^2k
    for k in N:Inf
)
>= sum(
    rising(1 + α, k) * rising(α - p, k) / rising(1 + α - p, k) / factorial(k) * y^k
    for k in 2N:Inf
)
```
We can then bound the last sum
[`Arblib.hypgeom_pfq_bound_factor!`](@ref) which returns `D` such that
```
abs(sum(
    rising(1 + α, k) * rising(α - p, k) / rising(1 + α - p, k) / factorial(k) * y^k
    for k in 2N:Inf
)) <= abs(D * rising(1 + α, 2N) * rising(α - p, 2N) / rising(1 + α - p, 2N) / factorial(2N) * y^2N)
```

# Integral from `1` to `a`
To compute this integral we use that close to `1` the factor
`log(inv(t))` is bounded and we can factor it out and integrate
explicitly. The integral in this case is given by
```
J1 = log(inv(Arb((1, a)))) * ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p dt
```
and we integrate explicitly using the primitive function given above.

# Integral from `a` to `b`
Done by numerically integrating
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p * log(inv(t)) dt
```

# Integral from `b` to `c`
This integral we split into several parts and factor out the
log-factor from each part separately. More precisely we split the
interval into 16 subintervals `subinterval[i] = [b[i], b[i+1]] and compute
```
sum(1:16) do i
    log(inv(Arb((b[i], b[i+1])))) * ∫_(b[i])^(b[i+1]) ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p dt
end
```
The subintervals are computed using
[`ArbExtras.bisect_interval_recursive`] by recursively bisecting at
the geometric midpoint.

# Integral from `c` to `π / x`
In this case we only compute an upper bound by noticing that the
integrand is negative and
```
log(inv(t)) < log(inv(c)) = -log(c)
```
Factoring out `-log(b)`. This gives us the integral
```
J4 <= -log(c) * ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p dt
```
integrated from `c` to `π / x`, which we integrate explicitly.
"""
function _T0_bhkdv_I3_M2(α, p, ϵ)
    C = gamma(1 + α) * sinpi(-α / 2)

    s = -α - 1

    integrand(t) =
        if isreal(t) && !iswide(t)
            rt = real(t)

            -t^p * log(t) * ArbExtras.enclosure_series(s, degree = 4) do s
                (rt - 1)^s + (1 + rt)^s - 2rt^s
            end
        else
            -((t - 1)^s + (1 + t)^s - 2t^s) * t^p * log(t)
        end

    # Compute
    # b^(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, b) + hypgeom_2f1(1 + α, α - p, 1 + α - p, -b) - 2) -
    # a^(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, a) + hypgeom_2f1(1 + α, α - p, 1 + α - p, -a) - 2)
    # Handling cancellations for small values of a and b
    hypgeom_2f1_helper(a, b) =
        if min(a, b) > 0.1
            # No need to do anything fancy
            b^(α - p) * (
                hypgeom_2f1(1 + α, α - p, 1 + α - p, b) +
                hypgeom_2f1(1 + α, α - p, 1 + α - p, -b) - 2
            ) -
            a^(α - p) * (
                hypgeom_2f1(1 + α, α - p, 1 + α - p, a) +
                hypgeom_2f1(1 + α, α - p, 1 + α - p, -a) - 2
            )
        else
            let N = 10
                main =
                    2(α - p) * sum(1:N-1) do k
                        rising(1 + α, 2k) / (2k + α - p) / factorial(2k) *
                        (b^(2k + α - p) - a^(2k + α - p))
                    end

                tail_lower = begin
                    @assert Arblib.isnonnegative(a) && Arblib.isnonnegative(b)
                    D = Arblib.Arblib.hypgeom_pfq_bound_factor!(
                        zero(Mag),
                        AcbVector([1 + α, α - p]),
                        2,
                        AcbVector([1 + α - p, 1]),
                        2,
                        Acb(max(a, b)),
                        UInt(2N),
                    )

                    -abs(
                        D * rising(1 + α, 2N) * rising(α - p, 2N) / rising(1 + α - p, 2N) / factorial(2N) *
                        (b^(2N + α - p) - a^(2N + α - p)),
                    )
                end

                # Tail is negative
                tail_upper = zero(main)

                return main + Arb((tail_lower, tail_upper))
            end
        end

    primitive_inv(a, b) = -inv(α - p) * hypgeom_2f1_helper(a, b)

    return x::Arb -> begin
        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        # J is decreasing in x so enough to compute with an upper
        # bound of `x
        xᵤ = ubound(Arb, x)

        a = Arb(1.001)
        if iszero(xᵤ)
            b = Arb(1e10)
            c = Arb(1e50)
        else
            b = ifelse(π / xᵤ < 1e10, π / xᵤ, Arb(1e10))
            c = ifelse(π / xᵤ < 1e50, π / xᵤ, Arb(1e50))
        end

        # Integration from 1 to a
        J1 = log(inv(Arb((1, a)))) * primitive_inv(one(a), inv(a))

        # Integration from a to b
        J2 = real(Arblib.integrate(integrand, a, b, rtol = 1e-5))

        # Integration from b to c
        J3 = if iszero(xᵤ) || b < π / xᵤ
            subintervals = ArbExtras.bisect_interval_recursive(b, c, 4, log_midpoint = true)

            sum(subintervals) do subinterval
                -Arb((log(subinterval[1]), log(subinterval[2]))) *
                primitive_inv(inv(subinterval[1]), inv(subinterval[2]))
            end
        else
            zero(J1)
        end

        # Integration from c to π / xᵤ
        J4 = if iszero(xᵤ) || c < π / xᵤ
            -log(c) * primitive_inv(inv(c), xᵤ / π)
        else
            zero(J1)
        end

        J = J1 + J2 + J3 + J4

        return C * invloginvx * J
    end
end

"""
    _T0_bhkdv_I3_R(α, p, ϵ)

Return a function for computing a bound of
```
2x^(1 + α) / log(inv(x)) * sum((-1)^m * zeta(-α - 2m) * x^2m / factorial(2m) * sum(binomial(2m, 2k) * ∫ t^(2k + p) * log(inv(t)) dt for k = 0:m-1) for m = 1:Inf)
```
where the integration is from `0` to `π / x` for `x < ϵ`.

# Implementation
We have
```
∫ t^(2k + p) * log(inv(t)) dt = (1 - (2k + 1 + p) * log(π / x)) / (2k + 1 + p)^2 * (π / x)^(2k + 1 + p)
```
Giving us
```
I3_R(x) = 2x^(1 + α) / log(inv(x)) * sum((-1)^m * zeta(-α - 2m) * x^2m / factorial(2m) * sum(binomial(2m, 2k) (1 - (2k + 1 + p) * log(π / x)) / (2k + 1 + p)^2 * (π / x)^(2k + 1 + p) for k = 0:m-1) for m = 1:Inf)
```
Similar simplifications as for when the weight is `x^p` allows us to
write it as
**IMPROVE**: Explain better
```
I3_R(x) = 2x^(2 + α - p) * π^(p - 1) / log(inv(x)) * sum((-1)^m * zeta(-α - 2m) * π^2m / factorial(2m) * sum(binomial(2m, 2k) * (1 - (2k + 1 + p) * log(π / x)) / (2k + 1 + p)^2 * (x / π)^((2(m - 1 - k))) for k = 0:m-1) for m = 1:Inf)
```

By splitting the inner sum into
```
sum(binomial(2m, 2k) / (2k + 1 + p)^2 * (x / π)^((2(m - 1 - k))) for k = 0:m-1) -
    -log(π / x) * sum(binomial(2m, 2k) / (2k + 1 + p) * (x / π)^((2(m - 1 - k))) for k = 0:m-1)
```
we can split `I3_R` into
```
I3_R1(x) = 2x^(2 + α - p) * π^(p - 1) / log(inv(x)) * sum((-1)^m * zeta(-α - 2m) * π^2m / factorial(2m) * sum(binomial(2m, 2k) / (2k + 1 + p)^2 * (x / π)^((2(m - 1 - k))) for k = 0:m-1) for m = 1:Inf)

I3_R2(x) = - 2x^(2 + α - p) * π^(p - 1) * log(π / x) / log(inv(x)) * sum((-1)^m * zeta(-α - 2m) * π^2m / factorial(2m) * sum(binomial(2m, 2k) / (2k + 1 + p) * (x / π)^((2(m - 1 - k))) for k = 0:m-1) for m = 1:Inf)
```
For `I3_R2` we can use that `log(π / x) / log(inv(x)) = 1 + log(π) /
log(inv(x))`. We can sum the first `N` terms explicitly, what remains
is to handle the tails
```
sum((-1)^m * zeta(-α - 2m) * π^2m / factorial(2m) * sum(binomial(2m, 2k) / (2k + 1 + p)^2 * (x / π)^((2(m - 1 - k))) for k = 0:m-1) for m = N:Inf)

sum((-1)^m * zeta(-α - 2m) * π^2m / factorial(2m) * sum(binomial(2m, 2k) / (2k + 1 + p) * (x / π)^((2(m - 1 - k))) for k = 0:m-1) for m = N:Inf)
```
For the first one we need to compute an upper bound, since the factor
in front is positive, whereas for the second one we need to compute a
lower bound since the factor in front is negative.

For the lower bound it is enough to notice that all terms in the sum
are positive. So we get the trivial lower bound zero. What remains is
thus to compute an upper bound for the first sum. Looking at the inner
sum we have
```
sum(binomial(2m, 2k) / (2k + 1 + p)^2 * (x / π)^((2(m - 1 - k))) for k = 0:m-1)
<= sum(binomial(2m, 2k) / (2k + 1 + p) * (x / π)^((2(m - 1 - k))) for k = 0:m-1)
<= sum(binomial(2m, 2k) / (2k + 1) * (x / π)^((2(m - 1 - k))) for k = 0:m-1)
<= sum(binomial(2m, 2k) / (2k + 1) * (ϵ / π)^((2(m - 1 - k))) for k = 0:m-1)
<= sum(binomial(2m, 2k) / (2k + 1) * (1 / 2)^((2(m - 1 - k))) for k = 0:m-1)
= (2^(-2m) * (1 + 3^(1 + 2m)) - 4) / (1 + 2m)
<= 3 * (3 / 2)^2m
```
Here we have used that `x < ϵ < π / 2` and explicitly computed the
last sum. Inserting this back into the outer sum we get
```
3 * sum((-1)^m * zeta(-α - 2m) * (3π / 2)^2m / factorial(2m) *  for m = N:Inf)
```
This is exactly `3clausenc_expansion_remainder(3Arb(π) / 2, -α, N)`.
"""
function _T0_bhkdv_I3_R(α, p, ϵ)
    # This is required for the bound of the tail of I3_R1
    @assert ϵ <= Arb(π) / 2

    return x -> begin
        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        # Enclosure of log(π / x) / log(inv(x)) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        N = 20

        # Sum in I3_R1 for m in 1:N-1
        I3_R1_sum = sum(1:N-1) do m
            term =
                (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) *
                sum(0:m-1) do k
                    binomial(2m, 2k) / (2k + 1 + p)^2 * (x / π)^(2(m - 1 - k))
                end
        end

        # Upper bound of tail for sum in I3_R1.
        I3_R1_sum_tail_upper = 3clausenc_expansion_remainder(3Arb(π) / 2, -α, N)

        I3_R1 =
            2abspow(x, 2 + α - p) *
            Arb(π)^(p - 1) *
            invloginvx *
            (I3_R1_sum + I3_R1_sum_tail_upper)


        # Sum in I3_R2 for m in 1:N-1
        I3_R2_sum = sum(1:N-1) do m
            term =
                (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) *
                sum(0:m-1) do k
                    binomial(2m, 2k) / (2k + 1 + p) * (x / π)^(2(m - 1 - k))
                end
        end

        # Lower bound of tail for sum in I3_R2. It is positive so a
        # trivial lower bound is zero
        I3_R2_sum_tail_lower = zero(α)

        # Note that log(π / x) / log(inv(x)) = 1 + log(π) / log(inv(x))
        I3_R2 =
            -2abspow(x, 2 + α - p) *
            Arb(π)^(p - 1) *
            (1 + log(Arb(π)) * invloginvx) *
            (I3_R2_sum + I3_R2_sum_tail_lower)


        return I3_R1 + I3_R2
    end
end
