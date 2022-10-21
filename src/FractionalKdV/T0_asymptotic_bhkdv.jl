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

For the integration from `1` to `a` we factor out the log-factor,
leaving us with
```
∫ _integrand_I_hat(x, t, α) * t^p dt
```
If we factor out `t^(1 - p)`, which is bounded, we get
```
Arb((1, a))^(1 - p) * ∫ _integrand_I_hat(x, t, α) * t dt
```
The integral can be computed explicitly using that the primitive
function is
```
(clausenc(x * (1 - t), 2 - α) / x^2 - t * clausens(x * (1 - t), 1 - α) / x) +
    (clausenc(x * (1 + t), 2 - α) / x^2 + t * clausens(x * (1 + t), 1 - α) / x) -
    2(clausenc(x * t, 2 - α) / x^2 + t * clausens(x * t, 1 - α) / x)
```
With our integration limits this gives us
```
(
    (clausenc(x * (1 - a), 2 - α) / x^2 - a * clausens(x * (1 - a), 1 - α) / x) +
    (clausenc(x * (1 + a), 2 - α) / x^2 + a * clausens(x * (1 + a), 1 - α) / x) -
    2(clausenc(x * a, 2 - α) / x^2 + a * clausens(x * a, 1 - α) / x)
) - (
    (clausenc(0, 2 - α) / x^2 - clausens(0, 1 - α) / x) +
    (clausenc(2x, 2 - α) / x^2 + clausens(2x, 1 - α) / x) -
    2(clausenc(x, 2 - α) / x^2 + clausens(x, 1 - α) / x)
)
```
"""
function _T0_bhkdv_I1(α, p, ϵ)
    I2 = _T0_bhkdv_I2(α, p, ϵ, return_parts = true)

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
        if x < 1e-5
            part2 = log(1 + 2Arb(ℯ) * Arb((x, π))) * I22
        else
            a = Arb(1.1)

            # Integral from 1 to a
            part21 =
                log(1 + 2Arb(ℯ) * x * Arb((1, a))) * (
                    (
                        (
                            clausenc(x * (1 - a), 2 - α) / x^2 -
                            a * clausens(x * (1 - a), 1 - α) / x
                        ) + (
                            clausenc(x * (1 + a), 2 - α) / x^2 +
                            a * clausens(x * (1 + a), 1 - α) / x
                        ) -
                        2(clausenc(x * a, 2 - α) / x^2 + a * clausens(x * a, 1 - α) / x)
                    ) - (
                        (clausenc(zero(x), 2 - α) / x^2 - clausens(zero(x), 1 - α) / x) + (clausenc(2x, 2 - α) / x^2 + clausens(2x, 1 - α) / x) -
                        2(clausenc(x, 2 - α) / x^2 + clausens(x, 1 - α) / x)
                    )
                )

            # Integral from a to π / x
            integrand(t) = _integrand_I_hat(x, t, α) * t^p * log(1 + 2Arb(ℯ) * x * t)
            part22 = real(Arblib.integrate(integrand, a, π / x, rtol = 1e-5))

            part2 = x^(1 + α) * (part21 + part22)
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
"""
function _T0_bhkdv_I2(α, p, ϵ; return_parts = false)
    # This is required for the bound of the tail of d2
    ϵ <= Arb(π) / 2 || throw(ArgumentError("we require that ϵ <= π / 2"))

    # This is a requirement in the lemma giving the bound. In practice
    # we would just get an indeterminate result if it doesn't hold.
    Arblib.overlaps(1 + α, p) && throw(ArgumentError("we require that 1 + α != p"))

    r = _integrand_compute_root(FractionalKdVAnsatz, zero(α), α)

    c1 =
        gamma(1 + α) *
        sinpi(-α / 2) *
        (
            2 / (α - p) +
            gamma(-α) * gamma(1 + p) / gamma(1 - α + p) +
            hypgeom_2f1(1 + α, 1 + p, 2 + p, -one(α)) / (1 + p) -
            2r^p * (
                2r^-α / (α - p) +
                r * hypgeom_2f1(1 + α, 1 + p, 2 + p, -r) / (1 + p) +
                r * hypgeom_2f1(1 + α, 1 + p, 2 + p, r) / (1 + p)
            )
        )

    c2 =
        gamma(1 + α) *
        sinpi(-α / 2) *
        (
            gamma(-α) * gamma(α - p) / gamma(-p) +
            (hypgeom_2f1(1 + α, α - p, 1 + α - p, -one(α)) - 2) / (α - p)
        )

    # c2 has a removable singularity for p = 1 and for some values of
    # α for p < 1. In practice we don't encounter these values because
    # those are covered by T0_p_one.
    isfinite(c2) || @error "non-finite enclosure for c in _T0_bhkdv_12" α p

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
        d2 +=
            -gamma(1 + α) * sinpi(-α / 2) * (1 + α) * (2 + α) / (2 + α - p) /
            (Arb(π)^(2 + α - p))

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
        if return_parts
            return (c1 + d1 * abspow(x, 3 + α)), (c2 + d2 * abspow(x, 2 + α - p))
        else
            return c1 + c2 + d1 * abspow(x, 3 + α) + d2 * abspow(x, 2 + α - p)
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
        part1 = real(Arblib.integrate(integrand, 0, a, check_analytic = true))

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

To compute the integral we integrate from `1` to `a` from some `a`
close to `1` using that `log(inv(t))` is bounded to factor it out and
integrate explicitly.

For the integral from `a` to `π / x` we integrate part of it
numerically and part of it explicitly. We integrate from `a` to `b =
min(π / x, 1e10)` numerically. If `b < π / x` we only compute an upper
bound of the integration from `b` to `π / x`. This is done by noticing
that the integrand is negative and we get an upper bound by using
```
log(inv(t)) < log(inv(b)) = -log(b)
```
and factoring out `-log(b)`. This gives us the integral
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p dt
```
from `b` to `π / x`, which we integrate explicitly.

The integral is given by
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p dt =
= inv(α - p) * inv(b)^(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, inv(b)) + hypgeom_2f1(1 + α, α - p, 1 + α - p, -inv(b)) - 2)
  -inv(α - p) * (x / π)^(α - p) * (hypgeom_2f1(1 + α, α - p, 1 + α - p, x / π) + hypgeom_2f1(1 + α, α - p, 1 + α - p, -x / π) - 2)
```
The `hypgeom_2f1` functions have large cancellations and need special
care.

"""
function _T0_bhkdv_I3_M2(α, p, ϵ)
    C = gamma(1 + α) * sinpi(-α / 2)

    integrand(t) = -((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p * log(t)
    integrand2(t) = ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^p

    a = Arb(1.001)

    # Integration from 1 to b
    J1 = zero(α) # TODO

    # Compute
    # hypgeom_2f1(1 + α, α - p, 1 + α - p, inv(b)) +
    # hypgeom_2f1(1 + α, α - p, 1 + α - p, -inv(b)) -
    # 2
    # Handling cancellations
    # FIXME: Bound tail
    hypgeom_2f1_helper(y) =
        2sum(1:10) do k
            rising(1 + α, 2k) * rising(α - p, 2k) / rising(1 + α - p, 2k) / factorial(2k) *
            y^2k
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

        b = min(π / x, Arb(1e10))

        # Integration from a to b
        # IMPROVE: Enclosure for wide x
        J2 = real(Arblib.integrate(integrand, a, b))

        if b < π / x
            J3 =
                -log(b) / (α - p) * (
                    inv(b)^(α - p) * hypgeom_2f1_helper(inv(b)) -
                    (x / π)^(α - p) * hypgeom_2f1_helper(x / π)
                )
        else
            J3 = zero(J1)
        end

        J = J1 + J2 + J3

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
