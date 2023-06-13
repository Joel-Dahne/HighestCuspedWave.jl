"""
    T0(u0::FractionalKdVAnsatz{Arb}, ::Ball, δ1, δ2, ϵ, skip_div_u0 = false)

Return a function such that `T0(u0)(x)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * ∫ (clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * u0.w(y) dy
```
where the integration is taken from `0` to `π`.

If the weight is given by `u0.w(x) = abs(x)` the compute the value
using [`T0_p_one`](@ref). Otherwise the integral is split into two
parts and computed using [`T01`](@ref) and [`T02`](@ref).

Th argument `δ1` is given to [`T01`](@ref) and `δ2` and `ϵ` are given
to [`T02`](@ref).

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.
"""
function T0(
    u0::FractionalKdVAnsatz{Arb},
    evaltype::Ball;
    δ1::Arb = Arb(1e-3),
    δ2::Arb = Arb(1e-2),
    ϵ::Arb = Arb(0.1),
    skip_div_u0 = false,
)
    # Use specialised implementation in the case the weight is abs(x)
    weightisx(u0) && return T0_p_one(u0, evaltype; skip_div_u0)

    f = T01(u0, evaltype, skip_div_u0 = true; δ1)
    g = T02(u0, evaltype, skip_div_u0 = true; δ2, ϵ)

    return x -> begin
        part1 = f(x) ## Integral on [0, x]

        # Short circuit on a non-finite result
        isfinite(part1) || return part1

        part2 = g(x) ## Integral on [x, π]

        isfinite(part2) || return part2

        if skip_div_u0
            return part1 + part2
        else
            return (part1 + part2) / u0(x)
        end
    end
end

"""
    T0(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic, M, ϵ, return_enclosure = false)

Return a function such that `T0(u0, Asymptotic())(x)` computes the
integral `T0(u0)(x)` using an evaluating strategy that works
asymptotically as `x` goes to zero.

In general it only computes an **upper bound** of the integral. If
`return_enclosure` is true then it returns an enclosure instead of an
upper bound by returning the interval between zero and the computed
upper bound.

If the weight is given by `u0.w(x) = abs(x)` the compute the value
using [`T0_p_one`](@ref). If `u0.use_bhkdv` is true then it uses
[`T0_bhkdv`](@ref). Otherwise the integral is split into two parts and
computed using [`T01`](@ref) and [`T02`](@ref).

The argument `M` gives the number of terms to use in the asymptotic
expansion of `u0` and `ϵ` determines the interval ``[0, ϵ]`` on which
the expansion is valid.
"""
function T0(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
    return_enclosure::Bool = false,
)
    # Use specialised implementation in the case the weight is x or if
    # use_bhkdv is true
    weightisx(u0) && return T0_p_one(u0, Asymptotic(); M, ϵ)
    u0.use_bhkdv && return T0_bhkdv(u0, Asymptotic(); M, ϵ)

    f = T01(u0, Asymptotic(); M, ϵ, return_enclosure)
    g = T02(u0, Asymptotic(); M, ϵ, return_enclosure)
    return x -> f(x) + g(x)
end

"""
    T0_p_one(u0, Ball())

Return a function such that `T0_p_one(u0)(x)` computes `T0(u0)(x)` for
`u0` with the weight `x`.

This is based on [`lemma_U0_primitive_weight_x`](@ref) but we give the
proof below as well. To make it easier to compute good enclosures the
terms are split slightly different from in the paper.

In this case the integral is given by
```
inv(π * x * u0(x)) * ∫abs(clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * y dy
```
with the integral taken from `0` to `π`. The switch of coordinates to
`t = y / x` gives us
```
x / (π * u0(x)) * ∫abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt
```
from `0` to `π / x`. Call the integral, without the factors in front,
`I`.

If we ignore the absolute value the integral can be computed
explicitly using that
```
∫(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt =
    (clausenc(x * (1 - t), 2 - α) / x^2 - t * clausens(x * (1 - t), 1 - α) / x) +
    (clausenc(x * (1 + t), 2 - α) / x^2 + t * clausens(x * (1 + t), 1 - α) / x) -
    2(clausenc(x * t, 2 - α) / x^2 + t * clausens(x * t, 1 - α) / x)
```
Call this function `primitive(t)`. The idea is to isolate the parts of
the interval where
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
is zero. We then split the interval into several parts where the sign
is constant, so we can remove the absolute value (taking the
appropriate sign) and integrate explicitly using `primitive(t)`.

To begin with we split it as `I = I1 + I2` where `I1` is the
integration on `[0, 1]` and `I2` on `[1, π / x]`.

On the interval `[0, 1]` the expression inside the absolute value has,
by [`lemma_I_hat_root`](@ref), a unique root which we can isolate with
[`_integrand_compute_root`](@ref), it is negative to the left of the
root and positive to the right. If we let `r` be the root then the
integral is given by
```
I1 = -I11 + I12
```
where `I11 = primitive(r) - primitive(0)` is the integral on `[0, r]`
and `I12 = primitive(1) - primitive(r)` is the integral on `[r, 1]`.

On the interval `[1, π / x]` the expression inside the absolute value
is positive, see [lemma_I_positive`](@ref), and we can just remove it.
This gives us the integral
```
I2 = primitive(π / x) - primitive(1)
```

Putting all of this together and simplifying we get
```
I1 + I2 = primitive(0) - 2primitive(r) + primitive(x / π)
```
We also get that
```
primitive(0) = 2clausencmzeta(x, 2 - α) / x^2
```
and
```
primitive(π / x) = 2(clausenc(x + π, 2 - α) - clausenc(Arb(π), 2 - α)) / x^2
```
where `clausenc(Arb(π), 2 - α)` can also be given as the negated
alternating zeta function, `-eta(2 - α)`.

We can notice that all terms in the result contains a division by `x`
and that we in the end multiply with `x` If we let
```
primitive_mul_x(t) = (clausenc(x * (1 - t), 2 - α) / x - t * clausens(x * (1 - t), 1 - α)) +
    (clausenc(x * (1 + t), 2 - α) / x + t * clausens(x * (1 + t), 1 - α)) -
    2(clausenc(x * t, 2 - α) / x + t * clausens(x * t, 1 - α))
```
we get
```
x * (I1 + I2) = primitive_mul_x(0) - 2primitive_mul_x(r) + primitive_mul_x(x / π)
```
"""
function T0_p_one(u0::FractionalKdVAnsatz, ::Ball = Ball(); skip_div_u0 = false)
    weightisx(u0) || error("only supports u0 with weight x")

    return x::Arb -> begin
        r = _integrand_compute_root(typeof(u0), x, u0.α)

        # primitive(0)
        primitive_mul_x_zero(x) = 2clausencmzeta(x, 2 - u0.α) / x

        # primitive(r)
        primitive_mul_x_r(x) =
            (
                clausenc(x * (1 - r), 2 - u0.α) + clausenc(x * (1 + r), 2 - u0.α) -
                2clausenc(x * r, 2 - u0.α)
            ) / x +
            r * (
                -clausens(x * (1 - r), 1 - u0.α) + clausens(x * (1 + r), 1 - u0.α) -
                2clausens(x * r, 1 - u0.α)
            )

        # primitive(π / x)
        # eta is only implemented for Acb in Arblib
        primitive_mul_x_pi_div_x(x) =
            2(clausenc(x + π, 2 - u0.α) + real(eta(Acb(2 - u0.α)))) / x

        # Compute a tighter enclosure by expanding in x
        I_mul_x = ArbExtras.enclosure_series(x) do x
            primitive_mul_x_zero(x) - 2primitive_mul_x_r(x) + primitive_mul_x_pi_div_x(x)
        end

        if skip_div_u0
            return I_mul_x / π
        else
            return I_mul_x / (π * u0(x))
        end
    end
end

"""
    T0_p_one(u0, Asymptotic())

Return a function such that `T0_p_one(u0, Asymptotic)(x)` computes
`T0(u0)(x)` for `u0` with the weight `x` using an evaluating strategy
that works asymptotically as `x` goes to zero.

# Implementation
It first splits the function as
```
T0(x) = inv(π) * inv(u0(x) / x^-α) * (U0(x) / x^(-α + 1))
```
where
```
U0(x) = x^2 * ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t
```
integrated from `0` to `π / x`. The factor `inv(u0(x) / x^-α)` is
computed using [`inv_u0_normalised`](@ref) so the remaining work is in
enclosing the `U0 / x^(-α + 1)` factor.

For enclosing `U0 / x^(-α + 1)` it uses the primitive function given
in the non-asymptotic version of this method. It expands all functions
at `x = 0` and explicitly cancels the division by `x^(-α + 1)`.
"""
function T0_p_one(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
)
    weightisx(u0) || error("only supports u0 with weight x")

    α = u0.α

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    r0 = _integrand_compute_root(typeof(u0), zero(α), α)

    # Asymptotic expansion of clausencmzeta(x, 2 - α) for x ∈ [0, ϵ * (1 + r0)]
    C1, e1, P1, E1 = clausenc_expansion(ϵ * (1 + r0), 2 - α, M, skip_constant = true)

    # Asymptotic expansion of clausens(x, 1 - α) for x ∈ [0, ϵ * (1 + r0)]
    C2, e2, P2, E2 = clausens_expansion(ϵ * (1 + r0), 1 - α, M)

    # Taylor expansion of clausenc(x + π, 2 - α) - clausenc(π, 2 - α) for x ∈ [0, ϵ]
    P3 = taylor_with_lagrange_remainder(
        x -> clausenc(x + π, 2 - α),
        Arb(0),
        Arb((0, ϵ)),
        degree = M,
        enclosure_degree = 2,
    )
    P3[0] = 0 # Zero constant term
    # Function is even around π so all odd terms are zero
    for i = 1:2:M
        P3[i] = 0
    end

    # Helper function for evaluating p(x * b) / x^exponent
    eval_poly(p, x, b, exponent) =
        sum(0:Arblib.degree(p)) do i
            p_i = p[i]
            if iszero(p_i)
                zero(p_i)
            else
                p_i * abspow(b, i) * abspow(x, (i + exponent))
            end
        end

    return x::Arb -> begin
        @assert x <= ϵ

        r = _integrand_compute_root(typeof(u0), x, u0.α)

        # clausencmzeta(x, 2 - α) / x^(-α + 1)
        term1(x) = C1 + eval_poly(P1, x, one(α), α - 1) + E1 * abspow(x, 2M + α - 1)

        # (clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x^(-α + 1)
        term2(x) = eval_poly(P3, x, one(α), α - 1)

        # (clausenc(x * (1 - r), 2 - α) + clausenc(x * (1 + r), 2 - α) - 2clausenc(x * r, 2 - α)) / x^(-α + 1)
        term3(x) = (
            C1 * (1 - r)^e1 +
            eval_poly(P1, x, 1 - r, α - 1) +
            E1 * abspow(x, 2M + α - 1) * (1 - r)^2M +
            C1 * (1 + r)^e1 +
            eval_poly(P1, x, 1 + r, α - 1) +
            E1 * abspow(x, 2M + α - 1) * (1 + r)^2M -
            2 * (
                C1 * r^e1 + eval_poly(P1, x, r, α - 1) + E1 * abspow(x, 2M + α - 1) * r^2M
            )
        )

        # r * (clausens(x * (1 - r), 1 - α) + clausens(x * (1 + r), 1 - α) - 2clausens(x * r, 1 - α)) / x^-α
        term4(x) =
            r * (
                -(
                    C2 * (1 - r)^e2 +
                    eval_poly(P2, x, 1 - r, α) +
                    E2 * abspow(x, 2M + 1 + α) * (1 - r)^2M
                ) +
                C2 * (1 + r)^e2 +
                eval_poly(P2, x, 1 + r, α) +
                E2 * abspow(x, 2M + 1 + α) * (1 + r)^2M -
                2 * (
                    C2 * r^e2 + eval_poly(P2, x, r, α) + E2 * abspow(x, 2M + 1 + α) * r^2M
                )
            )

        terms(x) = 2(term1(x) + term2(x) - term3(x) - term4(x))

        res = ArbExtras.enclosure_series(terms, x)

        return inv_u0(x) * res / π
    end
end
