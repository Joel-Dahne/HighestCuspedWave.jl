"""
    T0(u0::FractionalKdVAnsatz{Arb}, ::Ball)

**IMPROVE:** Tune the values for `δ0, δ1, δ2` depending on `u0.α` and
`u0.p`.
"""
function T0(
    u0::FractionalKdVAnsatz{Arb},
    evaltype::Ball;
    δ0::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    δ1::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    δ2::Arf = Arf(1e-2),
    ϵ::Arb = 1 + u0.α,
    skip_div_u0 = false,
)
    # Use specialised implementation in the case the weight is x
    isone(u0.p) && return T0_p_one(u0, evaltype; skip_div_u0)

    f = T01(u0, evaltype, skip_div_u0 = true; δ0, δ1)
    g = T02(u0, evaltype, skip_div_u0 = true; δ2, ϵ)

    return x -> begin
        ## Integral on [0, x]
        part1 = f(x)

        # Short circuit on a non-finite result
        isfinite(part1) || return part1

        ## Integral on [x, π]
        part2 = g(x)

        isfinite(part2) || return part2

        if skip_div_u0
            return part1 + part2
        else
            return (part1 + part2) / u0(x)
        end
    end
end

function T0(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M::Integer = 5, ϵ::Arb = Arb(1))
    # Use specialised implementation in the case the weight is x
    isone(u0.p) && return T0_p_one(u0, Asymptotic(); M, ϵ)

    f = T01(u0, Asymptotic(); M, ϵ)
    g = T02(u0, Asymptotic(); M, ϵ)
    return x -> f(x) + g(x)
end

"""
    T0_p_one(u0, Ball())

Compute the integral ``T_0`` for `u0` with `u0.p == 1`.

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
is zero. For the other parts of the interval we can remove the
absolute value (taking the appropriate sign) and integrate explicitly
using `primitive(t)`.

On the interval `[1, π / x]` the expression inside the absolute value
is positive, see [lemma_integrand_2`](@ref), and we can just remove
it. This gives us the integral
```
I2 = primitive(π / x) - primitive(1)
```

On the interval `[0, 1]` the expression inside the absolute value has
a unique root which we can isolate with
[`_integrand_compute_root`](@ref), it is negative to the left of the
root and positive to the right. If we let `r` be the root then the
integral is given by
```
I1 = -I11 + I13
```
where `I11 = primitive(root_lower) - primitive(0)` is the integral on
`[0, r]` and `I12 = primitive(1) - primitive(r)` is the integral on
`[r, 1]`.

Putting all of this together and simplify we get
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
function T0_p_one(u0::FractionalKdVAnsatz, evaltype::Ball = Ball(); skip_div_u0 = false)
    @assert isone(u0.p)

    return x::Arb -> begin
        r = _integrand_compute_root(u0, x)

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
        I_mul_x(x) =
            primitive_mul_x_zero(x) - 2primitive_mul_x_r(x) + primitive_mul_x_pi_div_x(x)

        res = ArbExtras.enclosure_series(I_mul_x, x)

        if skip_div_u0
            return res / π
        else
            return res / (π * u0(x))
        end
    end
end

"""
    T0_p_one(u0, Asymptotic())

Compute the integral ``T_0`` for `u0` with `u0.p == 1` using an
asymptotic approach that works for small values of `x`.

# Implementation
It first splits the function as
```
T0(x) = inv(π) * inv(u0(x) / x^-α) * (u0(x) / x^(-α + 1))
```
where `α = u0.α`, `p = u0.p` and
```
U0(x) = x^2 * ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t^u0.p
```
integrated from `0` to `π / x`. The factor `inv(u0(x) / x^-α)` is
computed using [`inv_u0_normalised`](@ref) so the remaining work is in
enclosing the `U0 / x^(-α + 1)` factor.

For enclosing `U0 / x^(-α + 1)` it uses the primitive function given
in the non-asymptotic version of this method. It expands all functions
at `x = 0` and explicitly cancels the division by `x^(-α + 1)`.

**IMPROVE:** Use `ArbSeries` to compute a tighter enclosure? Similar
  to the non-asymptotic version. Depends on how important it is for
  the timings.
**IMPROVE:** Better explain what it does.
"""
function T0_p_one(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
)
    @assert isone(u0.p)

    α = u0.α

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    r0 = _integrand_compute_root(u0, Arb(0))

    # Asymptotic expansion of clausencmzeta(x, 2 - α) for x ∈ [0, ϵ * (1 + r0)]
    C1, e1, P1, E1 = clausenc_expansion(ϵ * (1 + r0), 2 - α, M, skip_constant = true)

    # Asymptotic expansion of clausens(x, 1 - α) for x ∈ [0, ϵ * (1 + r0)]
    C2, e2, P2, E2 = clausens_expansion(ϵ * (1 + r0), 1 - α, M)

    # Taylor expansion of clausenc(x + π, 2 - α) - clausenc(π, 2 - α) for x ∈ [0, ϵ]
    P3 = taylor_with_remainder(
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
        x <= ϵ || @show x ϵ
        @assert x <= ϵ

        r = _integrand_compute_root(u0, x)

        # clausencmzeta(x, 2 - α) / x^(-α + 1)
        term1 = C1 + eval_poly(P1, x, one(x), α - 1) + E1 * abspow(x, 2M + α - 1)

        # (clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x^(-α + 1)
        term2 = eval_poly(P3, x, one(x), α - 1)

        # (clausenc(x * (1 - r), 2 - α) + clausenc(x * (1 + r), 2 - α) - 2clausenc(x * r, 2 - α)) / x^(-α + 1)
        term3 = begin
            (
                C1 * (1 - r)^e1 +
                eval_poly(P1, x, 1 - r, α - 1) +
                E1 * abspow(x, 2M + α - 1) * (1 - r)^2M +
                C1 * (1 + r)^e1 +
                eval_poly(P1, x, 1 + r, α - 1) +
                E1 * abspow(x, 2M + α - 1) * (1 + r)^2M -
                2 * (
                    C1 * r^e1 +
                    eval_poly(P1, x, r, α - 1) +
                    E1 * abspow(x, 2M + α - 1) * r^2M
                )
            )
        end

        # r * (clausens(x * (1 - r), 1 - α) + clausens(x * (1 + r), 1 - α) - 2clausens(x * r, 1 - α)) / x^-α
        term4 = begin
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
        end

        res = 2(term1 + term2 - term3 - term4)

        return inv_u0(x) * res / π
    end
end
