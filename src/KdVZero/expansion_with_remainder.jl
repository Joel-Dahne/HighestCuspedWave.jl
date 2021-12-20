"""
    truncate_with_remainder(p::ArbPoly, interval::Arb; degree::Integer)

Return a `ArbSeries` corresponding to a truncated version of `p` with
the last term being a remainder term which ensures that the truncated
version gives an enclosure of `p(x)` for all `x ∈ interval`.
"""
function truncate_with_remainder(p::ArbPoly, interval::Arb; degree::Integer)
    # Set the result to a truncated version of p
    res = ArbSeries(p; degree)

    # Set q to a p divided by x^degree and throwing away lower order
    # terms.
    q = Arblib.shift_right!(zero(p), p, degree)

    # Evaluate q on the given interval and set this as the remainder term
    res[degree] = q(interval)

    return res
end

"""
    compose_with_remainder(p::ArbSeries, q::ArbSeries, interval::Arb; degree)

Compute the composition `p(q(x))` and collapse all terms of degree
`degree` or higher into one remainder term of degree `degree` that is
valid for all `x ∈ interval`.

This method hence returns a polynomial `res` such that `p(q(x)) ∈
res(x)` for all `x ∈ interval`.
"""
function compose_with_remainder(
    p::ArbSeries,
    q::ArbSeries,
    interval::Arb;
    degree::Integer = Arblib._degree(p, q),
)
    # Compute a non-truncated composition as a ArbPoly
    poly = Arblib.compose!(ArbPoly(), p, q)

    # Truncate to the specified degree
    return truncate_with_remainder(poly, interval; degree)
end

"""
    compose_with_remainder(f, q::ArbSeries, interval::Arb; degree)

Compute an expansion of `f(q(x))` of the given degree. The last term
is a remainder term which ensures that the expansion gives an
enclosure of `f(q(x))` for all `x ∈ interval`.

It first computes an expansion `p` of `f` at the point `q[0]` with a
remainder term that makes it valid on `q(interval)`. That is, this
expansion satisfy
```
f(x) ∈ p(x - q[0])
```
for all `x ∈ q(interval)`.

It then computes the composition `p(q - q[0])` using
[`compose_with_remainder`](@ref) applied to `p` and `q - q[0]`.
"""
function compose_with_remainder(
    f,
    q::ArbSeries,
    interval::Arb;
    degree::Integer = Arblib.degree(q),
)
    # Compute expansion of f at q[0]
    p = f(ArbSeries((q[0], 1); degree))

    # Compute remainder term
    p[degree] = f(ArbSeries((q(interval), 1); degree))[degree]

    # q - q[0]
    qmq0 = ArbSeries(q)
    qmq0[0] = 0

    return compose_with_remainder(p, qmq0, interval; degree)
end

"""
    mul_with_remainder(p::ArbSeries, q::ArbSeries, interval::Arb; degree::Integer)

Compute an expansion of `p * q` with the last term being a remainder
term ensuring that it gives an enclosure of `p(x) * q(x)` for all `x ∈
interval`.
"""
function mul_with_remainder(
    p::ArbSeries,
    q::ArbSeries,
    interval::Arb;
    degree::Integer = Arblib._degree(p, q),
)
    # Compute non-truncated polynomial
    poly = Arblib.mul!(ArbPoly(), p, q)

    # Truncate to the specified degree
    return truncate_with_remainder(poly, interval; degree)
end

"""
    div_with_remainder(p::ArbSeries, q::ArbSeries, interval::Arb; degree::Integer)

Compute an expansion of `p / q` with the last term being a remainder
term ensuring that it gives an enclosure of `p(x) * q(x)` for all `x ∈
interval`.

This is done by first computing `inv(q)` with
[`compose_with_remainder`](@ref) and the `p * inv(q)` with
[`mul_with_remainder`](@ref).
"""
div_with_remainder(
    p::ArbSeries,
    q::ArbSeries,
    interval::Arb;
    degree::Integer = Arblib._degree(p, q),
) = mul_with_remainder(
    p,
    compose_with_remainder(inv, q, interval; degree),
    interval;
    degree,
)

"""
    clausencmzeta_with_remainder(x::Arb, s::ArbSeries, interval::Arb; degree)

Compute an expansion of `clausencmzeta(x, s)` in the parameter `s`
with the last term being a remainder term ensuring that it gives an
enclosure of `clausenc(x, s)` for all `s ∈ interval`.

The only difference to calling `clausencmzeta(x, s)` directly is that
the last term works as a remainder term.

Currently this is equivalent to `compose_with_remainder(s ->
clausencmzeta(x, s), s)`. The idea is that this will eventually
implement certain optimizations for this special case.
"""
function clausencmzeta_with_remainder(
    x::Arb,
    s::ArbSeries,
    interval::Arb;
    degree = Arblib.degree(s),
)
    # Compute expansion at s[0]
    p = clausencmzeta(x, ArbSeries((s[0], 1); degree))

    # Compute remainder term
    remainder_term = clausencmzeta(x, s(interval), degree) / factorial(degree)
    p[degree] = remainder_term

    # q - q[0]
    sms0 = ArbSeries(s)
    sms0[0] = 0

    return compose_with_remainder(p, sms0, interval; degree)
end
