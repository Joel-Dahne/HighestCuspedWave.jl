"""
    taylor_with_remainder(f, x0::Arb, interval::Arb; degree::Integer, enclosure_degree::Integer = 0)

Compute the Taylor expansion of `f` at the point `x0` with the last
term being a remainder term which ensures that the truncated version
gives an enclosure of `f(x)` for all `x ∈ interval`.

We require that `x0 ∈ interval`.

It computes a tighter enclosure of the remainder term using
[`ArbExtras.extrema_series`](@ref). The degree used for this can be
set with `enclosure_degree`. Setting it to a negative number makes it
compute it directly instead.
"""
function taylor_with_remainder(
    f,
    x0::Arb,
    interval::Arb;
    degree::Integer,
    enclosure_degree::Integer = 0,
)
    contains(interval, x0) || throw(
        ArgumentError(
            "expected x0 to be contained in interval, got x0 = $x0, interval = $interval",
        ),
    )

    # Compute expansion without remainder term
    res = f(ArbSeries((x0, 1), degree = degree - 1))

    # Make room for remainder term
    res = ArbSeries(res; degree)

    # Compute remainder term
    if enclosure_degree < 0
        res[degree] = f(ArbSeries((interval, 1); degree))
    else
        # We compute a tighter enclosure with the help of ArbExtras.extrema_series
        g(x::Arb) = f(ArbSeries((x, 1); degree))[degree] * factorial(degree)
        g(x::ArbSeries) =
            if iszero(Arblib.degree(x))
                ArbSeries(g(x[0]))
            else
                Arblib.derivative(
                    f(ArbSeries(x, degree = Arblib.degree(x) + degree)),
                    degree,
                )
            end
        res[degree] =
            Arb(
                ArbExtras.extrema_series(
                    g,
                    getinterval(interval)...,
                    degree = enclosure_degree,
                )[1:2],
            ) / factorial(degree)
    end

    return res
end

"""
    truncate_with_remainder(p::Union{ArbPoly,ArbSeries}, interval::Arb; degree::Integer)

Return a `ArbSeries` corresponding to a truncated version of `p` with
the last term being a remainder term which ensures that the truncated
version gives an enclosure of `p(x)` for all `x ∈ interval`.
"""
function truncate_with_remainder(
    p::Union{ArbPoly,ArbSeries},
    interval::Arb;
    degree::Integer,
)
    # Set the result to a truncated version of p
    res = ArbSeries(p; degree)

    Arblib.degree(p) <= degree && return res

    # Set q to a p divided by x^degree and throwing away lower order
    # terms.
    q = Arblib.shift_right!(ArbPoly(), p, degree)

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
    # Compute expansion of f at q[0] with degree - 1
    p = f(ArbSeries((q[0], 1), degree = degree - 1))

    # Increase the degree of p to make room for the remainder term
    p = ArbSeries(p; degree)

    # Compute remainder term
    # We compute a tighter enclosure with the help of ArbExtras.extrema_series
    g(x::Arb) = f(ArbSeries((x, 1); degree))[degree] * factorial(degree)
    g(x::ArbSeries) =
        if iszero(Arblib.degree(x))
            ArbSeries(g(x[0]))
        else
            Arblib.derivative(f(ArbSeries(x, degree = Arblib.degree(x) + degree)), degree)
        end
    p[degree] =
        Arb(ArbExtras.extrema_series(g, getinterval(q(interval))..., degree = 0)[1:2]) /
        factorial(degree)
    #p[degree] = f(ArbSeries((q(interval), 1); degree))[degree]

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
    clausenc_with_remainder(x::Arb, s::ArbSeries, interval::Arb; degree)

Compute an expansion of `clausenc(x, s)` in the parameter `s` with the
last term being a remainder term ensuring that it gives an enclosure
of `clausenc(x, s)` for all `s ∈ interval`.

The only difference to calling `clausenc(x, s)` directly is that the
last term works as a remainder term.

Currently this is equivalent to `compose_with_remainder(s ->
clausenc(x, s), s)`. The idea is that this will eventually implement
certain optimizations for this special case.
"""
function clausenc_with_remainder(
    x::Arb,
    s::ArbSeries,
    interval::Arb;
    degree = Arblib.degree(s),
)
    # Compute expansion at s[0] with degree - 1
    p = clausenc(x, ArbSeries((s[0], 1), degree = degree - 1))

    # Increase the degree of p to make room for the remainder term
    p = ArbSeries(p; degree)

    # Compute remainder term
    if s[0] == 0 || s[0] == 2
        @warn "non-rigorous remainder" maxlog = 1
        # FIXME: Properly implement this. For now it assumes
        # monotonicity for clausenc(x, s, degree) in s, which is not
        # true in general

        # We evaluate s like this to make the endpoint exact if they
        # are integers. Otherwise we typically get a value slightly
        # smaller or larger than the exact integer, which gives much
        # worse bounds.
        s_lower, s_upper =
            ArbExtras.extrema_polynomial(ArbPoly(s), getinterval(interval)...)

        remainder_term =
            union(clausenc(x, s_lower, degree), clausenc(x, s_upper, degree)) /
            factorial(degree)
    else
        remainder_term = clausenc(x, s(interval), degree) / factorial(degree)
    end

    p[degree] = remainder_term

    # s - s[0]
    sms0 = ArbSeries(s)
    sms0[0] = 0

    return compose_with_remainder(p, sms0, interval; degree)
end

"""
    clausens_with_remainder(x::Arb, s::ArbSeries, interval::Arb; degree)

Compute an expansion of `clausens(x, s)` in the parameter `s` with the
last term being a remainder term ensuring that it gives an enclosure
of `clausens(x, s)` for all `s ∈ interval`.

The only difference to calling `clausens(x, s)` directly is that the
last term works as a remainder term.

Currently this is equivalent to `compose_with_remainder(s ->
clausens(x, s), s)`. The idea is that this will eventually implement
certain optimizations for this special case.
"""
function clausens_with_remainder(
    x::Arb,
    s::ArbSeries,
    interval::Arb;
    degree = Arblib.degree(s),
)
    # Compute expansion at s[0] with degree - 1
    p = clausens(x, ArbSeries((s[0], 1), degree = degree - 1))

    # Increase the degree of p to make room for the remainder term
    p = ArbSeries(p; degree)

    # Compute remainder term
    if s[0] == 1
        @warn "non-rigorous remainder" maxlog = 1
        # FIXME: Properly implement this. For now it assumes
        # monotonicity for clausens(x, s, degree) in s, which is not
        # true in general

        # We evaluate s like this to make the endpoint exact if they
        # are integers. Otherwise we typically get a value slightly
        # smaller or larger than the exact integer, which gives much
        # worse bounds.
        s_lower, s_upper =
            ArbExtras.extrema_polynomial(ArbPoly(s), getinterval(interval)...)

        remainder_term =
            union(clausens(x, s_lower, degree), clausens(x, s_upper, degree)) /
            factorial(degree)
    else
        remainder_term = clausens(x, s(interval), degree) / factorial(degree)
    end

    p[degree] = remainder_term

    # s - s[0]
    sms0 = ArbSeries(s)
    sms0[0] = 0

    return compose_with_remainder(p, sms0, interval; degree)
end

"""
    clausencmzeta_with_remainder(x::Arb, s::ArbSeries, interval::Arb; degree)

Compute an expansion of `clausencmzeta(x, s)` in the parameter `s`
with the last term being a remainder term ensuring that it gives an
enclosure of `clausencmzeta(x, s)` for all `s ∈ interval`.

The only difference to calling `clausencmzeta(x, s)` directly is that
the last term works as a remainder term.

For wide values of `x` it computes each term in the expansion
separately, allowing it to use [`ArbExtras.extrema_series`](@ref) to
compute a tighter enclosure of the terms. While this could be done
also for the remainder term it currently isn't, it doesn't give much
improvement in the cases it could be relevant and comes with a big
performance cost.
"""
function clausencmzeta_with_remainder(
    x::Arb,
    s::ArbSeries,
    interval::Arb;
    degree = Arblib.degree(s),
)
    if iswide(x)
        p = ArbSeries(; degree)

        # For the constant term it already implements handling of
        # monotonicity and we don't have to use
        # ArbExtras.extrema_series
        p[0] = clausencmzeta(x, s[0])

        for β = 1:degree-1
            p[β] = Arb((ArbExtras.extrema_series(
                x -> clausencmzeta(x, s[0], β),
                getinterval(x)...,
                degree = 1,
            )[1:2]))
        end
    else
        # Compute expansion at s[0] with degree - 1
        p = clausencmzeta(x, ArbSeries((s[0], 1), degree = degree - 1))

        # Increase the degree of p to make room for the remainder term
        p = ArbSeries(p; degree)
    end

    # Compute remainder term
    if s[0] == 2
        @warn "non-rigorous remainder" maxlog = 1
        # FIXME: Properly implement this. For now it assumes
        # monotonicity for clausencmzeta(x, s, degree) in s, which is
        # not true in general

        # We evaluate s like this to make the endpoint exact if they
        # are integers. Otherwise we typically get a value slightly
        # smaller or larger than the exact integer, which gives much
        # worse bounds.
        s_lower, s_upper =
            ArbExtras.extrema_polynomial(ArbPoly(s), getinterval(interval)...)

        remainder_term =
            union(clausencmzeta(x, s_lower, degree), clausencmzeta(x, s_upper, degree)) /
            factorial(degree)
    else
        remainder_term = clausencmzeta(x, s(interval), degree) / factorial(degree)
    end

    p[degree] = remainder_term

    # s - s[0]
    sms0 = ArbSeries(s)
    sms0[0] = 0

    return compose_with_remainder(p, sms0, interval; degree)
end

"""
    abspow_with_remainder(x::Arb, y::ArbSeries, interval::Arb; degree)

Compute an expansion of `abspow(x, y)` in the parameter `y` with the
last term being a remainder term ensuring that it gives an enclosure
of `abspow(x, y)` for all `y ∈ interval`.

The only difference to calling `abspow(x, y)` directly is that the
last term works as a remainder term.

This is equivalent to `compose_with_remainder(s -> abspow(x, y), y,
interval)`, with one specialised optimizations. If `Arblib.is_x(-y)`
is true then instead of evaluating `y(interval)` it takes `-interval`
directly. These two are obviously equivalent, the reason for doing
this optimization is that in some cases we have an exponent which is
non-positive but whose left endpoint is zero, in that case
`y(interval)` is supposed to non-negative but naive evaluation doesn't
preserve this.
"""
function abspow_with_remainder(
    x::Arb,
    y::ArbSeries,
    interval::Arb;
    degree = Arblib.degree(y),
)
    # Compute expansion at y[0] with degree - 1
    p = abspow(x, ArbSeries((y[0], 1), degree = degree - 1))

    # Increase the degree of p to make room for the remainder term
    p = ArbSeries(p; degree)

    # Compute remainder term
    # We compute the interval in this way to avoid spurious negative
    # values of y in the special case Arblib.is_x(-y)
    yinterval = Arblib.is_x(-y) ? -interval : y(interval)

    p[degree] = abspow(x, (ArbSeries((yinterval, 1); degree)))[degree]

    # y - y[0]
    ymy0 = ArbSeries(y)
    ymy0[0] = 0

    return compose_with_remainder(p, ymy0, interval; degree)
end
