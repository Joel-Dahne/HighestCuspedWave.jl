"""
    <<(p::Union{ArbSeries,AcbSeries}, n::Integer)

Return `p` divided by `x^n`, updating the degree accordingly

It throws an error if the lower order coefficients are not all exactly
equal to zero.

Note that the naming is different from Arb where division by `x^n` is
referred to as right shift whereas here we call it a left shift.
"""
function Base.:(<<)(p::T, n::Integer) where {T<:Union{ArbSeries,AcbSeries}}
    n >= 0 || throw(ArgumentError("n needs to be non-negative, got $n"))
    for i = 0:n-1
        iszero(Arblib.ref(p, i)) ||
            throw(ArgumentError("coefficient $i not equal to zero, got $(p[i])"))
    end
    return Arblib.shift_right!(T(degree = Arblib.degree(p) - n, prec = precision(p)), p, n)
end

"""
    >>(p::Union{ArbSeries,AcbSeries}, n::Integer)

Return `p` multiplied by `x^n`, updating the degree accordingly.

Note that the naming is different from Arb where multiplication by
`x^n` is referred to as left shift whereas here we call it a right
shift.
"""
function Base.:(>>)(p::T, n::Integer) where {T<:Union{ArbSeries,AcbSeries}}
    n >= 0 || throw(ArgumentError("n needs to be non-negative, got $n"))
    return Arblib.shift_left!(T(degree = Arblib.degree(p) + n, prec = precision(p)), p, n)
end

"""
    taylor_with_error(f, x0::Arb, interval::Arb, degree::Integer)

Compute the Taylor expansion of `f` at the point `x0` to the given
`degree` and bound the error term on `interval`.

Returns a tuple `(p, E)` containing the Taylor expansion `p` and the
error term `E`. It satisfies that `f(x) ∈ p(x - x0) + E*(x -
x0)^(degree + 1)` for all `x ∈ interval`.

We required that `x0 ∈ interval`.
"""
function taylor_with_error(f, x0::Arb, interval::Arb, degree::Integer)
    contains(interval, x0) || throw(
        ArgumentError(
            "expected x0 to be contained in interval, got x0 = $x0, interval = $interval",
        ),
    )

    p = f(ArbSeries((x0, 1); degree))

    E = f(ArbSeries((interval, 1), degree = degree + 1))[degree+1]

    return p, E
end

"""
    taylor_with_remainder(f, x0::T, interval::T; degree::Integer, enclosure_degree::Integer = 0) where {T<:Union{Arb,Acb}}

Compute the Taylor expansion of `f` at the point `x0` with the last
term being a remainder term of degree `degree` which ensures that the
truncated version gives an enclosure of `f(x)` for all `x ∈ interval`.

We require that `x0 ∈ interval`.

For wide values of `interval` it computes a tighter enclosure of the
remainder term using [`ArbExtras.enclosure_series`](@ref). The degree
used for this can be set with `enclosure_degree`. Setting it to a
negative number makes it compute it directly instead. This argument is
only supported for `Arb` and not for `Acb`.
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

    if x0 == interval
        # In this case we can compute the full expansion directly
        return f(ArbSeries((x0, 1); degree))
    end

    # Compute expansion without remainder term
    res = f(ArbSeries((x0, 1), degree = degree - 1))

    # Make room for remainder term
    res = ArbSeries(res; degree)

    # Compute remainder term
    if enclosure_degree < 0 || !iswide(interval)
        res[degree] = f(ArbSeries((interval, 1); degree))[degree]
    else
        # We compute a tighter enclosure with the help of ArbExtras.enclosure_series
        res[degree] =
            ArbExtras.enclosure_series(
                ArbExtras.derivative_function(f, degree),
                interval,
                degree = enclosure_degree,
            ) / factorial(degree)
    end

    return res
end

function taylor_with_remainder(
    f,
    x0::Acb,
    interval::Acb;
    degree::Integer,
    enclosure_degree::Integer = 0,
)
    contains(interval, x0) || throw(
        ArgumentError(
            "expected x0 to be contained in interval, got x0 = $x0, interval = $interval",
        ),
    )

    if x0 == interval
        # In this case we can compute the full expansion directly
        return f(ArbSeries((x0, 1); degree))
    end

    # Compute expansion without remainder term
    res = f(AcbSeries((x0, 1), degree = degree - 1))

    # Make room for remainder term
    res = AcbSeries(res; degree)

    # Compute remainder term
    res[degree] = f(AcbSeries((interval, 1); degree))[degree]

    return res
end

"""
    fx_div_x(f, x::Union{Arb,Acb}[, order::Integer]; extra_degree::Integer = 0, enclosure_degree::Integer = 0, force = false)
    fx_div_x(f, x::ArbSeries[, order::Integer]; extra_degree::Integer = 0, enclosure_degree::Integer = 0, force = false)

Compute an enclosure of `f(x) / x` for a function `f` with a zero at
the origin.

If `order` is given it computes `f(x) / x^order`, assuming `f` has a
zero of the given order at zero.

Setting `extra_degree` to a value higher than `0` makes it use a
higher order expansion to enclose the value, this can give tighter
bounds in many cases. The argument `enclosure_degree` is passed to
[`taylor_with_remainder`](@ref).

If `force = false` it requires that the enclosure of `f` at zero is
exactly zero. If `f` is known to be exactly zero at zero but the
enclosure might be wider it can be forced to be zero by setting `force
= true`
"""
function fx_div_x(
    f,
    x::Union{Arb,Acb},
    order::Integer = 1;
    extra_degree::Integer = 0,
    enclosure_degree::Integer = 0,
    force = false,
)
    @assert Arblib.contains_zero(x)
    @assert order >= 1
    @assert extra_degree >= 0

    if iszero(x)
        # If x is exactly zero there is no need to compute extra terms
        # in the expansion
        extra_degree = 0
    end

    expansion = taylor_with_remainder(
        f,
        zero(x),
        x,
        degree = order + extra_degree;
        enclosure_degree,
    )

    isfinite(expansion) || return indeterminate(x)

    if force
        for i = 0:order-1
            @assert Arblib.contains_zero(Arblib.ref(expansion, i))
            expansion[i] = 0
        end
    end

    return (expansion << order)(x)
end

function fx_div_x(
    f,
    x::ArbSeries,
    order::Integer = 1;
    extra_degree::Integer = 0,
    enclosure_degree::Integer = 0,
    force = false,
)
    x0 = x[0]
    @assert Arblib.contains_zero(x0)
    @assert order >= 1
    @assert extra_degree >= 0

    expansion = taylor_with_remainder(
        f,
        zero(x0),
        x0,
        degree = Arblib.degree(x) + order + extra_degree;
        enclosure_degree,
    )

    isfinite(expansion) || return indeterminate(x)

    if force
        for i = 0:order-1
            @assert Arblib.contains_zero(Arblib.ref(expansion, i))
            expansion[i] = 0
        end
    end

    expansion_div_x = expansion << order

    # Set the result to the Taylor series of f(x) / x on interval
    res = ArbSeries(degree = Arblib.degree(x))
    for i = 0:Arblib.degree(res)
        res[i] = expansion_div_x(x0) / factorial(i)
        Arblib.derivative!(expansion_div_x, expansion_div_x)
    end

    return ArbExtras.compose_zero(res, x)
end
