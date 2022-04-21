export iswide

"""
    mince(x::Arb, n::Integer)

Return a vector with `n` balls covering the ball `x`.
"""
function mince(x::Arb, n::Integer)
    balls = Vector{Arb}(undef, n)
    xₗ, xᵤ = Arblib.getinterval(Arb, x)
    dx = (xᵤ - xₗ) / n
    for i in eachindex(balls)
        yₗ = xₗ + (i - 1) * dx
        yᵤ = xₗ + i * dx
        balls[i] = Arb((yₗ, yᵤ))
    end

    return balls
end

"""
    iswide(x; cutoff = 10)

Return true if `x` is wide in the meaning that the effective relative
accuracy of `x` measured in bits is more than `cutoff` lower than it's
precision. For `x` not of type `Arb` or `Acb` this always return
`false`.
"""
iswide(x::Union{Arb,Acb}; cutoff = 10) = Arblib.rel_accuracy_bits(x) < precision(x) - cutoff

iswide(::Number; cutoff = 10) = false

function stieltjes(::Type{Arb}, n::Integer)
    res = zero(Acb)
    a = one(Acb)

    ccall(
        (:acb_dirichlet_stieltjes, Arblib.libarb),
        Cvoid,
        (Ref{Arblib.acb_struct}, Ref{fmpz_struct}, Ref{Arblib.acb_struct}, Clong),
        res,
        fmpz_struct(convert(Int, n)),
        a,
        precision(res),
    )

    return real(res)
end

stieltjes(type, n::Integer) = convert(type, stieltjes(Arb, n))

"""
    unique_integer(x::Arb)

If `x` contains a unique integer return `true, n` where `n` is the
integer. Otherwise return `false, 0`

"""
function unique_integer(x::Arb)
    res = fmpz_struct()
    unique = ccall(
        Arblib.@libarb(arb_get_unique_fmpz),
        Int,
        (Ref{fmpz_struct}, Ref{Arblib.arb_struct}),
        res,
        x,
    )

    return !iszero(unique), Int(res)
end

"""
    abspow(x, y)

Compute `abs(x)^y `in a way that works if `x` overlaps with zero.
"""
function abspow(x::Arb, y::Arb)
    iszero(y) && return one(x)

    if iszero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(zero(x))
        Arblib.ispositive(y) && return zero(x)
        return Arblib.unit_interval!(zero(x))
    end

    if Arblib.contains_zero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(zero(x))
        x_upp = Arblib.abs_ubound(Arb, x)
        return Arb((zero(x), x_upp^y))
    end

    res = abs(x)
    return Arblib.pow!(res, res, y)
end

function abspow(x::ArbSeries, y::Arb)
    if Arblib.contains_zero(Arblib.ref(x, 0))
        # All non-constant terms are indeterminate, the constant term
        # is given by abs(x[0])^y
        res = ArbSeries(abspow(x[0], y), degree = Arblib.degree(x))
        for i = 1:Arblib.degree(x)
            Arblib.indeterminate!(Arblib.ref(res, i))
        end
        # Since we manually set the coefficients of the polynomial we
        # need to also manually set the degree.
        res.arb_poly.length = length(x)
        return res
    end

    res = abs(x)
    return Arblib.pow_arb_series!(res, res, y, length(res))
end

function abspow(x::Arb, y::ArbSeries)
    # This function is only partially implemented. In the case when x
    # overlaps zero it is currently based on differentiation and
    # isolation of extrema done by hand. If we want to support much
    # higher degrees we would need to do this algorithmically. In the
    # current version it is also not optimized at all, so it could be
    # much faster.

    # Helper function for creating indeterminate ArbSeries
    function nan(p::ArbSeries)
        indet = Arblib.indeterminate!(zero(Arb))
        res = zero(p)
        for i = 0:Arblib.degree(res)
            res[i] = indet
        end
        return res
    end

    iszero(y) && return one(y)

    if iszero(x)
        # If y[0] = 0 the constant term is zero but not the others. We
        # therefore need y[0] > 0 to get exactly zero.
        Arblib.contains_nonpositive(Arblib.ref(y, 0)) && return nan(y)

        return zero(y)
    end

    if Arblib.contains_zero(x)
        Arblib.contains_nonpositive(Arblib.ref(y, 0)) && return nan(y)

        # Differentiate with respect to the parameter of y manually
        # and enclose the terms

        deg = Arblib.degree(y)

        deg <= 3 || error("supports degree at most 3")

        res = zero(y)

        res[0] = abspow(x, y[0])
        if deg >= 1
            # res[1] = y[1] * log(x) * abspow(x, y[0])

            # Compute enclosure of log(x) * abspow(x, y[0])
            # Evaluate at x = 0, x = ubound(x) and possibly extrema
            f1(x) = log(x) * abspow(x, y[0])

            term = union(zero(x), f1(ubound(Arb, x)))

            extrema = exp(-1 / y[0])
            if Arblib.overlaps(x, extrema)
                term = union(term, f1(extrema))
            end

            res[1] = y[1] * term
        end
        if deg >= 2
            #res[2] = (2y[2] * log(x) + (y[1] * log(x))^2) / 2 * abspow(x, y[0])

            # Compute enclosure of (2y[2] * log(x) + (y[1] * log(x))^2) * abspow(x, y[0])
            f2(x) = (2y[2] * log(x) + (y[1] * log(x))^2) * abspow(x, y[0])

            term = union(zero(x), f2(ubound(Arb, x)))

            if iszero(y[1])
                # The extrema simplifies
                extrema = exp(-1 / y[0])
                if Arblib.overlaps(x, extrema)
                    term = union(term, f2(extrema))
                end
            else
                Δ = 4(y[0] * y[2] + y[1]^2)^2 - 8y[1]^2 * y[2]
                if Arblib.contains_nonnegative(Δ) # Otherwise there are no real roots
                    Arblib.sqrtpos(Δ)
                    extrema1 = exp((-2(y[0] * y[2] + y[1]^2) + Arblib.sqrtpos(Δ)) / 2y[1]^2)
                    extrema2 = exp((-2(y[0] * y[2] + y[1]^2) - Arblib.sqrtpos(Δ)) / 2y[1]^2)
                    if Arblib.overlaps(x, extrema1)
                        term = union(term, f2(extrema1))
                    end
                    if Arblib.overlaps(x, extrema2)
                        term = union(term, f2(extrema2))
                    end
                end

                res[2] = term / 2
            end
        end
        if deg >= 3
            # I DONT WANT TO DO THIS!!!
            @warn "not properly implemented for degree 3"
            res[3] =
                (6y[3] * log(x) + 6y[1] * y[2] * log(x)^2 + (y[1] * log(x))^3) / 6 * res[0]
        end

        return res
    end

    return abs(x)^y
end

abspow(x, y) = abs(x)^y

"""
    taylor_with_error(f, a::Arb, X::Arb, N::Integer)

Compute the Taylor expansion `P` of `f` of degree `N - 1` (i.e. `N`
terms) at the point `a` and bound the error term on `X`.

Returns a tuple `(P, E)` containing the Taylor expansion `P` and the
error term `E`. It satisfies that `f(x) ∈ P(x - a) + E*(x - a)^N` for
all `x ∈ X`.

We required that `a ∈ x`.
"""
function taylor_with_error(f, a::Arb, X::Arb, N::Integer)
    @assert contains(X, a)

    P = f(ArbSeries((a, 1), degree = N))

    E = f(ArbSeries((X, 1), degree = N + 1))[N]

    return P, E
end

"""
    abs(x::ArbSeries)

Compute the absolute value of `x`.

If `x[0]` contains zero then all non-constant terms are set to `NaN`,
otherwise either `-x` or `x` is returned depending on the sign of x.
"""
function Base.abs(x::ArbSeries)
    if Arblib.contains_zero(Arblib.ref(x, 0))
        res = zero(x)
        Arblib.abs!(Arblib.ref(res, 0), Arblib.ref(x, 0))
        for i = 1:Arblib.degree(x)
            Arblib.indeterminate!(Arblib.ref(res, i))
        end
        # Since we manually set the coefficients of the polynomial we
        # need to also manually set the degree.
        res.arb_poly.length = Arblib.degree(x) + 1
        return res
    elseif Arblib.isnegative(Arblib.ref(x, 0))
        return -x
    else
        return copy(x)
    end
end

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
            ArbExtras.enclosure_series(g, interval, degree = enclosure_degree) /
            factorial(degree)
    end

    return res
end

function taylor_with_remainder(f, x0::Acb, interval::Acb; degree::Integer)
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
    fx_div_x(f, x::Union{Arb,Acb}[, order::Integer]; extra_degree::Integer = 0, force = false)
    fx_div_x(f, x::ArbSeries[, order::Integer]; extra_degree::Integer = 0, force = false)

Compute an enclosure of `f(x) / x` for a function `f` with a zero at
the origin.

If `order` is given it computes `f(x) / x^order`, assuming `f` has a
zero of the given order at zero.

Setting `extra_degree` to a value higher than `0` makes it use a
higher order expansion to enclose the value, this can give tighter
bounds in many cases.

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

    expansion = taylor_with_remainder(f, zero(x), x, degree = order + extra_degree)

    isfinite(expansion) || return Arblib.indeterminate!(zero(x))

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
        degree = Arblib.degree(x) + order + extra_degree,
    )

    if !isfinite(expansion)
        # Return an indeterminate result
        res = zero(x)
        for i = 0:Arblib.degree(x)
            Arblib.indeterminate!(Arblib.ref(res, i))
        end
        # Since we manually set the coefficients of the polynomial we
        # need to also manually set the degree.
        res.arb_poly.length = length(x)
        return res
    end

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

    # Compose the Taylor series for the result with that of the input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end
