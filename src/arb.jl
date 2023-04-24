export iswide

# There is no version of this in Arblib
function Arblib.indeterminate!(x::Union{ArbSeries,AcbSeries})
    for i = 0:Arblib.degree(x)
        Arblib.indeterminate!(Arblib.ref(x, i))
    end
    # Since we manually set the coefficients of the polynomial we
    # need to also manually set the degree.
    Arblib.cstruct(x).length = Arblib.degree(x) + 1
    return x
end

"""
    indeterminate(x)

Construct an indeterminate version of `x`.
"""
indeterminate(x::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}) = Arblib.indeterminate!(zero(x))
indeterminate(::Type{T}) where {T<:Union{Arb,Acb}} = Arblib.indeterminate!(zero(T))
indeterminate(x::Union{ArbSeries,AcbSeries}) = Arblib.indeterminate!(zero(x))

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
`false`. For `x` of type `ArbSeries` or `AcbSeries` it checks the
first coefficient.
"""
iswide(x::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}; cutoff = 10) =
    Arblib.rel_accuracy_bits(x) < precision(x) - cutoff
iswide(x::Union{ArbSeries,AcbSeries}; cutoff = 10) = iswide(Arblib.ref(x, 0))
iswide(::Number; cutoff = 10) = false

"""
    is_approx_integer(x::Arb; tol = 0.01)

Return true if `x` is close to an integer but doesn't contain an
integer.

It compares the radius of `x` with the distance from the midpoint to
its nearest integer, if the quotient is greater than `tol` it returns
true.

This is useful when you want to take extra care with calculations when
the arguments are close to integers, for example in the presence of
removable singularities.
"""
function is_approx_integer(x::Arblib.ArbOrRef; tol = 0.01)
    Arblib.contains_int(x) && return false
    # Do the computations in Float64 since we don't have to be
    # rigorous anyway
    xf64 = Float64(x)
    return Float64(Arblib.radref(x)) / abs(xf64 - round(xf64)) > tol
end

"""
    stieltjes(T, n::Integer)

Compute the Stieltjes constant `γₙ` in type `T`.
"""
function stieltjes(::Type{Arb}, n::Integer)
    res = zero(Acb)
    a = one(Acb)

    ccall(
        Arblib.@libarb(acb_dirichlet_stieltjes),
        Cvoid,
        (Ref{Arblib.acb_struct}, Ref{fmpz_struct}, Ref{Arblib.acb_struct}, Clong),
        res,
        fmpz_struct(convert(Int, n)),
        a,
        precision(res),
    )

    return real(res)
end
stieltjes(T, n::Integer) = convert(T, stieltjes(Arb, n))

"""
    unique_integer(x::Arb)

If `x` contains a unique integer return `true, n` where `n` is the
integer. Otherwise return `false, 0`
"""
function unique_integer(x::Arblib.ArbOrRef)
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

# See documentation for abs(x::ArbSeries) for behaviour
function Arblib.abs!(res::ArbSeries, x::ArbSeries)
    Arblib.degree(res) == Arblib.degree(x) ||
        throw(ArgumentError("res and x should have the same degree"))
    sgn = Arblib.sgn_nonzero(Arblib.ref(x, 0))
    if sgn == 0
        Arblib.indeterminate!(res)
        Arblib.abs!(Arblib.ref(res, 0), Arblib.ref(x, 0))
        return Arblib.normalise!(res)
    elseif sgn < 0
        return Arblib.neg!(res, x)
    else
        return Arblib.set!(res, x)
    end
end

"""
    abs(x::ArbSeries)

Compute the absolute value of `x`.

If `x[0]` contains zero then all non-constant terms are set to `NaN`,
otherwise either `-x` or `x` is returned depending on the sign of
`x[0]`.
"""
Base.abs(x::ArbSeries) = Arblib.abs!(zero(x), x)
