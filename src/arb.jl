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

"""
    stieltjes(T, n::Integer)

Compute the Stieltjes constant `γₙ` in type `T`.
"""
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
stieltjes(T, n::Integer) = convert(T, stieltjes(Arb, n))

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
    abs(x::ArbSeries)

Compute the absolute value of `x`.

If `x[0]` contains zero then all non-constant terms are set to `NaN`,
otherwise either `-x` or `x` is returned depending on the sign of
`x[0]`.
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
