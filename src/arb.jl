export iswide

# For conversion from arb to Arb
function Arblib.set!(res::Arb, x::arb)
    ccall(Arblib.@libarb("arb_set"), Cvoid, (Ref{Arblib.arb_struct}, Ref{Nemo.arb}), res, x)
    return res
end

# For conversion from acb to Acb
Arblib.Acb(x::acb; prec::Integer = precision(parent(x))) = Arblib.set!(Acb(prec = prec), x)
function Arblib.set!(res::Acb, x::acb)
    ccall(Arblib.@libarb("acb_set"), Cvoid, (Ref{Arblib.acb_struct}, Ref{Nemo.acb}), res, x)
    return res
end

# For conversion from Arb to arb
function (r::ArbField)(x::Arb)
    z = arb(fmpz(0), r.prec)
    ccall(Arblib.@libarb("arb_set"), Cvoid, (Ref{Nemo.arb}, Ref{Arblib.arb_struct}), z, x)
    z.parent = r
    return z
end

# For conversion from Acb to acb
function (r::AcbField)(x::Acb)
    z = acb(fmpz(0), r.prec)
    ccall(Arblib.@libarb("acb_set"), Cvoid, (Ref{Nemo.acb}, Ref{Arblib.acb_struct}), z, x)
    z.parent = r
    return z
end

# Arb / fmpz
function Base.:/(x::Arblib.ArbOrRef, y::fmpz)
    res = zero(x)
    ccall(
        Arblib.@libarb("arb_div_fmpz"),
        Cvoid,
        (Ref{Arblib.arb_struct}, Ref{Arblib.arb_struct}, Ref{fmpz}, Int),
        res,
        x,
        y,
        precision(res),
    )
    return res
end

"""
    mince(xₗ::arb, xᵤ::arb, n::Integer; split = false)
Return a vector with `n` balls covering the interval `[xₗ, xᵤ]`.

If `split` is true then returns the intervals as tuples and not balls.
"""
function mince(xₗ::arb, xᵤ::arb, n::Integer; split = false)
    intervals = Vector{ifelse(split, NTuple{2,arb}, arb)}(undef, n)
    dx = (xᵤ - xₗ) / n
    for i in eachindex(intervals)
        yₗ = xₗ + (i - 1) * dx
        yᵤ = xₗ + i * dx
        if split
            intervals[i] = (yₗ, yᵤ)
        else
            intervals[i] = setunion(yₗ, yᵤ)
        end
    end
    return intervals
end

mince(x::arb, n::Integer; split = false) = mince(getinterval(x)..., n; split)

"""
    iswide(x; cutoff = 10)
Return true if `x` is wide in the meaning that the effective relative
accuracy of `x` measured in bits is more than `cutoff` lower than it's
precision.
"""
function iswide(x::arb; cutoff = 10)
    return ArbTools.rel_accuracy_bits(x) < precision(parent(x)) - cutoff
end

function iswide(x::acb; cutoff = 10)
    # TODO: Arb implements this as well, slightly different, but it
    # has to be implemented in ArbTools first
    return min(ArbTools.rel_accuracy_bits(real(x)), ArbTools.rel_accuracy_bits(imag(x))) <
           precision(parent(x)) - cutoff
end

iswide(x::Union{Arb,Acb}; cutoff = 10) = Arblib.rel_accuracy_bits(x) < precision(x) - cutoff

iswide(::Number; cutoff = 10) = false

function zeta(s::arb; d::Integer = 0)
    PP = ArbPolyRing(parent(s), :x)
    s_poly = PP([s, one(s)])
    a = one(s)
    res = PP()

    ccall(
        (:arb_poly_zeta_series, Nemo.libarb),
        Cvoid,
        (Ref{arb_poly}, Ref{arb_poly}, Ref{arb}, Cint, Clong, Clong),
        res,
        s_poly,
        a,
        0,
        d + 1,
        parent(s).prec,
    )

    return coeff(res, d)
end

function zeta(s::Arb; d::Integer = 0)
    s_series = ArbSeries([s, one(s)], degree = d)
    return Arblib.zeta_series!(zero(s_series), s_series, one(s), 0, d + 1)[d]
end

function zeta(s::T; d::Integer = 0) where {T}
    RR = RealField(precision(BigFloat))
    return convert(float(T), zeta(RR(s), d = d))
end

function stieltjes(type, n::Integer)
    CC = ComplexField(precision(BigFloat))
    res = zero(CC)
    a = one(CC)

    ccall(
        (:acb_dirichlet_stieltjes, Nemo.libarb),
        Cvoid,
        (Ref{acb}, Ref{fmpz}, Ref{acb}, Clong),
        res,
        fmpz(n),
        a,
        CC.prec,
    )

    convert(type, real(res))
end
function stieltjes(::Type{Arb}, n::Integer)
    res = zero(Acb)
    a = one(Acb)

    ccall(
        (:acb_dirichlet_stieltjes, Arblib.libarb),
        Cvoid,
        (Ref{Arblib.acb_struct}, Ref{fmpz}, Ref{Arblib.acb_struct}, Clong),
        res,
        fmpz(n),
        a,
        precision(res),
    )

    return real(res)
end
function stieltjes(n::Integer)
    stieltjes(Float64, n)
end

"""
    contains_pi(x1, x2)

Checks the interval `[x1, x2]` if it contains points of the form `kπ`.
Returns two booleans, the first one is false if it doesn't contain a
point on the form `2kπ` and the second one if it doesn't contain one
on the form `(2k+ 1)π`. If they are true it means they might contain
such a point. This is used to determine where the extrema of `Ci([x1,
x2])` can occur.

TODO: This is (hopefully) correct but not optimal.
"""
function contains_pi(x1::Arb, x2::Arb)
    @assert !(x1 > x2)

    # x1 or x2 equal to zero are the only cases when the division by π
    # can be exact, in which case it has to be handled differently.
    iszero(x1) && return (true, !(x2 < π))
    iszero(x2) && return (true, !(x1 > -Arb(π)))

    # We have k1ₗπ ≤ xₗ < (k1ᵤ + 1)π
    k1 = Arblib.floor!(zero(x1), x1 / π)
    # TODO: Implement unique_integer for Arb
    unique1ₗ, k1ₗ = unique_integer(
        ArbField(precision(k1))(Arblib.ceil!(zero(k1), Arblib.lbound(Arb, k1))),
    )
    unique1ᵤ, k1ᵤ = unique_integer(
        ArbField(precision(k1))(Arblib.floor!(zero(k1), Arblib.ubound(Arb, k1))),
    )
    k1ₗ, k1ᵤ = Int(k1ₗ), Int(k1ᵤ)
    @assert unique1ₗ && unique1ᵤ && k1ₗ * Arb(π) ≤ x1 < (k1ᵤ + 1) * Arb(π)
    # We have k2ₗπ ≤ xₗ < (k2ᵤ + 1)π
    k2 = Arblib.floor!(zero(x2), x2 / π)
    (unique2ₗ, k2ₗ) = unique_integer(
        ArbField(precision(k2))(Arblib.ceil!(zero(k2), Arblib.lbound(Arb, k2))),
    )
    (unique2ᵤ, k2ᵤ) = unique_integer(
        ArbField(precision(k2))(Arblib.floor!(zero(k2), Arblib.ubound(Arb, k2))),
    )
    k2ₗ, k2ᵤ = Int(k2ₗ), Int(k2ᵤ)
    @assert unique2ₗ && unique2ᵤ && k2ₗ * Arb(π) ≤ x2 < (k2ᵤ + 1) * Arb(π)

    if k1ₗ == k2ᵤ
        # No kπ
        return (false, false)
    elseif k1ₗ == k2ᵤ - 1
        # Might contain exactly one such point (or zero)
        if iseven(k1ₗ)
            return (false, true)
        else
            return (true, false)
        end
    else
        # Might contain two such points
        return (true, true)
    end
end

contains_pi(x1::arb, x2::arb) = contains_pi(Arb(x1), Arb(x2))

"""
    beta_inc(a, b, z)
Compute the (not regularised) incomplete beta function B(a, b; z).

Note that this method is different than
[`SpecialFunctions.beta_inc`](@ref) both in that it returns the
non-regularised value and that it only returns one value.
"""
beta_inc(a::Acb, b::Acb, z::Acb) = Arblib.hypgeom_beta_lower!(zero(z), a, b, z, 0)
beta_inc(a::Arb, b::Arb, z::Arb) = Arblib.hypgeom_beta_lower!(zero(z), a, b, z, 0)

beta_inc(a::acb, b::acb, z::acb) = parent(a)(beta_inc(Acb(a), Acb(b), Acb(z)))
beta_inc(a::arb, b::arb, z::arb) = parent(a)(beta_inc(Arb(a), Arb(b), Arb(z)))

"""
    beta_inc_zeroone(a, b, z)
Compute the (not regularised) incomplete beta function B(a, b; z)
assuming that `0 <= z <= 1`, discarding any other numbers in the
interval.

PROVE: That it's monotonically increasing in `z`.
"""
function beta_inc_zeroone(a::Arb, b::Arb, z::Arb)
    if 0 <= z <= 1
        return beta_inc(a, b, z)
    elseif z >= 0
        return setunion(beta_inc(a, b, Arblib.lbound(Arb, z)), beta_inc(a, b, one(z)))
    elseif z <= 1
        return setunion(beta_inc(a, b, zero(z)), beta_inc(a, b, Arblib.ubound(Arb, z)))
    else
        return setunion(beta_inc(a, b, zero(z)), beta_inc(a, b, one(z)))
    end
end

beta_inc_zeroone(a::arb, b::arb, z::arb) =
    parent(a)(beta_inc_zeroone(Arb(a), Arb(b), Arb(z)))

"""
    powpos(x, y)
Compute |x|^y in a way that works if x contains negative numbers.
"""
function abspow(x::arb, y::arb)
    if iszero(y)
        return one(x)
    end
    if contains_zero(x)
        contains_negative(y) && return parent(x)(NaN)
        x_upp = ArbTools.abs_ubound(x)
        return ArbTools.setunion(zero(x), x_upp^y)
    end

    return abs(x)^y
end

function abspow(x::Arb, y::Arb)
    if iszero(y)
        return one(x)
    end
    if Arblib.contains_zero(x)
        Arblib.contains_negative(y) && return Arb(NaN, prec = precision(x))
        # TODO: Replace with abs_ubound when implemented in Arblib
        x_upp = Arb(Arblib.get_abs_ubound!(Arf(prec = precision(x)), x))
        return Arb((zero(x), x_upp^y))
    end

    return abs(x)^y
end

abspow(x, y) = abs(x)^y

"""
    hypgeom_2f1(a, b, c, z)
Compute the (not regularised) Gauss hypergeometric function
₂F₁(a,b,c,z).
"""
function hypgeom_2f1(a::arb, b::arb, c::arb, z::arb)
    res = parent(z)()
    ccall(
        ("arb_hypgeom_2f1", Nemo.libarb),
        Cvoid,
        (Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Cint, Clong),
        res,
        a,
        b,
        c,
        z,
        0,
        precision(parent(z)),
    )
    return res
end

"""
    taylor_with_error(f, a::arb, X::arb, N::Integer)
Compute the Taylor expansion `P` of `f` of degree `N - 1` (i.e. `N`
terms) at the point `a` and bound the error term on `X`.

Returns a tuple `(P, E)` containing the Taylor expansion `P` and the
error term `E`. It satisfies that `f(x) ∈ P(x - a) + E*(x - a)^N` for
all `x ∈ X`.

We required that `a ∈ x`.
"""
function taylor_with_error(f, a::arb, X::arb, N::Integer)
    @assert contains(X, a)
    PP = ArbPolyRing(parent(a), :x)

    P = f(arb_series(PP([a, one(a)]), N))

    E = f(arb_series(PP([X, one(a)]), N + 1))[N]

    return (P, E)
end

"""
    abs(x::arb_series, n::Integer = length(x))
Compute the absolute value of `x`.

If `x[0]` contains zero then all non-constant terms are set to `NaN`,
otherwise either `-x` or `x` is returned depending on the sign of x.
"""
function Base.abs(x::arb_series, n::Integer = length(x))
    if contains_zero(x[0])
        return arb_series(
            parent(x.poly)([abs(x[0]); fill(base_ring(parent(x.poly))(NaN), n - 1)]),
        )
    elseif x[0] < 0
        return -x
    else
        return x
    end
end

"""
    abs(x::ArbSeries)
Compute the absolute value of `x`.

If `x[0]` contains zero then all non-constant terms are set to `NaN`,
otherwise either `-x` or `x` is returned depending on the sign of x.
"""
function Base.abs(x::ArbSeries)
    if Arblib.contains_zero(x[0])
        return ArbSeries([abs(x[0]); fill(Arb(NaN), Arblib.degree(x) - 1)])
    elseif x[0] < 0
        return -x
    else
        return x
    end
end

"""
    expint(s::acb, z::acb)
Compute the generalised exponential integral Eₛ(z).
"""
function expint(s::acb, z::acb)
    res = parent(s)()
    ccall(
        ("acb_hypgeom_expint", Nemo.libarb),
        Cvoid,
        (Ref{acb}, Ref{acb}, Ref{acb}, Int),
        res,
        s,
        z,
        precision(parent(z)),
    )
    return res
end

"""
    cosint(a::arb, z::arb)
Compute the generalised cosine integral Ciₐ(z).

It tries to give better enclosures by evaluating it at higher
precision if required.
"""
function cosint(a::arb, z::arb)
    RR = parent(z)
    CC = ComplexField(precision(RR))
    res = CC()
    ccall(
        ("acb_hypgeom_gamma_upper", Nemo.libarb),
        Cvoid,
        (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
        res,
        CC(a),
        CC(zero(z), z),
        0,
        precision(RR),
    )

    if iswide(res)
        res2 = CC()
        ccall(
            ("acb_hypgeom_gamma_upper", Nemo.libarb),
            Cvoid,
            (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
            res2,
            CC(a),
            CC(zero(z), z),
            0,
            4precision(RR),
        )
        res = CC(
            setintersection(real(res), real(res2)),
            setintersection(imag(res), imag(res2)),
        )
        if iswide(res, cutoff = 20)
            @warn "Wide enclosure when computing cosint res = $res"
        end
    end

    return real(exp(CC(zero(a), -RR(π) * a / 2)) * res)
end

"""
    cosintpi(a::arb, z)
Compute the generalised cosine integral Ciₐ(πz).

It tries to give better enclosures by evaluating it at higher
precision if required.
"""
function cosintpi(a::arb, z)
    RR = parent(a)
    CC = ComplexField(precision(RR))
    res = CC()
    ccall(
        ("acb_hypgeom_gamma_upper", Nemo.libarb),
        Cvoid,
        (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
        res,
        CC(a),
        CC(zero(z), RR(π) * z),
        0,
        precision(RR),
    )

    if iswide(res)
        RR2 = RealField(8precision(RR))
        CC2 = ComplexField(precision(RR2))
        res2 = CC2()
        ccall(
            ("acb_hypgeom_gamma_upper", Nemo.libarb),
            Cvoid,
            (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
            res2,
            CC2(a),
            CC2(zero(a), RR2(π) * RR2(z)),
            0,
            precision(RR2),
        )
        res = CC(
            setintersection(real(res), real(res2)),
            setintersection(imag(res), imag(res2)),
        )
        if iswide(res, cutoff = 20)
            @warn "Wide enclosure when computing cosintpi res = $res"
        end
    end

    return real(exp(CC(zero(a), -RR(π) * a / 2)) * res)
end

"""
    cosint_javi(x::arb)
Method that Javi used to compute cosint for some parameters that Arb
didn't handle well. To me it doesn't seem to perform much better
though.
"""
function cosint_javi(x::arb)
    # Good for some parameters
    res1 = -cosint(zero(x), x)

    s, c = sincos(x)
    lb = (s - (c + 1) / x) / x
    ub = (s - (c - 1) / x) / x
    res2 = setunion(lb, ub)

    if isfinite(res1)
        return setintersection(res1, res2)
    else
        return res2
    end
end

"""
    sinint(a::arb, z::arb)
Compute the generalised sine integral Siₐ(z).

It tries to give better enclosures by evaluating it at higher
precision if required.
"""
function sinint(a::arb, z::arb)
    RR = parent(z)
    CC = ComplexField(precision(RR))
    res = CC()
    ccall(
        ("acb_hypgeom_gamma_upper", Nemo.libarb),
        Cvoid,
        (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
        res,
        CC(a),
        CC(zero(z), z),
        0,
        precision(RR),
    )

    if iswide(res)
        res2 = CC()
        ccall(
            ("acb_hypgeom_gamma_upper", Nemo.libarb),
            Cvoid,
            (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
            res2,
            CC(a),
            CC(zero(z), z),
            0,
            4precision(RR),
        )
        res = CC(
            setintersection(real(res), real(res2)),
            setintersection(imag(res), imag(res2)),
        )
        if iswide(res, cutoff = 20)
            @warn "Wide enclosure when computing sinint res = $res"
        end
    end

    return -imag(exp(CC(zero(a), -RR(π) * a / 2)) * res)
end

"""
    sinintpi(a::arb, z)
Compute the generalised sinine integral Siₐ(πz).

It tries to give better enclosures by evaluating it at higher
precision if required.
"""
function sinintpi(a::arb, z)
    RR = parent(a)
    CC = ComplexField(precision(RR))
    res = CC()
    ccall(
        ("acb_hypgeom_gamma_upper", Nemo.libarb),
        Cvoid,
        (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
        res,
        CC(a),
        CC(zero(z), RR(π) * z),
        0,
        precision(RR),
    )

    if iswide(res)
        RR2 = RealField(8precision(RR))
        CC2 = ComplexField(precision(RR2))
        res2 = CC2()
        ccall(
            ("acb_hypgeom_gamma_upper", Nemo.libarb),
            Cvoid,
            (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
            res2,
            CC2(a),
            CC2(zero(a), RR2(π) * RR2(z)),
            0,
            precision(RR2),
        )
        res = CC(
            setintersection(real(res), real(res2)),
            setintersection(imag(res), imag(res2)),
        )
        if iswide(res, cutoff = 20)
            @warn "Wide enclosure when computing sinintpi res = $res"
        end
    end

    return -imag(exp(CC(zero(a), -RR(π) * a / 2)) * res)
end
