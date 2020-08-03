export iswide, Li, Ci, Si

"""
    iswide(x::arb)
Return true if x is wide in the meaning that the effective relative
accuracy of x measured in bits is more than 10 lower than it's parents
precision.
"""
function iswide(x::arb)
    # TODO: This might require some tuning, 10 might not be the
    # optimal number.
    return ArbTools.rel_accuracy_bits(x) < prec(parent(x)) - 10
end

function zeta(s::arb; d::Integer = 0)
    PP = ArbPolyRing(parent(s), :x)
    s_poly = PP([s, one(s)])
    a = one(s)
    res = PP()

    ccall((:arb_poly_zeta_series, Nemo.libarb), Cvoid,
          (Ref{arb_poly}, Ref{arb_poly}, Ref{arb}, Cint, Clong, Clong),
          res, s_poly, a, 0, d + 1, parent(s).prec)

    return coeff(res, d)
end

function zeta(s::T; d::Integer = 0) where {T}
    RealField(precision(BigFloat))
    return convert(float(T), zeta(RR(s), d = d))
end

function stieltjes(type, n::Integer)
    CC = ComplexField(precision(BigFloat))
    res = zero(CC)
    a = one(CC)

    ccall((:acb_dirichlet_stieltjes, Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{fmpz}, Ref{acb}, Clong),
          res, fmpz(n), a, CC.prec)

    convert(type, real(res))
end
function stieltjes(n::Integer)
    stieltjes(Float64, n)
end

"""
    Li(z, s)
Compute the polylogarithm Liₛ(z)
"""
function Li(z::acb, s::acb)
    res = parent(z)()
    ccall(("acb_polylog", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Clong), res, s, z, prec(parent(z)))
    return res
end

function Li(z::acb, s::Integer)
    res = parent(z)()
    ccall(("acb_polylog_si", Nemo.libarb), Cvoid,
          (Ref{acb}, Clong, Ref{acb}, Clong), res, s, z, prec(parent(z)))
    return res
end

# Not required
#function Li(z::arb, s::Integer)
#    CC = ComplexField(prec(parent(z)))
#    return real(Li(CC(z), s))
#end
#
#function Li(z::arb, s)
#    CC = ComplexField(prec(parent(z)))
#    return real(Li(CC(z), CC(s)))
#end
#
#function Li(z::T, s) where {T <: Real}
#    CC = ComplexField(precision(BigFloat))
#    res = real(Li(CC(z), CC(s)))
#    return convert(float(T), res)
#end

"""
    Ci(x, s)
Compute the Clausian function Ciₛ(x).

If x is a wide (real) ball (as determined by iswide(x)) it computes a
tighter enclosure by using that Ci 2π periodic, monotonic for x ∈ [0,
π] and even, so that it's enough to evaluate on the endpoints.
"""
function Ci(x::acb, s)
    im = x.parent(0, 1)
    return (Li(exp(im*x), s) + Li(exp(-im*x), s))/2
end

function Ci(x::arb, s::arb)
    if iswide(x)
        x_lower = max(ArbTools.abs_lbound(x), parent(x)(0))
        x_upper = min(ArbTools.abs_ubound(x), parent(x)(π))
        return setunion(Ci(x_upper, s), Ci(x_lower, s))
    end
    CC = ComplexField(prec(parent(x)))
    return real(Li(exp(CC(zero(x), x)), CC(s)))
end

function Ci(x::arb, s::Integer)
    if iswide(x)
        x_upper = min(ArbTools.abs_ubound(x), parent(x)(π))
        x_lower = max(ArbTools.abs_lbound(x), parent(x)(0))
        return setunion(Ci(x_lower, s), Ci(x_upper, s))
    end
    CC = ComplexField(prec(parent(x)))
    return real(Li(exp(CC(zero(x), x)), s))
end

function Ci(x::T, s) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = real(Li(CC(real(z), imag(z)), CC(s)))
    return convert(float(T), res)
end

function Ci(x::T, s::Integer) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = real(Li(CC(real(z), imag(z)), s))
    return convert(float(T), res)
end

"""
    Si(x, s)
Compute the Clausian function Siₛ(x). Assumes that x ∈ (-π, π).
"""
function Si(x::acb, s::acb)
    im = x.parent(0, 1)
    return (Li(exp(im*x), s) - Li(exp(-im*x), s))/2
end

function Si(x::arb, s::arb)
    if iswide(x)
        # TODO: This seems to perform worse
        return ball(Si(midpoint(x), s), (x - midpoint(x))*Ci(x, s - 1))
    end
    CC = ComplexField(prec(parent(x)))
    return imag(Li(exp(CC(zero(x), x)), CC(s)))
end

function Si(x::arb, s::Integer)
    if iswide(x)
        # TODO: This seems to perform worse
        return ball(Si(midpoint(x), s), (x - midpoint(x))*Ci(x, s - 1))
    end
    CC = ComplexField(prec(parent(x)))
    return imag(Li(exp(CC(zero(x), x)), s))
end

function Si(x::T, s) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = imag(Li(CC(real(z), imag(z)), CC(s)))
    return convert(float(T), res)
end

function Si(x::T, s::Integer) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = imag(Li(CC(real(z), imag(z)), s))
    return convert(float(T), res)
end

"""
    hypbeta(a, b, z)
Compute the incomplete beta function B(a, b; z)
"""
function hypbeta(a::acb, b::acb, z::acb)
    res = parent(z)()
    ccall(("acb_hypgeom_beta_lower", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Cint, Clong), res, a, b, z, 0, prec(parent(z)))
    return res
end

function hypbeta(a::arb, b::arb, z::arb)
    res = parent(z)()
    ccall(("arb_hypgeom_beta_lower", Nemo.libarb), Cvoid,
          (Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Cint, Clong), res, a, b, z, 0, prec(parent(z)))
    return res
end

function hypbeta(a, b, z::T) where {T}
    CC = ComplexField(precision(BigFloat))
    res = hypbeta(CC(a), CC(b), CC(z))
    accuracy = min(
        ArbTools.rel_accuracy_bits(real(res)),
        ArbTools.rel_accuracy_bits(imag(res))
    )
    if accuracy < 53
        @warn "Unsufficient precision $accuracy in hypbeta for z = $z"
        @show res
    end
    return convert(complex(float(T)), res)
end

function hypbeta(a, b, z::T) where {T <: Real}
    RR = RealField(precision(BigFloat))
    res = hypbeta(R(a), RR(b), RR(z))
    accuracy = ArbTools.rel_accuracy_bits(res)
    if accuracy < 53
        @warn "Unsufficient precision $accuracy in hypbeta for z = $z"
        @show res
    end
    return convert(float(T), res)
end

"""
    powpos(x, y)
Compute |x|^y in a way that works if x contains negative numbers.
"""
function abspow(x::arb, y::arb)
    if contains_zero(x)
        x_upp = ArbTools.abs_ubound(x)
        return ArbTools.setinterval(zero(x), x_upp^y)
    end

    return abs(x)^y
end

abspow(x, y) = abs(x)^y
