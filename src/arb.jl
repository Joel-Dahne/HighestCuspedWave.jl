export iswide, Li, Ci, Si

"""
    iswide(x)
Return true if x is wide in the meaning that the effective relative
accuracy of x measured in bits is more than 10 lower than it's parents
precision.
"""
function iswide(x::arb)
    # TODO: This might require some tuning, 10 might not be the
    # optimal number.
    return ArbTools.rel_accuracy_bits(x) < prec(parent(x)) - 10
end

function iswide(x::acb)
    # TODO: Arb implements this as well, slightly different, but it
    # has to be implemented in ArbTools first
    return min(
        ArbTools.rel_accuracy_bits(real(x)),
        ArbTools.rel_accuracy_bits(imag(x)),
    ) < prec(parent(x)) - 10
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
    RR = RealField(precision(BigFloat))
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

TODO: Optimize it for wide s values
"""
function Li(z::acb, s::acb)
    if iswide(s)
        # TODO: Check that everything here is correct

        # TODO: Tune this
        n = 3 # Degree of Taylor expansion

        s_mid = parent(s)(midpoint(real(s)), midpoint(imag(s)))
        PP = AcbPolyRing(parent(s), :x)
        w = PP()

        # Compute the rest term of the Taylor expansion
        s_poly = PP([s, one(s)])
        ccall(("acb_poly_polylog_series", Nemo.libarb), Cvoid,
              (Ref{acb_poly}, Ref{acb_poly}, Ref{acb}, Int, Int),
              w, s_poly, z, n + 2, prec(parent(z)))
        restterm = (s - s_mid)^n*coeff(w, n)

        # Compute the Taylor polynomial at the midpoint of x
        s_poly = PP([s_mid, one(s)])
        ccall(("acb_poly_polylog_series", Nemo.libarb), Cvoid,
              (Ref{acb_poly}, Ref{acb_poly}, Ref{acb}, Int, Int),
              w, s_poly, z, n + 1, prec(parent(z)))

        # Evaluate the Taylor polynomial on s - s_mid and add the rest
        # term
        return evaluate(w, s - s_mid) + restterm
    end
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
    Ci(x::arb_series, s, n::Integer = length(x))
Compute `n` terms of the Taylor series of Ciₛ(x).

It's computed by directly computing the Taylor coefficients by
differentiating Ciₛ and then composing with `x`.
"""
function Ci(x::arb_series, s, n::Integer = length(x))
    res = arb_series(parent(x.poly)(), n)
    x₀ = x[0]

    for i in 0:n-1
        if i%2 == 0
            res[i] = (-1)^(div(i, 2))*Ci(x₀, s - i)/factorial(i)
        else
            res[i] = -(-1)^(div(i, 2))*Si(x₀, s - i)/factorial(i)
        end
    end

    # Compose the Taylor series for the Clausian with that of the
    # input
    x_tmp = arb_series(deepcopy(x.poly))
    x_tmp[0] = base_ring(parent(x.poly))(0)

    return Nemo.compose(res, x_tmp, n)
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
    beta_inc(a, b, z)
Compute the (not regularised) incomplete beta function B(a, b; z)
"""
function beta_inc(a::acb, b::acb, z::acb)
    res = parent(z)()
    ccall(("acb_hypgeom_beta_lower", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Cint, Clong), res, a, b, z, 0, prec(parent(z)))
    return res
end

function beta_inc(a::arb, b::arb, z::arb)
    res = parent(z)()
    ccall(("arb_hypgeom_beta_lower", Nemo.libarb), Cvoid,
          (Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Cint, Clong), res, a, b, z, 0, prec(parent(z)))
    return res
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
