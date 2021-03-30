# Contains methods related to Ci and Si
export Li, Ci, Si

"""
    Li(z, s)
Compute the polylogarithm ``Li_s(z)``
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
              w, s_poly, z, n + 1, prec(parent(z)))
        restterm = (s - s_mid)^n*coeff(w, n)
        # Compute the Taylor polynomial at the midpoint of x
        s_poly = PP([s_mid, one(s)])
        ccall(("acb_poly_polylog_series", Nemo.libarb), Cvoid,
              (Ref{acb_poly}, Ref{acb_poly}, Ref{acb}, Int, Int),
              w, s_poly, z, n, prec(parent(z)))

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

function Li(z::Acb, s::Union{Acb,Integer})
    if iswide(s)
        # TODO: Check that everything here is correct

        # TODO: Tune this
        n = 3 # Degree of Taylor expansion

        s_mid = Acb(Arblib.midref(real(s)), Arblib.midref(imag(s)))

        # Compute the rest term of the Taylor expansion
        w = Arblib.polylog_series!(
            AcbSeries(degree = n, prec = precision(z)),
            AcbSeries([s, 1]),
            z,
            n + 1,
        )

        restterm = (s - s_mid)^n * w[n]

        # Compute the Taylor polynomial at the midpoint of x
        w = Arblib.polylog_series!(
            AcbSeries(degree = n, prec = precision(z)),
            AcbSeries([s_mid, 1]),
            z,
            n + 1,
        )

        # Evaluate the Taylor polynomial on s - s_mid and add the rest
        # term
        return w(s - s_mid) + restterm
    end

    return Arblib.polylog!(zero(z), s, z)
end

"""
    Li(z, s, β)
Compute the polylogarithm ``Li_s^{(β)}(z)``.

That is, `Li(z, s)` differentiated `β` times w.r.t. `s` evaluated at `z`.
"""
function Li(z::acb, s::acb, β::Integer)
    res = parent(z)()
    PP = AcbPolyRing(parent(z), :x)

    w = PP()
    s_poly = PP([s, one(s)])

    ccall(("acb_poly_polylog_series", Nemo.libarb), Cvoid,
          (Ref{acb_poly}, Ref{acb_poly}, Ref{acb}, Clong, Clong),
          w, s_poly, z, β + 1, prec(parent(z)))

    return coeff(w, β)*factorial(β)
end

function Li(z::Acb, s::Acb, β::Integer)
    s_poly = AcbSeries([s, 1])
    w = Arblib.polylog_series!(AcbSeries(degree = β, prec = precision(z)), s_poly, z, β + 1)

    return w[β]*factorial(β)
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
tighter enclosure by using that `Ci` is 2π periodic, monotonic for `x ∈
[0, π]` and even, so that it's enough to evaluate on the endpoints and
possibly at zero or π if `x` contains points on the form `2kπ` or (2k
+ 1)π` respectively.
"""
function Ci(x::acb, s)
    im = x.parent(0, 1)
    return (Li(exp(im*x), s) + Li(exp(-im*x), s))/2
end

function Ci(x::arb, s::arb)
    if iswide(x)
        xₗ = ArbTools.lbound(x)
        xᵤ = ArbTools.ubound(x)
        (include_zero, include_pi) = contains_pi(xₗ, xᵤ)
        res = setunion(Ci(xₗ, s), Ci(xᵤ, s))
        if include_zero
            res = setunion(res, Ci(zero(x), s))
        end
        if include_pi
            res = setunion(res, Ci(parent(x)(π), s))
        end
        return res
    end
    CC = ComplexField(prec(parent(x)))
    return real(Li(exp(CC(zero(x), x)), CC(s)))
end

function Ci(x::arb, s::Integer)
    if iswide(x)
        xₗ = ArbTools.lbound(x)
        xᵤ = ArbTools.ubound(x)
        (include_zero, include_pi) = contains_pi(xₗ, xᵤ)
        res = setunion(Ci(xₗ, s), Ci(xᵤ, s))
        if include_zero
            res = setunion(res, Ci(zero(x), s))
        end
        if include_pi
            res = setunion(res, Ci(parent(x)(π), s))
        end
        return res
    end
    CC = ComplexField(prec(parent(x)))
    return real(Li(exp(CC(zero(x), x)), s))
end

# TODO: Optimize for wide x
function Ci(x::Arb, s::Arb)
    real(Li(exp(Acb(0, x)), Acb(s)))
end

# TODO: Optimize for wide x
function Ci(x::Arb, s::Integer)
    real(Li(exp(Acb(0, x)), s))
end

Ci(x::Real, s::Arb) = Ci(Arb(x), s)

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
    Ci(x, s, β)
Compute the Clausian function Ciₛ^(β)(x).

That is, `Ciₛ` differentiated `β` times w.r.t. `s` evaluated at `x`.

TODO: If x is a wide (real) ball (as determined by iswide(x)) it
computes a tighter enclosure by using that Ci 2π periodic, monotonic
for x ∈ [0, π] and even, so that it's enough to evaluate on the
endpoints and possibly at zero or π if `x` contains points on the form
`2kπ` or (2k + 1)π` respectively.
"""
function Ci(x::acb, s, β::Integer)
    im = x.parent(0, 1)
    return (Li(exp(im*x), s, β) + Li(exp(-im*x), s, β))/2
end

function Ci(x::arb, s, β::Integer)
    CC = ComplexField(prec(parent(x)))
    return real(Li(exp(CC(zero(x), x)), CC(s), β))
end

# TODO: Optimize for wide x
function Ci(x::Acb, s, β::Integer)
    s = convert(Acb, s)
    return (Li(exp(im*x), s, β) + Li(exp(-im*x), s, β))/2
end

# TODO: Optimize for wide x
function Ci(x::Arb, s, β::Integer)
    return real(Li(exp(Acb(0, x)), convert(Acb, s), β))
end

function Ci(x::T, s, β::Integer) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = real(Li(CC(real(z), imag(z)), CC(s), β))
    return convert(float(T), res)
end

# For non-integer values of β we don't have an Arb-implementation and
# fall back to a finite sum. This is extremely inefficient and
# NON-RIGOROUS.
# TODO: Why does this have the wrong sign???
function Ci(x, s, β; N = 1000)
    res = zero(x)

    for k in 1:N
        res += cos(k * x) * log(k)^β / k^s
    end

    return res
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
    Ci_expansion(x, s, M::Integer)
Compute the asymptotic expansion of `Ci` at zero up to order `2M - 2`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms in an `arb_series` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
enciso18:convex_whith.

It satisfies that `Ci(y, s) ∈ C*abs(y)^e + P(y) + E*y^(2M)` for all `|y|
<= |x|`.
"""
function Ci_expansion(x::arb, s::arb, M::Integer)
    π = parent(x)(pi)

    # Non-analytic term
    C = Nemo.gamma(1 - s)*sinpi(s/2)
    e = s - 1

    # Analytic term
    P = arb_series(ArbPolyRing(parent(x), :x)(), 2M - 1)
    for m = 0:M-1
        P[2m] = (-1)^m*zeta(s - 2m)/factorial(fmpz(2m))
    end

    # Error term
    E = ball(
        zero(x),
        2(2π)^(1 + s - 2M)*zeta(2M + 1 - s)/(4π^2 - x^2)
    )

    return (C, e, P, E)
end

"""
    Si(x, s)
Compute the Clausian function Siₛ(x).

If x is a wide (real) ball (as determined by iswide(x)) it computes a
tighter enclosure by first checking if the derivative doesn't contains
zero, if not it uses monotonicity to only evaluate at endpoints. If
the derivative does contain zero it uses a zero order approximation
instead.
"""
function Si(x::acb, s::acb)
    im = x.parent(0, 1)
    return (Li(exp(im*x), s) - Li(exp(-im*x), s))/2
end

function Si(x::arb, s::arb)
    if iswide(x)
        # Compute derivative
        dSi = Ci(x, s - 1)
        if contains_zero(dSi)
            # Use a zero order approximation
            ball(Si(midpoint(x), s), (x - midpoint(x))*dSi)
        else
            # Use that it's monotone
            xₗ, xᵤ = getinterval(x)
            return setunion(Si(xₗ, s), Si(xᵤ, s))
        end
    end
    CC = ComplexField(prec(parent(x)))
    return imag(Li(exp(CC(zero(x), x)), CC(s)))
end

function Si(x::arb, s::Integer)
    if iswide(x)
        # Compute derivative
        dSi = Ci(x, s - 1)
        if contains_zero(dSi)
            # Use a zero order approximation
            ball(Si(midpoint(x), s), (x - midpoint(x))*dSi)
        else
            # Use that it's monotone
            xₗ, xᵤ = getinterval(x)
            return setunion(Si(xₗ, s), Si(xᵤ, s))
        end
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
    Si_expansion(x, s, M::Integer)
Compute the asymptotic expansion of `Si` at zero up to order `2M - 1`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms in an `arb_series` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
enciso18:convex_whith.

It satisfies that `Si(y, s) ∈ C*sign(y)*abs(y)^e + P(y) + E*abs(y)^(2M + 1)`
for all `|y| <= |x|`.
"""
function Si_expansion(x::arb, s::arb, M::Integer)
    π = parent(x)(pi)

    # Non-analytic term
    C = Nemo.gamma(1 - s)*cospi(s/2)
    e = s - 1

    # Analytic term
    P = arb_series(ArbPolyRing(parent(x), :x)(), 2M)
    for m = 0:M-1
        P[2m + 1] = (-1)^m*zeta(s - 2m - 1)/factorial(fmpz(2m + 1))
    end

    # Error term
    E = ball(
        zero(x),
        2(2π)^(s - 2M)*zeta(2M + 2 - s)/(4π^2 - x^2)
    )

    return (C, e, P, E)
end

# For non-integer values of β we don't have an Arb-implementation and
# fall back to a finite sum. This is extremely inefficient and
# NON-RIGOROUS.
# TODO: Why does this have the wrong sign???
function Si(x, s, β; N = 1000)
    res = zero(x)

    for k in 1:N
        res += sin(k * x) * log(k)^β / k^s
    end

    return res
end
