export iswide, Li, Ci, Si

# For conversion from arb to Arb
function Arblib.set!(res::Arb, x::arb)
    ccall(Arblib.@libarb("arb_set"), Cvoid, (Ref{Arblib.arb_struct}, Ref{Nemo.arb}), res, x)
    return res
end

# For conversion from Arb to arb
function (r::ArbField)(x::Arb)
    z = arb(fmpz(0), r.prec)
    ccall(Arblib.@libarb("arb_set"), Cvoid, (Ref{Nemo.arb}, Ref{Arblib.arb_struct}), z, x)
    z.parent = r
    return z
end

"""
    mince(xₗ::arb, xᵤ::arb, n::Integer; split = false)
Return a vector with `n` balls covering the interval `[xₗ, xᵤ]`.

If `split` is true then returns the intervals as tuples and not balls.
"""
function mince(xₗ::arb, xᵤ::arb, n::Integer; split = false)
    intervals = Vector{ifelse(split, NTuple{2,arb}, arb)}(undef, n)
    dx = (xᵤ - xₗ)/n
    for i in eachindex(intervals)
        yₗ = xₗ + (i - 1)*dx
        yᵤ = xₗ + i*dx
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
    iswide(x)
Return true if x is wide in the meaning that the effective relative
accuracy of x measured in bits is more than 10 lower than it's parents
precision.
"""
function iswide(x::arb; cutoff = 10)
    # TODO: This might require some tuning, 10 might not be the
    # optimal number.
    return ArbTools.rel_accuracy_bits(x) < prec(parent(x)) - cutoff
end

function iswide(x::acb; cutoff = 10)
    # TODO: Arb implements this as well, slightly different, but it
    # has to be implemented in ArbTools first
    return min(
        ArbTools.rel_accuracy_bits(real(x)),
        ArbTools.rel_accuracy_bits(imag(x)),
    ) < prec(parent(x)) - cutoff
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

    ccall((:acb_dirichlet_stieltjes, Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{fmpz}, Ref{acb}, Clong),
          res, fmpz(n), a, CC.prec)

    convert(type, real(res))
end
function stieltjes(::Type{Arb}, n::Integer)
    res = zero(Acb)
    a = one(Acb)

    ccall((:acb_dirichlet_stieltjes, Arblib.libarb), Cvoid,
          (Ref{Arblib.acb_struct}, Ref{fmpz}, Ref{Arblib.acb_struct}, Clong),
          res, fmpz(n), a, precision(res))

    return real(res)
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

Li(z::Acb, s::Union{Acb,Integer}) = Arblib.polylog!(zero(z), s, z)

"""
    Li(z, s, β)
Compute the polylogarithm `Liₛ^(β)(z)`.

That is, `Liₛ` differentiated `β` times w.r.t. `s` evaluated at `z`.
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
    contains_pi(x1::arb, x2::arb)
Checks the interval `[x1, x2]` if it contains points of the form `kπ`.
Returns two booleans, the first one is false if it doesn't contain a
point on the form `2kπ` and the second one if it doesn't contain one
on the form `(2k+ 1)π`. If they are true it means they might contain
such a point. This is used to determine where the extrema of `Ci([x1,
x2])` can occur.

TODO: This is (hopefully) correct but not optimal.
"""
function contains_pi(x1::arb, x2::arb)
    @assert !(x1 > x2)

    # x1 or x2 equal to zero are the only cases when the division by π
    # can be exact, in which case it has to be handled differently.
    if iszero(x1)
        return (true, !(x2 < parent(x2)(π)))
    end
    if iszero(x2)
        return (true, !(x1 > -parent(x1)(π)))
    end

    # We have k1ₗπ ≤ xₗ < (k1ᵤ + 1)π
    k1 = floor(x1/parent(x1)(π))
    (unique1ₗ, k1ₗ) = unique_integer(ceil(ArbTools.lbound(k1)))
    (unique1ᵤ, k1ᵤ) = unique_integer(floor(ArbTools.ubound(k1)))
    @assert unique1ₗ && unique1ᵤ && k1ₗ*parent(x1)(π) ≤ x1 < (k1ᵤ + 1)*parent(x1)(π)
    # We have k2ₗπ ≤ xₗ < (k2ᵤ + 1)π
    k2 = floor(x2/parent(x2)(π))
    (unique2ₗ, k2ₗ) = unique_integer(ceil(ArbTools.lbound(k2)))
    (unique2ᵤ, k2ᵤ) = unique_integer(floor(ArbTools.ubound(k2)))
    @assert unique2ₗ && unique2ᵤ && k2ₗ*parent(x2)(π) ≤ x2 < (k2ᵤ + 1)*parent(x2)(π)

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


"""
    Ci(x, s)
Compute the Clausian function Ciₛ(x).

If x is a wide (real) ball (as determined by iswide(x)) it computes a
tighter enclosure by using that Ci 2π periodic, monotonic for x ∈ [0,
π] and even, so that it's enough to evaluate on the endpoints and
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

"""
    beta_inc(a, b, z)
Compute the (not regularised) incomplete beta function B(a, b; z).
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
    beta_inc_zeroone(a::arb, b::arb, z::arb)
Compute the (not regularised) incomplete beta function B(a, b; z)
assuming that `0 <= z <= 1`, discarding any other numbers in the
interval.

# PROVE: That it's monotonically increasing in `z`.
"""
function beta_inc_zeroone(a::arb, b::arb, z::arb)
    if 0 <= z <= 1
        return beta_inc(a, b, z)
    elseif z >= 0
        return setunion(beta_inc(a, b, ArbTools.lbound(z)), beta_inc(a, b, one(z)))
    elseif z <= 1
        return setunion(beta_inc(a, b, zero(z)), beta_inc(a, b, ArbTools.ubound(z)))
    else
        return setunion(beta_inc(a, b, zero(z)), beta_inc(a, b, one(z)))
    end
end


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
    ccall(("arb_hypgeom_2f1", Nemo.libarb), Cvoid,
          (Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Cint, Clong),
          res, a, b, c, z, 0, prec(parent(z)))
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
        return arb_series(parent(x.poly)(
            [abs(x[0]); fill(base_ring(parent(x.poly))(NaN), n - 1)]
        ))
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
    ccall(("acb_hypgeom_expint", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Int), res, s, z, prec(parent(z)))
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
    CC = ComplexField(prec(RR))
    res = CC()
    ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int), res, CC(a), CC(zero(z), z), 0, prec(RR))

    if iswide(res)
        res2 = CC()
        ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
              (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int), res2, CC(a), CC(zero(z), z), 0, 4prec(RR))
        res = CC(setintersection(real(res), real(res2)),
                 setintersection(imag(res), imag(res2)))
        if iswide(res, cutoff = 20)
            @warn "Wide enclosure when computing cosint res = $res"
        end
    end

    return real(exp(CC(zero(a), -RR(π)*a/2))*res)
end

"""
    cosintpi(a::arb, z)
Compute the generalised cosine integral Ciₐ(πz).

It tries to give better enclosures by evaluating it at higher
precision if required.
"""
function cosintpi(a::arb, z)
    RR = parent(a)
    CC = ComplexField(prec(RR))
    res = CC()
    ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int), res, CC(a), CC(zero(z), RR(π)*z), 0, prec(RR))

    if iswide(res)
        RR2 = RealField(8prec(RR))
        CC2 = ComplexField(prec(RR2))
        res2 = CC2()
        ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
              (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
              res2, CC2(a), CC2(zero(a), RR2(π)*RR2(z)), 0, prec(RR2))
        res = CC(setintersection(real(res), real(res2)),
                 setintersection(imag(res), imag(res2)))
        if iswide(res, cutoff = 20)
            @warn "Wide enclosure when computing cosintpi res = $res"
        end
    end

    return real(exp(CC(zero(a), -RR(π)*a/2))*res)
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
    lb = (s - (c + 1)/x)/x
    ub = (s - (c - 1)/x)/x
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
    CC = ComplexField(prec(RR))
    res = CC()
    ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int), res, CC(a), CC(zero(z), z), 0, prec(RR))

    if iswide(res)
        res2 = CC()
        ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
              (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int), res2, CC(a), CC(zero(z), z), 0, 4prec(RR))
        res = CC(setintersection(real(res), real(res2)),
                 setintersection(imag(res), imag(res2)))
        if iswide(res, cutoff = 20)
            @warn "Wide enclosure when computing sinint res = $res"
        end
    end

    return -imag(exp(CC(zero(a), -RR(π)*a/2))*res)
end

"""
    sinintpi(a::arb, z)
Compute the generalised sinine integral Siₐ(πz).

It tries to give better enclosures by evaluating it at higher
precision if required.
"""
function sinintpi(a::arb, z)
    RR = parent(a)
    CC = ComplexField(prec(RR))
    res = CC()
    ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int), res, CC(a), CC(zero(z), RR(π)*z), 0, prec(RR))

    if iswide(res)
        RR2 = RealField(8prec(RR))
        CC2 = ComplexField(prec(RR2))
        res2 = CC2()
        ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
              (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int),
              res2, CC2(a), CC2(zero(a), RR2(π)*RR2(z)), 0, prec(RR2))
        res = CC(setintersection(real(res), real(res2)),
                 setintersection(imag(res), imag(res2)))
        if iswide(res, cutoff = 20)
            @warn "Wide enclosure when computing sinintpi res = $res"
        end
    end

    return -imag(exp(CC(zero(a), -RR(π)*a/2))*res)
end
