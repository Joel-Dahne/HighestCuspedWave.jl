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
    Ci_expansion(x, s, M::Integer)
Compute the asymptotic expansion of `Ci` at zero up to order `2M - 2`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms in an `arb_series` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
enciso18:convex_whith.

It satisfies that `Ci(y, s) ∈ C*abs(y)^s + P(y) + E*y^(2M)` for all `|y|
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
    Si_expansion(x, s, M::Integer)
Compute the asymptotic expansion of `Si` at zero up to order `2M - 1`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms in an `arb_series` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
enciso18:convex_whith.

It satisfies that `Si(y, s) ∈ C*sign(y)*abs(y)^s + P(abs(y)) + E*y^(2M + 1)` for all `|y|
<= |x|`.
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
    if contains_zero(x)
        x_upp = ArbTools.abs_ubound(x)
        return ArbTools.setinterval(zero(x), x_upp^y)
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
Compute the generalised cosine integral Cₐ(z).

TODO: Depending on the choice of precision this gives very poor error
bounds. For example with 256 bits of precision, `a = 1.7` and `z =
175` it gives extremely bad bounds. Both 64 bits and 512 bits works
well.
"""
function cosint(a::arb, z::arb)
    RR = parent(z)
    CC = ComplexField(prec(RR))
    res = CC()
    ccall(("acb_hypgeom_gamma_upper", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Int, Int), res, CC(a), CC(zero(z), z), 0, prec(RR))
    return real(exp(CC(zero(a), -RR(π)*a/2))*res)
end
