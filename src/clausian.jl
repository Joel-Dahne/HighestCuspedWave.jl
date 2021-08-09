# Contains methods related to Ci and Si
export Li, Ci, Si, Ci_tilde

###
### Li
###

"""
    Li(z, s)
Compute the polylogarithm ``Li_s(z)``
"""
function Li(z::Acb, s::Union{Acb,Integer})
    if iswide(s) # If this is true then s is always an Acb
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

Li(z::acb, s::acb) = parent(z)(Li(Acb(z), Acb(s)))

Li(z::acb, s::Integer) = parent(z)(Li(Acb(z), s))

"""
    Li(z, s, β)
Compute the polylogarithm ``Li_s^{(β)}(z)``.

That is, `Li(z, s)` differentiated `β` times w.r.t. `s` evaluated at `z`.
"""
function Li(z::Acb, s::Acb, β::Integer)
    s_poly = AcbSeries([s, 1])
    w = Arblib.polylog_series!(AcbSeries(degree = β, prec = precision(z)), s_poly, z, β + 1)

    return w[β] * factorial(β)
end

Li(z::acb, s::acb, β::Integer) = parent(z)(Li(Acb(z), Acb(s), β))

###
### Ci
###

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
    return (Li(exp(im * x), s) + Li(exp(-im * x), s)) / 2
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
    CC = ComplexField(precision(parent(x)))
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
    CC = ComplexField(precision(parent(x)))
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

function Ci(x::T, s) where {T<:Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im * x)
    res = real(Li(CC(real(z), imag(z)), CC(s)))
    return convert(float(T), res)
end

function Ci(x::T, s::Integer) where {T<:Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im * x)
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
    return (Li(exp(im * x), s, β) + Li(exp(-im * x), s, β)) / 2
end

function Ci(x::arb, s, β::Integer)
    CC = ComplexField(precision(parent(x)))
    return real(Li(exp(CC(zero(x), x)), CC(s), β))
end

# TODO: Optimize for wide x
function Ci(x::Acb, s, β::Integer)
    s = convert(Acb, s)
    return (Li(exp(im * x), s, β) + Li(exp(-im * x), s, β)) / 2
end

# TODO: Optimize for wide x
function Ci(x::Arb, s, β::Integer)
    return real(Li(exp(Acb(0, x)), convert(Acb, s), β))
end

function Ci(x::T, s, β::Integer) where {T<:Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im * x)
    res = real(Li(CC(real(z), imag(z)), CC(s), β))
    return convert(float(T), res)
end

# For non-integer values of β we don't have an Arb-implementation and
# fall back to a finite sum. This is extremely inefficient and
# NON-RIGOROUS.
# The sign here does not agree with the implementation above for
# integer values of β, there is a factor (-1)^β missing. I'm still not
# sure how to handle this in a good way.
function Ci(x, s, β; N = 10000)
    res = zero(x)

    for k = ifelse(iszero(β), 1, 2):N
        res += cos(k * x) * log(oftype(x, k))^β / k^s
    end

    return res
end

function Ci_tilde(x, s, β; N = 100000)
    res = zero(x)

    for k = ifelse(iszero(β), 1, 2):N
        res += (cos(k * x) - 1) * log(oftype(x, k))^β / k^s
    end

    return res
end

function Ci_tilde(x::Arb, s::Arb, β::Arb; N = 100000)
    res = zero(x)
    tmp1, tmp2 = Arb(prec = 256), Arb(prec = 256)

    for k = ifelse(iszero(β), 1, 2):N
        #res += (cos(k * x) - 1) * log(oftype(x, k))^β / k^s
        #Arblib.mul!(tmp1, x, k)
        #Arblib.cos!(tmp1, tmp1)
        #Arblib.sub!(tmp1, tmp1, 1)

        #Arblib.set!(tmp2, k)
        #Arblib.log!(tmp2, tmp2)
        #Arblib.pow!(tmp2, tmp2, β)

        #Arblib.mul!(tmp1, tmp1, tmp2)

        #Arblib.set!(tmp2, k)
        #Arblib.pow!(tmp2, tmp2, s)
        #Arblib.div!(tmp1, tmp1, tmp2)

        #Arblib.add!(res, res, tmp1)
    end
    res += let x = Float64(x), s = Float64(s), β = Float64(β)
        res = 0

        for k = ifelse(iszero(β), 1, 2):N
            res += (cos(k * x) - 1) * log(k)^β / k^s
        end

        res
    end

    if Arblib.isnonpositive(β)
        # The rest term is between -2log(N + 1)^β * ζ(s, N + 1) and 0
        restterm = Arb((-2log(Arb(N + 1))^β * SpecialFunctions.zeta(s, Arb(N + 1)), 0))
        res += restterm
    else
        # Add rest term is bounded in absolute value by 2ζ(s - ϵ, N +
        # 1) with ϵ = β * log(log(N + 1)) / log(N + 1)
        # TODO: Potentially use the derivatives of ζ(s, N + 1) with
        # respect to s?
        ϵ = β * log(log(Arb(N + 1))) / log(Arb(N + 1))
        @show ϵ
        restterm = 2SpecialFunctions.zeta(s - ϵ, Arb(N + 1))
        Arblib.add_error!(res, restterm)
    end

    return res
end

"""
    Ci_alternative(x, s, β; N = 1000)

Alternative implementation of `Ci` for non-integer `β` using the
Euler-Maclaurin formula.
"""
function Ci_alternative(x::Arb, s::Arb, β::Arb; N = 1000, p = 2)
    @assert !iszero(β) # The implementations assumes β is non-zero

    res = zero(x)

    # Parameters for the Euler-Maclaurin formula using the notation
    # from
    # https://en.wikipedia.org/wiki/Euler%E2%80%93Maclaurin_formula
    f(j) = cos(j * x) / j^s * log(j)^β
    m = Arb(2, prec = precision(res))
    n = Arb(N, prec = precision(res)) # TODO: tune this

    # Integral
    integral_part = real(Arblib.integrate(f, m, n))
    @show integral_part
    res += integral_part

    # Endpoint values
    endpoint_part = zero(res)
    h = p ÷ 2

    fm = f(ArbSeries([m, 1], degree = 2h - 1))
    fn = f(ArbSeries([n, 1], degree = 2h - 1))

    res += (fn[0] + fm[0]) / 2

    for k = 1:h
        endpoint_part +=
            Arblib.bernoulli!(zero(res), unsigned(2k)) / 2k * (fn[2k-1] - fm[2k-1])
    end

    @show endpoint_part
    res += endpoint_part

    # Remainder term
    integrand(j; analytic::Bool) = begin
        j_floor = Arblib.real_floor!(zero(j), j, analytic)
        isnan(j_floor) && return j_floor
        # TODO: Check analyticity of f
        B = Arblib.bernoulli_poly!(zero(j), unsigned(p), j - j_floor)
        return f(AcbSeries([j, 1], degree = p))[p] * B
    end

    atol = Mag(Arblib.midref(res))
    atol = Arblib.mul_2exp!(atol, atol, Arblib.rel_error_bits(res))

    remainder =
        (-1)^(p + 1) * real(Arblib.integrate(integrand, m, n, check_analytic = true; atol))

    #remainder = zero(res)
    #for start = 2:N-1
    #    remainder +=
    #        (-1)^(p + 1) * real(
    #            Arblib.integrate(integrand, start, start + 1, check_analytic = true; atol)
    #        )
    #end

    @show remainder
    res += remainder

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

    for i = 0:n-1
        if i % 2 == 0
            res[i] = (-1)^(div(i, 2)) * Ci(x₀, s - i) / factorial(i)
        else
            res[i] = -(-1)^(div(i, 2)) * Si(x₀, s - i) / factorial(i)
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
Compute the asymptotic expansion of `Ciₛ(x)` at zero up to order `2M -
2`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a series `P` (`arb_series` or
`ArbSeries`) and the error term `E`. The `M` is the same as in Lemma
2.1 in enciso18:convex_whith.

It satisfies that `Ci(y, s) ∈ C*abs(y)^e + P(y) + E*y^(2M)` for all
`|y| <= |x|`.
"""
function Ci_expansion(x::arb, s::arb, M::Integer)
    π = parent(x)(pi)

    # Non-analytic term
    C = Nemo.gamma(1 - s) * sinpi(s / 2)
    e = s - 1

    # Analytic term
    P = arb_series(ArbPolyRing(parent(x), :x)(), 2M - 1)
    for m = 0:M-1
        P[2m] = (-1)^m * zeta(s - 2m) / factorial(fmpz(2m))
    end

    # Error term
    E = ball(zero(x), 2(2π)^(1 + s - 2M) * zeta(2M + 1 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

function Ci_expansion(x::Arb, s::Arb, M::Integer)
    Arblib.ispositive(s) || throw(ArgumentError("s must be positive"))
    # TODO: Check this
    M > (s + 1) / 2 || throw(ArgumentError("M must be larger that (s + 1) / 2"))

    π = oftype(x, pi)

    # Non-analytic term
    if s == 2
        C = -π / 2
    else
        C = SpecialFunctions.gamma(1 - s) * sinpi(s / 2)
    end
    e = s - 1

    # Analytic term
    P = ArbSeries(degree = 2M - 2, prec = precision(x))
    for m = 0:M-1
        P[2m] = (-1)^m * zeta(s - 2m) / factorial(oftype(x, 2m))
    end

    # Error term
    E = Arblib.add_error!(zero(x), 2(2π)^(1 + s - 2M) * zeta(2M + 1 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

"""
    Ci_expansion(x, s, β, M::Integer)
Compute the asymptotic expansion of `Ciₛ`(β)(x)` at zero up to order
`2M - 2`.

It currently only supports `β == 1`, for `β == 0` use `Ci_expansion(x,
s)` instead.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a series `P` (`arb_series` or
`ArbSeries`) and the error term `E`. The `M` is the same as in Lemma
2.1 in enciso18:convex_whith.

It satisfies that `Ci(y, s, β) ∈ C*abs(y)^e + P(y) + E*y^(2M)` for all
`|y| <= |x|`.

TODO: Implement these
"""
function Ci_expansion(x::Arb, s::Arb, β::Integer, M::Integer)
    β == 1 || throw(ArgumentError("only implemented for β = 1, got $β"))

end

###
### Si
###

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
    return (Li(exp(im * x), s) - Li(exp(-im * x), s)) / 2
end

function Si(x::arb, s::arb)
    if iswide(x)
        # Compute derivative
        dSi = Ci(x, s - 1)
        if contains_zero(dSi)
            # Use a zero order approximation
            ball(Si(midpoint(x), s), (x - midpoint(x)) * dSi)
        else
            # Use that it's monotone
            xₗ, xᵤ = getinterval(x)
            return setunion(Si(xₗ, s), Si(xᵤ, s))
        end
    end
    CC = ComplexField(precision(parent(x)))
    return imag(Li(exp(CC(zero(x), x)), CC(s)))
end

function Si(x::arb, s::Integer)
    if iswide(x)
        # Compute derivative
        dSi = Ci(x, s - 1)
        if contains_zero(dSi)
            # Use a zero order approximation
            ball(Si(midpoint(x), s), (x - midpoint(x)) * dSi)
        else
            # Use that it's monotone
            xₗ, xᵤ = getinterval(x)
            return setunion(Si(xₗ, s), Si(xᵤ, s))
        end
    end
    CC = ComplexField(precision(parent(x)))
    return imag(Li(exp(CC(zero(x), x)), s))
end

function Si(x::T, s) where {T<:Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im * x)
    res = imag(Li(CC(real(z), imag(z)), CC(s)))
    return convert(float(T), res)
end

function Si(x::T, s::Integer) where {T<:Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im * x)
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
    C = Nemo.gamma(1 - s) * cospi(s / 2)
    e = s - 1

    # Analytic term
    P = arb_series(ArbPolyRing(parent(x), :x)(), 2M)
    for m = 0:M-1
        P[2m+1] = (-1)^m * zeta(s - 2m - 1) / factorial(fmpz(2m + 1))
    end

    # Error term
    E = ball(zero(x), 2(2π)^(s - 2M) * zeta(2M + 2 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

# For non-integer values of β we don't have an Arb-implementation and
# fall back to a finite sum. This is extremely inefficient and
# NON-RIGOROUS.
# TODO: Why does this have the wrong sign???
function Si(x, s, β; N = 1000)
    res = zero(x)

    for k = ifelse(iszero(β), 1, 2):N
        res += sin(k * x) * log(k)^β / k^s
    end

    return res
end
