export polylog, clausenc, clausencmzeta, clausens

###
### polylog
###

"""
    polylog(s, z)

Compute the polylogarithm \$Li_s(z)\$.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.
"""
function polylog(s::Union{Acb,Integer}, z::Acb)
    if iswide(s) # If this is true then s is always an Acb
        # Degree of Taylor expansion, could possibly be tuned
        degree = 2

        s_mid = Acb(Arblib.midref(Arblib.realref(s)), Arblib.midref(Arblib.imagref(s)))

        # Compute the rest term of the Taylor expansion
        w = Arblib.polylog_series!(
            AcbSeries(degree = degree + 1, prec = precision(z)),
            AcbSeries([s, 1]),
            z,
            degree + 2,
        )

        restterm = (s - s_mid)^(degree + 1) * w[degree+1]

        # Compute the Taylor polynomial at the midpoint of x
        w_mid = Arblib.polylog_series!(
            AcbSeries(prec = precision(z); degree),
            AcbSeries([s_mid, 1]),
            z,
            degree + 1,
        )

        # Evaluate the Taylor polynomial on s - s_mid and add the rest
        # term
        res = w_mid(s - s_mid) + restterm

        # If the resulting enclosure is not contained in the enclosure
        # coming from w[0] then take their intersection. Notice that
        # they will always intersect, so taking the intersection is
        # always okay.
        if !Arblib.contains(w[0], res)
            Arblib.intersection!(
                Arblib.realref(res),
                Arblib.realref(res),
                Arblib.realref(w[0]),
            )
            Arblib.intersection!(
                Arblib.imagref(res),
                Arblib.imagref(res),
                Arblib.imagref(w[0]),
            )
        end

        return res
    end

    return Arblib.polylog!(zero(z), s, z)
end

###
### clausenc
###

"""
    _clausenc_polylog(x::Arb, s::Union{Arb,Integer})

Evaluation of the `clausenc` function through the polylog function.

It uses the formula
```
clausenc(x, s) = real(polylog(exp(im * x), s))
```
"""
function _clausenc_polylog(x::Arb, s::Union{Arb,Integer})
    z = exp(Acb(0, x, prec = precision(x)))
    s = s isa Integer ? s : Acb(s, prec = precision(x))
    return real(polylog(s, z))
end

"""
    _clausenc_zeta(x::Arb, s::Arb)

Evaluation of the `clausenc` function through the zeta function.

It uses the formula
```
clausenc(x, s) = let v = 1 - s
    gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
end
```
Based on formula [25.13.2](https://dlmf.nist.gov/25.13) for the
periodic zeta function and then taking the real part.

It currently only handles `0 < x < 2π` and `s` not overlapping any
integer.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.
"""
function _clausenc_zeta(x::Arb, s::Arb)
    # Check that x > 0
    Arblib.ispositive(x) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    inv2pi = inv(2Arb(π, prec = precision(x)))
    xinv2pi = x * inv2pi
    onemxinv2pi = let onemxinv2pi = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(onemxinv2pi, xinv2pi)
        Arblib.add!(onemxinv2pi, onemxinv2pi, 1)
    end

    # Check that 1 - x / 2π > 0, i.e. x < 2π
    Arblib.ispositive(onemxinv2pi) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    v = let v = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(v, s)
        Arblib.add!(v, v, 1)
    end

    f(v) =
        SpecialFunctions.gamma(v) *
        inv2pi^v *
        cospi(v / 2) *
        (SpecialFunctions.zeta(v, xinv2pi) + SpecialFunctions.zeta(v, onemxinv2pi))

    if iswide(s)
        res = ArbExtras.extrema_series(f, Arblib.getinterval(v)..., degree = 2)[1:2]
        return Arb(res)
    end

    return f(v)
end

"""
    _clausenc_zeta(x::Acb, s::Arb)

Evaluation of the `clausenc` function through the zeta function.

This uses the same formula as the method with `x::Arb`.

It currently only handles `0 < real(x) < 2π` and `s` not overlapping
any integer.

Not that this does **not** handle wide values of `s` in any special
way. This method is currently only used in the integration for
bounding the error term and in that case getting the most accurate
bound is not important. It could be something to consider later on.

- *TODO:* Check if this formula holds for complex `x`, it seems to
   give correct results at least.
"""
function _clausenc_zeta(x::Acb, s::Arb)
    # Check that real(x) > 0
    Arblib.ispositive(Arblib.realref(x)) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    inv2pi = inv(2Arb(π, prec = precision(x)))
    xinv2pi = x * inv2pi
    onemxinv2pi = let onemxinv2pi = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(onemxinv2pi, xinv2pi)
        Arblib.add!(onemxinv2pi, onemxinv2pi, 1)
    end

    # Check that real(1 - x / 2π) > 0, i.e. real(x) < 2π
    Arblib.ispositive(Arblib.realref(onemxinv2pi)) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    v = let v = Arb(prec = precision(x)) # We do it like this to preserve the precision
        Arblib.neg!(v, s)
        Arblib.add!(v, v, 1)
    end

    f(v) =
        SpecialFunctions.gamma(v) *
        inv2pi^v *
        cospi(v / 2) *
        (
            SpecialFunctions.zeta(Acb(v), xinv2pi) +
            SpecialFunctions.zeta(Acb(v), onemxinv2pi)
        )

    return f(v)
end

"""
    clausenc(x, s)

Compute the Clausen function \$C_s(x)\$.

If `x` is a wide (real) ball, as determined by `iswide(x)`, it
computes a tighter enclosure by using that the function is
`2π`-periodic, monotonic for `x ∈ [0, π]` and even, so that it's
enough to evaluate on the endpoints and possibly at zero or `π` if `x`
contains points on the form `2kπ` or (2k + 1)π` respectively. In the
wide case it computes the endpoints at a reduced precision given by
```
prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
```

**TODO:** Figure out how to compute accurate bounds for `s`
overlapping an integer or being close to an integer. The current
approach doesn't work if `s overlaps an integer and gives terrible
bounds for wide `s` close to an integer.
"""
function clausenc(x::Arb, s::Union{Arb,Integer})
    if iswide(x)
        prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
        xₗ, xᵤ = Arblib.getinterval(Arb, setprecision(x, prec))
        (include_zero, include_pi) = contains_pi(xₗ, xᵤ)
        res = union(clausenc(xₗ, s), clausenc(xᵤ, s))
        if include_zero
            res = union(res, clausenc(zero(xₗ), s))
        end
        if include_pi
            res = union(res, clausenc(Arb(π; prec), s))
        end
        return setprecision(res, precision(x))
    end

    # If s is not an integer and 0 < x < 2π call _clausenc_zeta(x, s)
    if s isa Arb && !Arblib.contains_int(s)
        if Arblib.ispositive(x) && x < 2Arb(π)
            return _clausenc_zeta(x, s)
        end
    end

    return _clausenc_polylog(x, s)
end

function clausenc(x::Acb, s::Arb)
    # If s is not an integer and 0 < x < 2π call _clausenc_zeta(x, s)
    if s isa Arb && !Arblib.contains_int(s)
        if Arblib.ispositive(Arblib.realref(x)) && Arblib.realref(x) < 2Arb(π)
            return _clausenc_zeta(x, s)
        end
    end

    s = Acb(s)
    return (polylog(s, exp(im * x)) + polylog(s, exp(-im * x))) / 2
end

function clausenc(x::Acb, s::Union{Acb,Integer})
    if s isa Acb && isreal(s)
        clausenc(x, real(s))
    else
        (polylog(s, exp(im * x)) + polylog(s, exp(-im * x))) / 2
    end
end

clausenc(x::S, s::T) where {S<:Real,T<:Real} = convert(
    float(promote_type(S, T)),
    clausenc(convert(Arb, x), s isa Integer ? s : convert(Arb, s)),
)

"""
    clausenc(x::ArbSeries, s)

Compute the Taylor series of the Clausen function \$C_s(x)\$.

It's computed by directly computing the Taylor coefficients by
differentiating `clausenc` and then composing with `x`.
"""
function clausenc(x::ArbSeries, s)
    res = zero(x)
    x₀ = x[0]

    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausenc(x₀, s - i) / factorial(i)
        else
            res[i] = -(-1)^(i ÷ 2) * clausens(x₀, s - i) / factorial(i)
        end
    end

    # Compose the Taylor series for the result with that of the input
    x[0] = 0
    Arblib.compose_series!(res, res, x, Arblib.degree(res) + 1)
    x[0] = x₀

    return res
end

"""
    clausenc_expansion(x, s, M::Integer)

Compute the asymptotic expansion of `clausenc(x, s)` at zero up to
order `2M - 2`, meaning that the error term is of order `2M`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a `ArbSeries` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
enciso18:convex_whith.

It satisfies that `clausenc(y, s) ∈ C*abs(y)^e + P(y) + E*y^(2M)` for
all `|y| <= |x|`.

If `skip_constant = true` it doesn't compute the constant term in the
expansion. This is useful if you want to compute the expansion for
`clausenc(x, s) - clausenc(0, s)`.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure of the coefficients using a Taylor expansion in `s`.
"""
function clausenc_expansion(x::Arb, s::Arb, M::Integer; skip_constant = false)
    Arblib.ispositive(s) || throw(ArgumentError("s must be positive"))
    # TODO: Check this
    M > (s + 1) / 2 || throw(ArgumentError("M must be larger that (s + 1) / 2"))

    π = oftype(x, pi)

    # Non-analytic term
    if s == 2
        C = -π / 2
    else
        if iswide(s)
            C = Arb(
                ArbExtras.extrema_series(
                    s -> SpecialFunctions.gamma(1 - s) * sinpi(s / 2),
                    Arblib.getinterval(s)...,
                    degree = 2,
                )[1:2],
            )
        else
            C = SpecialFunctions.gamma(1 - s) * sinpi(s / 2)
        end
    end
    e = s - 1

    # Analytic term
    P = ArbSeries(degree = 2M - 2, prec = precision(x))
    start = skip_constant ? 1 : 0
    for m = start:M-1
        if iswide(s)
            z = Arb(
                ArbExtras.extrema_series(s -> zeta(s - 2m), Arblib.getinterval(s)...)[1:2],
            )
        else
            z = zeta(s - 2m)
        end
        P[2m] = (-1)^m * z / factorial(2m)
    end

    # Error term
    E = Arblib.add_error!(zero(x), 2(2π)^(1 + s - 2M) * zeta(2M + 1 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

###
### clausencmzeta
###

"""
    _clausencmzeta_zeta(x::Arb, s::Arb)

Evaluation of the `clausencmzeta` function through the zeta function.

This method is similar to [`_clausenc_zeta`](@ref) with the addition
that for wide `s` it also takes into account the subtraction of
`zeta(s)` when computing the enclosure.
"""
function _clausencmzeta_zeta(x::Arb, s::Arb)
    # Implements y -> 1 - y in a way that preserves the precision
    onem(y::Arb) =
        let res = zero(y)
            Arblib.neg!(res, y)
            Arblib.add!(res, res, 1)
        end
    onem(y::ArbSeries) = 1 - y # This one already preserves precision

    # Check that x > 0
    Arblib.ispositive(x) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    inv2pi = inv(2Arb(π, prec = precision(x)))
    xinv2pi = x * inv2pi
    onemxinv2pi = onem(xinv2pi)

    # Check that 1 - x / 2π > 0, i.e. x < 2π
    Arblib.ispositive(onemxinv2pi) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    v = onem(s)

    f(v) =
        SpecialFunctions.gamma(v) *
        inv2pi^v *
        cospi(v / 2) *
        (SpecialFunctions.zeta(v, xinv2pi) + SpecialFunctions.zeta(v, onemxinv2pi)) -
        zeta(onem(v))

    if iswide(s)
        res = ArbExtras.extrema_series(f, Arblib.getinterval(v)..., degree = 2)[1:2]
        return Arb(res)
    end

    return f(v)
end

"""
    clausencmzeta(x, s)

Compute `clausenc(x, s) - zeta(s)`. Notice that `clausenc(0, s) =
zeta(s)` so this is `clausenc` normalized to be zero at `x = 0`.

Typically this method just calls [`clausenc`](@ref) and [`zeta`](@ref)
directly and gives no direct benefit.

However, if `s` is a wide ball it can better handle the cancellations
between the two terms by computing them correctly. This is handled by
the method [`_clausencmzeta_zeta`](@ref). Since this method only
supports `0 < x < π` and `s` not overlapping an integer we fall back
to calling the methods directly if this is not satisfied.

Otherwise it behaves like `clausenc(x, s) - zeta(s)` with the only
difference being that it converts `x` and `s` to the same type to
begin with.
"""
function clausencmzeta(x::Arb, s::Arb)
    if !iswide(s) || Arblib.contains_int(s) || !Arblib.ispositive(x) || !(x < 2Arb(π))
        # If s is not wide or we are in a case which
        # _clausencmzeta_zeta doesn't support, call the clausenc and
        # zeta functions directly.
        return clausenc(x, s) - zeta(s)
    end

    if iswide(x)
        prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
        xₗ, xᵤ = Arblib.getinterval(Arb, setprecision(x, prec))
        # We know that 0 < x < π so it is always monotonically
        # decreasing
        res = Arb((clausencmzeta(xᵤ, s), clausencmzeta(xₗ, s)))

        return setprecision(res, precision(x))
    end

    return _clausencmzeta_zeta(x, s)
end

function clausencmzeta(x::ArbSeries, s::Arb)
    res = zero(x)
    x₀ = x[0]

    res[0] = clausencmzeta(x₀, s)
    for i = 1:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausenc(x₀, s - i) / factorial(i)
        else
            res[i] = -(-1)^(i ÷ 2) * clausens(x₀, s - i) / factorial(i)
        end
    end

    # Compose the Taylor series for the result with that of the input
    x[0] = 0
    Arblib.compose_series!(res, res, x, Arblib.degree(res) + 1)
    x[0] = x₀

    return res
end

clausencmzeta(x::ArbSeries, s) = clausenc(x, s) - zeta(Arb(s, prec = precision(x)))

function clausencmzeta(x, s)
    x, s = promote(x, s)
    return clausenc(x, s) - zeta(s)
end

###
### clausens
###

"""
    _clausens_polylog(x::Arb, s::Union{Arb,Integer})

Evaluation of the `clausens` function through the polylog function.

It uses the formula
```
clausens(x, s) = imag(polylog(exp(im * x), s))
```
"""
function _clausens_polylog(x::Arb, s::Union{Arb,Integer})
    z = exp(Acb(0, x, prec = precision(x)))
    s = s isa Integer ? s : Acb(s, prec = precision(x))
    return imag(polylog(s, z))
end

"""
    _clausens_zeta(x::Arb, s::Arb)

Evaluation of the `clausens` function through the zeta function.

It uses the formula
```
clausens(x, s) = let v = 1 - s
    gamma(v) * inv(2π)^v * sinpi(v / 2) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
end
```
Based on formula [25.13.2](https://dlmf.nist.gov/25.13) for the
periodic zeta function and then taking the imaginary part.

It currently only handles `0 < x < 2π` and `s` not overlapping any
integer.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.
"""
function _clausens_zeta(x::Arb, s::Arb)
    # Check that x > 0
    Arblib.ispositive(x) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    inv2pi = inv(2Arb(π, prec = precision(x)))
    xinv2pi = x * inv2pi
    onemxinv2pi = let onemxinv2pi = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(onemxinv2pi, xinv2pi)
        Arblib.add!(onemxinv2pi, onemxinv2pi, 1)
    end

    # Check that 1 - x / 2π > 0, i.e. x < 2π
    Arblib.ispositive(onemxinv2pi) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    v = let v = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(v, s)
        Arblib.add!(v, v, 1)
    end

    f(v) =
        SpecialFunctions.gamma(v) *
        inv2pi^v *
        sinpi(v / 2) *
        (SpecialFunctions.zeta(v, xinv2pi) - SpecialFunctions.zeta(v, onemxinv2pi))

    if iswide(s)
        res = ArbExtras.extrema_series(f, Arblib.getinterval(v)..., degree = 2)[1:2]
        return Arb(res)
    end

    return f(v)
end

"""
    clausens(x, s)

Compute the Clausen function \$S_s(x)\$.

If `x` is a wide (real) ball, as determined by `iswide(x)`, it
computes a tighter enclosure by first checking if the derivative
doesn't contains zero, if not it uses monotonicity to only evaluate at
endpoints. If the derivative does contain zero it uses a zeroth order
approximation instead. In the wide case it computes the endpoints at a
reduced precision given by
```
prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
```
"""
function clausens(x::Arb, s::Union{Arb,Integer})
    if iswide(x)
        orig_prec = precision(x)
        prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
        x = setprecision(x, prec)
        # Compute derivative
        dclausens = clausenc(x, s - 1)
        if Arblib.contains_zero(dclausens)
            # Use a zero order approximation
            res = Arblib.add_error!(
                clausens(Arblib.midpoint(Arb, x), s),
                (x - Arblib.midpoint(Arb, x)) * dclausens,
            )
        else
            # Use that it's monotone
            xₗ, xᵤ = Arblib.getinterval(Arb, x)
            res = union(clausens(xₗ, s), clausens(xᵤ, s))
        end
        return setprecision(res, orig_prec)
    end

    # If s is not an integer and 0 < x < 2π call _clausenc_zeta(x, s)
    if s isa Arb && !Arblib.contains_int(s)
        if Arblib.ispositive(x) && x < 2Arb(π)
            return _clausens_zeta(x, s)
        end
    end

    return _clausens_polylog(x, s)
end

clausens(x::Acb, s) = (polylog(s, exp(im * x)) - polylog(s, exp(-im * x))) / 2

clausens(x::S, s::T) where {S<:Real,T<:Real} = convert(
    float(promote_type(S, T)),
    clausens(convert(Arb, x), s isa Integer ? s : convert(Arb, s)),
)
