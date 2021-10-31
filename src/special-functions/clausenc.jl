export clausenc, clausencmzeta

"""
    _clausenc_polylog(x::Arb, s::Union{Arb,Integer})

Evaluation of the `clausenc` function through the polylog function.

It uses the formula
```
clausenc(x, s) = real(polylog(s, exp(im * x)))
```
"""
function _clausenc_polylog(x::Arb, s::Union{Arb,Integer})
    z = exp(Acb(0, x, prec = precision(x)))
    s = s isa Integer ? s : Acb(s, prec = precision(x))
    return real(polylog(s, z))
end

"""
    _clausenc_polylog(x::Arb, s::Arb, β::Integer)

Evaluation of the `clausenc(x, s, β)` function through the polylog
function.

It uses the formula
```
clausenc(x, s) = real(polylog(s, exp(im * x)))
```
and computes the derivative with respect to `s` using `AcbSeries`.
"""
function _clausenc_polylog(x::Arb, s::Arb, β::Integer)
    z = exp(Acb(0, x, prec = precision(x)))
    return real(polylog(AcbSeries([s, 1], degree = β), z)[β]) * factorial(β)
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

It currently only handles `0 < x < 2π`. If `s` overlaps with any
non-negative integer the result will be indeterminate.

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

    f(v) = gamma(v) * inv2pi^v * cospi(v / 2) * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))

    if iswide(s)
        res = ArbExtras.extrema_series(f, getinterval(v)..., degree = 2)[1:2]
        return Arb(res)
    end

    return f(v)
end

"""
    _clausenc_zeta(x::Acb, s::Arb)

Evaluation of the `clausenc` function through the zeta function.

This uses the same formula as the method with `x::Arb`.

It currently only handles `0 < real(x) < 2π`. If `s` overlaps with any
non-negative integer the result will be indeterminate.

Not that this does **not** handle wide values of `s` in any special
way. This method is currently only used in the integration for
bounding the error term and in that case getting the most accurate
bound is not important. It could be something to consider later on.

- **TODO:** Check if this formula holds for complex `x`, it seems to
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
        gamma(v) *
        inv2pi^v *
        cospi(v / 2) *
        (zeta(Acb(v), xinv2pi) + zeta(Acb(v), onemxinv2pi))

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

The case when `s` is a wide ball is in general handled by the
underlying methods [`_clausenc_polylog`](@ref) and
[`_clausenc_zeta`](@ref). The exception is when `s` overlaps with a
non-negative integer, in which case both the above methods give
indeterminate results. In that case we compute at the midpoint of `s`
and bound the error by using a bound for the derivative in `s`. For s
> 1 the derivative in `s` is bounded by `dzeta(s)`, this can be
seen by looking at the Fourier series for `clausenc(x, s, 1)` and
noticing that it attains it maximum at `x = 0` where it precisely
equals `dzeta(s)`.
- **TODO:** Figure out how to bound this for `s = 1` and `s = 0`. In
  this case the derivative in `s` blows up at `x = 0` so we can't use
  a uniform bound. For now we compute a bound assuming monotonicity in
  `s`, which is not true.
"""
function clausenc(x::Arb, s::Union{Arb,Integer})
    if iswide(x)
        prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
        xₗ, xᵤ = getinterval(Arb, setprecision(x, prec))
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

    contains_nonnegative_int = s isa Arb && Arblib.contains_int(s) && !Arblib.isnegative(s)

    # Handle the case when s contains a non-negative integer but it
    # not exactly an integer
    if contains_nonnegative_int && !iszero(Arblib.radref(s))
        if s > 1
            # Evaluate at midpoint
            smid = midpoint(Arb, s)
            res = clausenc(x, smid)

            # Bound derivative in s using derivative of zeta function
            derivative = dzeta(s)
            error = (s - smid) * abs(derivative)

            return res + error
        end

        # FIXME: For now we bound it assuming monotonicity in s, this
        # is not true in practice. We only do this for s close to 0
        # and 1
        @warn "Incorrect bound for clausenc with s = $s, it assumes monotonicity" maxlog =
            100
        if Arblib.radref(s) < 1e-2
            # s is not negative and not greater than 1, since the
            # radius is small it must therefore contain either 0 or 1.
            sₗ, sᵤ = getinterval(Arb, s)
            return union(clausenc(x, sₗ), clausenc(x, sᵤ))
        end
    end

    # If s is not a non-negative integer and 0 < x < 2π call
    # _clausenc_zeta(x, s)
    if s isa Arb && !contains_nonnegative_int
        if Arblib.ispositive(x) && x < 2Arb(π)
            return _clausenc_zeta(x, s)
        end
    end

    return _clausenc_polylog(x, s)
end

function clausenc(x::Acb, s::Arb)
    # If s is not a non-negative integer and 0 < real(x) < 2π call
    # _clausenc_zeta(x, s)
    if Arblib.contains_int(s) && !Arblib.isnegative(s)
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
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

"""
    clausenc(x, s, β)

Compute \$C_s^{(β)}(x)\$, that is `clausenc(x, s)` differentiated `β`
times w.r.t. `s`.

If `x` is a wide (real) ball, as determined by `iswide(x)`, it
computes a tighter enclosure by using the monotonicity properties of
the function. Currently this is only implemented for `β = 1`, `0 < x <
2π` and `s = 2` or `s = 3`. The function has a critical point at `x =
π`, one on the interval `0 < x < π` and one (due to being even around
`π`) on `π < x < 2π`. The extrema can occur either on one of these
critical points or on the endpoints of the ball. For efficiency
reasons the critical point on `0 < x < π` is precomputed for `s = 2`
and `s = 3` (the one on `π < x < 2π` is given by symmetry). In the
wide case it computes the endpoints at a reduced precision given by
```
prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
```
- **PROVE:** That there is only once critical point.
"""
function clausenc(x::Arb, s::Arb, β::Integer)
    if iszero(x) && s > 1
        isone(β) && return dzeta(s)
        s_series = ArbSeries([s, 1], degree = β)
        return zeta(s_series, one(s))[β] * factorial(β)
    end

    if iswide(x) && β == 1 && 0 < x < 2Arb(π) && (s == 2 || s == 3)
        prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
        xₗ, xᵤ = getinterval(Arb, setprecision(x, prec))
        res = union(clausenc(xₗ, s, β), clausenc(xᵤ, s, β))
        if s == 2
            critical_point =
                Arb("[1.010782703526315549251222370194235400 +/- 7.10e-37]"; prec)
        elseif s == 3
            critical_point = Arb("[1.219556773337345811161114646108970 +/- 5.13e-34]"; prec)
        end
        if Arblib.overlaps(x, critical_point) ||
           Arblib.overlaps(x, 2Arb(π) - critical_point)
            # Depending on the precision critical_point might count as
            # wide so we explicitly call _clausenc_zeta to avoid
            # infinite recursion.
            res = union(res, _clausenc_polylog(critical_point, s, β))
        end

        if Arblib.overlaps(x, Arb(π))
            res = union(res, clausenc(Arb(π; prec), s, β))
        end

        return setprecision(res, precision(x))
    end

    return _clausenc_polylog(x, s, β)
end

clausenc(x::S, s::T, β::Integer) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausenc(convert(Arb, x), convert(Arb, s), β))

"""
    clausenc(x::ArbSeries, s, β)

Compute the Taylor series with respect to `x` of \$C_s^{(β)}(x)\$,
that is `clausenc(x, s)` differentiated `β` times w.r.t. `s`.

It's computed by directly computing the Taylor coefficients by
differentiating \$C_s^{(\beta)}\$ and then composing with `x`.
"""
function clausenc(x::ArbSeries, s, β::Integer)
    res = zero(x)
    x₀ = x[0]

    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausenc(x₀, s - i, β) / factorial(i)
        else
            res[i] = -(-1)^(i ÷ 2) * clausens(x₀, s - i, β) / factorial(i)
        end
    end

    # Compose the Taylor series for the Clausian with that of the
    # input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

"""
    clausenc_expansion(x, s, M::Integer)

Compute the asymptotic expansion of `clausenc(x, s)` at zero up to
order `2M - 2`, meaning that the error term is of order `2M`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a `ArbSeries` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
arXiv:1810.10935.

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
    M > (s + 1) / 2 || throw(ArgumentError("M must be larger that (s + 1) / 2"))

    π = oftype(x, pi)

    # Non-analytic term
    if s == 2
        C = -π / 2
    else
        if iswide(s)
            C = Arb(
                ArbExtras.extrema_series(
                    s -> gamma(1 - s) * sinpi(s / 2),
                    getinterval(s)...,
                    degree = 2,
                )[1:2],
            )
        else
            C = gamma(1 - s) * sinpi(s / 2)
        end
    end
    e = s - 1

    # Analytic term
    P = ArbSeries(degree = 2M - 2, prec = precision(x))
    start = skip_constant ? 1 : 0
    for m = start:M-1
        if iswide(s)
            z = Arb(ArbExtras.extrema_series(s -> zeta(s - 2m), getinterval(s)...)[1:2])
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
        gamma(v) * inv2pi^v * cospi(v / 2) * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi)) -
        zeta(onem(v))

    if iswide(s)
        res = ArbExtras.extrema_series(f, getinterval(v)..., degree = 2)[1:2]
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
between the two terms by computing them together.

If `s > 1` it uses that the function is non-decreasing in `s`. That
this is the case can be checked by looking at the Fourier series of
the `clausenc` and the series for `zeta` and check that all terms in
their difference are non-negative.

For `s <= 1` this is no longer the case. Instead it uses Taylor
expansions to compute accurate bounds, this is implemented in
[`_clausencmzeta_zeta`](@ref). Since this method only supports `0 < x
< 2π` and `s` not overlapping a non-negative integer we fall back to
calling the methods directly if this is not satisfied.
- **IMPROVE:** Better handle this case? This case might not occur at
  all in the code though since the case `s < 1` mostly appears as
  derivatives w.r.t. `x` and then the `mzeta` part is not relevant.

Otherwise it behaves like `clausenc(x, s) - zeta(s)` with the only
difference being that it converts `x` and `s` to the same type to
begin with.
"""
function clausencmzeta(x::Arb, s::Arb)
    if iswide(s) && s > 1
        # Use that the function is non-decreasing in s to compute at
        # lower and upper bound of s
        sₗ, sᵤ = getinterval(Arb, s)
        return Arb((clausencmzeta(x, sₗ), clausencmzeta(x, sᵤ)))
    end

    if !iswide(s) ||
       (Arblib.contains_int(s) && !Arblib.isnegative(s)) ||
       !(Arblib.ispositive(x) && x < 2Arb(π))
        # If s is not wide or we are in a case which
        # _clausencmzeta_zeta doesn't support, call the clausenc and
        # zeta functions directly.
        return clausenc(x, s) - zeta(s)
    end

    if iswide(x)
        prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
        xₗ, xᵤ = getinterval(Arb, setprecision(x, prec))

        res = Arb((clausencmzeta(xᵤ, s), clausencmzeta(xₗ, s)))

        # We know that 0 < x < 2π so the only possible critical point
        # is π
        if Arblib.overlaps(x, Arb(π))
            res = union(res, clausencmzeta(Arb(π; prec), s))
        end

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
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

clausencmzeta(x::ArbSeries, s) = clausenc(x, s) - zeta(Arb(s, prec = precision(x)))

function clausencmzeta(x, s)
    x, s = promote(x, s)
    return clausenc(x, s) - zeta(s)
end

"""
    clausencmzeta(x, s, β)

Compute `clausenc(x, s, β) - clausenc(0, s, β)`. Notice that
`clausenc(0, s, β)` is given by `zeta(s)` differentiated `β` times
with respect to `s`.

This method just calls [`clausenc`](@ref) and [`zeta`](@ref) directly
and gives no direct benefit other than converting `x` and `s` to the
same type to begin with.
"""
clausencmzeta(x::Union{Arb,ArbSeries}, s, β::Integer) =
    if isone(β)
        clausenc(x, s, β) - dzeta(s)
    else
        clausenc(x, s, β) - zeta(ArbSeries([s, 1], degree = β))[β] * factorial(β)
    end

function clausencmzeta(x, s, β)
    x, s = promote(x, s)
    res = clausenc(x, s, β)
    if isone(β)
        return res - dzeta(s)
    else
        res - zeta(ArbSeries([s, 1], degree = β))[β] * factorial(β)
    end
end
