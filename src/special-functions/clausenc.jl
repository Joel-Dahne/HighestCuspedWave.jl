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
    _clausenc_polylog(x::Arb, s::Union{Arb,Integer})

Evaluation of the `clausenc` function through the polylog function as
a power series in `s`.
"""
function _clausenc_polylog(x::Arb, s::ArbSeries)
    z = exp(Acb(0, x, prec = precision(x)))
    return ArbSeries(real.(Arblib.coeffs(polylog(AcbSeries(s), z))))
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

The formula is well defined as long as `s` doesn't overlap with any
non-negative integer. See further down for how those cases are
handled.

It currently only handles `0 < x < 2π`.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.

# Handling `s = 0`
If `s` is zero then `zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` diverges
but `cospi(v / 2)` goes to zero. To compute their product we use the
following approach. Let `zeta_deflated(s, a) = zeta(s, a) + 1 / (1 -
s)` be the deflated zeta function. We can then write the product as
```
cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) =
    cospi(v / 2) * (zeta_deflated(v, x / 2π) + zeta_deflated(v, 1 - x / 2π))
    - 2cospi(v / 2) / (1 - v)
```
The first term is zero since `zeta_deflated(1, a)` is finite and
`cospi(1 / 2)` is zero. For the second term we get, using L'Hoptial,
```
-2cospi(v / 2) / (1 - v) = -π * sinpi(v / 2)
```
Multiplying it with `gamma(v) * inv(2π)^v` and setting `v = 1` gives
us `-1 / 2`.

**TODO:** Handle `s` overlapping zero but not being exactly zero. We
need to be able to compute an enclosure of `cospi(v / 2) / (1 - v)`
and
```
cospi(v / 2) * (zeta_deflated(v, x / 2π) + zeta_deflated(v, 1 - x / 2π))
```
The former should be straight forward. The latter might be harder,
both factors are well defined but Arb doesn't compute `zeta_deflated`
for `s` overlapping `1` but not exactly equal to `1`.

# Handling `s` being a positive integer
If `s` is a positive integer then `gamma(v)` diverges, if the
integer is even then `cospi(v / 2)` is zero and if the integer is odd
then `zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` is zero. To see that
`zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` is zero when `s` is an odd
non-negative integer, i.e. when `v` is an even non-positive integer,
we can use formula [25.11.14](https://dlmf.nist.gov/25.11.E14)
together with [24.4.1](https://dlmf.nist.gov/24.4.E1).

For even `s` we thus get, using L'Hopital,
```
gamma(v) * cospi(v / 2) = cospi(v / 2) / rgamma(v) = -π / 2 * sinpi(v / 2) / drgamma(v)
```
where `rgamma` is the reciprocal gamma function and `drgamma(v)` its
derivative. This gives us the formula
```
clausenc(x, s) = let v = 1 - s
    -π / 2 * sinpi(v / 2) / drgamma(v) * inv(2π)^v * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
end
```
For odd `s` we instead have
```
gamma(v) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) = (dzeta(v, x / 2π) + dzeta(v, 1 - x / 2π)) / drgamma(v)
```
and we get
```
clausenc(x, s) = let v = 1 - s
    inv(2π)^v * cospi(v / 2) * (dzeta(v, x / 2π) + dzeta(v, 1 - x / 2π)) / drgamma(v)
end
```

**TODO:** Handle `s` overlapping a non-negative integer but not an
exact integer. We need to be able to compute an enclosure of
```
gamma(v) * cospi(v / 2)
```
in the even case and
```
gamma(v) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
```
in the odd case.
"""
function _clausenc_zeta(x::Arb, s::Arb)
    # For s = 0 we don't have to check 0 < x < 2π and can return
    # directly
    iszero(s) && return Arb(-1 // 2)

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

    # Handle the case when s is exactly a non-negative integer.
    if isinteger(s) && Arblib.ispositive(s)
        # Note that in this case s is never wide.

        # Compute drgamma(v)
        drgamma = rgamma(ArbSeries((v, 1)))[1]

        # Get the unique integer
        unique, s_integer = unique_integer(s)
        @assert unique

        if iseven(s_integer)
            return -Arb(π, prec = precision(x)) / 2 * sinpi(v / 2) / drgamma *
                   inv2pi^v *
                   (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))
        else
            # Compute dzeta(v, x / 2π) + dzeta(v, 1 - x / 2π)
            dzeta = let v_series = ArbSeries((v, 1))
                zeta(v_series, xinv2pi)[1] + zeta(v_series, onemxinv2pi)[1]
            end

            return inv2pi^v * cospi(v / 2) * dzeta / drgamma
        end
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
    _clausenc_zeta(x::Arb, s::ArbSeries)

Evaluation of the `clausenc` function through the zeta function as a
power series in `s`.

It currently only handles `0 < x < 2π`.

It supports non-negative integer values of `s` in a similar way as
`_clausenc_zeta(x::Arb, s::Arb)` does.

**TODO:** Handle `s` overlapping a non-negative integer but not being
exact.

# Handle `s = 0`
In this case we want to compute
```
gamma(v) * inv(2π)^v * cospi(v / 2) * (
    zeta_deflated(v, x / 2π) + zeta_deflated(v, 1 - x / 2π) - 2 / (1 - v)
)
```
The only problematic part is `-2cospi(v / 2) / (1 - v) = 2cospi(v / 2)
/ (v - 1)`. But since `v[0] = 1` we can compute this by simply
shifting the expansion of `2cospi(v / 2)

# Handling `s` being a positive integer
In this case we want to compute the series of
```
cospi(v / 2) / rgamma(v)
```
and
```
(zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / rgamma(v)
```
This we do by computing them to a degree one higher than the input and
then explicitly cancelling in the division.
"""
function _clausenc_zeta(x::Arb, s::ArbSeries)
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

    v = 1 - s

    s0 = s[0]

    if iszero(s0)
        # This implementation doesn't work if the degree of s is zero.
        # In this case we just call _clausenc_zeta(x::Arb, s::Arb)
        # directly
        if iszero(Arblib.degree(s))
            return ArbSeries(_clausenc_zeta(x, s0))
        end

        zeta_deflated(s, a) = Arblib.zeta_series!(zero(s), s, a, 1, length(s))

        # Compute 2cospi(v / 2) / (v - 1) to the same degree as v
        cosdivv = let w = ArbSeries(v, degree = Arblib.degree(v) + 1)
            2(cospi(w / 2) << 1) / ((w - 1) << 1)
        end

        return gamma(v) *
               inv2pi^v *
               (
                   cospi(v / 2) *
                   (zeta_deflated(v, xinv2pi) + zeta_deflated(v, onemxinv2pi)) + cosdivv
               )
    end

    # Handle the case when s is exactly a non-negative integer.
    if isinteger(s0) && Arblib.ispositive(s0)
        w = ArbSeries(v, degree = Arblib.degree(v) + 1)

        # Compute rgamma(v) / (v - v[0]) to the same degree as v
        denominator = rgamma(w) << 1

        # Get the unique integer
        unique, s_integer = unique_integer(s0)
        @assert unique

        if iseven(s_integer)
            # Compute cospi(v / 2) / (v - v[0]) to the same degree as v
            numerator = cospi(w / 2) << 1

            return numerator / denominator *
                   inv2pi^v *
                   (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))
        else
            # Compute (zeta(v, xinv2pi) + zeta(v, onemxinv2pi)) / (v -
            # v[0]) to the same degree as v
            numerator = let tmp = (zeta(w, xinv2pi) + zeta(w, onemxinv2pi))
                # The constant term is typically not exactly zero
                @assert Arblib.contains_zero(Arblib.ref(tmp, 0))
                tmp[0] = 0
                tmp << 1
            end

            return numerator / denominator * inv2pi^v * cospi(v / 2)
        end
    end

    return gamma(v) * inv2pi^v * cospi(v / 2) * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))
end

"""
    _clausenc_zeta(x::Arb, s::Arb, β::Integer)

Evaluation of the `clausenc(x, s, β)` function through the zeta
function.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.
**IMPROVE:** The degree used in this expansion could be tuned more. At
the moment it is set to a quite high value, but that seems to be
beneficial for `KdVZeroansatz` with a wide interval for `α`. Possibly
we want to use a higher value for higher values of `s` and also wider
values of `s`.
"""
function _clausenc_zeta(x::Arb, s::Arb, β::Integer)
    if iswide(s)
        f(s::Arb) = _clausenc_zeta(x, ArbSeries((s, 1), degree = β))[β] * factorial(β)
        f(s::ArbSeries) = begin
            @assert !iszero(Arblib.degree(s)) # Doesn't work in this case
            Arblib.derivative(
                _clausenc_zeta(x, ArbSeries(s, degree = Arblib.degree(s) + β)),
                β,
            )
        end

        res = Arb(ArbExtras.extrema_series(f, getinterval(s)..., degree = 10)[1:2])
    else
        res = _clausenc_zeta(x, ArbSeries((s, 1), degree = β))[β] * factorial(β)
    end

    return res
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
prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))
```
where `min_prec` is `32` in general but `64` if `s` is close to an
integer, determined by checking if the midpoint withing `1e-2` of an
integer, in which case higher precision is typically needed.
- **IMPROVE:** This could be tuned more, but is probably not needed.

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
        f64_s = Float64(s)
        if abs(f64_s - round(f64_s)) < 1e-2
            min_prec = 64
        else
            min_prec = 32
        end
        prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))

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

    # Handle the case when s contains a non-negative integer but it
    # not exactly an integer
    contains_nonnegative_int = s isa Arb && Arblib.contains_int(s) && !Arblib.isnegative(s)
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
        @warn "Incorrect bound for clausenc with s = $s, it assumes monotonicity" maxlog = 1
        if Arblib.radref(s) < 1e-2
            # s is not negative and not greater than 1, since the
            # radius is small it must therefore contain either 0 or 1.
            sₗ, sᵤ = getinterval(Arb, s)
            return union(clausenc(x, sₗ), clausenc(x, sᵤ))
        end
    end

    # _clausenc(x, s) is only defined for 0 < x < 2π
    if Arblib.ispositive(x) && x < 2Arb(π)
        return _clausenc_zeta(x, convert(Arb, s))
    else
        return _clausenc_polylog(x, s)
    end
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
    clausenc(x::Arb, s::ArbSeries)

Compute the Taylor series of the Clausen function \$C_s(x)\$ in the
parameter `s`.

It uses [`_clausenc_zeta`](@ref) when `0 < x < 2π and
[`_clausenc_polylog`](@ref) otherwise. It currently doesn't support
`s[0]` overlapping an integer but not being exactly an integer.

- **TODO:** Implement support for `s` overlapping integers. This will
  be needed to enclose remainder terms.
"""
function clausenc(x::Arb, s::ArbSeries)
    # _clausenc_zeta(x, s) is only defined for 0 < x < 2π
    if Arblib.ispositive(x) && x < 2Arb(π)
        return _clausenc_zeta(x, s)
    else
        return _clausenc_polylog(x, s)
    end
end

"""
    clausenc(x, s, β)

Compute \$C_s^{(β)}(x)\$, that is `clausenc(x, s)` differentiated `β`
times w.r.t. `s`.

If `x` is a wide (real) ball, as determined by `iswide(x)`, it
computes a tighter enclosure by using the monotonicity properties of
the function. Currently this is only implemented for `β = 1`, `0 < x <
2π` and `s = 0`, `s = 2` or `s = 3`.

In all of the above cases the function has a critical point at `x =
π`, since it is even around this point. For `s = 0` this is the only
critical point, whereas for `s = 2` and `s = 3` it also has one on the
interval `0 < x < π` and one (mirrored, due to being even around `π`)
on `π < x < 2π`.
- **PROVE:** That there are other critical points than those mentioned
  here.

For efficiency reasons the critical point on `0 < x < π` is
precomputed for `s = 2` and `s = 3` (the one on `π < x < 2π` is given
by symmetry).

In the wide case it computes the endpoints at a reduced precision
given by

```
prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
```
"""
function clausenc(x::Arb, s::Arb, β::Integer)
    if iszero(x) && s > 1
        isone(β) && return dzeta(s)
        s_series = ArbSeries([s, 1], degree = β)
        return zeta(s_series, one(s))[β] * factorial(β)
    end

    if iswide(x) && β == 1 && 0 < x < 2Arb(π) && (s == 0 || s == 2 || s == 3)
        prec = min(max(Arblib.rel_accuracy_bits(x) + 32, 32), precision(x))
        xₗ, xᵤ = getinterval(Arb, setprecision(x, prec))
        res = union(clausenc(xₗ, s, β), clausenc(xᵤ, s, β))
        if s == 2 || s == 3
            if s == 2
                critical_point =
                    Arb("[1.010782703526315549251222370194235400 +/- 7.10e-37]"; prec)
            elseif s == 3
                critical_point =
                    Arb("[1.219556773337345811161114646108970 +/- 5.13e-34]"; prec)
            end
            if Arblib.overlaps(x, critical_point) ||
               Arblib.overlaps(x, 2Arb(π) - critical_point)
                # Depending on the precision critical_point might count as
                # wide so we explicitly call _clausenc_zeta to avoid
                # infinite recursion.
                res = union(res, _clausenc_zeta(critical_point, s, β))
            end
        end

        if Arblib.overlaps(x, Arb(π))
            res = union(res, clausenc(Arb(π; prec), s, β))
        end

        return setprecision(res, precision(x))
    end

    # _clausenc_zeta(x, s, β) is only defined for 0 < x < 2π
    if Arblib.ispositive(x) && x < 2Arb(π)
        return _clausenc_zeta(x, s, β)
    else
        return _clausenc_polylog(x, s, β)
    end
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

- **TODO:** Handle the case when `s` overlaps an even integer. In that
  `gamma(1 - s) * sinpi(s / 2)` has a removable singularity.
- **TODO:** Handle the case when `s` overlaps an odd integer. In that
  case the non-analytic term coincides with one of the analytic
  terms. Their coefficients blow up in different directions. We
  might just not handle this case?
"""
function clausenc_expansion(x::Arb, s::Arb, M::Integer; skip_constant = false)
    Arblib.ispositive(s) || throw(ArgumentError("s must be positive, got s = $s"))
    M > (s + 1) / 2 ||
        throw(ArgumentError("M must be larger that (s + 1) / 2, got M = $M, s = $s"))

    π = oftype(x, pi)

    # Non-analytic term
    if s == 2
        C = -π / 2
    else
        contains_int, n = unique_integer(s)

        if !contains_int
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
        else
            if iseven(n)
                # FIXME: This assumes monotonicity with respect to
                # `s`. In practice this is true close to even
                # integers, but it is not always true
                @warn "Incorrect bound for leading term for Clausen with s = $s, it assumes monotonicity" maxlog =
                    1
                sₗ, sᵤ = getinterval(Arb, s)
                C = union(gamma(1 - sₗ) * sinpi(sₗ / 2), gamma(1 - sᵤ) * sinpi(sᵤ / 2))
            else
                # TODO: How should we handle this case? Just return NaN?
                C = Arblib.indeterminate!(zero(s))
            end
        end
    end
    e = s - 1

    # Analytic term
    P = ArbSeries(degree = 2M - 2, prec = precision(x))
    start = skip_constant ? 1 : 0
    for m = start:M-1
        if iswide(s)
            z = Arb(
                ArbExtras.extrema_series(s -> zeta(s - 2m), getinterval(s)..., degree = 2)[1:2],
            )
            if !isfinite(z)
                # TODO: In some cases, when s overlaps zero (but not
                # always), the above returns NaN but the evaluation
                # below works. Take a look at this.
                z = zeta(s - 2m)
            end
        else
            z = zeta(s - 2m)
        end
        P[2m] = (-1)^m * z / factorial(2m)
    end

    # Error term
    E = Arblib.add_error!(zero(x), 2(2π)^(1 + s - 2M) * zeta(2M + 1 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

"""
    clausenc_expansion_remainder(x::Arb, s::Arb, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausenc(x, s)` at zero up to order `2M - 2`, meaning that the
remainder is of order `2M`.

This is the `E` occurring in [`clausenc_expansion`](@ref).

An upper bound for the absolute value of the remainder is given by
```
2(2π)^(1 + s - 2M) * zeta(2M + 1 - s) / (4π^2 - x^2)
```
and this functions returns a ball centered at zero with this radius.
"""
clausenc_expansion_remainder(x::Arb, s::Arb, M::Integer) =
    let π = Arb(π, prec = precision(x))
        Arblib.add_error!(zero(x), 2(2π)^(1 + s - 2M) * zeta(2M + 1 - s) / (4π^2 - x^2))
    end

"""
    clausenc_expansion_remainder(x::Arb, s::ArbSeries, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausenc(x, s)` at zero up to order `2M - 2`, meaning that the
remainder is of order `2M`.

The remainder term is given by
```
x^2M * sum((-1)^m * zeta(s - 2m) * x^(2(m - M)) / factorial(2m) for m = M:Inf)
```
where we enclose the sum. We want to compute an enclosure of each term
in the expansion in `s`.

**FIXME:** Properly implement this. For now we compute a finite number
of terms in the expansion and sum them. For the constant term we do
use the rigorous enclosure.
"""
function clausenc_expansion_remainder(x::Arb, s::ArbSeries, M::Integer)
    @warn "remainder not rigorously bounded" maxlog = 1

    res = zero(s)
    for m = M:10
        term = (-1)^m * zeta(s - 2m) * x^(2(m - M)) / factorial(2m)
        res += term
    end

    res[0] = clausenc_expansion_remainder(x, s[0], M)

    return res
end

"""
    clausenc_expansion_odd_s_singular(ϵ::Arb, s::Arb, e::Arb)

For `s` overlapping an odd positive integer `2k + 1` the exponents for
the two terms `x^(s - 1)` and `x^2k` in the expansion overlap and
their coefficients blow up. This method returns `C` such that `C *
x^e` gives an enclosure for the sum of these two terms for all `x` in
the interval `[0, ϵ]`.

It requires that that `k` is at least `1` and that `0 < e < s - 1`. It
also assumes that `0 <= ϵ < π`, any negative parts of `ϵ` are ignored.

For now it only implements `s` overlapping `3`. This seems to be the
only case we actually need. Everything below assumes that `s` overlaps
with `3`, i.e. `k = 1`.

In this case the sum of the two terms we are interested in are
```
gamma(1 - s) * sinpi(s / 2) * x^(s - 1) - zeta(s - 2) / 2 * x^2
```
Expanding at `s = 3` we have
```
(-3 // 4 + log(x) / 2) * x^2 + O(s - 3)
```
Ignoring the `O(s - 3)` term for now we can factor out `x^e`, giving
us
```
x^e * ((-3 // 4 + log(x) / 2) * x^(2 - e))
```
We are therefore left computing an enclosure for
```
((-3 // 4 + log(x) / 2) * x^(2 - e))
```
on the interval `[0, ϵ]`. Splitting it into two terms we have
```
-3 // 4 * x^(2 - e)
```
Since `e < s - 1` and `s` overlaps `3` we have `2 - e > 0` so this is
zero at `x = 0` and decreasing in `x`, allowing us to compute an
enclosure. The term
```
log(x) / 2 * x^(2 - e)
```
is also zero for `x = 0`. The derivative is given by
```
x^(1 - e) / 2 + (2 - e) * log(x) / 2 * x^(1 - e) = (1 + (2 - e) * log(x)) * x^(1 - e) / 2
```
which has the unique zero
```
x = exp(1 / (e - 2))
```
so we can enclose it by evaluating it at `x = 0`, `x = ϵ` and also `x
= exp(1 / (e - 2))` if this falls in the interval `[0, ϵ]`.
- **TODO:** Handle the remainder terms from `O(s - 3)` which we ignore
  for now.
"""
function clausenc_expansion_odd_s_singular(ϵ::Arb, s::Arb, e::Arb)
    # Check requirement on e
    0 < e < s - 1 ||
        throw(ArgumentError("expected e < s - 1, got e = $e, s - 1 = $(s - 1)"))

    # Compute the integer k
    contains_int, n = unique_integer(s)

    contains_int || throw(ArgumentError("expected s overlapping integer, got s = $s"))

    isodd(n) || throw(ArgumentError("expected s overlapping odd integer, got s = $s"))

    n == 3 ||
        throw(ArgumentError("method currently only supports s overlapping 3, got s = $s"))

    iszero(ϵ) && return zero(ϵ)

    # It's enough to work with the upper bound of ϵ
    ϵ = ubound(Arb, ϵ)

    0 < ϵ < π || throw(DomainError(ϵ, "method only supports ϵ on the interval [0, π)"))

    # We assume n = 3 from here

    # Enclosure of -3 // 4 * x^(2 - e)
    term1_zero = zero(ϵ)
    term1_ϵ = -3 // 4 * ϵ^(2 - e)

    term1 = union(term1_zero, term1_ϵ)

    # Enclosure of log(x) / 2 * x^(2 - e)
    term2_zero = zero(ϵ)
    term2_ϵ = log(ϵ) / 2 * ϵ^(2 - e)

    term2 = union(term2_zero, term2_ϵ)

    critical_point = exp(1 / (e - 2))
    if !(ϵ < critical_point)
        term2_critical_point = log(critical_point) / 2 * critical_point^(e - 2)
        term2 = union(term2, term2_critical_point)
    end

    C = term1 + term2

    return C
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

clausencmzeta(x::Arb, s::ArbSeries) = clausenc(x, s) - zeta(s)

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
function clausencmzeta(x::Union{Arb,ArbSeries}, s, β::Integer)
    if isone(β)
        return clausenc(x, s, β) - dzeta(Arb(s))
    else
        return clausenc(x, s, β) - zeta(ArbSeries([s, 1], degree = β))[β] * factorial(β)
    end
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
