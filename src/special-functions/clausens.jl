export clausens

"""
    _clausens_polylog(x::Arb, s::Union{Arb,Integer})

Evaluation of the `clausens` function through the polylog function.

It uses the formula
```
clausens(x, s) = imag(polylog(s, exp(im * x)))
```
"""
function _clausens_polylog(x::Arb, s::Union{Arb,Integer})
    z = exp(Acb(0, x, prec = precision(x)))
    s = s isa Integer ? s : Acb(s, prec = precision(x))
    return imag(polylog(s, z))
end

"""
    _clausens_polylog(x::Arb, s::Union{Arb,Integer})

Evaluation of the `clausens` function through the polylog function as
a power series in `s`.
"""
function _clausens_polylog(x::Arb, s::ArbSeries)
    z = exp(Acb(0, x, prec = precision(x)))
    return ArbSeries(imag.(Arblib.coeffs(polylog(AcbSeries(s), z))))
end

"""
    _clausens_polylog(x::Arb, s::Arb, β::Integer)

Evaluation of the `clausens(x, s, β)` function through the polylog
function.

It uses the formula
```
clausens(x, s) = imag(polylog(s, exp(im * x)))
```
and computes the derivative with respect to `s` using `AcbSeries`.
"""
function _clausens_polylog(x::Arb, s::Arb, β::Integer)
    z = exp(Acb(0, x, prec = precision(x)))
    return imag(polylog(AcbSeries([s, 1], degree = β), z)[β]) * factorial(β)
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

The formula is well defined as long as `s` doesn't overlap with any
non-negative integer. See further down for how those cases are
handled.

It currently only handles `0 < x < 2π`.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.

# Handling `s = 0`
If `s` is zero then both `zeta(v, x / 2π)` and `zeta(v, 1 - x / 2π)`
blow up, however their difference is bounded. To see this let
`zeta_deflated(s, a) = zeta(s, a) + 1 / (1 - s)` be the deflated zeta
function. Looking at the difference the term `1 / (1 - v)` is
cancelled and we get
```
zeta(v, x / 2π) - zeta(v, 1 - x / 2π) =
    zeta_deflated(v, x / 2π) - zeta_deflated(v, 1 - x / 2π)
```

We can also note that for `s = 0` we have
```
gamma(v) * inv(2π)^v * sinpi(v / 2) = inv(2π)
```

**TODO:** Handle `s` overlapping zero but not being exactly zero. We
need to be able to compute an enclosure of
```
zeta_deflated(v, x / 2π) - zeta_deflated(v, 1 - x / 2π)
```
Both terms are well defined but Arb doesn't compute `zeta_deflated`
for `s` overlapping `1` but not exactly equal to `1`.

# Handling `s` being a positive integer
If `s` is a positive integer then `gamma(v)` diverges, if the
integer is even then `zeta(v, x / 2π) - zeta(v, 1 - x / 2π)` is zero
and if the integer is odd then `sinpi(v / 2)` is zero. To see that
`zeta(v, x / 2π) - zeta(v, 1 - x / 2π)` is zero when `s` is an even
non-negative integer, i.e. when `v` is an odd non-positive integer, we
can use formula [25.11.14](https://dlmf.nist.gov/25.11.E14) together
with [24.4.1](https://dlmf.nist.gov/24.4.E1).

For even `s` we thus get, using L'Hopital,
```
gamma(v) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π)) =
    (zeta(v, x / 2π) - zeta(v, 1 - x / 2π)) / rgamma(v)=
    (dzeta(v, x / 2π) - dzeta(v, 1 - x / 2π)) / drgamma(v)
```
where `rgamma` is the reciprocal gamma function and `drgamma(v)` its
derivative. This gives us the formula
```
clausens(x, s) = let v = 1 - s
    inv(2π)^v * sinpi(v / 2) * (dzeta(v, x / 2π) - dzeta(v, 1 - x / 2π)) / drgamma(v)
end
```
For odd `s` we instead have
```
gamma(v) * sinpi(v / 2) = π / 2 * cospi(v / 2) / drgamma(v)
```
and we get
```
clausens(x, s) = let v = 1 - s
    π / 2 * cospi(v / 2) / drgamma(v) * inv(2π)^v * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
end
```

**TODO:** Handle `s` overlapping a non-negative integer but not an
exact integer. We need to be able to compute an enclosure of
```
gamma(v) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
```
in the even case and
```
gamma(v) * sinpi(v / 2)
```
in the odd case.
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

    if iszero(s)
        # Function for enclosing zeta_deflated(s, a)
        zeta_deflated(s, a) = Arblib.zeta_series!(ArbSeries(0), ArbSeries(v), a, 1, 1)[0]

        return inv2pi * (zeta_deflated(v, xinv2pi) - zeta_deflated(v, onemxinv2pi))
    end

    if isinteger(s) && Arblib.ispositive(s)
        # Handle the case when s is exactly a non-negative integer.
        # Note that in this case s is never wide.

        # Compute drgamma(v)
        drgamma = rgamma(ArbSeries((v, 1)))[1]

        # Get the unique integer
        unique, s_integer = unique_integer(s)
        @assert unique

        if iseven(s_integer)
            # Compute dzeta(v, x / 2π) + dzeta(v, 1 - x / 2π)
            dzeta = let v_series = ArbSeries((v, 1))
                zeta(v_series, xinv2pi)[1] - zeta(v_series, onemxinv2pi)[1]
            end

            return inv2pi^v * sinpi(v / 2) * dzeta / drgamma
        else
            return Arb(π, prec = precision(x)) / 2 * cospi(v / 2) / drgamma *
                   inv2pi^v *
                   (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))
        end
    end

    f(v) = gamma(v) * inv2pi^v * sinpi(v / 2) * (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))

    if iswide(s)
        return ArbExtras.enclosure_series(f, v, degree = 2)
    end

    return f(v)
end

"""
    _clausens_zeta(x::Arb, s::ArbSeries)

Evaluation of the `clausens` function through the zeta function as a
power series in `s`.

It currently only handles `0 < x < 2π`.

It supports non-negative integer values of `s` in a similar way as
`_clausens_zeta(x::Arb, s::Arb)` does.

**TODO:** Handle `s` overlapping a non-negative integer but not being
exact.

# Handle `s = 0`
In this case we want to compute
```
gamma(v) * inv(2π)^v * sinpi(v / 2) * (
    zeta_deflated(v, x / 2π) - zeta_deflated(v, 1 - x / 2π)
)
```
for which all factors are well defined.

# Handling `s` being a positive integer
In this case we want to compute the series of
```
(zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / rgamma(v)
```
and
```
sinpi(v / 2) / rgamma(v)
```
This we do by computing them to a degree one higher than the input and
then explicitly cancelling in the division.
"""
function _clausens_zeta(x::Arb, s::ArbSeries)
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
        zeta_deflated(s, a) = Arblib.zeta_series!(zero(s), s, a, 1, length(s))

        return gamma(v) *
               inv2pi^v *
               sinpi(v / 2) *
               (zeta_deflated(v, xinv2pi) - zeta_deflated(v, onemxinv2pi))
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
            # Compute (zeta(v, xinv2pi) - zeta(v, onemxinv2pi)) / (v -
            # v[0]) to the same degree as v
            numerator = let tmp = (zeta(w, xinv2pi) - zeta(w, onemxinv2pi))
                # The constant term is typically not exactly zero
                @assert Arblib.contains_zero(Arblib.ref(tmp, 0))
                tmp[0] = 0
                tmp << 1
            end

            return numerator / denominator * inv2pi^v * sinpi(v / 2)
        else
            # Compute sinpi(v / 2) / (v - v[0]) to the same degree as v
            numerator = sinpi(w / 2) << 1

            return numerator / denominator *
                   inv2pi^v *
                   (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))
        end
    end

    return gamma(v) * inv2pi^v * sinpi(v / 2) * (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))
end

"""
    _clausens_zeta(x::Arb, s::Arb, β::Integer)

Evaluation of the `clausens(x, s, β)` function through the zeta
function.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.
**IMPROVE:** The degree used in this expansion could be tuned more. At
the moment it is set to a quite high value, but that seems to be
beneficial for `KdVZeroansatz` with a wide interval for `α`. Possibly
we want to use a higher value for higher values of `s` and also wider
values of `s`.
"""
function _clausens_zeta(x::Arb, s::Arb, β::Integer)
    if iswide(s)
        f(s::Arb) = _clausens_zeta(x, ArbSeries((s, 1), degree = β))[β] * factorial(β)
        f(s::ArbSeries) = begin
            @assert !iszero(Arblib.degree(s)) # Doesn't work in this case
            Arblib.derivative(
                _clausens_zeta(x, ArbSeries(s, degree = Arblib.degree(s) + β)),
                β,
            )
        end

        res = ArbExtras.enclosure_series(f, s, degree = 10)
    else
        res = _clausens_zeta(x, ArbSeries((s, 1), degree = β))[β] * factorial(β)
    end

    return res
end

"""
    clausens(x, s)

Compute the Clausen function ``S_s(x)``.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

If `x` contains zero and we don't have `s > 1` it returns an
indeterminate result, except in the case that `s` is exactly an odd
negative integer in which case the result is exactly zero.

If `x` contains zero and `s > 1` it looks at the lower and upper bound
of `x` to check if `x` includes the critical points on ``(0, π)`` or
``(-π, 0)``. This is done by checking that the points have not passed
`-π` or `π` respectively and that derivative in `x` at the points is
positive, in this case it doens't include the critical point.
- **PROVE:** That the above is correct. More precisely it is enough to
    prove that there is a unique critical point on ``(0, π)``.
If it doesn't include any of the critical points it evaluates at the
endpoints of `x`. Otherwise it uses the trivial upper bound of the
absolute value given by `zeta(s)`.

If `x` is a wide ball (not containing zero), as determined by
`iswide(x)`, it computes a tighter enclosure by first checking if the
derivative doesn't contains zero, if not it uses monotonicity to only
evaluate at endpoints. If the derivative does contain zero it uses a
zeroth order approximation instead. In the wide case it computes the
endpoints at a reduced precision given by
```
prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))
```
where `min_prec` is `32` in general but `64` if `s` is close to an
integer, determined by checking if the midpoint withing `1e-2` of an
integer, in which case higher precision is typically needed.
- **IMPROVE:** This could be tuned more, but is probably not needed.

The case when `s` is a wide ball is in general handled by the
underlying method [`_clausens_zeta`](@ref). The exception is when `s`
overlaps with a non-negative integer, in which case it gives an
indeterminate result. In that case we compute at the midpoint of `s`
and bound the error by using a bound for the derivative in `s`. For s
> 1 the derivative in `s` is bounded by `dzeta(s)`, this can be seen
by looking at the Fourier series for `clausens(x, s, 1)` and noticing
that it is always bounded by `dzeta(s)`.
- **TODO:** Figure out how to bound this for `s = 1` and `s = 0`. In
  this case the derivative in `s` blows up at `x = 0` so we can't use
  a uniform bound. For now we compute a bound assuming monotonicity in
  `s`, which is not true.
- **TODO:** Move the logic for handling `s` overlapping integers to
  [`_clausenc_zeta`](@ref).
"""
function clausens(x::Arb, s::Arb)
    x, haszero, haspi, has2pi = _reduce_argument_clausen(x)

    @assert !(has2pi && !haszero)

    # Handle the special case when x contains zero
    if haszero
        if !(s > 1)
            # Identically equal to zero for odd negative integers,
            # otherwise not finite
            if isinteger(s) && isodd(unique_integer(s)[2])
                return zero(x)
            else
                return Arblib.indeterminate!(zero(x))
            end
        elseif haspi
            # We could give a better bound by checking if we should
            # include use the positive or negative version. But this
            # is likely not so important.
            z = zeta(s) # Trivial upper bound
            return union(-z, z)
        elseif iszero(x)
            return zero(x)
        else
            xₗ, xᵤ = getinterval(Arb, x)
            # Check for inclusion of critical point on (-π, 0)
            if iszero(xₗ)
                res_left = zero(x)
            elseif Arblib.ispositive(clausenc(xₗ, s - 1))
                res_left = clausens(xₗ, s)
            else
                res_left = -zeta(s) # Trivial lower bound
            end

            # Check for inclusion of critical point on (0, π)
            if iszero(xᵤ)
                res_right = zero(x)
            elseif Arblib.ispositive(clausenc(xᵤ, s - 1))
                res_right = clausens(xᵤ, s)
            else
                res_right = zeta(s) # Trivial upper bound
            end

            return Arb((res_left, res_right))
        end
    end

    # We can now assume that 0 < x < 2π

    if iswide(x)
        orig_prec = precision(x)
        s_f64 = Float64(s)
        min_prec = abs(s_f64 - round(s_f64)) < 1e-2 ? 64 : 32
        prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))
        x = setprecision(x, prec)

        # Compute derivative
        dclausens = clausenc(x, s - 1)
        if Arblib.contains_zero(dclausens)
            # Use a zero order approximation
            mid = Arblib.midpoint(Arb, x)
            res = Arblib.add_error!(clausens(mid, s), (x - mid) * dclausens)
        else
            # Use that it's monotone
            xₗ, xᵤ = Arblib.getinterval(Arb, x)
            res = union(clausens(xₗ, s), clausens(xᵤ, s))
        end
        return setprecision(res, orig_prec)
    end

    # Handle the case when s contains a non-negative integer but it
    # not exactly an integer
    # TODO: This should optimally be handled internally by _clausenc_zeta
    contains_nonnegative_int = Arblib.contains_int(s) && !Arblib.isnegative(s)
    if contains_nonnegative_int && !iszero(Arblib.radref(s))
        if s > 1
            # Evaluate at midpoint
            smid = Arblib.midpoint(Arb, s)
            res = clausens(x, smid)

            # Bound derivative in s using derivative of zeta function
            derivative = dzeta(s)
            error = (s - smid) * abs(derivative)

            return res + error
        end

        # FIXME: For now we bound it assuming monotonicity in s, this
        # is not true in practice. We only do this for s close to 0
        # and 1
        @warn "Incorrect bound for clausens with s = $s, it assumes monotonicity" maxlog =
            100
        if Arblib.radref(s) < 1e-2
            # s is not negative and not greater than 1, since the
            # radius is small it must therefore contain either 0 or 1.
            sₗ, sᵤ = getinterval(Arb, s)
            return union(clausens(x, sₗ), clausens(x, sᵤ))
        end
    end

    return _clausens_zeta(x, s)
end

clausens(x::Acb, s) = (polylog(s, exp(im * x)) - polylog(s, exp(-im * x))) / 2

clausens(x::S, s::T) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausens(convert(Arb, x), convert(Arb, s)))

"""
    clausens(x::ArbSeries, s)

Compute the Taylor series of the Clausen function \$S_s(x)\$.

It's computed by directly computing the Taylor coefficients by
differentiating `clausens` and then composing with `x`.
"""
function clausens(x::ArbSeries, s)
    res = zero(x)
    x₀ = x[0]

    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausens(x₀, s - i) / factorial(i)
        else
            res[i] = (-1)^(i ÷ 2) * clausenc(x₀, s - i) / factorial(i)
        end
    end

    # Compose the Taylor series for the result with that of the input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

"""
    clausens(x::Arb, s::ArbSeries)

Compute the Taylor series of the Clausen function \$S_s(x)\$ in the
parameter `s`.

It uses [`_clausens_zeta`](@ref) when `0 < x < 2π and
[`_clausens_polylog`](@ref) otherwise. It currently doesn't support
`s[0]` overlapping an integer but not being exactly an integer.

- **TODO:** Implement support for `s` overlapping integers. This will
  be needed to enclose remainder terms.
"""
function clausens(x::Arb, s::ArbSeries)
    # _clausenc_zeta(x, s) is only defined for 0 < x < 2π
    if Arblib.ispositive(x) && x < 2Arb(π)
        return _clausens_zeta(x, s)
    else
        return _clausens_polylog(x, s)
    end
end

"""
    clausens(x, s, β)

Compute \$S_s^{(β)}(x)\$, that is `clausens(x, s)` differentiated `β`
times w.r.t. `s`.

- **IMPROVE**: Handle wide (real) balls better, similar to how
  `clausens(x, s)` does it.
"""
function clausens(x::Arb, s::Arb, β::Integer)
    # _clausenc_zeta(x, s, β) is only defined for 0 < x < 2π
    if Arblib.ispositive(x) && x < 2Arb(π)
        return _clausens_zeta(x, s, β)
    else
        return _clausens_polylog(x, s, β)
    end
end

clausens(x::S, s::T, β::Integer) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausens(convert(Arb, x), convert(Arb, s), β))

"""
    clausens_expansion(x, s, M::Integer)

Compute the asymptotic expansion of `clausens(x, s)` at zero up to
order `2M - 1`, meaning that the error term is of order `2M + 1`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a `ArbSeries` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
arXiv:1810.10935.

It satisfies that `clausens(y, s) ∈ C*abs(y)^e + P(y) + E*y^(2M + 1)` for
all `|y| <= |x|`.

Note that this method doesn't handle wide values of `s` in any special
way. This has not been needed anywhere so far.
"""
function clausens_expansion(x::Arb, s::Arb, M::Integer)
    Arblib.ispositive(s) || throw(ArgumentError("s must be positive, got s = $s"))
    M > (s + 1) / 2 ||
        throw(ArgumentError("M must be larger that (s + 1) / 2, got M = $M, s = $s"))

    π = oftype(x, pi)

    # Non-analytic term
    C = gamma(1 - s) * cospi(s / 2)
    e = s - 1

    # Analytic term
    P = ArbSeries(degree = 2M - 1, prec = precision(x))
    for m = 0:M-1
        P[2m+1] = (-1)^m * zeta(s - 2m - 1) / factorial(2m + 1)
    end

    # Error term
    E = Arblib.add_error!(zero(x), 2(2π)^(s - 2M) * zeta(2M + 2 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

"""
    clausens_expansion_remainder(x::Arb, s::Arb, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausens(x, s)` at zero up to order `2M - 1`, meaning that the
remainder is of order `2M + 1`.

This is the `E` occurring in [`clausens_expansion`](@ref).

An upper bound for the absolute value of the remainder is given by
```
2(2π)^(s - 2M) * zeta(2M + 2 - s) / (4π^2 - x^2)
```
and this functions returns a ball centered at zero with this radius.
"""
clausens_expansion_remainder(x::Arb, s::Arb, M::Integer) =
    let π = Arb(π, prec = precision(x))
        Arblib.add_error!(zero(x), 2(2π)^(s - 2M) * zeta(2M + 2 - s) / (4π^2 - x^2))
    end

"""
    clausens_expansion_remainder(x::Arb, s::ArbSeries, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausens(x, s)` at zero up to order `2M - 1`, meaning that the
remainder is of order `2M + 1`.

The remainder term is given by
```
x^(2M + 1) * sum((-1)^m * zeta(s - 2m - 1) * x^(2(m - M)) / factorial(2m + 1) for m = M:Inf)
```
where we enclose the sum. We want to compute an enclosure of each term
in the expansion in `s`.

**FIXME:** Properly implement this. For now we compute a finite number
of terms in the expansion and sum them. For the constant term we do
use the rigorous enclosure.
"""
function clausens_expansion_remainder(x::Arb, s::ArbSeries, M::Integer)
    @warn "remainder not rigorously bounded" maxlog = 1

    res = zero(s)
    for m = M:9
        term = (-1)^m * zeta(s - 2m - 1) * x^(2(m - M)) / factorial(2m + 1)
        res += term
    end

    res[0] = clausens_expansion_remainder(x, s[0], M)

    return res
end
