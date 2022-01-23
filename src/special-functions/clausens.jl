export clausens

"""
    _clausens_polylog(x::Arb, s::Union{Arb,Integer})
    _clausens_polylog(x::Arb, s::ArbSeries)
    _clausens_polylog(x::Arb, s::Arb, β::Integer)

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

function _clausens_polylog(x::Arb, s::ArbSeries)
    z = exp(Acb(0, x, prec = precision(x)))
    return ArbSeries(imag.(Arblib.coeffs(polylog(AcbSeries(s), z))))
end

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
enclosure using [`enclosure_series`](@ref).

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

# Handling `s` overlapping a unique positive integer
If `s` is a positive integer then `gamma(v)` diverges, if the
integer is even then `zeta(v, x / 2π) - zeta(v, 1 - x / 2π)` is zero
and if the integer is odd then `sinpi(v / 2)` is zero. To see that
`zeta(v, x / 2π) - zeta(v, 1 - x / 2π)` is zero when `s` is an even
non-negative integer, i.e. when `v` is an odd non-positive integer, we
can use formula [25.11.14](https://dlmf.nist.gov/25.11.E14) together
with [24.4.1](https://dlmf.nist.gov/24.4.E1).

For even `s` we thus want to compute an enclosure of
```
gamma(v) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
```
and for odd `s` we want to compute an enclosure of
```
gamma(v) * sinpi(v / 2)
```
In both cases we rewrite them using the reciprocal gamma function
[`rgamma`](@ref), giving us
```
(zeta(v, x / 2π) - zeta(v, 1 - x / 2π)) / rgamma(v)

sinpi(v / 2) / rgamma(v)
```
To be able to use [`fx_div_x`](@ref) we let `t = v - v_integer`, where
`v_integer` is the integer that `v` overlaps with, and write them as
```
((zeta(t + v_integer, x / 2π) - zeta(t + v_integer, 1 - x / 2π)) / t) / (rgamma(t + v_integer) / t)

(sinpi((t + v_integer) / 2) / t) / (rgamma(t + v_integer) / t)
```
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

    unique, s_integer = unique_integer(s)

    if unique && s_integer == 0
        # Enclosure of zeta_deflated(v, x / 2π) - zeta_deflated(v, 1 - x / 2π)
        z = zeta_deflated(v, xinv2pi) - zeta_deflated(v, onemxinv2pi)

        if iswide(s)
            rest = ArbExtras.enclosure_series(
                v -> gamma(v) * inv2pi^v * sinpi(v / 2),
                v,
                degree = 2,
            )
        else
            rest = gamma(v) * inv2pi^v * sinpi(v / 2)
        end

        return rest * z
    elseif unique && s_integer > 0
        v_integer = 1 - s_integer

        # Enclosure of rgamma(v) / (v - v_integer)
        rgamma_div_v = fx_div_x(t -> rgamma(t + v_integer), v - v_integer, extra_degree = 2)

        if iseven(s_integer)
            # Enclosure of (zeta(v, x / 2π) - zeta(v, 1 - x / 2π)) / (v - v_integer)
            zeta_div_v = fx_div_x(
                t -> zeta(t + v_integer, xinv2pi) - zeta(t + v_integer, onemxinv2pi),
                v - v_integer,
                extra_degree = 2,
                force = true,
            )

            # Enclosure of gamma(v) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
            gammazeta = zeta_div_v / rgamma_div_v

            if iswide(s)
                rest =
                    ArbExtras.enclosure_series(v -> inv2pi^v * sinpi(v / 2), v, degree = 2)
            else
                rest = inv2pi^v * sinpi(v / 2)
            end

            return gammazeta * rest
        else
            # Enclosure of sinpi(v / 2) / (v - v_integer)
            sin_div_v =
                fx_div_x(t -> sinpi((t + v_integer) / 2), v - v_integer, extra_degree = 2)

            # Enclosure of gamma(v) * sinpi(v / 2)
            gammasin = sin_div_v / rgamma_div_v

            if iswide(s)
                rest = ArbExtras.enclosure_series(
                    v -> inv2pi^v * (zeta(v, xinv2pi) - zeta(v, onemxinv2pi)),
                    v,
                    degree = 2,
                )
            else
                rest = inv2pi^v * (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))
            end

            return gammasin * rest
        end
    else

        f(v) =
            gamma(v) * inv2pi^v * sinpi(v / 2) * (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))

        if iswide(s)
            return ArbExtras.enclosure_series(f, v, degree = 2)
        else
            return f(v)
        end
    end
end

"""
    _clausens_zeta(x::Arb, s::ArbSeries)

Evaluation of the `clausens` function through the zeta function as a
power series in `s`.

It supports `s` overlapping non-negative integers in the same way
`_clausenc_zeta(x::Arb, s::Arb)` does, using [`fx_div_x`](@ref) to
handle the removable singularities.
"""
function _clausens_zeta(x::Arb, s::ArbSeries)
    # Handle the case when s has degree 0, so that we can assume that
    # the degree is positive from here on.
    iszero(Arblib.degree(s)) && return ArbSeries(_clausens_zeta(x, s[0]))

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

    unique, s0_integer = unique_integer(s0)

    if unique && s0_integer == 0
        # Enclosure of zeta_deflated(v, x / 2π) - zeta_deflated(v, 1 - x / 2π)
        z = zeta_deflated(v, xinv2pi) - zeta_deflated(v, onemxinv2pi)

        return gamma(v) * inv2pi^v * sinpi(v / 2) * z
    elseif unique && s0_integer > 0
        v0_integer = 1 - s0_integer

        # Enclosure of rgamma(v) / (v - v_integer)
        rgamma_div_v =
            fx_div_x(t -> rgamma(t + v0_integer), v - v0_integer, extra_degree = 2)

        if iseven(s0_integer)
            # Enclosure of (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / (v - v_integer)
            zeta_div_v = fx_div_x(
                t -> zeta(t + v0_integer, xinv2pi) - zeta(t + v0_integer, onemxinv2pi),
                v - v0_integer,
                extra_degree = 2,
                force = true,
            )

            # Enclosure of gamma(v) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
            gammazeta = zeta_div_v / rgamma_div_v

            rest = inv2pi^v * sinpi(v / 2)

            return gammazeta * rest
        else
            # Enclosure of sinpi(v / 2) / (v - v_integer)
            sin_div_v =
                fx_div_x(t -> sinpi((t + v0_integer) / 2), v - v0_integer, extra_degree = 2)

            # Enclosure of gamma(v) * sinpi(v / 2)
            gammasin = sin_div_v / rgamma_div_v

            rest = inv2pi^v * (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))

            return gammasin * rest
        end
    else
        return gamma(v) *
               inv2pi^v *
               sinpi(v / 2) *
               (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))
    end
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
indeterminate result.
- **IMPROVE:** Use that for odd negative integers it is exactly zero.

If `x` contains zero and `s > 1` it uses the asymptotic expansion at
`x = 0` from [`clausens_expansion`](@ref). This is unless `x` is wide
enough to also include `π`, in which case it uses the trivial upper
bound of the absolute value given by `zeta(s)`.

If `x` doesn't contain zero we are assured by
[`_reduce_argument_clausen`](@ref) that `0 < x < 2π`.

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
"""
function clausens(x::Arb, s::Arb)
    x, haszero, haspi, has2pi = _reduce_argument_clausen(x)

    @assert !(has2pi && !haszero)

    # Handle the special case when x contains zero
    if haszero
        if !(s > 1)
            return Arblib.indeterminate!(zero(x))
        elseif haspi
            # We could give a better bound by checking if we should
            # include use the positive or negative version. But this
            # is likely not so important.
            z = zeta(s) # Trivial upper bound
            return union(-z, z)
        elseif iszero(x)
            return zero(x)
        else
            # Compute asymptotic expansion
            M = 4
            C, e, P, E = clausens_expansion(x, s, M)

            # IMPROVE: We could compute tighter enclosures by
            # evaluating these more carefully, but this is likely not
            # needed.
            # Evaluate asymptotic expansion on [lbound(x), 0]
            res_left = let x = Arb((lbound(x), 0))
                -C * abspow(x, e) + P(x) + E * x^(2M + 1)
            end
            # Evaluate asymptotic expansion on [0, ubound(x)]
            res_right = let x = Arb((0, ubound(x)))
                C * abspow(x, e) + P(x) + E * x^(2M + 1)
            end

            return union(res_left, res_right)
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

    return _clausens_zeta(x, s)
end

# Rarely used - only naive implementation
function clausens(x::Acb, s::Union{Acb,Integer})
    @warn "Apparently this method is used!" x s maxlog = 1
    (polylog(s, exp(im * x)) - polylog(s, exp(-im * x))) / 2
end

clausens(x::S, s::T) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausens(convert(Arb, x), convert(Arb, s)))

"""
    clausens(x::ArbSeries, s)

Compute the Taylor series of the Clausen function ``S_s(x)``.

It's computed by directly computing the Taylor coefficients by
differentiating `clausens` and then composing with `x`.
"""
function clausens(x::ArbSeries, s)
    x₀ = _reduce_argument_clausen(x[0])[1]

    res = zero(x)
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

Compute the Taylor series of the Clausen function ``S_s(x)`` in the
parameter `s`.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

If `x` contains zero it returns an indeterminate result. For `s > 1`
it would be possible to compute a finite result, but this has not been
needed.
"""
function clausens(x::Arb, s::ArbSeries)
    x, haszero, _, _ = _reduce_argument_clausen(x)

    if haszero
        res = zero(s)
        for i = 0:Arblib.degree(res)
            res[0] = Arblib.indeterminate!(zero(x))
        end
        return res
    else
        return _clausens_zeta(x, s)
    end
end

"""
    clausens(x, s, β)

Compute ``S_s^{(β)}(x)``, that is `clausens(x, s)` differentiated `β`
times w.r.t. `s`.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

If `x` contains zero and we don't have `s > 1` it returns an
indeterminate result. If `s > 1` and `x`is not exactly zero it
computes an extremely naive enclosure using that an upper bound of the
absolute value is given by `zeta(s)` differentiated `β` times. If `x =
0` the result is zero.
- **IMPROVE:** Compute a tighter enclosure if `x` is not exactly zero.
  Either checking the derivative at the endpoint, similarly to
  [`clausens`](@ref), or expanding at `x = 0`.

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
"""
function clausens(x::Arb, s::Arb, β::Integer)
    x, haszero, haspi, has2pi = _reduce_argument_clausen(x)

    @assert !(has2pi && !haszero)

    # Handle the special case when x contains zero
    if haszero
        if !(s > 1)
            # Only finite for s > 1
            return Arblib.indeterminate!(zero(x))
        elseif haspi
            # Absolute value upper bounded by corresponding derivative
            # of zeta function
            r = isone(β) ? dzeta(s) : zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
            return union(-r, r)
        elseif iszero(x)
            return zero(x)
        else
            # IMPROVE: Expand at x = 0? Check derivative at abs_ubound(x)?
            # Absolute value upper bounded by corresponding derivative
            # of zeta function
            r = isone(β) ? dzeta(s) : zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
            return union(-r, r)
        end
    end

    # We can now assume that 0 < x < 2π

    if iswide(x)
        orig_prec = precision(x)
        s_f64 = Float64(s)
        min_prec = abs(s_f64 - round(s_f64)) < 1e-2 ? 64 : 32
        prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))
        x = setprecision(x, prec)

        # Compute derivative. We here call _clausenc_zeta directly to
        # avoid potentially infinite recursion. This is okay since 0 <
        # x < 2π.
        dclausens = _clausenc_zeta(x, s - 1, β)
        if Arblib.contains_zero(dclausens)
            # Use a zero order approximation
            mid = Arblib.midpoint(Arb, x)
            res = Arblib.add_error!(clausens(mid, s, β), (x - mid) * dclausens)
        else
            # Use that it's monotone
            xₗ, xᵤ = Arblib.getinterval(Arb, x)
            res = union(clausens(xₗ, s, β), clausens(xᵤ, s, β))
        end
        return setprecision(res, orig_prec)
    end

    return _clausens_zeta(x, s, β)
end

clausens(x::S, s::T, β::Integer) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausens(convert(Arb, x), convert(Arb, s), β))

"""
    clausens_expansion(x, s, M::Integer)

Compute the asymptotic expansion of `clausens(x, s)` at zero up to
order `2M - 1`, meaning that the error term is of order `2M + 1`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a `ArbSeries` `P` and the
remainder term `E`. The `M` is the same as in Lemma 2.1 in
arXiv:1810.10935.

It satisfies that
```
clausens(y, s) ∈ C * sign(y) * abs(y)^e + P(y) + E * y^(2M + 1)
```
for all `abs(y) <= abs(x)`.

Note that this method doesn't handle wide values of `s` in any special
way and doesn't support `s` overlapping positive integers, this has
not been needed.
"""
function clausens_expansion(x::Arb, s::Arb, M::Integer)
    # Non-analytic term
    C = gamma(1 - s) * cospi(s / 2)
    e = s - 1

    # Analytic terms
    P = ArbSeries(degree = 2M - 1, prec = precision(x))
    for m = 0:M-1
        P[2m+1] = (-1)^m * zeta(s - 2m - 1) / factorial(2m + 1)
    end

    # Error term
    E = clausens_expansion_remainder(x, s, M)

    return (C, e, P, E)
end

"""
    clausens_expansion_remainder(x::Arb, s::Arb, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausens(x, s)` at zero up to order `2M - 1`, meaning that the
remainder is of order `2M + 1`.

This is the `E` occurring in [`clausens_expansion`](@ref).

It requires that `abs(x) < 2π` and `2M >= s + 1`. In this case an
upper bound for the absolute value of the remainder is given by
```
2(2π)^(s - 2M) * abs(cospi(s / 2)) * zeta(2M + 2 - s) / (4π^2 - x^2)
```
and this functions returns a ball centered at zero with this radius.
"""
function clausens_expansion_remainder(x::Arb, s::Arb, M::Integer)
    pi = Arb(π)

    abs(x) < 2pi || throw(DomainError(x, "x must be less than 2π"))
    2M >= s + 1 || throw(DomainError(M, "must have 2M >= s + 1, got s = $s"))

    return Arblib.add_error!(
        zero(x),
        2(2pi)^(s - 2M) * abs(cospi(s / 2)) * zeta(2M + 2 - s) / (4pi^2 - x^2),
    )
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
    for m = M:M+10
        term = (-1)^m * zeta(s - 2m - 1) * x^(2(m - M)) / factorial(big(2m + 1))
        res += term
    end

    res[0] = clausens_expansion_remainder(x, s[0], M)

    return res
end
