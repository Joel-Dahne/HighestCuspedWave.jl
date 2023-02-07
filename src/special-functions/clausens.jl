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
    _clausens_zeta_f(v::Arb, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)
    _clausens_zeta_f(v::ArbSeries, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)

Compute
```
gamma(v) * inv(2π)^v * sinpi(v / 2) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
```
It assumes that the arguments after `v` are set to `2π`, `x / 2π` and
`1 - x / 2π` respectively. This function is used internally by
[`_clausens_zeta`](@ref).
"""
function _clausens_zeta_f(v::Arb, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)
    res = zeta(v, xinv2pi)
    tmp = zeta(v, onemxinv2pi)
    Arblib.sub!(res, res, tmp)

    Arblib.gamma!(tmp, v)
    Arblib.mul!(res, res, tmp)

    Arblib.pow!(tmp, inv2pi, v)
    Arblib.mul!(res, res, tmp)

    Arblib.mul_2exp!(tmp, v, -1)
    Arblib.sin_pi!(tmp, tmp)
    Arblib.mul!(res, res, tmp)

    return res
end

function _clausens_zeta_f(v::ArbSeries, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)
    res = zeta(v, xinv2pi)
    tmp = zeta(v, onemxinv2pi)
    Arblib.sub!(res, res, tmp)

    Arblib.gamma_series!(tmp, v, length(tmp))
    Arblib.mullow!(res, res, tmp, length(res))

    Arblib.zero!(tmp)
    tmp[0] = inv2pi
    Arblib.pow_series!(tmp, tmp, v, length(tmp))
    Arblib.mullow!(res, res, tmp, length(res))

    Arblib.mul_2exp!(tmp, v, -1)
    Arblib.sin_pi_series!(tmp, tmp, length(tmp))
    Arblib.mullow!(res, res, tmp, length(res))

    return res
end

"""
    _clausens_zeta_f_one(v::Union{Arb,ArbSeries}, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)

Compute
```
gamma(v) * inv(2π)^v * sinpi(v / 2) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
```
in a way that accounts for the removable singularity at `v = 1`.

If `v = 1` then both `zeta(v, x / 2π)` and `zeta(v, 1 - x / 2π)`
blow up, however their difference is bounded. To see this let
`zeta_deflated(s, a) = zeta(s, a) + 1 / (1 - s)` be the deflated zeta
function. Looking at the difference the term `1 / (1 - v)` is
cancelled and we get
```
zeta(v, x / 2π) - zeta(v, 1 - x / 2π) =
    zeta_deflated(v, x / 2π) - zeta_deflated(v, 1 - x / 2π)
```

It assumes that the arguments after `v` are set to `inv(2π)`, `x / 2π`
and `1 - x / 2π` respectively. This function is used internally by
[`_clausens_zeta`](@ref).
"""
_clausens_zeta_f_one(v::Union{Arb,ArbSeries}, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb) =
    gamma(v) *
    inv2pi^v *
    sinpi(v / 2) *
    (zeta_deflated(v, xinv2pi) - zeta_deflated(v, onemxinv2pi))

"""
    _clausens_zeta_f_even(v::Union{Arb,ArbSeries}, v0::Integer, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)

Compute
```
gamma(v) * inv(2π)^v * sinpi(v / 2) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
```
in a way that accounts for the removable singularity at even `v <= 0`.
The argument `v0` should be set to the integer where the the removable
singularity is located.

In this case `gamma(v)` diverges but `cospi(v / 2)` is zero. We
rewrite it as
```
(sinpi(v / 2) / (v - v0)) / (rgamma(v) / (v - v0))
```
and use [`fx_div_x`](@ref) to compute the factors.

It assumes that the arguments after `v0` are set to `2π`, `x / 2π` and
`1 - x / 2π` respectively. This function is used internally by
[`_clausens_zeta`](@ref).
"""
function _clausens_zeta_f_even(
    v::Union{Arb,ArbSeries},
    v0::Integer,
    inv2pi::Arb,
    xinv2pi::Arb,
    onemxinv2pi::Arb,
)
    # Enclosure of rgamma(v) / (v - v_integer)
    rgamma_div_v = fx_div_x(t -> rgamma(t + v0), v - v0, extra_degree = 2)

    # Enclosure of sinpi(v / 2) / (v - v_integer)
    sin_div_v = fx_div_x(t -> sinpi((t + v0) / 2), v - v0, extra_degree = 2)

    return inv2pi^v * sin_div_v / rgamma_div_v * (zeta(v, xinv2pi) - zeta(v, onemxinv2pi))
end

"""
    _clausens_zeta_f_odd(v::Union{Arb,ArbSeries}, v0::Integer, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)

Compute
```
gamma(v) * inv(2π)^v * sinpi(v / 2) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
```
in a way that accounts for the removable singularity at odd `v <= 0`.
The argument `v0` should be set to the integer where the the removable
singularity is located.

In this case `gamma(v)` diverges but `zeta(v, x / 2π) - zeta(v, 1 - x
/ 2π)` is zero. To see that `zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` is
zero when `v` is an odd non-positive integer we can use formula
[25.11.14](https://dlmf.nist.gov/25.11.E14) together with
[24.4.1](https://dlmf.nist.gov/24.4.E1). We rewrite it as
```
((zeta(v, x / 2π) - zeta(v, 1 - x / 2π)) / (v - v0)) / (rgamma(v) / (v - v0))
```
and use [`fx_div_x`](@ref) to compute the factors.

It assumes that the arguments after `v0` are set to `2π`, `x / 2π` and
`1 - x / 2π` respectively. This function is used internally by
[`_clausens_zeta`](@ref).
"""
function _clausens_zeta_f_odd(
    v::Union{Arb,ArbSeries},
    v0::Integer,
    inv2pi::Arb,
    xinv2pi::Arb,
    onemxinv2pi::Arb,
)
    # Enclosure of rgamma(v) / (v - v0)
    rgamma_div_v = fx_div_x(t -> rgamma(t + v0), v - v0, extra_degree = 2)

    # Enclosure of (zeta(v, x / 2π) - zeta(v, 1 - x / 2π)) / (v - v0)
    zeta_div_v = fx_div_x(
        t -> zeta(t + v0, xinv2pi) - zeta(t + v0, onemxinv2pi),
        v - v0,
        extra_degree = 2,
        force = true,
    )

    return inv2pi^v * sinpi(v / 2) * zeta_div_v / rgamma_div_v
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

The formula is well defined as long as `v` doesn't overlap with any
integer `<= 1` (meaning `s` doesn't overlap with any non-negative
integer). For integer `v <= 1` there is a removable singularity that
has to be handled. This is done using the functions
[`_clausens_zeta_f_one`](@ref), [`_clausens_zeta_f_even`](@ref) and
[`_clausens_zeta_f_odd`](@ref), handling the cases when `v` is one and
`v <= 0` and even or odd respectively.

It currently only handles `0 < x < 2π`. If `s` is wide, as determined
by `iswide(s)` it computes a tighter enclosure of the coefficients
using a Taylor expansion in `s`.
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

    # Function for computing
    # gamma(v) * inv(2π)^v * sinpi(v / 2) * (zeta(v, x / 2π) - zeta(v, 1 - x / 2π))
    # It takes care of picking the right function to use near the
    # removable singularities.
    f(v) = begin
        v0 = v isa Arb ? v : Arblib.ref(v, 0)
        rounded_mid_v0 = round(Int, Float64(v0))

        if rounded_mid_v0 <= 1 && is_approx_integer(v0, tol = 0.001)
            # v0 is never 0 in this case so no need to update degree
            # in case v::ArbSeries
            Arblib.union!(v0, v0, Arb(rounded_mid_v0))
        end

        unique, v0_integer = unique_integer(v0)

        if unique && v0_integer == 1
            _clausens_zeta_f_one(v, inv2pi, xinv2pi, onemxinv2pi)
        elseif unique && iseven(v0_integer) && v0_integer <= 0
            _clausens_zeta_f_even(v, v0_integer, inv2pi, xinv2pi, onemxinv2pi)
        elseif unique && isodd(v0_integer) && v0_integer <= 0
            _clausens_zeta_f_odd(v, v0_integer, inv2pi, xinv2pi, onemxinv2pi)
        else
            _clausens_zeta_f(v, inv2pi, xinv2pi, onemxinv2pi)
        end
    end

    if iswide(s)
        return ArbExtras.enclosure_series(f, v, degree = 2)
    else
        return f(v)
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
        res = ArbExtras.enclosure_series(
            ArbExtras.derivative_function(s -> _clausens_zeta(x, s), β),
            s,
            degree = 10,
        )
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

If `x` contains zero and `s > 1` it first checks if `x` contains `π`,
in which case it returns a trivial bound. Otherwise it checks if it is
increasing by checking if `clausenc(abs_ubound(x), s - 1) > 0`, this
uses the fact that `clausenc(x, s - 1)` is decreasing on ``[0, π]``.
If it's not monotone it returns a trivial bound.

If `x` doesn't contain zero we are assured by
[`_reduce_argument_clausen`](@ref) that `0 < x < 2π`.

If `x` is a wide ball (not containing zero), as determined by
`iswide(x)`, it computes a tighter enclosure by first checking if the
derivative doesn't contains zero, if not it uses monotonicity to only
evaluate at endpoints. If the derivative does contain zero it uses a
zeroth order approximation instead.

It accepts the keyword argument `deriv_x` which can be set to a
precomputed value of the derivative, more precisely it should be an
enclosure of `clausenc(x, s - 1)`. If deriv_x` is not set it computes
the derivative automatically. This argument is useful when computing
expansions in `x` of Clausen functions since the derivative then
corresponds to the next coefficient in the expansion and we can avoid
computing it multiple times this way.
"""
function clausens(x::Arb, s::Arb; deriv_x::Union{Nothing,Arb} = nothing)
    x, haszero, haspi, has2pi = _reduce_argument_clausen(x)

    @assert !(has2pi && !haszero)

    prec, original_prec = _choose_precision_clausen(x, s)
    if prec != original_prec
        x = setprecision(x, prec)
        s = setprecision(s, prec)
    end

    if haszero # Handle the special case when x contains zero
        if !(s > 1)
            res = indeterminate(x)
        elseif haspi
            # We could give a better bound by checking if we should
            # include use the positive or negative version. But this
            # is likely not so important.
            z = zeta(s) # Trivial upper bound
            res = union(-z, z)
        elseif iszero(x)
            res = zero(x)
        elseif Arblib.ispositive(clausenc(abs_ubound(Arb, x), s - 1))
            # Monotone on the interval
            xₗ, xᵤ = ArbExtras.enclosure_getinterval(x)
            res = Arb((-clausens(-xₗ, s), clausens(xᵤ, s)))
        else
            z = zeta(s) # Trivial upper bound
            res = union(-z, z)
        end
    elseif iswide(x) # We can now assume that 0 < x < 2π
        # Compute derivative
        if isnothing(deriv_x)
            deriv_x = clausenc(x, s - 1)
        end
        if Arblib.contains_zero(deriv_x)
            # Use a zero order approximation
            mid = midpoint(Arb, x)
            res = Arblib.add_error!(_clausens_zeta(mid, s), (x - mid) * deriv_x)
        else
            # Use that it's monotone
            xₗ, xᵤ = ArbExtras.enclosure_getinterval(x)
            res = _clausens_zeta(xₗ, s)
            Arblib.union!(res, res, _clausens_zeta(xᵤ, s))
        end
    else
        res = _clausens_zeta(x, s)
    end

    if prec != original_prec
        res = setprecision(res, original_prec)
    end

    return res
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
    s = convert(Arb, s)

    x₀ = _reduce_argument_clausen(x[0])[1]

    res = zero(x)

    # Precompute clausenc functions for i = 1:2:Arblib.degree(x)+1.
    # They are used both for the i-th coefficients and also in the
    # computation of clausens for the (i-1)-th coefficient.
    clausencs = [clausenc(x₀, s - i) for i = 1:2:Arblib.degree(x)+1]
    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] =
                (-1)^(i ÷ 2) * clausens(x₀, s - i, deriv_x = clausencs[i÷2+1]) /
                factorial(i)
        else
            res[i] = (-1)^(i ÷ 2) * clausencs[(i+1)÷2] / factorial(i)
        end
    end

    return ArbExtras.compose_zero(res, x)
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
            res[0] = indeterminate(x)
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
zeroth order approximation instead.
"""
function clausens(x::Arb, s::Arb, β::Integer)
    x, haszero, haspi, has2pi = _reduce_argument_clausen(x)

    @assert !(has2pi && !haszero)

    prec, original_prec = _choose_precision_clausen(x, s)
    if prec != original_prec
        x = setprecision(x, prec)
        s = setprecision(s, prec)
    end

    if haszero # Handle the special case when x contains zero
        if !(s > 1)
            # Only finite for s > 1
            res = indeterminate(x)
        elseif haspi
            # Absolute value upper bounded by corresponding derivative
            # of zeta function
            r = isone(β) ? dzeta(s) : zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
            res = union(-r, r)
        elseif iszero(x)
            res = zero(x)
        else
            # IMPROVE: Expand at x = 0? Check derivative at abs_ubound(x)?
            # Absolute value upper bounded by corresponding derivative
            # of zeta function
            r = isone(β) ? dzeta(s) : zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
            res = union(-r, r)
        end
    elseif iswide(x) # We can now assume that 0 < x < 2π
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
            xₗ, xᵤ = ArbExtras.enclosure_getinterval(x)
            res = union(clausens(xₗ, s, β), clausens(xᵤ, s, β))
        end
    else
        res = _clausens_zeta(x, s, β)
    end

    if prec != original_prec
        res = setprecision(res, original_prec)
    end

    return res
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

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure of the coefficients using a Taylor expansion in `s`.
"""
function clausens_expansion(x::Arb, s::Arb, M::Integer)
    # Non-analytic term
    if iswide(s)
        C = ArbExtras.enclosure_series(s -> gamma(1 - s) * cospi(s / 2), s, degree = 2)
    else
        C = gamma(1 - s) * cospi(s / 2)
    end
    e = s - 1

    # Analytic terms
    P = ArbSeries(degree = 2M - 1, prec = precision(x))
    for m = 0:M-1
        if iswide(s)
            z = ArbExtras.enclosure_series(s -> zeta(s - 2m - 1), s, degree = 2)

            if !isfinite(z)
                # In some cases, when s overlaps zero (but not
                # always), the above returns NaN but the evaluation
                # below works.
                z = zeta(s - 2m - 1)
            end
        else
            z = zeta(s - 2m - 1)
        end
        P[2m+1] = (-1)^m * z / factorial(2m + 1)
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

This is the `E` occurring in [`clausens_expansion`](@ref) and is given
by
```
sum((-1)^m * zeta(s - 2m - 1) * x^(2m + 1) / factorial(2m) for m = M:Inf) / x^(2M + 1)
```

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
    clausens_expansion_remainder(x::Arb, s::Arb, β::Integer, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausens(x, s, β)` at zero up to order `2M - 1`, meaning that the
remainder is of order `2M + 1`.

This is the tail in the expansion at `x = 0` and is given by
```
sum((-1)^m * dzeta(s - 2m - 1, β) * x^(2m + 1) / factorial(2m + 1) for m = M:Inf) / x^(2M + 1)
```
where we by `dzeta(s - 2m, β)` mean the zeta-function differentiated
`β` times.

It requires that `abs(x) < 2π` and `2M >= s + 1`.

The upper bound of the absolute value of the remainder is given by a
somewhat awkward expression involving a multitude of special
functions. See the paper for details.

This functions returns a ball centered at zero with the upper bound as
radius.
"""
function clausens_expansion_remainder(x::Arb, s::Arb, β::Integer, M::Integer)
    β == 0 && return clausens_expansion_remainder(x, s, M)

    abs(x) < 2pi || throw(DomainError(x, "x must be less than 2π"))
    2M >= s + 1 || throw(DomainError(M, "must have 2M >= s + 1, got s = $s"))

    twopi = 2Arb(π)

    # Function for computing multinomial
    # Denominator is smaller than numerator, so overflow
    # would always be caught by factorial(β)
    multinomial(β, j1, j2, j3) =
        Arb(factorial(β) // (factorial(j1) * factorial(j2) * factorial(j3)))

    # Upper bound of absolute value of polygamma functions for
    # 1 <= k <= β
    polygamma_bounds = [abs(real(polygamma(Acb(k), Acb(2 + 2M - s)))) for k = 1:β]

    # Upper bounds of zeta function and derivatives up to β
    zeta_expansion = zeta(ArbSeries((2 + 2M - s, 1), degree = β))
    zeta_bounds(j3) = abs(zeta_expansion[j3]) * factorial(j3)

    # p_k for 0 <= k <= β
    ps = _clausen_expansion_remainder_ps(β)

    res = zero(x)

    for j1 = 0:β
        for j2 = 0:β-j1
            j3 = β - j1 - j2

            # Compute upper bound of
            # sum(abs(p_j2(2 + 2m - s)) * (x / 2π)^(2m + 1) for m = M:Inf)
            term = zero(x)
            for (exponents, coefficient) in ps[j2]
                q0 = Arb(isempty(exponents) ? 0 : exponents[1])

                # Enclosure of
                # sum((m + M + 1)^(q0 / 2) * (x / 2π)^2m for m = 0:Inf)
                S = lerch_phi((x / twopi)^2, -q0 / 2, Arb(M + 1))

                # Add factor to get an upper bound of
                # sum((1 + 2m)^(q0 / 2) * (x / 2π)^2m for m = M:Inf) / x^2M
                S *= twopi^(-2M - 1) / 2^(q0 / 2)

                # Upper bound of polygamma-factors
                polygamma_factor = one(x)
                for i = 2:length(exponents)
                    polygamma_factor *= polygamma_bounds[i-1]^exponents[i]
                end

                term += coefficient * polygamma_factor * S
            end

            # Multiply with bounds for remaining factors

            # Upper bound of part of f factor
            F = 2(log(twopi) + Arb(π) / 2)^j1 * twopi^(s - 1)

            term *= multinomial(β, j1, j2, j3) * F * zeta_bounds(j3)

            res += term
        end
    end

    return Arblib.add_error!(zero(res), res)
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

We want to compute an enclosure of each term in the expansion in `s`.

We use that `clausenc_expansion_remainder(x, s[0], β, M)` gives an
enclosure of the `β` derivative of the sum. The term in the expansion
is thus given by dividing this by `factorial(β)`.
"""
function clausens_expansion_remainder(x::Arb, s::ArbSeries, M::Integer)
    res = zero(s)

    # Compute series for s[0]
    for β = 0:Arblib.degree(s)
        res[β] = clausens_expansion_remainder(x, s[0], β, M) / factorial(β)
    end

    return ArbExtras.compose_zero(res, s)
end
