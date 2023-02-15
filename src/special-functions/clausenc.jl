export clausenc, clausencmzeta

"""
    _reduce_argument_clausen(x::Arb)

Perform argument reduction for the Clausen functions. It returns `(y,
haszero, haspi, has2pi)` where `y` is the reduced argument, `haszero`,
and `has2pi` are true if `y` overlaps with zero, and `2π`
respectively, `haspi` is true if `x` overlaps with `π` **or** `-π`.

The possible results are divided into three cases
1. `x` covers a full period, i.e. the radius is at least `π`. In this
  case it returns the full interval ``[0, 2π]`` and `haszero, haspi,
  has2pi` are all true.
2. `x` doesn't contain a multiple of `2π`. In this case the result
  satisfies `0 < y < 2π`, `haszero` and `has2pi` are false and `haspi`
  is true if `y` contains `π`.
3. `x` contains precisely one multiple of `2π`. In this case `-2π < y
  < 2π` and `haszero` is true while `has2pi` is false, `haspi` is true
  if `y` contains `-π` or `π`.

Note that if `haszero` is false then `y` is guaranteed to satisfy `0 <
y < 2π`. This is used by many of the functions.

It has issues with very small but negative values of `x`. In this case
`x + 2π < 2π` in theory but due to overestimations we might get that
they overlap. For this reason it can be beneficial to pass the
absolute value of `x` to this function and handle the sign outside of
it.
"""
function _reduce_argument_clausen(x::Arb)
    twopi = let tmp = Arb(π, prec = precision(x))
        Arblib.mul_2exp!(tmp, tmp, 1)
    end

    if Arblib.ispositive(x) && x < twopi
        # Happy path, when 0 < x < 2π
        # Optimize this case for performance so reuse twopi as much as
        # possible
        haszero = has2pi = false
        pi = Arblib.mul_2exp!(twopi, twopi, -1)
        haspi = Arblib.overlaps(x, pi)

        y = Arblib.set!(twopi, x)
    else
        pi = Arb(π, prec = precision(x))

        xdiv2pi = x / twopi

        if Arblib.cmp_2exp(Arblib.radref(xdiv2pi), -1) >= 0
            # If the radius of x / 2π is not less than 1 / 2 then y
            # spans the full interval [0, 2π]
            y = Arb((0, twopi))
            haszero = haspi = has2pi = true
        else
            # Add/subtract a multiple of 2π so that the midpoint is on
            # the interval [0, 2π]
            k = Arblib.floor!(zero(Arf), Arblib.midref(xdiv2pi))
            y = x - k * twopi

            # We want to avoid the case when y overlaps 2π
            if Arblib.overlaps(y, twopi)
                y -= twopi
            end

            haszero = Arblib.contains_zero(y)
            haspi = Arblib.overlaps(y, pi) || Arblib.overlaps(y, -pi)
            has2pi = Arblib.overlaps(y, twopi)

            # Neither of these checks should be required
            # mathematically, though overestimations might lead to
            # issues here?

            # We guarantee this
            @assert haszero || 0 < y < twopi "$x $y"
            # It should never contain -2π or 2π in this case
            @assert !(Arblib.overlaps(y, -twopi) || has2pi)
        end
    end

    return y, haszero, haspi, has2pi
end

"""
    _choose_precision_clausen(x::Union{Arb,Acb}, s::Arb)

Compute a precision to use when computing Clausen functions with
parameter `s` and argument `x`. It returns `prec, original_prec` where
`prec` is the precision to use for the internal computations and
`original_prec` is the precision to use for the return value.

The precision to use is based on the relative accuracy of `x` and `s`
as well as how close `s` is to an integer.

The original precision is given by `Arblib._precision(x, s)` and the
precision to use for the internal computations is never higher than
this. The minimum precise is `32` if `s` is not close to an integer
and `64` if it is, close is this case is determined by the midpoint of
`s` being close than `1e-2` to an integer. In general the precision is
then given by
```
max(
    min(Arblib.rel_accuracy_bits(x), Arblib.rel_accuracy_bits(s)) + min_prec,
    min_prec,
)
```
"""
function _choose_precision_clausen(x::Union{Arb,Acb}, s::Arb)
    original_prec = Arblib._precision(x, s)

    prec = if Arblib.isexact(x) && Arblib.isexact(s)
        # In this case
        # min(Arblib.rel_accuracy_bits(x), Arblib.rel_accuracy_bits(s)) + min_prec
        # overflows
        original_prec
    else
        s_f64 = Float64(s)
        min_prec = abs(s_f64 - round(s_f64)) < 1e-2 ? 64 : 32

        min(
            max(
                min(Arblib.rel_accuracy_bits(x), Arblib.rel_accuracy_bits(s)) + min_prec,
                min_prec,
            ),
            original_prec,
        )
    end

    return prec, original_prec
end

"""
    _clausenc_polylog(x::Arb, s::Union{Arb,Integer})
    _clausenc_polylog(x::Arb, s::ArbSeries)
    _clausenc_polylog(x::Arb, s::Arb, β::Integer)

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

function _clausenc_polylog(x::Arb, s::ArbSeries)
    z = exp(Acb(0, x, prec = precision(x)))
    return ArbSeries(real.(Arblib.coeffs(polylog(AcbSeries(s), z))))
end

function _clausenc_polylog(x::Arb, s::Arb, β::Integer)
    z = exp(Acb(0, x, prec = precision(x)))
    return real(polylog(AcbSeries([s, 1], degree = β), z)[β]) * factorial(β)
end

"""
    _clausenc_zeta_f(v::Arb, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)
    _clausenc_zeta_f(v::ArbSeries, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)

Compute
```
gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
```
It assumes that the arguments after `v` are set to `2π`, `x / 2π` and
`1 - x / 2π` respectively. This function is used internally by
[`_clausenc_zeta`](@ref).
"""
function _clausenc_zeta_f(v::Arb, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)
    res = zeta(v, xinv2pi)
    tmp = zeta(v, onemxinv2pi)
    Arblib.add!(res, res, tmp)

    Arblib.gamma!(tmp, v)
    Arblib.mul!(res, res, tmp)

    Arblib.pow!(tmp, inv2pi, v)
    Arblib.mul!(res, res, tmp)

    Arblib.mul_2exp!(tmp, v, -1)
    Arblib.cos_pi!(tmp, tmp)
    Arblib.mul!(res, res, tmp)

    return res
end

function _clausenc_zeta_f(v::ArbSeries, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)
    res = zeta(v, xinv2pi)
    tmp = zeta(v, onemxinv2pi)
    Arblib.add!(res, res, tmp)

    Arblib.gamma_series!(tmp, v, length(tmp))
    Arblib.mullow!(res, res, tmp, length(res))

    Arblib.zero!(tmp)
    tmp[0] = inv2pi
    Arblib.pow_series!(tmp, tmp, v, length(tmp))
    Arblib.mullow!(res, res, tmp, length(res))

    Arblib.mul_2exp!(tmp, v, -1)
    Arblib.cos_pi_series!(tmp, tmp, length(tmp))
    Arblib.mullow!(res, res, tmp, length(res))

    return res
end

"""
    _clausenc_zeta_f_one(v::Union{Arb,ArbSeries}, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)

Compute
```
gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
```
in a way that accounts for the removable singularity at `v = 1`.

If `v = 1` then `zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` diverges but
`cospi(v / 2)` goes to zero. To compute their product we use the deflated zeta function
```
zeta_deflated(s, a) = zeta(s, a) + 1 / (1 - s)
```
We can then write the product as
```
cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) =
    cospi(v / 2) * (zeta_deflated(v, x / 2π) + zeta_deflated(v, 1 - x / 2π))
    + 2cospi(v / 2) / (v - 1)
```
The first term is well defined and can be evaluated directly. The
second term has a removable singularity and is handled using
[`fx_div_x`](@ref).

It assumes that the arguments after `v` are set to `inv(2π)`, `x / 2π`
and `1 - x / 2π` respectively. This function is used internally by
[`_clausenc_zeta`](@ref).
"""
function _clausenc_zeta_f_one(
    v::Union{Arb,ArbSeries},
    inv2pi::Arb,
    xinv2pi::Arb,
    onemxinv2pi::Arb,
)
    # Enclosure of cospi(v / 2) / (v - 1)
    cos_div_v = fx_div_x(t -> cospi((t + 1) / 2), v - 1, extra_degree = 2)

    # Enclosure of zeta_deflated(v, x / 2π) + zeta_deflated(v, 1 - x / 2π)
    z = zeta_deflated(v, xinv2pi) + zeta_deflated(v, onemxinv2pi)

    return gamma(v) * inv2pi^v * (cospi(v / 2) * z + 2cos_div_v)
end

"""
    _clausenc_zeta_f_even(v::Union{Arb,ArbSeries}, v0::Integer, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)

Compute
```
gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
```
in a way that accounts for the removable singularity at even `v <= 0`.
The argument `v0` should be set to the integer where the the removable
singularity is located.

In this case `gamma(v)` diverges but `zeta(v, x / 2π) + zeta(v, 1 - x
/ 2π)` is zero. To see that `zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` is
zero when `v` is an even non-positive integer we can use formula
[25.11.14](https://dlmf.nist.gov/25.11.E14) together with
[24.4.1](https://dlmf.nist.gov/24.4.E1). We rewrite it as
```
((zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / (v - v0)) / (rgamma(v) / (v - v0))
```
and use [`fx_div_x`](@ref) to compute the factors.

It assumes that the arguments after `v0` are set to `2π`, `x / 2π` and
`1 - x / 2π` respectively. This function is used internally by
[`_clausenc_zeta`](@ref).
"""
function _clausenc_zeta_f_even(
    v::Union{Arb,ArbSeries},
    v0::Integer,
    inv2pi::Arb,
    xinv2pi::Arb,
    onemxinv2pi::Arb,
)
    # Enclosure of rgamma(v) / (v - v0)
    rgamma_div_v = fx_div_x(t -> rgamma(t + v0), v - v0, extra_degree = 2)

    # Enclosure of (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / (v - v0)
    zeta_div_v = fx_div_x(
        t -> zeta(t + v0, xinv2pi) + zeta(t + v0, onemxinv2pi),
        v - v0,
        extra_degree = 2,
        force = true,
    )

    return inv2pi^v * cospi(v / 2) * zeta_div_v / rgamma_div_v
end

"""
    _clausenc_zeta_f_odd(v::Union{Arb,ArbSeries}, v0::Integer, inv2pi::Arb, xinv2pi::Arb, onemxinv2pi::Arb)

Compute
```
gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
```
in a way that accounts for the removable singularity at odd `v <= 0`.
The argument `v0` should be set to the integer where the the removable
singularity is located.

In this case `gamma(v)` diverges but `cospi(v / 2)` is zero. We
rewrite it as
```
(cospi(v / 2) / (v - v0)) / (rgamma(v) / (v - v0))
```
and use [`fx_div_x`](@ref) to compute the factors.

It assumes that the arguments after `v0` are set to `2π`, `x / 2π` and
`1 - x / 2π` respectively. This function is used internally by
[`_clausenc_zeta`](@ref).
"""
function _clausenc_zeta_f_odd(
    v::Union{Arb,ArbSeries},
    v0::Integer,
    inv2pi::Arb,
    xinv2pi::Arb,
    onemxinv2pi::Arb,
)
    # Enclosure of rgamma(v) / (v - v_integer)
    rgamma_div_v = fx_div_x(t -> rgamma(t + v0), v - v0, extra_degree = 2)

    # Enclosure of cospi(v / 2) / (v - v_integer)
    cos_div_v = fx_div_x(t -> cospi((t + v0) / 2), v - v0, extra_degree = 2)

    return inv2pi^v * cos_div_v / rgamma_div_v * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))
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

The formula is well defined as long as `v` doesn't overlap with any
integer `<= 1` (meaning `s` doesn't overlap with any non-negative
integer). For integer `v <= 1` there is a removable singularity that
has to be handled. This is done using the functions
[`_clausenc_zeta_f_one`](@ref), [`_clausenc_zeta_f_even`](@ref) and
[`_clausenc_zeta_f_odd`](@ref), handling the cases when `v` is one and
`v <= 0` and even or odd respectively.

It currently only handles `0 < x < 2π`. If `s` is wide, as determined
by `iswide(s)` it computes a tighter enclosure of the coefficients
using a Taylor expansion in `s`.
"""
function _clausenc_zeta(x::Arb, s::Arb)
    # Check that x > 0
    Arblib.ispositive(x) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    inv2pi = let tmp = Arb(π, prec = precision(x))
        Arblib.mul_2exp!(tmp, tmp, 1)
        Arblib.inv!(tmp, tmp)
    end
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
    # gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
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
            _clausenc_zeta_f_one(v, inv2pi, xinv2pi, onemxinv2pi)
        elseif unique && iseven(v0_integer) && v0_integer <= 0
            _clausenc_zeta_f_even(v, v0_integer, inv2pi, xinv2pi, onemxinv2pi)
        elseif unique && isodd(v0_integer) && v0_integer <= 0
            _clausenc_zeta_f_odd(v, v0_integer, inv2pi, xinv2pi, onemxinv2pi)
        else
            _clausenc_zeta_f(v, inv2pi, xinv2pi, onemxinv2pi)
        end
    end

    if iswide(s)
        return ArbExtras.enclosure_series(f, v, degree = 2)
    else
        return f(v)
    end
end

"""
    _clausenc_zeta(x::Acb, s::Arb)

Evaluation of the `clausenc` function through the zeta function.

This uses the same formula as the method with `x::Arb`, valid for `0 <
real(x) < 2π`. If `s` overlaps with any non-negative integer the
result will be indeterminate.

This method is only used in the integration for bounding the error
term. It is therefore not as optimized as many of the other methods.
"""
function _clausenc_zeta(x::Acb, s::Arb)
    twopi = 2Arb(π, prec = Arblib._precision(x, s))

    Arblib.ispositive(Arblib.realref(x)) || throw(
        DomainError(x, "method only supports x with real part on the interval (0, 2π)"),
    )
    Arblib.realref(x) < twopi || throw(
        DomainError(x, "method only supports x with real part on the interval (0, 2π)"),
    )

    inv2pi = inv(twopi)
    xinv2pi = x * inv2pi

    v = let v = zero(s) # We do it like this to preserve the precision
        Arblib.neg!(v, s)
        Arblib.add!(v, v, 1)
    end

    onemxinv2pi = let onemxinv2pi = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(onemxinv2pi, xinv2pi)
        Arblib.add!(onemxinv2pi, onemxinv2pi, 1)
    end

    return ArbExtras.enclosure_series(
        v -> gamma(v) * inv2pi^v * cospi(v / 2),
        v,
        degree = 2,
    ) * (zeta(Acb(v), xinv2pi) + zeta(Acb(v), onemxinv2pi))
end

"""
    _clausenc_zeta(x::Arb, s::ArbSeries)

Evaluation of the `clausenc` function through the zeta function as a
power series in `s`.

It supports `s` overlapping non-negative integers in the same way
`_clausenc_zeta(x::Arb, s::Arb)` does, using [`fx_div_x`](@ref) to
handle the removable singularities.
"""
function _clausenc_zeta(x::Arb, s::ArbSeries)
    # Handle the case when s has degree 0, so that we can assume that
    # the degree is positive from here on.
    iszero(Arblib.degree(s)) && return ArbSeries(_clausenc_zeta(x, s[0]))

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

    if s[0] > 0 && is_approx_integer(s[0])
        s = copy(s)
        s[0] = union(s[0], Arb(round(Float64(s[0]))))
    end

    v = 1 - s

    s0 = s[0]

    unique, s0_integer = unique_integer(s0)

    if unique && s0_integer == 0
        cos_div_v = fx_div_x(t -> cospi((t + 1) / 2), v - 1, extra_degree = 1)

        return gamma(v) *
               inv2pi^v *
               (
                   cospi(v / 2) *
                   (zeta_deflated(v, xinv2pi) + zeta_deflated(v, onemxinv2pi)) + 2cos_div_v
               )
    elseif unique && s0_integer > 0
        v0_integer = 1 - s0_integer

        # Enclosure of rgamma(v) / (v - v_integer)
        rgamma_div_v =
            fx_div_x(t -> rgamma(t + v0_integer), v - v0_integer, extra_degree = 2)

        if iseven(s0_integer)
            # Enclosure of cospi(v / 2) / (v - v_integer)
            cos_div_v =
                fx_div_x(t -> cospi((t + v0_integer) / 2), v - v0_integer, extra_degree = 2)

            # Enclosure of gamma(v) * cospi(v / 2)
            gammacos = cos_div_v / rgamma_div_v

            rest = inv2pi^v * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))

            return gammacos * rest
        else
            # Enclosure of (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / (v - v_integer)
            zeta_div_v = fx_div_x(
                t -> zeta(t + v0_integer, xinv2pi) + zeta(t + v0_integer, onemxinv2pi),
                v - v0_integer,
                extra_degree = 2,
                force = true,
            )

            # Enclosure of gamma(v) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
            gammazeta = zeta_div_v / rgamma_div_v

            rest = inv2pi^v * cospi(v / 2)

            return gammazeta * rest
        end
    else
        return gamma(v) *
               inv2pi^v *
               cospi(v / 2) *
               (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))
    end
end

"""
    _clausenc_zeta(x::Float64, s::Float64)

Evaluation of the `clausenc` function through the zeta function.

This uses the same formula as the method with `x::Arb`. But it doesn't
any special handling when `s` is an integer, it will return `NaN`.
"""
function _clausenc_zeta(x::Float64, s::Float64)
    x = mod2pi(x)
    v = 1 - s
    return gamma(v) / (2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
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
        res = ArbExtras.enclosure_series(
            ArbExtras.derivative_function(s -> _clausenc_zeta(x, s), β),
            s,
            degree = 10,
        )
    else
        res = _clausenc_zeta(x, ArbSeries((s, 1), degree = β))[β] * factorial(β)
    end

    return res
end

"""
    clausenc(x, s)

Compute the Clausen function ``C_s(x)``.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

If `x` contains zero and we don't have `s > 1` it returns an
indeterminate result.
- **IMPROVE:** Use that for even negative integers it is exactly zero.

If `x` contains zero and `s > 1` it computes an enclosure using that
the extrema is attained at `x = 0`, `abs_ubound(x)` or `π`, where we
have `clausenc(0, s) = zeta(s)` and `clausenc(π, s) = -eta(s)`.

If `x` doesn't contain zero we are assured by
[`_reduce_argument_clausen`](@ref) that `0 < x < 2π`.

If `x` is a wide ball (not containing zero), as determined by
`iswide(x)`, it computes a tighter enclosure by using that the
function is `2π`-periodic, monotone for `x ∈ [0, π]` and even, so that
it's enough to evaluate on the endpoints and possibly at `π` if it is
contained in `x`.
"""
function clausenc(x::Arb, s::Arb)
    # The function is even and this avoids issues with very small
    # negative x
    x = abs(x)

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
            # Extrema at x = 0 and x = π
            res = union(zeta(s), -Arblib.realref(eta(Acb(s))))
        elseif iszero(x)
            res = zeta(s)
        else
            # Extrema at x = 0 and upper bound of abs(x)
            res = union(zeta(s), clausenc(ArbExtras.enclosure_ubound(abs(x)), s))
        end
    elseif iswide(x) # We can now assume that 0 < x < 2π
        xₗ, xᵤ = ArbExtras.enclosure_getinterval(x)

        res = _clausenc_zeta(xₗ, s)
        Arblib.union!(res, res, _clausenc_zeta(xᵤ, s))
        if haspi
            Arblib.union!(res, res, clausenc(Arb(π; prec), s))
        end
    else
        res = _clausenc_zeta(x, s)
    end

    if prec != original_prec
        res = setprecision(res, original_prec)
    end

    return res
end

function clausenc(x::Acb, s::Arb)
    x_real, haszero, _, _ = _reduce_argument_clausen(real(x))
    x = Acb(x_real, Arblib.imagref(x))

    prec, original_prec = _choose_precision_clausen(x, s)
    if prec != original_prec
        x = setprecision(x, prec)
        s = setprecision(s, prec)
    end

    # If x overlaps zero or s is a non-negative integer use polylog formulation
    if haszero || (isinteger(s) && Arblib.isnonnegative(s))
        res = (polylog(Acb(s), exp(im * x)) + polylog(Acb(s), exp(-im * x))) / 2
    else
        res = _clausenc_zeta(x, s)
    end

    if prec != original_prec
        res = setprecision(res, original_prec)
    end

    return res
end

function clausenc(x::Float64, s::Float64)
    if isinteger(s) && s >= 0
        return Float64(clausenc(Arb(mod2pi(x)), Arb(s)))
    else
        return _clausenc_zeta(x, s)
    end
end

clausenc(x::S, s::T) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausenc(convert(Arb, x), convert(Arb, s)))

"""
    clausenc(x::ArbSeries, s)

Compute the Taylor series of the Clausen function ``C_s(x)``.

It's computed by directly computing the Taylor coefficients by
differentiating `clausenc` and then composing with `x`.
"""
function clausenc(x::ArbSeries, s)
    s = convert(Arb, s)

    x₀ = _reduce_argument_clausen(x[0])[1]

    res = zero(x)

    # Precompute clausenc functions for i = 0:2:Arblib.degree(x)+1.
    # They are used both for the i-th coefficients and also in the
    # computation of clausens for the (i-1)-th coefficient.
    clausencs = [clausenc(x₀, s - i) for i = 0:2:Arblib.degree(x)+1]
    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausencs[i÷2+1] / factorial(i)
        else
            res[i] =
                -(-1)^(i ÷ 2) * clausens(x₀, s - i, deriv_x = clausencs[(i+1)÷2+1]) /
                factorial(i)
        end
    end

    return ArbExtras.compose_zero!(res, res, x)
end

"""
    clausenc(x::Arb, s::ArbSeries)

Compute the Taylor series of the Clausen function ``C_s(x)`` in the
parameter `s`.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

If `x` contains zero it returns an indeterminate result. For `s > 1`
it would be possible to compute a finite result, but this has not been
needed.
"""
function clausenc(x::Arb, s::ArbSeries)
    # The function is even and this avoids issues with very small
    # negative x
    x = abs(x)

    x, haszero, _, _ = _reduce_argument_clausen(x)

    if haszero
        return indeterminate(s)
    else
        return _clausenc_zeta(x, s)
    end
end

"""
    clausenc(x, s, β)

Compute ``C_s^{(β)}(x)``, that is `clausenc(x, s)` differentiated `β`
times w.r.t. `s`.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

If `x` contains zero and we don't have `s > 1` it returns an
indeterminate result. If `s > 1` and `x` is not exactly zero it
computes an extremely naive enclosure using that `abs(clausenc(x, s,
β)) < abs(clausenc(0, s, β))`
- **IMPROVE:** Compute a tighter enclosure in this case. Either
  checking the derivative at the endpoint, similarly to
  [`clausens`](@ref), or expanding at `x = 0`.

If `x` is a wide ball (not containing zero), as determined by
`iswide(x)`, it computes a tighter enclosure by first checking if the
derivative doesn't contains zero, if not it uses monotonicity to only
evaluate at endpoints. If the derivative does contain zero it uses a
zeroth order approximation instead.
"""
function clausenc(x::Arb, s::Arb, β::Integer)
    # The function is even and this avoids issues with very small
    # negative x
    x = abs(x)

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
        else
            # Value at x = 0
            z = isone(β) ? dzeta(s) : zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
            if haspi
                # Absolute value upper bounded by value at x = 0
                res = union(-z, z)
            elseif iszero(x)
                res = z
            else
                # IMPROVE: Expand at x = 0? Check derivative at abs_ubound(x)?
                # Absolute value upper bounded by value at x = 0
                res = union(-z, z)
            end
        end
    elseif iswide(x) # We can now assume that 0 < x < 2π
        # Compute derivative. We here call _clausens_zeta directly to
        # avoid potentially infinite recursion. This is okay since 0 <
        # x < 2π.
        dclausenc = -_clausens_zeta(x, s - 1, β)
        if Arblib.contains_zero(dclausenc)
            # Use a zero order approximation
            mid = Arblib.midpoint(Arb, x)
            res = Arblib.add_error!(clausenc(mid, s, β), (x - mid) * dclausenc)
        else
            # Use that it's monotone
            xₗ, xᵤ = ArbExtras.enclosure_getinterval(x)
            res = union(clausenc(xₗ, s, β), clausenc(xᵤ, s, β))
        end
    else
        res = _clausenc_zeta(x, s, β)
    end

    if prec != original_prec
        res = setprecision(res, original_prec)
    end

    return res
end

clausenc(x::S, s::T, β::Integer) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausenc(convert(Arb, x), convert(Arb, s), β))

"""
    clausenc(x::ArbSeries, s, β)

Compute the Taylor series with respect to `x` of ``C_s^{(β)}(x)``,
that is `clausenc(x, s)` differentiated `β` times w.r.t. `s`.

It's computed by directly computing the Taylor coefficients by
differentiating ``C_s^{(\beta)}`` and then composing with `x`.
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

    return ArbExtras.compose_zero!(res, res, x)
end

"""
    clausenc_expansion_main(s::Arb)

Compute
```
gamma(1 - s) * sinpi(s / 2)
```
Which is the coefficient for the non-analytic term in the expansion of
`clausenc(x, s)` at zero.

In general the value is computed directly using the formula above. The
exception is when `s` is close to a positive even integer, in which
case it has a removable singularity that has to be dealt with.

We consider `s` to be close to a positive even integer if
```
s > 1 && iseven(s0) && !(abs(s - s0) > 0.1)
```
where `s0 = round(Int, Float64(s))`. In that case we rewrite the
expression as
```
gamma(1 - s) * sinpi(s / 2) = gamma(1 - s) / rising(1 - s, s0) * sinpi(s / 2)
    = gamma(1 - s) / rising(1 - s, s0 - 1) * sinpi(s / 2) / (s0 - s)
    = -(-1)^(s0 ÷ 2) * π / 2 * gamma(1 - s) / rising(1 - s, s0 - 1) * sinc((s0 - s) / 2)
```
"""
function clausenc_expansion_main(s::Arb)
    s0 = round(Int, Float64(s))

    if s > 1 && iseven(s0) && !(abs(s - s0) > 0.1)
        # Treat the removable singularity
        C =
            -(-1)^(s0 ÷ 2) * Arb(π) / 2 * ArbExtras.enclosure_series(s, degree = 2) do s
                sinc_s0ms = _sinc((s0 - s) / 2)

                if s isa ArbSeries
                    # _sinc performs poorly for arguments close to zero, it is
                    # often better to slightly widen the argument to include
                    # zero
                    t = (s0 - s) / 2
                    t[0] = union(t[0], Arb(0))
                    sinc_s0ms_wide = _sinc(t)
                    for i = 0:Arblib.degree(s)
                        sinc_s0ms[i] = intersect(sinc_s0ms[i], sinc_s0ms_wide[i])
                    end
                end

                gamma(1 + s0 - s) / rising(1 - s, s0 - 1) * sinc_s0ms
            end
    elseif iswide(s)
        # When s is close to a positive odd integer the series
        # enclosure is not very good, a direct evaluation gives a
        # better enclosure. We can therefore improve the situation by
        # taking the intersection of the two enclosures. Since it is a
        # cheap computation we always do this.
        C = intersect(
            ArbExtras.enclosure_series(s -> gamma(1 - s) * sinpi(s / 2), s, degree = 2),
            gamma(1 - s) * sinpi(s / 2),
        )
    else
        C = gamma(1 - s) * sinpi(s / 2)
    end

    return C
end

"""
    clausenc_expansion(x, s, M::Integer)

Compute the asymptotic expansion of `clausenc(x, s)` at zero up to
order `2M - 2`, meaning that the error term is of order `2M`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a `ArbSeries` `P` and the
remainder term `E`. The `M` is the same as in Lemma 2.1 in
arXiv:1810.10935.

It satisfies that
```
clausenc(y, s) ∈ C*abs(y)^e + P(y) + E*y^(2M)
```
for all `abs(y) <= abs(x)`.

If `skip_constant = true` it doesn't compute the constant term in the
expansion. This is useful if you want to compute the expansion for
`clausenc(x, s) - clausenc(0, s)`.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure of the coefficients using a Taylor expansion in `s`.

The coefficient `C` is given by `gamma(1 - s) * sinpi(s / 2)` has a
singularity for positive integer values of `s` where the `gamma`
function blows up. For even values this singularity is removable and
can be handled using [`fx_div_x`](@ref).

For odd values of `s` the singularity is not removable, instead the
non-analytic term coincides with one of the analytic terms, with their
coefficients blowing up in different directions. This case is
currently not handled in any special way and it just returns an
indeterminate result. Note that it is not clear how to handle this
case if we want to allow `s` overlapping odd integers.
"""
function clausenc_expansion(x::Arb, s::Arb, M::Integer; skip_constant = false)
    # Coefficient and exponent for non-analytic term
    C = clausenc_expansion_main(s)
    e = s - 1

    # Analytic terms
    P = ArbSeries(degree = 2M - 2, prec = precision(x))
    for m = (skip_constant ? 1 : 0):M-1
        if iswide(s)
            term = ArbExtras.enclosure_series(s, degree = 2) do s
                z = zeta(s - 2m)
                if s isa ArbSeries && abs(Float64(s[0] - 2m)) <= 0.1 && iswide(s)
                    # zeta doesn't handle wide arguments that are
                    # close to zero but not centered around zero very
                    # well. In this case it can be better to force the
                    # argument to be centered at zero instead.
                    t = s - 2m
                    t[0] = union(t[0], -t[0])
                    z2 = zeta(t)

                    if isfinite(z) && isfinite(z2)
                        z = ArbSeries(intersect.(Arblib.coeffs(z), Arblib.coeffs(z2)))
                    elseif isfinite(z2)
                        z = z2
                    end
                end
                z
            end
        else
            term = zeta(s - 2m)
        end
        P[2m] = (-1)^m * term / factorial(2m)
    end

    # Error term
    E = clausenc_expansion_remainder(x, s, M)

    return (C, e, P, E)
end

"""
    _clausen_expansion_remainder_ps(β::Integer)

Helper function for computing a representation of `p_k` for `k = 0:β`

Here `p_k` is defined by the recurrence ``p_{k + 1}(s) =
ψ^{(0)}(s)p_k(s) + p_k'(s)`` and ``p_0 = 1``.

It returns an `OffsetVector` with indices `0:β`. Each element is a
dictionary which maps a list of exponents to a coefficient. More
precisely the element `[q0, q1, ..., qk] => c` in the dictionary
represents the term
```
c * ψ(s, 0)^q0 * ψ(s, 1)^q1 * ... *ψ(s, k)^qk
```
where `ψ(s, i)` here means the function ``ψ^{(i)}(s)``. Note that this
doesn't actually evaluate the `ψ` functions in any way, it only stores
the representation of the terms.
"""
function _clausen_expansion_remainder_ps(β::Integer)
    p_0 = OrderedDict(Int[] => 1)
    ps = OffsetVector([p_0], 0:0)

    for k = 1:β
        p_k = empty(ps[k-1])

        # The term ψ^(0) * p_k
        for (exponents, coefficient) in ps[k-1]
            if isempty(exponents)
                new_exponents = [1] # ψ^(0) * 1 = ψ^(0)
            else
                new_exponents = copy(exponents)
                new_exponents[1] += 1 # Multiply by ψ^(0)
            end
            p_k[new_exponents] = coefficient
        end

        # The term p_k'
        for (exponents, coefficient) in ps[k-1]
            for l = 1:length(exponents)
                if !iszero(exponents[l])
                    # Differentiate ψ^(l+1)^exponents[l]
                    new_exponents = copy(exponents)

                    new_exponents[l] -= 1
                    if l == length(exponents)
                        push!(new_exponents, 1)
                    else
                        new_exponents[l+1] += 1
                    end

                    p_k[new_exponents] =
                        get(p_k, new_exponents, 0) + coefficient * exponents[l]
                end
            end
        end

        push!(ps, p_k)
    end

    return ps
end


"""
    clausenc_expansion_remainder(x::Arb, s::Arb, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausenc(x, s)` at zero up to order `2M - 2`, meaning that the
remainder is of order `2M`.

This is the `E` occurring in [`clausenc_expansion`](@ref) and is given
by
```
sum((-1)^m * zeta(s - 2m) * x^2m / factorial(2m) for m = M:Inf) / x^2M
```

It requires that `abs(x) < 2π`, `s >= 0` and `2M >= s + 1`. In this
case an upper bound for the absolute value of the remainder is given
by
```
2(2π)^(1 + s - 2M) * abs(sinpi(s / 2)) * zeta(2M + 1 - s) / (4π^2 - x^2)
```
and this functions returns a ball centered at zero with this radius.
"""
function clausenc_expansion_remainder(x::Arb, s::Arb, M::Integer)
    pi = Arb(π)

    abs(x) < 2pi || throw(DomainError(x, "x must be less than 2π"))
    s >= 0 || throw(DomainError(M, "s must be positive, got s = $s"))
    2M >= s + 1 || throw(DomainError(M, "must have 2M >= s + 1, got s = $s"))

    return Arblib.add_error!(
        zero(x),
        2(2pi)^(1 + s - 2M) * abs(sinpi(s / 2)) * zeta(2M + 1 - s) / (4pi^2 - x^2),
    )
end

"""
    clausenc_expansion_remainder(x::Arb, s::Arb, β::Integer, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausenc(x, s, β)` at zero up to order `2M - 2`, meaning that the
remainder is of order `2M`.

This is the tail in the expansion at `x = 0` and is given by
```
sum((-1)^m * dzeta(s - 2m, β) * x^2m / factorial(2m) for m = M:Inf) / x^2M
```
where we by `dzeta(s - 2m, β)` mean the zeta-function differentiated
`β` times.

It requires that `abs(x) < 2π`, `s >= 0` and `2M >= s + 1`.

The upper bound of the absolute value of the remainder is given by a
somewhat awkward expression involving a multitude of special
functions. See the paper for details.

This functions returns a ball centered at zero with the upper bound as
radius.
"""
function clausenc_expansion_remainder(x::Arb, s::Arb, β::Integer, M::Integer)
    β == 0 && return clausenc_expansion_remainder(x, s, M)

    abs(x) < 2pi || throw(DomainError(x, "x must be less than 2π"))
    s >= 0 || throw(DomainError(M, "s must be positive, got s = $s"))
    2M >= s + 1 || throw(DomainError(M, "must have 2M >= s + 1, got s = $s"))

    twopi = 2Arb(π)

    # Function for computing multinomial
    # Denominator is smaller than numerator, so overflow
    # would always be caught by factorial(β)
    multinomial(β, j1, j2, j3) =
        Arb(factorial(β) // (factorial(j1) * factorial(j2) * factorial(j3)))

    # Upper bound of absolute value of polygamma functions for
    # 1 <= k <= β
    polygamma_bounds = [abs(real(polygamma(Acb(k), Acb(1 + 2M - s)))) for k = 1:β]

    # Upper bounds of zeta function and derivatives up to β
    zeta_expansion = zeta(ArbSeries((1 + 2M - s, 1), degree = β))
    zeta_bounds(j3) = abs(zeta_expansion[j3]) * factorial(j3)

    # p_k for 0 <= k <= β
    ps = _clausen_expansion_remainder_ps(β)

    res = zero(x)

    for j1 = 0:β
        for j2 = 0:β-j1
            j3 = β - j1 - j2

            # Compute upper bound of
            # sum(abs(p_j2(1 + 2m - s)) * (x / 2π)^2m for m = M:Inf)
            term = zero(x)
            for (exponents, coefficient) in ps[j2]
                q0 = Arb(isempty(exponents) ? 0 : exponents[1])

                # Enclosure of
                # sum((m + M + 1 / 2)^(q0 / 2) * (x / 2π)^2m for m = 0:Inf)
                S = lerch_phi((x / twopi)^2, -q0 / 2, M + Arb(1 // 2))

                # Add factor to get an upper bound of
                # sum((1 + 2m)^(q0 / 2) * (x / 2π)^2m for m = M:Inf) / x^2M
                S *= twopi^(-2M) / 2^(q0 / 2)

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
    clausenc_expansion_remainder(x::Arb, s::ArbSeries, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausenc(x, s)` at zero up to order `2M - 2`, meaning that the
remainder is of order `2M`.

The remainder term is given by
```
x^2M * sum((-1)^m * zeta(s - 2m) * x^(2(m - M)) / factorial(2m) for m = M:Inf)
```

We want to compute an enclosure of each term in the expansion in `s`.

We use that `clausenc_expansion_remainder(x, s[0], β, M)` gives an
enclosure of the `β` derivative of the sum. The term in the expansion
is thus given by dividing this by `factorial(β)`.
"""
function clausenc_expansion_remainder(x::Arb, s::ArbSeries, M::Integer)
    res = zero(s)

    # Compute series for s[0]
    for β = 0:Arblib.degree(s)
        res[β] = clausenc_expansion_remainder(x, s[0], β, M) / factorial(β)
    end

    return ArbExtras.compose_zero(res, s)
end

"""
    clausenc_expansion_odd_s_singular_K1_K2(s::Arb, m::Integer)

Compute `K1, K2`, which together with `K3` from
[`clausenc_expansion_odd_s_singular_K1_K2`](@ref) satisfy that
```
K1 * abs(x)^(s - 1) + K2 * x^2m + K3 * x^2m * (abs(x)^(s - (2m + 1)) - 1) / (s - (2m + 1))
```
gives an enclosure of
```
gamma(1 - s) * sinpi(s / 2) * abs(x)^(s - 1) + (-1)^m * zeta(s - 2m) * x^2m / factorial(2m)
```
for any `x`. It is done in a way that works well for `s` overlapping
the removable singularity at `2m + 1` and wide `s` close to `2m + 1`.

The factor
```
x^2m * (x^(s - (2m + 1)) - 1) / (2m + 1 - s)
```
can be computed using [`x_pow_s_x_pow_t_m1_div_t`](@ref) with
```
x_pow_s_x_pow_t_m1_div_t(x, 2m, s - (2m + 1))
```

# Derivation
Using the deflated zeta function and adding and subtracting `(-1)^m *
x^(s - 1) / ((2m + 1 - s) * factorial(2m))` we can split the function
into three terms
```
(gamma(1 - s) * sinpi(s / 2) - (-1)^m / ((2m + 1 - s) * factorial(2m))) * x^(s - 1)
```
```
-(-1)^m * zeta_deflated(s - 2m) / factorial(2m) * x^2m
```
and
```
-(-1)^m / factorial(2m) * x^2m * (x^(s - (2m + 1)) - 1) / (2m + 1 - s)
```
Where the factors in front of them are `K1`, `K2` and `K3`
respectively.

# Computing `K1`
For `K1` we use that
```
gamma(1 - s) = gamma(2m + 2 - s) / rising(1 - s, 2m + 1)
             = gamma(2m + 2 - s) / ((2m + 1 - s) * rising(1 - s, 2m))
```
to rewrite it as
```
K1 = (
    gamma(2m + 2 - s) / rising(1 - s, 2m) * sinpi(s / 2)
    - (-1)^m / factorial(2m)
) / (2m + 1 - s)
```
which has a removable singularity at `s = 2m + 1`.
"""
function clausenc_expansion_odd_s_singular_K1_K2(s::Arb, m::Integer)
    m >= 1 || throw(ArgumentError("m should be positive"))

    K1_f(s) =
        if (s isa Arb && contains(s, 2m + 1)) || (s isa ArbSeries && contains(s[0], 2m + 1))
            # Handle the removable singularity
            fx_div_x(2m + 1 - s, force = true) do t
                gamma(1 + t) / rising(t - 2m, 2m) * sinpi((2m + 1 - t) / 2) -
                (-1)^m // factorial(2m)
            end
        else
            (
                gamma(2m + 2 - s) / rising(1 - s, 2m) * sinpi(s / 2) -
                (-1)^m // factorial(2m)
            ) / (2m + 1 - s)
        end

    K1 = ArbExtras.enclosure_series(s, degree = 2) do s
        if iswide(s)
            # Include the removable singularity to get better enclosures
            @assert s isa ArbSeries
            s = copy(s)
            s[0] = union(s[0], Arb(2m + 1))
        end
        K1_f(s)
    end
    K2 = ArbExtras.enclosure_series(s, degree = 2) do s
        if iswide(s)
            # Include the removable singularity to get better enclosures
            @assert s isa ArbSeries
            s = copy(s)
            s[0] = union(s[0], Arb(2m + 1))
        end
        (-1)^m * zeta_deflated(s - 2m, Arb(1)) / factorial(2m)
    end

    return K1, K2
end

"""
    clausenc_expansion_odd_s_singular_K3(m::Integer)

Compute `K3 = -(-1)^m / factorial(2m)` converted to `Arb`, as
described in [`clausenc_expansion_odd_s_singular_K1_K2`](@ref).
"""
clausenc_expansion_odd_s_singular_K3(m::Integer) = convert(Arb, -(-1)^m // factorial(2m))

"""
    clausencmzeta(x, s)

Compute `clausenc(x, s) - zeta(s)` for `s > 1`.

Note that for `s > 1` we have `clausenc(0, s) = zeta(s)` so this is
`clausenc` normalized to be zero at `x = 0`. When `s <= 1` we no
longer have `clausenc(0, s) = zeta(s)` and instead `clausenc(0, s)`
has a singularity, the expression `clausenc(x, s) - zeta(s)` is still
well defined however.

Typically this method just calls [`clausenc`](@ref) and [`zeta`](@ref)
directly and gives no direct benefit. However, if `s > 1` and is a
wide ball it can better handle the cancellations between the two terms
by computing them together.

For `s > 1` it uses that the function is non-decreasing in `s`. That
this is the case can be checked by looking at the Fourier series of
the `clausenc` and the series for `zeta` and check that all terms in
their difference are non-negative.

We could compute a better enclosure for `s <= 1` by treating the terms
together and expanding in `s`. However this is not relevant for our
use case.
"""
function clausencmzeta(x::Arb, s::Arb)
    if s > 1 && iswide(s)
        prec, original_prec = _choose_precision_clausen(x, s)
        if prec != original_prec
            x = setprecision(x, prec)
            s = setprecision(s, prec)
        end

        sₗ, sᵤ = getinterval(Arb, s)
        sₗ > 1 || return indeterminate(setprecision(x, original_prec))
        res_lower = clausenc(x, sₗ) - zeta(sₗ)
        res_upper = clausenc(x, sᵤ) - zeta(sᵤ)
        return Arb((res_lower, res_upper), prec = original_prec)
    else
        return clausenc(x, s) - zeta(s)
    end
end

function clausencmzeta(x::ArbSeries, s::Arb)
    res = zero(x)
    x₀ = x[0]

    res[0] = clausencmzeta(x₀, s)

    # Precompute clausenc functions for i = 2:2:Arblib.degree(x)+1.
    # They are used both for the i-th coefficients and also in the
    # computation of clausens for the (i-1)-th coefficient.
    clausencs = [clausenc(x₀, s - i) for i = 2:2:Arblib.degree(x)+1]
    for i = 1:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausencs[i÷2] / factorial(i)
        else
            res[i] =
                -(-1)^(i ÷ 2) * clausens(x₀, s - i, deriv_x = clausencs[(i+1)÷2]) /
                factorial(i)
        end
    end

    return ArbExtras.compose_zero!(res, res, x)
end

clausencmzeta(x::ArbSeries, s) = clausenc(x, s) - zeta(Arb(s, prec = precision(x)))

clausencmzeta(x::Arb, s::ArbSeries) = clausenc(x, s) - zeta(s)

function clausencmzeta(x, s)
    x, s = promote(x, s)
    return clausenc(x, s) - zeta(s)
end

"""
    clausencmzeta(x, s, β::Integer)

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
        return clausenc(x, s, β) - zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
    end
end

function clausencmzeta(x, s, β::Integer)
    x, s = promote(x, s)
    if isone(β)
        return clausenc(x, s, β) - dzeta(s)
    else
        return clausenc(x, s, β) - zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
    end
end
