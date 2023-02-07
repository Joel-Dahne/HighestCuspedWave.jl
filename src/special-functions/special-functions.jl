"""
    _sinc(x)

The same as `sinc(x)` but for `x::ArbSeries` it allows evaluation
around zero.
"""
_sinc(x) = sinc(x)

function _sinc(x::ArbSeries)
    if Arblib.degree(x) >= 1 && Arblib.contains_zero(Arblib.ref(x, 0))
        return fx_div_x(sin, π * x)
    else
        return sinc(x)
    end
end

"""
    dzeta(s)

Compute the Zeta function differentiated once with respect to `s`.
"""
dzeta(s::Arb) = zeta(ArbSeries((s, 1)))[1]
dzeta(s) = convert(float(typeof(s)), dzeta(Arb(s)))

"""
    beta_inc(a, b, z)

Compute the (not regularised) incomplete beta function ``B(a, b; z)``.

Note that this method is different than
[`SpecialFunctions.beta_inc`](@ref) both in that it returns the
non-regularised value and that it only returns one value.
"""
beta_inc(a::Acb, b::Acb, z::Acb) = Arblib.hypgeom_beta_lower!(zero(z), a, b, z, 0)
beta_inc(a::Arb, b::Arb, z::Arb) = Arblib.hypgeom_beta_lower!(zero(z), a, b, z, 0)

"""
    beta_inc_zeroone(a, b, z)

Compute the (not regularised) incomplete beta function ``B(a, b; z)``
assuming that `0 <= z <= 1`, discarding any other numbers in the
interval.

It uses the monotonicity in `z` to only have to compute the values at
the endpoint. That the function is increasing in `z` on `[0, 1]` is
easily seen from the integral representation
```
∫_0^z t^(a - 1) (1 - t)^(b - 1) dt
```
"""
function beta_inc_zeroone(a::Arb, b::Arb, z::Arb)
    z_lower = max(lbound(Arb, z), zero(z))
    z_upper = min(ubound(Arb, z), one(z))
    return Arb((beta_inc(a, b, z_lower), beta_inc(a, b, z_upper)))
end

"""
    hypgeom_2f1(a, b, c, z)

Compute the hypergeometric ``₂F₁(a, b, c, z)`` function.
"""
hypgeom_2f1(a::Arb, b::Arb, c::Arb, z::Arb) = Arblib.hypgeom_2f1!(zero(z), a, b, c, z, 0)
hypgeom_2f1(a::Acb, b::Acb, c::Acb, z::Acb) =
    Arblib.hypgeom_2f1!(zero(z), a, b, c, z, flags = 0)

"""
    rgamma(x)

Compute the reciprocal gamma function, defined by `rgamma(x) = 1 /
gamma(x)`.
"""
rgamma(x::Arb) = Arblib.rgamma!(zero(x), x)
rgamma(x::ArbSeries) = Arblib.rgamma_series!(zero(x), x, length(x))

"""
    rising(x, n)

Compute the rising factorial ``(x)_n``.
"""
rising(x::Arb, n::Arb) = Arblib.rising!(zero(x), x, n)
rising(x::Arb, n::Integer) = Arblib.rising!(zero(x), x, convert(UInt, n))
rising(x::ArbSeries, n::Integer) =
    Arblib.rising_ui_series!(zero(x), x, convert(UInt, n), length(x))

"""
    _zeta_deflated(s::ArbSeries, a::Arb)

Compute the deflated zeta function.

The implementation in Arb doesn't support evaluation for `s`
overlapping one but not exactly one. This method is a reimplementation
of the code from Arb but adjuste to handle the removable singularity.

Most of it follows the code in `_acb_poly_zeta_cpx_series`. The part
which requires adjustments is in `_acb_poly_zeta_em_sum`, so this
implementation is adjusted.
"""
function _zeta_deflated(s::ArbSeries, a::Arb)
    Arblib.ispositive(a) || return indeterminate(s)

    # _acb_poly_acb_invpow_cpx is not documented and hence not
    # included in Arblib
    function _acb_poly_acb_invpow_cpx!(
        res::AcbRefVector,
        n::Acb,
        c::Acb,
        trunc::Clong,
        prec::Clong,
    )
        ccall(
            (:_acb_poly_acb_invpow_cpx, Arblib.libarb),
            Cvoid,
            (
                Ptr{Arblib.acb_struct},
                Ref{Arblib.acb_struct},
                Ref{Arblib.acb_struct},
                Clong,
                Clong,
            ),
            res,
            n,
            c,
            trunc,
            prec,
        )
    end

    # Corresponding to FLINT_BIT_COUNT
    bit_count(n::Union{UInt64,Int64}) = 64 - leading_zeros(n)

    # Compute series for constant part of s
    d = length(s)
    z = AcbRefVector(d)
    s0 = Acb(s[0])

    # From _acb_poly_zeta_cpx_series
    bound = Mag()
    vb = ArbRefVector(d)
    prec = precision(s)
    bound_prec = 40 + prec ÷ 20
    # N and M are given as pointers
    N_ptr = UInt[0]
    M_ptr = UInt[0]

    Arblib.zeta_em_choose_param!(
        bound,
        N_ptr,
        M_ptr,
        s0,
        Acb(a),
        min(d, 2),
        prec,
        bound_prec,
    )
    N = only(N_ptr)
    M = only(M_ptr)

    Arblib.zeta_em_bound!(vb, s0, Acb(a), N, M, d, bound_prec)


    # Implementation of _acb_poly_zeta_em_sum
    let prec = prec + 2 * (bit_count(N) + bit_count(d))
        t = AcbRefVector(d + 1; prec)
        u = AcbRefVector(d; prec)
        v = AcbRefVector(d; prec)
        term = AcbRefVector(d; prec)
        sum = AcbRefVector(d; prec)
        Na = Acb(a + N; prec)

        # sum 1/(k+a)^(s+x)
        # IMPROVE: Do not only use naive sum?
        Arblib.powsum_series_naive!(sum, s0, Acb(a), one(s0), N, d, prec)

        # t = 1/(N+a)^(s+x); take one extra term for deflation
        _acb_poly_acb_invpow_cpx!(t, Na, s0, d + 1, prec)

        # This is the updated part compared to Arb
        begin
            # Compute coefficients of ((N + a) * t - 1) / (s - 1)
            w = let s = ArbSeries((s[0], 1), degree = d - 1)
                # (N + a) * inv(N + a)^s - 1 = inv(N + a)^(s - 1) - 1
                num(v) = inv(real(Na))^v - 1

                if contains(real(s0), 1)
                    w = fx_div_x(num, s - 1)
                else
                    w = num(s - 1) / (s - 1)
                end
                Arblib.coeffs(w)
            end

            # Compute sum += ((N + a) * t - 1) / (s - 1)
            sum .+= w

            _acb_poly_acb_invpow_cpx!(t, Na, s0, d, prec)
        end

        # sum += u = 1/2 * t
        Arblib.mul_2exp!(u, t, d, -1)
        Arblib.add!(sum, sum, u)

        # Euler-Maclaurin formula tail
        if d < 5 || d < M ÷ 10
            Arblib.zeta_em_tail_naive!(u, s0, Na, t, M, d)
        else
            Arblib.zeta_em_tail_bsplit!(u, s0, Na, t, M, d)
        end

        Arblib.add!(z, sum, u)
    end

    for i = 1:d
        Arblib.add_error!(Arblib.realref(z[i]), vb[i])
    end

    # This corresponds to _acb_poly_zeta_series
    @assert all(isreal, z)
    res = ArbSeries(real.(z))

    return ArbExtras.compose_zero(res, s)
end

_zeta_deflated(s::Arb, a::Arb) = _zeta_deflated(ArbSeries(s), a)[0]

"""
    zeta_deflated(s, a)

Compute the deflated zeta function.

The deflated zeta function is defined by
```
zeta_deflated(s, a) = zeta(s, a) + 1 / (1 - s)
```

Arb doesn't support evaluation for `s` overlapping one. For this we
use the implementation in [`_zeta_deflated`](@ref). While it would be
possible to use that implementation for all values of `s` we opt to
only use it when `s` overlaps one but is not exactly one.
"""
function zeta_deflated(s::Arb, a::Arb)
    contains(s, 1) && !isone(s) && return _zeta_deflated(s, a)

    res = ArbSeries(s)
    Arblib.zeta_series!(res, res, a, 1, length(res))
    return res[0]
end

function zeta_deflated(s::ArbSeries, a::Arb)
    s0 = Arblib.ref(s, 0)
    contains(s0, 1) && !isone(s0) && return _zeta_deflated(s, a)

    return Arblib.zeta_series!(zero(s), s, a, 1, length(s))
end

"""
    lerch_phi(z::Arb, s::Arb, a::Arb)

Compute the Lerch transcendent
```
lerch_phi(z, s, a) = sum(z^n / (n + a)^s for n = 0:Inf)
```
Returns an indeterminate value if the result is not real.
"""
function lerch_phi(z::Arb, s::Arb, a::Arb)
    res = Arblib.dirichlet_lerch_phi!(Acb(prec = precision(z)), Acb(z), Acb(s), Acb(a))

    if isreal(res)
        return real(res)
    else
        return indeterminate(z)
    end
end

function abspow!(res::Arb, x::Arblib.ArbOrRef, y::Arb)
    iszero(y) && return Arblib.one!(res)

    if iszero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(res)
        Arblib.ispositive(y) && return Arblib.zero!(res)
        return Arblib.unit_interval!(res)
    end

    if Arblib.contains_zero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(res)
        upper = abs_ubound(Arb, x) # One extra allocation
        Arblib.pow!(upper, upper, y)
        Arblib.zero!(res)
        return Arblib.union!(res, res, upper)
    end

    if res === y
        # In this case we need an extra allocation to not overwrite y
        y = copy(y)
    end
    Arblib.abs!(res, x)
    return Arblib.pow!(res, res, y)
end

function abspow!(res::ArbSeries, x::ArbSeries, y::Arb)
    Arblib.degree(res) == Arblib.degree(x) ||
        throw(ArgumentError("res and x should have the same degree"))

    sgn = Arblib.sgn_nonzero(Arblib.ref(x, 0))

    if sgn == 0
        # We don't have to be that careful with allocations here.

        # All non-constant terms are indeterminate, the constant term
        # is given by abs(x[0])^y
        res[0] = abspow(x[0], y)
        for i = 1:Arblib.degree(res)
            res[i] = indeterminate(Arb)
        end

        return res
    elseif sgn < 0
        Arblib.neg!(res, x)
        Arblib.pow_arb_series!(res, res, y, length(res))
    else
        Arblib.pow_arb_series!(res, x, y, length(res))
    end

    return res
end


"""
    abspow(x, y)

Compute `abs(x)^y `in a way that works if `x` overlaps with zero.

For complex `x` is complex it doesn't use `abs(x)` but `+x` in the
right half plane and `-x` in the left half plane. For zero we are
conservative and return an indeterminate result. This means that for
non-zero `x` it represents an analytic continuation of `abs(x)^y` and
can thus be used in [`Arblib.integrate`](@ref).
"""
function abspow(x::Arb, y::Arb)
    iszero(y) && return one(x)

    if iszero(x)
        Arblib.contains_negative(y) && return indeterminate(x)
        Arblib.ispositive(y) && return zero(x)
        return Arblib.unit_interval!(zero(x))
    end

    if Arblib.contains_zero(x)
        Arblib.contains_negative(y) && return indeterminate(x)
        x_upp = Arblib.abs_ubound(Arb, x)
        return Arb((zero(x), x_upp^y))
    end

    res = abs(x)
    return Arblib.pow!(res, res, y)
end

function abspow(x::ArbSeries, y::Arb)
    if Arblib.contains_zero(Arblib.ref(x, 0))
        # All non-constant terms are indeterminate, the constant term
        # is given by abs(x[0])^y
        res = indeterminate(x)
        res[0] = abspow(x[0], y)
        return res
    end

    res = abs(x)
    return Arblib.pow_arb_series!(res, res, y, length(res))
end

function abspow(x::Arb, y::ArbSeries)
    # This function is only partially implemented. In the case when x
    # overlaps zero it is currently based on differentiation and
    # isolation of extrema done by hand. If we want to support much
    # higher degrees we would need to do this algorithmically. In the
    # current version it is also not optimized at all, so it could be
    # much faster.

    iszero(y) && return one(y)

    if iszero(x)
        # If y[0] = 0 the constant term is zero but not the others. We
        # therefore need y[0] > 0 to get exactly zero.
        Arblib.contains_nonpositive(Arblib.ref(y, 0)) && return indeterminate(y)

        return zero(y)
    end

    if Arblib.contains_zero(x)
        Arblib.contains_nonpositive(Arblib.ref(y, 0)) && return indeterminate(y)

        # Differentiate with respect to the parameter of y manually
        # and enclose the terms

        deg = Arblib.degree(y)

        deg <= 3 || error("supports degree at most 3")

        y0 = y[0]

        res = zero(y)

        res[0] = abspow(x, y[0])
        if deg >= 1
            # res[1] = y[1] * log(x) * abspow(x, y[0])

            # Compute enclosure of log(x) * abspow(x, y[0])
            f1(x) = logabspow(x, 1, y0)

            # Evaluate at endpoints
            term = union(f1(zero(x)), f1(ubound(Arb, x)))

            # Evaluate at possible extrema
            # IMPROVE: If y[0] is extremely small and somewhat wide
            # then the enclosure of extrema might overlap zero due to
            # radius having low precision.
            extrema = exp(-1 / y0)
            if Arblib.overlaps(x, extrema)
                Arblib.union!(term, term, f1(extrema))
            end

            res[1] = y[1] * term
        end
        if deg >= 2
            # res[2] = (2y[2] * log(x) + (y[1] * log(x))^2) * abspow(x, y[0])
            #        = 2y[2] * logabspow(x, 1, y[0]) + y[1]^2 * logabspow(x, 2, y[0])
            y1, y2 = y[1], y[2]

            f2(x) = 2y2 * logabspow(x, 1, y0) + y1^2 * logabspow(x, 2, y0)

            # Evaluate at endpoints
            term = union(f2(zero(x)), f2(ubound(Arb, x)))

            # Evaluate at possible extrema
            if iszero(y[1])
                # The extrema simplifies
                extrema = exp(-1 / y0)
                if Arblib.overlaps(x, extrema)
                    term = union(term, f2(extrema))
                end
            else
                Δ = (2y1 + y0 * 2y1)^2 - 4y0 * y1 * 2y2
                if Arblib.contains_nonnegative(Δ) # Otherwise there are no real roots
                    A = -(2y1^2 + y0 * 2y2)
                    B = 2y0 * y1^2
                    sqrt_Δ = Arblib.sqrtpos(Δ)
                    extrema1 = exp((A + sqrt_Δ) / B)
                    extrema2 = exp((A - sqrt_Δ) / B)
                    if Arblib.overlaps(x, extrema1)
                        Arblib.union!(term, term, f2(extrema1))
                    end
                    if Arblib.overlaps(x, extrema2)
                        Arblib.union!(term, term, f2(extrema2))
                    end
                end
            end

            res[2] = term / 2
        end

        return res
    end

    return abs(x)^y
end

abspow(x::Complex, y) =
    if real(x) > 0
        return x^y
    elseif real(x) < 0
        return (-x)^y
    else
        # Simple way to return NaN + im * NaN of correct type
        return x * (NaN + im * NaN)
    end

abspow(x::Acb, y) =
    if Arblib.ispositive(real(x))
        return x^y
    elseif Arblib.isnegative(real(x))
        return (-x)^y
    else
        return indeterminate(x)
    end

abspow(x, y) = abs(x)^y

"""
    logabspow(x, i::Integer, y)

Compute `log(abs(x))^i * abs(x)^y` in a way that works for `x`
overlapping zero.

For `x` overlapping zero it uses that for `x > 0` and `y != 0` the
unique critical point of `log(x)^i * x^y` is given by `exp(-i / y)`.
It is hence enough to evaluate at the endpoints of the interval as
well as possibly this point.
"""
logabspow(x, i::Integer, y) = iszero(i) ? abspow(x, y) : log(abs(x))^i * abspow(x, y)

function logabspow(x::Arb, i::Integer, y::Arb)
    iszero(i) && return abspow(x, y)

    if Arblib.contains_zero(x)
        if Arblib.ispositive(y)
            iszero(x) && return zero(x)

            # Evaluate at endpoints of x
            xᵤ = abs_ubound(Arb, x)
            res = union(zero(x), log(xᵤ)^i * xᵤ^y)

            # Check if critical point is contained in x, if so evaluate on
            # it
            critical_point = exp(-i / y)
            if Arblib.overlaps(x, critical_point)
                res = union(res, log(critical_point)^i * critical_point^y)
            end

            return res
        elseif iszero(y) && i < 0
            iszero(x) && return zero(x)

            xᵤ = abs_ubound(Arb, x)

            # log(1) = 0 so we get an indeterminate value there
            xᵤ < 1 || return indeterminate(x)

            # Monotone for 0 < x < 1, evaluate on endpoints
            return union(zero(x), log(xᵤ)^i)
        else
            # Non-zero at x = 0
            return indeterminate(x)
        end
    end

    return log(abs(x))^i * abspow(x, y)
end

function logabspow(x::ArbSeries, i::Integer, y::Arb)
    iszero(i) && return abspow(x, y)

    if Arblib.contains_zero(x[0])
        # Return an indeterminate result except for the constant term
        res = indeterminate(x)
        res[0] = logabspow(x[0], i, y)

        return res
    end

    return log(abs(x))^i * abspow(x, y)
end

"""
    x_pow_s_x_pow_t_m1_div_t(x::Arb, s::Arb, t::Arb)

Compute an enclosure of
```
abs(x)^s * (abs(x)^t - 1) / t
```
in a way that works well for `t` overlapping zero.

To get a good enclosure it uses monotonicity in `t`. The derivative in
`t` can be written as
```
abs(x)^s * (1 + (log(abs(x)) * t - 1) * exp(log(abs(x)) * t)) / t^2
```
This is either zero or the sign depends only on `1 + (log(abs(x)) * t
- 1) * exp(log(abs(x)) * t)´. If we let `v = log(abs(x)) * t` we have
to study `1 + (v - 1) * exp(v)` which has the unique root `v = 0` and
is positive for all other values of `v`. It is hence non-decreasing in
`t` and we can evaluate on the endpoints.

For `t = 0` it reduces to `abs(x)^s * log(abs(x))`

For `x` overlapping zero we write it as
```
(abspow(x, s + t) - abspow(x, s)) / t
```
Which we can evaluate directly if `t` doesn't overlap zero, otherwise
using [`fx_div_x`](@ref) to handle the removable singularity.

**IMPROVE:** It doesn't handle very small `s` and `x` overlapping zero
very well. This is however hard to do in general since the function
grows extremely quickly in `x` near `x = 0` in that case. It would be
possibly to improve the current version by checking for monotonicity
in `x`, however we might not need to do that.
"""
function x_pow_s_x_pow_t_m1_div_t(x::Arblib.ArbOrRef, s::Arb, t::Arb)
    iszero(t) && return logabspow(x, 1, s)

    if Arblib.contains_zero(t) || iswide(t)
        # Use that it is non-decreasing in t
        if Arblib.contains_zero(x)
            # Avoid reimplementing this special case
            tₗ, tᵤ = getinterval(Arb, t)
            return Arb((
                x_pow_s_x_pow_t_m1_div_t(x, s, tₗ),
                x_pow_s_x_pow_t_m1_div_t(x, s, tᵤ),
            ))
        else
            # Manually implement this case for performance reasons

            # resₗ = resᵤ = log(abs(x))
            resₗ = abs(x)
            Arblib.log!(resₗ, resₗ)
            resᵤ = copy(resₗ)

            t0 = lbound(Arb, t)
            Arblib.mul!(resₗ, resₗ, t0)
            Arblib.expm1!(resₗ, resₗ)
            Arblib.div!(resₗ, resₗ, t0)

            Arblib.get_ubound!(Arblib.midref(t0), t)
            Arblib.mul!(resᵤ, resᵤ, t0)
            Arblib.expm1!(resᵤ, resᵤ)
            Arblib.div!(resᵤ, resᵤ, t0)

            Arblib.union!(resₗ, resₗ, resᵤ)
            abspow!(resᵤ, x, s)
            return Arblib.mul!(resₗ, resₗ, resᵤ)
        end
    end

    if Arblib.contains_zero(x)
        x = convert(Arb, x) # This case doesn't handle ArbRef
        # We don't work as hard to get a good enclosure or good
        # performance in this case
        if Arblib.contains_zero(t)
            # Handle removable singularity
            res = fx_div_x(t, force = true) do t
                abspow(x, s + t) - abspow(x, s)
            end

            return res
        else
            # Evaluate directly and with the removable zero added and
            # take the intersection
            res1 = (abspow(x, s + t) - abspow(x, s)) / t
            res2 = fx_div_x(union(t, zero(t)), force = true) do t
                abspow(x, s + t) - abspow(x, s)
            end
            if !isfinite(res1)
                return res2
            elseif !isfinite(res2)
                return res1
            else
                return intersect(res1, res2)
            end
        end
    end

    res = abs(x)
    Arblib.log!(res, res)
    Arblib.mul!(res, res, t)
    Arblib.expm1!(res, res)
    Arblib.div!(res, res, t)
    return Arblib.mul!(res, res, abspow(x, s))
end

"""
    x_pow_s_x_pow_t_m1_div_t(x::ArbSeries, s::Arb, t::Arb)

Compute an enclosure of
```
abs(x)^s * (abs(x)^t - 1) / t
```
in a way that works well for `t` overlapping zero.

In case `x` overlaps the constant term is computed using
`x_pow_t_div_t(x[0], s, t)` and all higher order terms are set to an
indeterminate value, even if they in some cases could be finite.

To compute the series for non-zero `x` we first compute the series for
`(abs(x)^t - 1) / t` and then multiply with the one for `abs(x)^s`. To
get the series for `(abs(x)^t - 1) / t` we first compute the
derivative, then integrate and set the constant according to
`x_pow_t_div_t(x[0], 0, t)`. The derivative is given by
```
abs(x)' * abs(x)^(t - 1)
```
where `x'` denotes the derivative of `x`.
"""
function x_pow_s_x_pow_t_m1_div_t(x::ArbSeries, s::Arb, t::Arb)
    if Arblib.contains_zero(Arblib.ref(x, 0))
        res = indeterminate(x)
        res[0] = x_pow_s_x_pow_t_m1_div_t(Arblib.ref(x, 0), s, t)
        return res
    end

    # Only constant term
    iszero(Arblib.degree(x)) &&
        return ArbSeries(x_pow_s_x_pow_t_m1_div_t(Arblib.ref(x, 0), s, t))

    res = let absx = abs(x)
        # Compute x' * x^(t - 1) with one degree lower
        dres = let tmp1 = ArbSeries(absx, degree = Arblib.degree(x) - 1)
            xtm1 = abspow(tmp1, t - 1)
            Arblib.derivative!(tmp1, absx)
            Arblib.mullow!(tmp1, tmp1, xtm1, length(tmp1))
        end

        # Integrate
        res = Arblib.integral!(absx, dres) # Reuse absx
    end

    # Set the constant of integration correctly
    res[0] = x_pow_s_x_pow_t_m1_div_t(Arblib.ref(x, 0), zero(s), t)

    # Multiply with series for abs(x)^s
    return Arblib.mullow!(res, res, abspow(x, s), length(res))
end
