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
hypgeom_2f1(a::Acb, b::Acb, c::Acb, z::Acb) = Arblib.hypgeom_2f1!(zero(z), a, b, c, z, 0)

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
rising(x::Arb, n::Integer) = Arblib.rising!(zero(x), x, unsigned(n))
rising(x::ArbSeries, n::Integer) =
    Arblib.rising_ui_series!(zero(x), x, unsigned(n), length(x))

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
    if !Arblib.ispositive(a)
        # Return an indeterminate result
        res = zero(s)
        for i = 0:Arblib.degree(s)
            Arblib.indeterminate!(Arblib.ref(res, i))
        end
        # Since we manually set the coefficients of the polynomial we
        # need to also manually set the degree.
        res.arb_poly.length = length(s)
        return res
    end

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

    # Compose with non-constant part
    s_tmp = copy(s)
    s_tmp[0] = 0

    return Arblib.compose(res, s_tmp)
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
    lerch_phi(z, s, a)

Compute a naive enclosure of the Lerch transendent
```
lerch_phi(z, s, a) = sum(z^n / (n + a)^s for n = 0:Inf)
```

It only supports `abs(z) < 1`, `s < 0` and `a > 0`. We treat the whole
series as a tail, for bounding it see
[https://fredrikj.net/blog/2022/02/computing-the-lerch-transcendent/](Fredrik
Johanssons blog).

Note that this is a very naive implementation and only intended for
use in [`clausenc_expansion_remainder`](@ref). Once the version of Arb
that implements this function is released we should switch to using
that one.
"""
function lerch_phi(z::Arb, s::Arb, a::Arb)
    @assert abs(z) < 1
    @assert Arblib.isnonpositive(s)
    @assert Arblib.ispositive(a)

    # Find C
    C = exp(-s / a)

    abs(C * z) < 1 || throw(ArgumentError("C to large to satisfy C * z < 1"))

    return a^-s / (1 - C * abs(z))
end
