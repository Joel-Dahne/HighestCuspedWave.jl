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
    zeta_deflated(s, a)

Compute the deflated zeta function.

The deflated zeta function is defined by
```
zeta_deflated(s, a) = zeta(s, a) + 1 / (1 - s)
```

**TODO:** Arb doesn't support evaluation for `s` overlapping one.
Since this is the main case where we need this function we currently
implement a non-rigorous enclosure for that case. This should be
updated with a rigorous enclosure.
"""
function zeta_deflated(s::Arb, a::Arb)
    if contains(s, 1) && !isone(s)
        @warn "Non-rigorous enclosure of zeta_deflated" s maxlog = 1

        # Assume monotonicity
        sₗ, sᵤ = getinterval(Arb, s)
        return Arb((zeta_deflated(sₗ, a), zeta_deflated(sᵤ, a)))
    end

    res = ArbSeries(s)
    Arblib.zeta_series!(res, res, a, 1, length(res))
    return res[0]
end

function zeta_deflated(s::ArbSeries, a::Arb)
    s0 = Arblib.ref(s, 0)
    if contains(s0, 1) && !isone(s0)
        @warn "Non-rigorous enclosure of zeta_deflated" s maxlog = 1

        # Assume monotonicity
        sₗ, sᵤ = copy(s), copy(s)
        s0ₗ, s0ᵤ = getinterval(Arb, s0)
        # If s0ₗ or s0ᵤ is very close to 1 we get very bad enclosures.
        # If one of them is very close take it to have the same
        # distance as the other one, this still would given an
        # enclosure but is obviously larger than it needs to be.
        if 1 - s0ₗ < Arb("1e-10")
            s0ₗ = 1 - (s0ᵤ - 1)
        end
        if s0ᵤ - 1 < Arb("1e-10")
            s0ᵤ = 1 + (1 - s0ₗ)
        end
        sₗ[0], sᵤ[0] = s0ₗ, s0ᵤ

        res1 = zeta_deflated(sₗ, a)
        res2 = zeta_deflated(sᵤ, a)
        return ArbSeries(union.(Arblib.coeffs(res1), Arblib.coeffs(res2)))
    end

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
