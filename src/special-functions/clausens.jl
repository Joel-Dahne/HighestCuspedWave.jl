export clausens

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

"""
    clausens_expansion(x, s, M::Integer)

Compute the asymptotic expansion of `clausens(x, s)` at zero up to
order `2M - 1`, meaning that the error term is of order `2M`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a `ArbSeries` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
enciso18:convex_whith.

It satisfies that `clausens(y, s) ∈ C*abs(y)^e + P(y) + E*y^(2M)` for
all `|y| <= |x|`.

Note that this method doesn't handle wide values of `s` in any special
way. This has not been needed anywhere so far.
"""
function clausens_expansion(x::Arb, s::Arb, M::Integer)
    Arblib.ispositive(s) || throw(ArgumentError("s must be positive"))
    # TODO: Check this
    M > (s + 1) / 2 || throw(ArgumentError("M must be larger that (s + 1) / 2"))

    π = oftype(x, pi)

    # Non-analytic term
    C = SpecialFunctions.gamma(1 - s) * cospi(s / 2)
    e = s - 1

    # Analytic term
    P = ArbSeries(degree = 2M - 1, prec = precision(x))
    for m = 0:M-1
        P[2m + 1] = (-1)^m * zeta(s - 2m - 1) / factorial(2m + 1)
    end

    # Error term
    E = Arblib.add_error!(zero(x), 2(2π)^(s - 2M) * zeta(2M + 2 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end
