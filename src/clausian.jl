# Contains methods related to Ci and Si
export Li, Ci, Si, Ci_tilde

###
### Li
###

"""
    Li(z, s)
Compute the polylogarithm ``Li_s(z)``
"""
function Li(z::Acb, s::Union{Acb,Integer})
    if iswide(s) # If this is true then s is always an Acb
        # TODO: Check that everything here is correct

        # TODO: Tune this
        n = 3 # Degree of Taylor expansion

        s_mid = Acb(Arblib.midref(real(s)), Arblib.midref(imag(s)))

        # Compute the rest term of the Taylor expansion
        w = Arblib.polylog_series!(
            AcbSeries(degree = n, prec = precision(z)),
            AcbSeries([s, 1]),
            z,
            n + 1,
        )

        restterm = (s - s_mid)^n * w[n]

        # Compute the Taylor polynomial at the midpoint of x
        w = Arblib.polylog_series!(
            AcbSeries(degree = n, prec = precision(z)),
            AcbSeries([s_mid, 1]),
            z,
            n + 1,
        )

        # Evaluate the Taylor polynomial on s - s_mid and add the rest
        # term
        return w(s - s_mid) + restterm
    end

    return Arblib.polylog!(zero(z), s, z)
end

Li(z::acb, s::acb) = parent(z)(Li(Acb(z), Acb(s)))

Li(z::acb, s::Integer) = parent(z)(Li(Acb(z), s))

"""
    Li(z, s, β)
Compute the polylogarithm ``Li_s^{(β)}(z)``.

That is, `Li(z, s)` differentiated `β` times w.r.t. `s` evaluated at `z`.
"""
function Li(z::Acb, s::Acb, β::Integer)
    s_poly = AcbSeries([s, 1])
    w = Arblib.polylog_series!(AcbSeries(degree = β, prec = precision(z)), s_poly, z, β + 1)

    return w[β] * factorial(β)
end

Li(z::acb, s::acb, β::Integer) = parent(z)(Li(Acb(z), Acb(s), β))

###
### Ci
###

"""
    Ci(x, s)
Compute the Clausian function Ciₛ(x).

If x is a wide (real) ball (as determined by iswide(x)) it computes a
tighter enclosure by using that `Ci` is 2π periodic, monotonic for `x ∈
[0, π]` and even, so that it's enough to evaluate on the endpoints and
possibly at zero or π if `x` contains points on the form `2kπ` or (2k
+ 1)π` respectively.
"""
Ci(x::Acb, s) = (Li(exp(im * x), s) + Li(exp(-im * x), s)) / 2
Ci(x::acb, s::acb) = parent(x)(Ci(Acb(x), Acb(s)))
Ci(x::acb, s::Integer) = parent(x)(Ci(Acb(x), s))

function Ci(x::Arb, s::Union{Arb,Integer})
    if iswide(x)
        xₗ, xᵤ = Arblib.getinterval(Arb, x)
        (include_zero, include_pi) = contains_pi(xₗ, xᵤ)
        res = union(Ci(xₗ, s), Ci(xᵤ, s))
        if include_zero
            res = union(res, Ci(zero(x), s))
        end
        if include_pi
            res = union(res, Ci(Arb(π), s))
        end
        return res
    end
    s = s isa Integer ? s : Acb(s)
    return real(Li(exp(Acb(0, x)), s))
end

Ci(x::arb, s::arb) = parent(x)(Ci(Arb(x), Arb(s)))
Ci(x::arb, s::Integer) = parent(x)(Ci(Arb(x), s))
Ci(x::S, s::T) where {S<:Real,T<:Real} = convert(
    float(promote_type(S, T)),
    Ci(convert(Arb, x), s isa Integer ? s : convert(Arb, s)),
)

"""
    Ci(x, s, β)
Compute the Clausian function Ciₛ^(β)(x).

That is, `Ciₛ` differentiated `β` times w.r.t. `s` evaluated at `x`.

If x is a wide (real) ball (as determined by iswide(x)) it computes a
tighter enclosure by using the monotonicity properties of `Ci`.
Currently this is only implemented for `β = 1`, `0 < x < 2π` and `s =
2` or `s = 3`. The function has a critical point at `x = π`, one on
the interval `0 < x < π` and one (due to being even around `π`) on `π
< x < 2π`. The extrema can occur either on one of these critical
points or on the endpoints of the ball. For efficiency reasons the
critical point on `0 < x < π` is precomputed for `s = 2` and `s = 3`
(the one on `π < x < 2π` is given by symmetry).

PROVE: That there is only once critical point.
TODO: Use that `Ci` is 2π periodic to allow for `x` outside `[0, 2π]`.
This might not be needed in the end though.
"""
Ci(x::Acb, s::Acb, β::Integer) = (Li(exp(im * x), s, β) + Li(exp(-im * x), s, β)) / 2
Ci(x::acb, s::acb, β::Integer) = parent(x)(Ci(Acb(x), Acb(s), β))

function Ci(x::Arb, s::Arb, β::Integer)
    iszero(x) && s > 1 && return zeta(s, d = β)

    if false && iswide(x) && β == 1 && 0 < x < 2Arb(π) && (s == 2 || s == 3)
        xₗ, xᵤ = Arblib.getinterval(Arb, x)

        res = union(Ci(xₗ, s, β), Ci(xᵤ, s, β))

        if s == 2
            critical_point = Arb("[1.010782703526315549251222370194235400 +/- 7.10e-37]")
        elseif s == 3
            critical_point = Arb("[1.219556773337345811161114646108970 +/- 5.13e-34]")
        end
        if Arblib.overlaps(x, critical_point) ||
           Arblib.overlaps(x, 2Arb(π) - critical_point)
            # Depending on the precision critical_point might count as
            # wide so we explicitly call Li to avoid infinite
            # recursion.
            res = union(res, real(Li(exp(Acb(0, critical_point)), convert(Acb, s), β)))
        end

        if Arblib.overlaps(x, Arb(π))
            res = union(res, Ci(Arb(π), s, β))
        end

        return res
    end

    return real(Li(exp(Acb(0, x)), convert(Acb, s), β))
end
Ci(x::arb, s::arb, β::Integer) = parent(x)(Ci(Arb(x), Arb(s), β))
Ci(x::S, s::T, β::Integer) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), Ci(convert(Arb, x), convert(Arb, s), β))

"""
    Ci(x::Union{ArbSeries,arb_series}, s)

Compute the Taylor series of Ciₛ(x).

It's computed by directly computing the Taylor coefficients by
differentiating Ciₛ and then composing with `x`.
"""
function Ci(x::ArbSeries, s)
    res = zero(x)
    x₀ = x[0]

    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * Ci(x₀, s - i) / factorial(i)
        else
            res[i] = -(-1)^(i ÷ 2) * Si(x₀, s - i) / factorial(i)
        end
    end

    # Compose the Taylor series for the Clausian with that of the
    # input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

function Ci(x::arb_series, s)
    res = arb_series(parent(x.poly)(), length(x))
    x₀ = x[0]

    for i = 0:length(x)-1
        if i % 2 == 0
            res[i] = (-1)^(div(i, 2)) * Ci(x₀, s - i) / factorial(i)
        else
            res[i] = -(-1)^(div(i, 2)) * Si(x₀, s - i) / factorial(i)
        end
    end

    # Compose the Taylor series for the Clausian with that of the
    # input
    x_tmp = arb_series(deepcopy(x.poly))
    x_tmp[0] = base_ring(parent(x.poly))(0)

    return Nemo.compose(res, x_tmp)
end

"""
    Ci(x::Union{ArbSeries,arb_series}, s, β)

Compute the Taylor series of Ciₛ(x).

It's computed by directly computing the Taylor coefficients by
differentiating Ciₛ and then composing with `x`.
"""
function Ci(x::ArbSeries, s, β::Integer)
    res = zero(x)
    x₀ = x[0]

    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * Ci(x₀, s - i, β) / factorial(i)
        else
            res[i] = -(-1)^(i ÷ 2) * Si(x₀, s - i, β) / factorial(i)
        end
    end

    # Compose the Taylor series for the Clausian with that of the
    # input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

function Ci(x::arb_series, s, β::Integer)
    res = arb_series(parent(x.poly)(), length(x))
    x₀ = x[0]

    for i = 0:length(x)-1
        if i % 2 == 0
            res[i] = (-1)^(div(i, 2)) * Ci(x₀, s - i, β) / factorial(i)
        else
            res[i] = -(-1)^(div(i, 2)) * Si(x₀, s - i, β) / factorial(i)
        end
    end

    # Compose the Taylor series for the Clausian with that of the
    # input
    x_tmp = arb_series(deepcopy(x.poly))
    x_tmp[0] = base_ring(parent(x.poly))(0)

    return Nemo.compose(res, x_tmp)
end

"""
    Ci_expansion(x, s, M::Integer)
Compute the asymptotic expansion of `Ciₛ(x)` at zero up to order `2M -
2`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a series `P` (`arb_series` or
`ArbSeries`) and the error term `E`. The `M` is the same as in Lemma
2.1 in enciso18:convex_whith.

It satisfies that `Ci(y, s) ∈ C*abs(y)^e + P(y) + E*y^(2M)` for all
`|y| <= |x|`.
"""
function Ci_expansion(x::arb, s::arb, M::Integer)
    π = parent(x)(pi)

    # Non-analytic term
    C = Nemo.gamma(1 - s) * sinpi(s / 2)
    e = s - 1

    # Analytic term
    P = arb_series(ArbPolyRing(parent(x), :x)(), 2M - 1)
    for m = 0:M-1
        P[2m] = (-1)^m * zeta(s - 2m) / factorial(fmpz(2m))
    end

    # Error term
    E = ball(zero(x), 2(2π)^(1 + s - 2M) * zeta(2M + 1 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

function Ci_expansion(x::Arb, s::Arb, M::Integer)
    Arblib.ispositive(s) || throw(ArgumentError("s must be positive"))
    # TODO: Check this
    M > (s + 1) / 2 || throw(ArgumentError("M must be larger that (s + 1) / 2"))

    π = oftype(x, pi)

    # Non-analytic term
    if s == 2
        C = -π / 2
    else
        C = SpecialFunctions.gamma(1 - s) * sinpi(s / 2)
    end
    e = s - 1

    # Analytic term
    P = ArbSeries(degree = 2M - 2, prec = precision(x))
    for m = 0:M-1
        P[2m] = (-1)^m * zeta(s - 2m) / factorial(oftype(x, 2m))
    end

    # Error term
    E = Arblib.add_error!(zero(x), 2(2π)^(1 + s - 2M) * zeta(2M + 1 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

"""
    Ci_expansion(x, s, β, M::Integer)
Compute the asymptotic expansion of `Ciₛ`(β)(x)` at zero up to order
`2M - 2`.

It currently only supports `β == 1`, for `β == 0` use `Ci_expansion(x,
s)` instead.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a series `P` (`arb_series` or
`ArbSeries`) and the error term `E`. The `M` is the same as in Lemma
2.1 in enciso18:convex_whith.

It satisfies that `Ci(y, s, β) ∈ C*abs(y)^e + P(y) + E*y^(2M)` for all
`|y| <= |x|`.

TODO: Implement these
"""
function Ci_expansion(x::Arb, s::Arb, β::Integer, M::Integer)
    β == 1 || throw(ArgumentError("only implemented for β = 1, got $β"))

end

###
### Si
###

"""
    Si(x, s)
Compute the Clausian function Siₛ(x).

If x is a wide (real) ball (as determined by iswide(x)) it computes a
tighter enclosure by first checking if the derivative doesn't contains
zero, if not it uses monotonicity to only evaluate at endpoints. If
the derivative does contain zero it uses a zero order approximation
instead.
"""
Si(x::Acb, s) = (Li(exp(im * x), s) - Li(exp(-im * x), s)) / 2
Si(x::acb, s::acb) = parent(x)(Ci(Acb(x), Acb(s)))
Si(x::acb, s::Integer) = parent(x)(Ci(Acb(x), s))

function Si(x::Arb, s::Union{Arb,Integer})
    if iswide(x) # If this is true then s is always an Arb
        # Compute derivative
        dSi = Ci(x, s - 1)
        if Arblib.contains_zero(dSi)
            # Use a zero order approximation
            return Arblib.add_error!(
                Si(Arblib.midpoint(Arb, x), s),
                (x - Arblib.midpoint(Arb, x)) * dSi,
            )
        else
            # Use that it's monotone
            xₗ, xᵤ = Arblib.getinterval(Arb, x)
            return union(Si(xₗ, s), Si(xᵤ, s))
        end
    end
    s = s isa Integer ? s : Acb(s)
    return imag(Li(exp(Acb(0, x)), s))
end

Si(x::arb, s::arb) = parent(x)(Si(Arb(x), Arb(s)))
Si(x::arb, s::Integer) = parent(x)(Si(Arb(x), s))
Si(x::S, s::T) where {S<:Real,T<:Real} = convert(
    float(promote_type(S, T)),
    Si(convert(Arb, x), s isa Integer ? s : convert(Arb, s)),
)

"""
    Si(x, s, β)
Compute the Clausian function Siₛ^(β)(x).

That is, `Siₛ` differentiated `β` times w.r.t. `s` evaluated at `x`.

TODO: Handle wide (real) balls better similar to how `Si(x, s)` does
it.
"""
Si(x::Acb, s::Acb, β::Integer) = (Li(exp(im * x), s, β) - Li(exp(-im * x), s, β)) / 2
Si(x::acb, s::acb, β::Integer) = parent(x)(Si(Acb(x), Acb(s), β))

Si(x::Arb, s::Arb, β::Integer) = imag(Li(exp(Acb(0, x)), convert(Acb, s), β))
Si(x::arb, s::arb, β::Integer) = parent(x)(Si(Arb(x), Arb(s), β))
Si(x::S, s::T, β::Integer) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), Si(convert(Arb, x), convert(Arb, s), β))

"""
    Si_expansion(x, s, M::Integer)
Compute the asymptotic expansion of `Si` at zero up to order `2M - 1`.

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms in an `arb_series` `P` and the
error term `E`. The `M` is the same as in Lemma 2.1 in
enciso18:convex_whith.

It satisfies that `Si(y, s) ∈ C*sign(y)*abs(y)^e + P(y) + E*abs(y)^(2M + 1)`
for all `|y| <= |x|`.
"""
function Si_expansion(x::arb, s::arb, M::Integer)
    π = parent(x)(pi)

    # Non-analytic term
    C = Nemo.gamma(1 - s) * cospi(s / 2)
    e = s - 1

    # Analytic term
    P = arb_series(ArbPolyRing(parent(x), :x)(), 2M)
    for m = 0:M-1
        P[2m+1] = (-1)^m * zeta(s - 2m - 1) / factorial(fmpz(2m + 1))
    end

    # Error term
    E = ball(zero(x), 2(2π)^(s - 2M) * zeta(2M + 2 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end

function Si_expansion(x::Arb, s::Arb, M::Integer)
    π = Arb(pi)

    # Non-analytic term
    C = SpecialFunctions.gamma(1 - s) * cospi(s / 2)
    e = s - 1

    # Analytic term
    P = ArbSeries(degree = 2M - 1)
    for m = 0:M-1
        P[2m+1] = (-1)^m * zeta(s - 2m - 1) / factorial(2m + 1)
    end

    # Error term
    E = Arblib.add_error!(zero(x), 2(2π)^(s - 2M) * zeta(2M + 2 - s) / (4π^2 - x^2))

    return (C, e, P, E)
end
