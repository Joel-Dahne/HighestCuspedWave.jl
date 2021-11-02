export iswide

"""
    mince(x::Arb, n::Integer)

Return a vector with `n` balls covering the ball `x`.
"""
function mince(x::Arb, n::Integer)
    balls = Vector{Arb}(undef, n)
    xₗ, xᵤ = Arblib.getinterval(Arb, x)
    dx = (xᵤ - xₗ) / n
    for i in eachindex(balls)
        yₗ = xₗ + (i - 1) * dx
        yᵤ = xₗ + i * dx
        balls[i] = Arb((yₗ, yᵤ))
    end

    return balls
end

"""
    iswide(x; cutoff = 10)

Return true if `x` is wide in the meaning that the effective relative
accuracy of `x` measured in bits is more than `cutoff` lower than it's
precision. For `x` not of type `Arb` or `Acb` this always return
`false`.
"""
iswide(x::Union{Arb,Acb}; cutoff = 10) = Arblib.rel_accuracy_bits(x) < precision(x) - cutoff

iswide(::Number; cutoff = 10) = false

"""
    dzeta(s)

Compute the Zeta function differentiated once with respect to `s`.
"""
dzeta(s::Arb) = zeta(ArbSeries([s, 1]))[1]
dzeta(s) = convert(float(typeof(s)), dzeta(Arb(s)))

function stieltjes(::Type{Arb}, n::Integer)
    res = zero(Acb)
    a = one(Acb)

    # This call uses fmpz from Nemo
    ccall(
        (:acb_dirichlet_stieltjes, Arblib.libarb),
        Cvoid,
        (Ref{Arblib.acb_struct}, Ref{Nemo.fmpz}, Ref{Arblib.acb_struct}, Clong),
        res,
        convert(Nemo.fmpz, n),
        a,
        precision(res),
    )

    return real(res)
end

stieltjes(type, n::Integer) = convert(type, stieltjes(Arb, n))

"""
    unique_integer(x::Arb)

If `x` contains a unique integer return `true, n` where `n` is the
integer. Otherwise return `false, 0`

"""
function unique_integer(x::Arb)
    # This call uses fmpz from Nemo
    res = Nemo.fmpz()
    unique = ccall(
        Arblib.@libarb(arb_get_unique_fmpz),
        Int,
        (Ref{Nemo.fmpz}, Ref{Arblib.arb_struct}),
        res,
        x,
    )

    return !iszero(unique), Int(res)
end

"""
    contains_pi(x1, x2)

Checks the interval `[x1, x2]` if it contains points of the form `kπ`.
Returns two booleans, the first one is false if it doesn't contain a
point on the form `2kπ` and the second one if it doesn't contain one
on the form `(2k+ 1)π`. If they are true it means they might contain
such a point. This is used to determine where the extrema of `Ci([x1,
x2])` can occur.

**TODO:** This is (hopefully) correct but not optimal.
"""
function contains_pi(x1::Arb, x2::Arb)
    @assert !(x1 > x2)

    # x1 or x2 equal to zero are the only cases when the division by π
    # can be exact, in which case it has to be handled differently.
    iszero(x1) && return (true, !(x2 < π))
    iszero(x2) && return (true, !(x1 > -Arb(π)))

    # We have k1ₗπ ≤ xₗ < (k1ᵤ + 1)π
    k1 = Arblib.floor!(zero(x1), x1 / π)

    unique1ₗ, k1ₗ = unique_integer(Arblib.ceil!(zero(k1), Arblib.lbound(Arb, k1)))
    unique1ᵤ, k1ᵤ = unique_integer(Arblib.floor!(zero(k1), Arblib.ubound(Arb, k1)))
    k1ₗ, k1ᵤ = Int(k1ₗ), Int(k1ᵤ)
    @assert unique1ₗ && unique1ᵤ && k1ₗ * Arb(π) ≤ x1 < (k1ᵤ + 1) * Arb(π)

    # We have k2ₗπ ≤ xₗ < (k2ᵤ + 1)π
    k2 = Arblib.floor!(zero(x2), x2 / π)
    unique2ₗ, k2ₗ = unique_integer(Arblib.ceil!(zero(k2), Arblib.lbound(Arb, k2)))
    unique2ᵤ, k2ᵤ = unique_integer(Arblib.floor!(zero(k2), Arblib.ubound(Arb, k2)))
    k2ₗ, k2ᵤ = Int(k2ₗ), Int(k2ᵤ)
    @assert unique2ₗ && unique2ᵤ && k2ₗ * Arb(π) ≤ x2 < (k2ᵤ + 1) * Arb(π)

    if k1ₗ == k2ᵤ
        # No kπ
        return (false, false)
    elseif k1ₗ == k2ᵤ - 1
        # Might contain exactly one such point (or zero)
        if iseven(k1ₗ)
            return (false, true)
        else
            return (true, false)
        end
    else
        # Might contain two such points
        return (true, true)
    end
end

"""
    beta_inc(a, b, z)

Compute the (not regularised) incomplete beta function \$B(a, b; z)\$.

Note that this method is different than
[`SpecialFunctions.beta_inc`](@ref) both in that it returns the
non-regularised value and that it only returns one value.
"""
beta_inc(a::Acb, b::Acb, z::Acb) = Arblib.hypgeom_beta_lower!(zero(z), a, b, z, 0)
beta_inc(a::Arb, b::Arb, z::Arb) = Arblib.hypgeom_beta_lower!(zero(z), a, b, z, 0)

"""
    beta_inc_zeroone(a, b, z)

Compute the (not regularised) incomplete beta function \$B(a, b; z)\$
assuming that `0 <= z <= 1`, discarding any other numbers in the
interval.

**PROVE:** That it's monotonically increasing in `z`.
"""
function beta_inc_zeroone(a::Arb, b::Arb, z::Arb)
    z_lower = Arblib.ispositive(z) ? lbound(Arb, z) : zero(z)
    z_upper = z < 1 ? ubound(Arb, z) : one(z)
    return Arb((beta_inc(a, b, z_lower), beta_inc(a, b, z_upper)))
end

"""
    abspow(x, y)

Compute `|x|^y `in a way that works if `x` contains negative numbers.
"""
function abspow(x::Arb, y::Arb)
    iszero(y) && return one(x)

    if iszero(x)
        Arblib.contains_negative(y) && return Arb(NaN, prec = precision(x))
        Arblib.ispositive(y) && return zero(x)
        return Arblib.unit_interval!(zero(x))
    end

    if Arblib.contains_zero(x)
        Arblib.contains_negative(y) && return Arb(NaN, prec = precision(x))
        x_upp = Arblib.abs_ubound(Arb, x)
        return Arb((zero(x), x_upp^y))
    end

    return abs(x)^y
end

abspow(x, y) = abs(x)^y

hypgeom_2f1(a::Arb, b::Arb, c::Arb, z::Arb) = Arblib.hypgeom_2f1!(zero(z), a, b, c, z, 0)

hypgeom_2f1(a::Acb, b::Acb, c::Acb, z::Acb) = Arblib.hypgeom_2f1!(zero(z), a, b, c, z, 0)

"""
    taylor_with_error(f, a::Arb, X::Arb, N::Integer)

Compute the Taylor expansion `P` of `f` of degree `N - 1` (i.e. `N`
terms) at the point `a` and bound the error term on `X`.

Returns a tuple `(P, E)` containing the Taylor expansion `P` and the
error term `E`. It satisfies that `f(x) ∈ P(x - a) + E*(x - a)^N` for
all `x ∈ X`.

We required that `a ∈ x`.
"""
function taylor_with_error(f, a::Arb, X::Arb, N::Integer)
    @assert contains(X, a)

    P = f(ArbSeries([a, one(a)], degree = N))

    E = f(ArbSeries([X, one(a)], degree = N + 1))[N]

    return P, E
end

"""
    abs(x::ArbSeries)

Compute the absolute value of `x`.

If `x[0]` contains zero then all non-constant terms are set to `NaN`,
otherwise either `-x` or `x` is returned depending on the sign of x.
"""
function Base.abs(x::ArbSeries)
    if Arblib.contains_zero(x[0])
        return ArbSeries([abs(x[0]); fill(Arb(NaN), Arblib.degree(x))])
    elseif x[0] < 0
        return -x
    else
        return x
    end
end
