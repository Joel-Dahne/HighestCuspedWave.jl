export TaylorModel

"""
    TaylorModel(p, I, x0)

Struct representing a Taylor model. The remainder term `Δ` is stored
as the last term in `p`.

The type and most of the methods are based on
- Joldes, M. M. (2011). Rigorous polynomial approximations and
  applications (Doctoral dissertation). Ecole normale supérieure
  de lyon - ENS LYON.
More precisely we implement what in that thesis is referred to as
Taylor models with relative remainder. Most of the methods we define
on Taylor models are straight forward adaptions of the algorithms
given by Joldes.
"""
struct TaylorModel
    p::ArbSeries
    I::Arb
    x0::Arb
end

"""
    TaylorModel(f, I::Arb, x0::Arb; degree::Integer, enclosure_degree::Integer = -1)

Construct a Taylor model of `f` on the interval `I` centered at `x0`
with the given degree.

For wide values of `I` it computes a tighter enclosure of the
remainder term using [`ArbExtras.enclosure_series`](@ref). The degree
used for this can be set with `enclosure_degree`. Setting it to a
negative number makes it compute it directly instead.
"""
function TaylorModel(f, I::Arb, x0::Arb; degree::Integer, enclosure_degree::Integer = -1)
    contains(I, x0) || throw(
        ArgumentError("expected x0 to be contained in interval, got x0 = $x0, I = $I"),
    )

    if x0 == I
        p = f(ArbSeries((x0, 1), degree = degree + 1))
    else
        p = f(ArbSeries((x0, 1); degree))

        # Make room for remainder term
        p = ArbSeries(p, degree = degree + 1)

        # Compute remainder term
        if enclosure_degree < 0 || !iswide(I)
            p[degree+1] = f(ArbSeries((I, 1), degree = degree + 1))[degree+1]
        else
            # We compute a tighter enclosure with the help of ArbExtras.enclosure_series
            p[degree+1] =
                ArbExtras.enclosure_series(
                    ArbExtras.derivative_function(f, degree + 1),
                    I,
                    degree = enclosure_degree,
                ) / factorial(degree + 1)
        end
    end

    return TaylorModel(p, I, x0)
end

Base.zero(M::TaylorModel) = TaylorModel(zero(M.p), M.I, M.x0)
Base.one(M::TaylorModel) = TaylorModel(one(M.p), M.I, M.x0)
Base.iszero(M::TaylorModel) = iszero(M.p)
Base.isone(M::TaylorModel) = isone(M.p)

Arblib.degree(M::TaylorModel) = Arblib.degree(M.p) - 1

function Base.show(io::IO, ::MIME"text/plain", M::TaylorModel)
    println(io, "Taylor model of degree $(Arblib.degree(M)) centered at x0 = $(M.x0)")
    println(io, "I = $(M.I), Δ = $(M.p[end])")
    print(io, "p = $(ArbPoly(M.p[0:end-1]))")
end

"""
    checkcompatible(::Type{Bool}, M1::TaylorModel, M2::TaylorModel)

Return `true` if `M1` and `M2` have the same degree, center and
interval.

**IMPROVE:** At the moment the midpoints have to be exactly equal,
meaning that only point-intervals are allowed. This seems to be enough
for what we need.
"""
checkcompatible(::Type{Bool}, M1::TaylorModel, M2::TaylorModel) =
    Arblib.degree(M1) == Arblib.degree(M2) && isequal(M1.I, M2.I) && M1.x0 == M2.x0

"""
    checkcompatible(::Type{Bool}, M1::TaylorModel, M2::TaylorModel)

Throw an error if `M1` and `M2` are not compatible according to
`checkcompatible(Bool, M1, M2)`.

**IMPROVE:** Add info about what failed.
"""
checkcompatible(M1::TaylorModel, M2::TaylorModel) =
    if !checkcompatible(Bool, M1, M2)
        d1 = Arblib.degree(M1)
        d2 = Arblib.degree(M2)
        Arblib.degree(M1) == Arblib.degree(M2) || error(
            "Taylor models with non-compatible degrees, degree(M1) = $d1, degree(M2) = $d2",
        )
        error("Taylor models with non-compatible intervals")
    end

"""
    overlaps(M1::TaylorModel, M2::TaylorModel; require_compatible::Bool = true)

Return true if the coefficients and remainder term of `M1` and `M2`
overlaps.

By default it throws an error if `M1` and `M2` are not compatible,
according to [`checkcompatible`](@ref). Setting `require_compatible`
to false removes the requirement that `M1.I` and `M2.I` should be
equal, it still requires that they have the same degree and midpoint.
"""
function Arblib.overlaps(M1::TaylorModel, M2::TaylorModel; require_compatible::Bool = true)
    if require_compatible
        checkcompatible(M1, M2)
    else
        Arblib.degree(M1) == Arblib.degree(M2) || error(
            "Taylor models with non-compatible degrees, degree(M1) = $(Arblib.degree(M1)), degree(M2) = $(Arblib.degree(M2))",
        )
        M1.x0 == M2.x0 || error("Taylor models with different midpoints")
    end

    return Arblib.overlaps(M1.p, M2.p)
end

(M::TaylorModel)(x::Arb) = M.p(x - M.x0)

"""
    truncate_with_remainder(p::ArbPoly, I::Arb, x0::Arb; degree::Integer)

Compute a polynomial `q` corresponding to a truncated version of `p`
with the last term being a remainder term which ensures `p(x - x0) ⊂
res(x - x0)` for all `x ∈ I`.
"""
function truncate_with_remainder(p::ArbPoly, I::Arb, x0::Arb; degree::Integer)
    # Set q to a truncated version of p
    q = Arblib.truncate!(copy(p), degree + 1)

    Arblib.degree(p) <= degree && return q

    # Set p_div_x to a p divided by x^degree, throwing away lower order
    # terms.
    p_div_x = Arblib.shift_right!(zero(p), p, degree)

    # Evaluate q on the given interval and set this as the remainder term
    q[degree] = p_div_x(I - x0)

    return q
end

"""
    truncate(M::TaylorModel; degree::Integer)

Compute a Taylor model enclosing `M` but of a lower `degree`.
"""
function truncate(M::TaylorModel; degree::Integer)
    degree <= Arblib.degree(M) || throw(
        ArgumentError(
            "can't truncate degree $(Arblib.degree(M)) TaylorModel to degree $degree",
        ),
    )

    degree == Arblib.degree(M) && return M

    # Set the polynomial to a truncated version of M.p
    return TaylorModel(
        ArbSeries(
            truncate_with_remainder(M.p.poly, M.I, M.x0, degree = degree + 1),
            degree = degree + 1,
        ),
        M.I,
        M.x0,
    )
end

"""
    compose(f, M::TaylorModel)

Compute a Taylor model of `f` applied to `M`. If `M` is a Taylor model
of the function `g` then this gives a Taylor model of the function `f
∘ g`.
"""
function compose(f, M::TaylorModel)
    degree = Arblib.degree(M)

    # Compute expansion of f at M.p[0]
    p = f(ArbSeries((M.p[0], 1); degree))

    # Increase the degree of p to make room for the remainder term
    p = ArbSeries(p, degree = degree + 1)

    # Compute remainder term
    J = M(M.I) # Interval to compute remainder on
    if isfinite(J)
        # We compute a tighter enclosure with the help of ArbExtras.enclosure_series
        remainder_term =
            ArbExtras.enclosure_series(ArbExtras.derivative_function(f, degree + 1), J) /
            factorial(degree + 1)
    else
        remainder_term = indeterminate(J)
    end

    p[degree+1] = remainder_term

    # Compute a non-truncated composition
    q = ArbExtras.compose_zero(p.poly, M.p.poly)

    # Truncate the polynomial to the specified degree
    return TaylorModel(
        ArbSeries(
            truncate_with_remainder(q, M.I, M.x0, degree = degree + 1),
            degree = degree + 1,
        ),
        M.I,
        M.x0,
    )
end

function Base.:+(M1::TaylorModel, M2::TaylorModel)
    checkcompatible(M1, M2)
    return TaylorModel(M1.p + M2.p, M1.I, M1.x0)
end

function Base.:-(M1::TaylorModel, M2::TaylorModel)
    checkcompatible(M1, M2)
    TaylorModel(M1.p - M2.p, M1.I, M1.x0)
end

function Base.:*(M1::TaylorModel, M2::TaylorModel)
    checkcompatible(M1, M2)

    # Compute non-truncated polynomial
    p = M1.p.poly * M2.p.poly

    # Truncate to the specified degree
    return TaylorModel(
        ArbSeries(
            truncate_with_remainder(p, M1.I, M1.x0, degree = Arblib.degree(M1) + 1),
            degree = Arblib.degree(M1) + 1,
        ),
        M1.I,
        M1.x0,
    )
end

Base.:/(M1::TaylorModel, M2::TaylorModel) = M1 * compose(inv, M2)

# Scalar functions

Base.:+(M::TaylorModel, c::Union{Arb,Integer}) = TaylorModel(M.p + c, M.I, M.x0)
Base.:-(M::TaylorModel, c::Union{Arb,Integer}) = TaylorModel(M.p - c, M.I, M.x0)
Base.:*(M::TaylorModel, c::Union{Arb,Integer}) = TaylorModel(M.p * c, M.I, M.x0)
Base.:/(M::TaylorModel, c::Union{Arb,Integer}) = TaylorModel(M.p / c, M.I, M.x0)
Base.:+(c::Union{Arb,Integer}, M::TaylorModel) = M + c
Base.:-(c::Union{Arb,Integer}, M::TaylorModel) = TaylorModel(c - M.p, M.I, M.x0)
Base.:*(c::Union{Arb,Integer}, M::TaylorModel) = M * c
Base.:/(c::Union{Arb,Integer}, M::TaylorModel) = c * compose(inv, M)

Base.:-(M::TaylorModel) = TaylorModel(-M.p, M.I, M.x0)


function Base.:(<<)(M::TaylorModel, n::Integer)
    n <= Arblib.degree(M) || error("shift must be less than degree of TaylorModel")
    TaylorModel(M.p << n, M.I, M.x0)
end
function Base.:(>>)(M::TaylorModel, n::Integer)
    TaylorModel(M.p >> n, M.I, M.x0)
end

"""
    div_removable(M1::TaylorModel, M2::TaylorModel, order::Integer = 1; force = true)

Compute `M1 / M2` in the case of a removable singularity of the given
`order`. The degree of the output is the degree of the input minus
`order`.
"""
function div_removable(M1::TaylorModel, M2::TaylorModel, order::Integer = 1; force = true)
    checkcompatible(M1, M2)
    order >= 1 || error("order must be positive")
    order <= Arblib.degree(M1) || error("order must be lower than degree of input")

    if force
        # Optimize in case all things happen to be zero?
        M1 = TaylorModel(copy(M1.p), M1.I, M1.x0)
        M2 = TaylorModel(copy(M2.p), M2.I, M2.x0)
        for i = 0:order-1
            @assert Arblib.contains_zero(Arblib.ref(M1.p, i))
            @assert Arblib.contains_zero(Arblib.ref(M2.p, i))
            M1.p[i] = 0
            M2.p[i] = 0
        end
    end

    return (M1 << order) / (M2 << order)
end

# TaylorModel implementations of some functions that are used a lot

"""
    clausenc(x::Arb, s::TaylorModel)

Compute a Taylor model of `clausenc(x, s)` in the parameter `s`.

The only difference to calling `compose(s -> clausenc(x, s), s)` is
that it computes the remainder term by direct evaluation instead of
through `ArbSeries`.
"""
function clausenc(x::Arb, s::TaylorModel)
    degree = Arblib.degree(s)

    # Compute expansion of f at s.p[0]
    p = clausenc(x, ArbSeries((s.p[0], 1); degree))

    # Increase the degree of p to make room for the remainder term
    p = ArbSeries(p, degree = degree + 1)

    # Compute remainder term
    J = s(s.I) # Interval to compute remainder on
    remainder_term = clausenc(x, J, degree + 1) / factorial(degree + 1)
    p[degree+1] = remainder_term

    # Compute a non-truncated composition
    q = ArbExtras.compose_zero(p.poly, s.p.poly)

    # Truncate the polynomial to the specified degree
    return TaylorModel(
        ArbSeries(
            truncate_with_remainder(q, s.I, s.x0, degree = degree + 1),
            degree = degree + 1,
        ),
        s.I,
        s.x0,
    )
end

"""
    clausens(x::Arb, s::TaylorModel)

Compute a Taylor model of `clausens(x, s)` in the parameter `s`.

The only difference to calling `compose(s -> clausens(x, s), s)` is
that it computes the remainder term by direct evaluation instead of
through `ArbSeries`.
"""
function clausens(x::Arb, s::TaylorModel)
    degree = Arblib.degree(s)

    # Compute expansion of f at s.p[0]
    p = clausens(x, ArbSeries((s.p[0], 1); degree))

    # Increase the degree of p to make room for the remainder term
    p = ArbSeries(p, degree = degree + 1)

    # Compute remainder term
    J = s(s.I) # Interval to compute remainder on
    remainder_term = clausens(x, J, degree + 1) / factorial(degree + 1)
    p[degree+1] = remainder_term

    # Compute a non-truncated composition
    q = ArbExtras.compose_zero(p.poly, s.p.poly)

    # Truncate the polynomial to the specified degree
    return TaylorModel(
        ArbSeries(
            truncate_with_remainder(q, s.I, s.x0, degree = degree + 1),
            degree = degree + 1,
        ),
        s.I,
        s.x0,
    )
end

"""
    clausencmzeta(x::Arb, s::ArbSeries, interval::Arb; degree)

Compute a Taylor model of `clausencmzeta(x, s)` in the parameter `s`.

For wide values of `x` it computes each term in the expansion
separately, allowing it to use [`ArbExtras.enclosure_series`](@ref) to
compute a tighter enclosure of the terms. While this could be done
also for the remainder term it currently isn't, it doesn't give much
improvement in the cases it could be relevant and comes with a big
performance cost.
"""
function clausencmzeta(x::Arb, s::TaylorModel)
    degree = Arblib.degree(s)

    if iswide(x)
        p = ArbSeries(degree = degree + 1)

        # For the constant term it already implements handling of
        # monotonicity and we don't have to use
        # ArbExtras.enclosure_series
        p[0] = clausencmzeta(x, s.p[0])

        for β = 1:degree
            p[β] =
                ArbExtras.enclosure_series(x -> clausencmzeta(x, s.p[0], β), x, degree = 1)
        end
    else
        # Compute expansion at s[0] with degree - 1
        p = clausencmzeta(x, ArbSeries((s.p[0], 1), degree = degree))

        # Increase the degree of p to make room for the remainder term
        p = ArbSeries(p, degree = degree + 1)
    end

    # Compute remainder term
    remainder_term = clausencmzeta(x, s(s.I), degree + 1) / factorial(degree + 1)

    p[degree+1] = remainder_term

    # Compute a non-truncated composition
    q = ArbExtras.compose_zero(p.poly, s.p.poly)

    # Truncate to the specified degree
    return TaylorModel(
        ArbSeries(
            truncate_with_remainder(q, s.I, s.x0, degree = degree + 1),
            degree = degree + 1,
        ),
        s.I,
        s.x0,
    )
end

"""
    abspow(x::Arb, y::TaylorModel)

Compute a Taylor model of `abspow(x, y)` in the argument `y`.

This is equivalent to `compose_with_remainder(s -> abspow(x, y), y,
interval)`, with one specialised optimizations. If `Arblib.is_x(-y)`
is true then instead of evaluating `y(interval)` it takes `-interval`
directly. These two are obviously equivalent, the reason for doing
this optimization is that in some cases we have an exponent which is
non-positive but whose left endpoint is zero, in that case
`y(interval)` is supposed to non-negative but naive evaluation doesn't
preserve this.
"""
function abspow(x::Arb, y::TaylorModel)
    degree = Arblib.degree(y)

    # Compute expansion at y[0] with degree - 1
    p = abspow(x, ArbSeries((y.p[0], 1); degree))

    # Increase the degree of p to make room for the remainder term
    p = ArbSeries(p, degree = degree + 1)

    # Compute remainder term
    # We compute the interval in this way to avoid spurious negative
    # values of y in the special case Arblib.is_x(-y)
    yinterval = Arblib.is_x(-y.p) ? -y.I : y(y.I)

    p[degree+1] = abspow(x, (ArbSeries((yinterval, 1), degree = degree + 1)))[degree+1]

    # Compute a non-truncated composition
    q = ArbExtras.compose_zero(p.poly, y.p.poly)

    # Truncate to the specified degree
    return TaylorModel(
        ArbSeries(
            truncate_with_remainder(q, y.I, y.x0, degree = degree + 1),
            degree = degree + 1,
        ),
        y.I,
        y.x0,
    )
end
