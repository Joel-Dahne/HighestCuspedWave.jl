export KdVZeroAnsatz

"""
    KdVZeroAnsatz(α::Arb, α0::Arb; degree::Integer = 2) <: AbstractAnsatz{Arb}

Represents an ansatz for an approximate solution for the fractional
KdV equations as an expansion around `α0`.

The computed expansions are valid on the interval given by `α` and we
require that `α0` is contained in `α`.

The ansatz is given by
```
a[0] * clausencmzeta(x, 1 - α) + a[1] * clausencmzeta(x, 1 - α + p0) + a[2] * clausencmzeta(x, 1 - α + 2p0)
```
where `a0`, `a1`, `a2` and `p0` all give as expansions in around `α0`.

Any expansions computed with this approximation are typically computed
to the given degree, though some methods compute to a lower or higher
degree.

The value of `α` is assumed to be non-positive. However the enclosure
for `α` might contain spurious positive parts. We explicitly remove
these so that the right endpoint of the enclosure is non-negative.

For the case `α0 = 0` there are typically several removable
singularities to handle and for that reason this case often have a
special implementation.

It stores precomputed expansions for `a` and `p0`. The expansions for
`a` and `p0` are computed to the degree `degree` with the exception of
`a[0]` which is computed to one degree higher, the reason for this is
that for `α0 = 0` the constant term in `a[0]` is zero and in many
cases we divide `a[0]` by `α` and still want to have sufficiently high
degree. The last term in each expansion is a remainder term which
ensures that that evaluating the expansion gives an enclosure of the
value for the full interval in `α`, i.e. `a[i] ∈ u0.a[i](α - α0)` for
each `α ∈ u0.α` and similarly for `p0`.

!!! Note that this way of having the last term in the expansion
    represent a remainder term is a non-standard way of using
    `ArbSeries`. The normal way is to have each term represent the
    corresponding derivative and any remainder term is stored outside
    of the expansion. For that reason we have to take care when
    computing with these expansions. Most of the functionality for
    this is implemented in several `_with_remainder` methods, e.g.
    [`compose_with_remainder`](@ref).
"""
struct KdVZeroAnsatz <: AbstractAnsatz{Arb}
    α::Arb
    α0::Arb
    a::OffsetVector{ArbSeries,Vector{ArbSeries}}
    p0::ArbSeries
    degree::Int

    function KdVZeroAnsatz(α::Arb, α0::Arb = Arb(0); degree::Integer = 2)
        contains(α, α0) || throw(ArgumentError("expected α0 to be contained in α"))

        # α is non-positive, this removes any spurious positive parts
        α = -Arblib.nonnegative_part!(zero(α), -α)

        a = expansion_as(KdVZeroAnsatz, α, α0; degree)
        p0 = expansion_p0(KdVZeroAnsatz, α, α0; degree)
        u0 = new(α, α0, a, p0, degree)

        return u0
    end
end

function Base.getproperty(u0::KdVZeroAnsatz, name::Symbol)
    if name == :N0
        return 2
    elseif name == :w
        return abs
    else
        return getfield(u0, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u0::KdVZeroAnsatz)
    print(io, "KdVZeroAnsatz α = $(u0.α), α0 = $(u0.α0), degree = $(u0.degree)")
end
