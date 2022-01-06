export KdVZeroAnsatz

"""
    KdVZeroAnsatz(α::Arb; degree::Integer = 2) <: AbstractAnsatz{Arb}

Represents an ansatz for an approximate solution for the fractional
KdV equations as an expansion around `α = 0`. The computed expansions
are valid on the interval given by `α`.

The ansatz is given by
```
a[0] * clausencmzeta(x, 1 - α) + a[1] * clausencmzeta(x, 1 - α + p0) + a[2] * clausencmzeta(x, 1 - α + 2p0)
```
where `a0`, `a1`, `a2` and `p0` all are computed through expansions in
`α` around zero.

Any expansions computed with this approximation are typically computed
to the given degree, though some methods compute to a lower or higher
degree.

The value of `α` is assumed to be non-positive. However the enclosure
for `α` might contain spurious positive parts. We explicitly remove
these so that the right endpoint of the enclosure is non-negative.

It stores precomputed expansions for `a` and `p0`. The expansions for
`a` and `p0` are computed to the degree `degree` with the exception of
`a[0]` which is computed to one degree higher, the reason for this is
that the constant term in `a[0]` is zero and in many cases we divide
`a[0]` by `α` and still want to have sufficiently high degree. The
last term in each expansion is a remainder term which ensures that
that evaluating the expansion gives an enclosure of the value for the
full interval in `α`, i.e. `a[i] ∈ u0.a[i](α)` for each `α ∈ u0.α` and
similarly for `p0`.

!!! Note that this way of having the last term in the expansion
    represent a remainder term is a non-standard way of using
    `ArbSeries`, typically each term represents the corresponding
    derivative and any remainder term is stored outside of the
    expansion, and for that reason we have to take care when computing
    with these expansions.
"""
struct KdVZeroAnsatz <: AbstractAnsatz{Arb}
    α::Arb
    a::OffsetVector{ArbSeries,Vector{ArbSeries}}
    p0::ArbSeries
    degree::Int

    function KdVZeroAnsatz(α::Arb; degree::Integer = 2)
        # α is non-positive, this removes any spurious positive parts
        α = -Arblib.nonnegative_part!(zero(α), -α)

        a = expansion_as(KdVZeroAnsatz, α; degree)
        p0 = expansion_p0(KdVZeroAnsatz, α; degree)
        u0 = new(α, a, p0, degree)

        return u0
    end
end

function Base.getproperty(u0::KdVZeroAnsatz, name::Symbol)
    if name == :N0
        return 2
    elseif name == :w
        return x -> abs(x)
    else
        return getfield(u0, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u0::KdVZeroAnsatz)
    print(io, "KdVZeroAnsatz α = $(u0.α), degree = $(u0.degree)")
end
