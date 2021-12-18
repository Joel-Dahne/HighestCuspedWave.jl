export KdVZeroAnsatz

"""
    KdVZeroAnsatz(α::Arb; degree::Integer = 2) <: AbstractAnsatz{Arb}

Represents an ansatz for an approximate solution for the fractional
KdV equations as an expansion around `α = 0`.

The ansatz is given by
```
a[0] * clausencmzeta(x, 1 - α) + a[1] * clausencmzeta(x, 1 - α + p0) + a[2] * clausencmzeta(x, 1 - α + 2p0)
```
where `a0`, `a1`, `a2` and `p0` all are computed through expansions in
`α` around zero.

It stores the value of `α` as well as precomputed expansions for `a`
and `p0`. The expansions for `a` and `p0` are computed to the degree
`degree`. The last term in each expansion is a remainder term which
ensures that that evaluating the expansion gives an enclosure of the
value for the full interval in `α`, i.e. `a[i] ∈ u0.a[i](α)` for each
`α ∈ u0.α` and similarly for `p0`.

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
