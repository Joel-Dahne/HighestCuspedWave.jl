export KdVZeroAnsatz

"""
    KdVZeroAnsatz<:AbstractAnsatz{Arb}

Represents an ansatz for an approximate solution for the fractional
KdV equations as an expansion around `α = 0`.

The ansatz is given by
```
a[0] * clausencmzeta(x, 1 - α) + a[1] * clausencmzeta(x, 1 - α + p0) + a[2] * clausencmzeta(x, 1 - α + 2p)
```
where `a0`, `a1`, `a2` and `p0` all are computed through expansions in
`α` around zero.

It stores the value of `α` as well as precomputed `a` and `p0`.
"""
struct KdVZeroAnsatz <: AbstractAnsatz{Arb}
    α::Arb
    a::OffsetVector{ArbSeries,Vector{ArbSeries}}
    p0::ArbSeries

    function KdVZeroAnsatz(α::Arb)
        a = expansion_as(KdVZeroAnsatz)
        p0 = expansion_p0(KdVZeroAnsatz)
        u0 = new(α, a, p0)

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
    print(io, "KdVZeroAnsatz α = $(u0.α)")
end
