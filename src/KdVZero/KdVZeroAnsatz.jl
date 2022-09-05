export KdVZeroAnsatz

"""
    KdVZeroAnsatz(α::Arb, α0::Arb; degree::Integer = 2) <: AbstractAnsatz{Arb}

Represents an ansatz for an approximate solution for the fractional
KdV equations as an expansion centered at `α0` and valid on `α`.

The ansatz is given by
```
a[0] * clausencmzeta(x, 1 - α) + a[1] * clausencmzeta(x, 1 - α + p0) + a[2] * clausencmzeta(x, 1 - α + 2p0)
```
where `a0`, `a1`, `a2` and `p0` are given by expansions around `α0`.

The value of `α` is assumed to be non-positive. However the enclosure
for `α` might contain spurious positive parts. We explicitly remove
these so that the right endpoint of the enclosure is non-negative.

The expansions are handled using Taylor models through the
[`TaylorModel`](@ref) type. In general the models are computed to the
given degree, though in some cases lower or higher degrees are used.

It stores precomputed Taylor models for `a` and `p0`. The Taylor
models for `a` and `p0` are computed to the degree `degree` with the
exception of `a[0]` which is computed to one degree higher, the reason
for this is that for `α0 = 0` the constant term in `a[0]` is zero and
in many cases we divide `a[0]` by `α` and still want to have
sufficiently high degree.

For the case `α0 = 0` there are typically several removable
singularities to handle and for that reason this case often have a
special implementation.
"""
struct KdVZeroAnsatz <: AbstractAnsatz{Arb}
    α::Arb
    α0::Arb
    a::OffsetVector{TaylorModel,Vector{TaylorModel}}
    p0::TaylorModel
    degree::Int

    function KdVZeroAnsatz(α::Arb, α0::Arb = Arb(0); degree::Integer = 2)
        contains(α, α0) || throw(ArgumentError("expected α0 to be contained in α"))

        # α is non-positive, this removes any spurious positive parts
        α = -Arblib.nonnegative_part!(zero(α), -α)

        # TODO: These return Taylor models which we for now have to
        # convert to the raw underlying polynomial
        p0 = expansion_p0(KdVZeroAnsatz, α0, α, degree = degree - 1)
        a = expansion_as(KdVZeroAnsatz, α0, α, degree = degree - 1; p0)

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
