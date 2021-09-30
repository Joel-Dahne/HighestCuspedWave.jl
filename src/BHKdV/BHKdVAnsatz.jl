export BHKdVAnsatz

"""
    BHKdVAnsatz{T}(ϵ, v0, b)

Represents an ansatz used for the interval `α ∈ (-1, -1 + ϵ]`. The
ansatz consists of one main term depending on `α` and one fixed tail
which do not depend on `α`.

The main term, which depends on `α`, is given by
```
a0 * (Ci(x, 1 - α) - Ci(x, 1 - α + p0) - (zeta(1 - α) - zeta(1 - α + p0)))
```
where both `a0` and `p0` depend on `α`. The expression for `a0` is
that given by `finda0(α)`. The expression for `p0` is **not completely
determined yet**, one candidate is
```
p0 = 1 - α + (1 - α)^2 / 2
```
It should have the property that it converges to `0` as `α -> -1` but
also that `1 - α + p0 > 2`.
- **TODO:** Determine the value for `p0`.

The tail, which is fixed with respect to `α`, is the same as that used
for `BHAnsatz`. It consists of a sum of Clausen functions together
with a sum of Fourier terms. The Clausen functions used are determined
by the `v0`, which works similar to `BHAnsatz`, and the Fourier terms
by the vector `b`, again this works similar to `BHAnsatz`.

The weight is fixed to be given by
```
u0.w(x) = abs(x) * log(10 + inv(abs(x)))
```
The reason for the `10` inside the `log` is to make sure that the
non-asymptotic value of the norm is sufficiently low.
"""
struct BHKdVAnsatz{T} <: AbstractAnsatz{T}
    ϵ::T
    v0::FractionalKdVAnsatz{T}
    b::Vector{T}
end

"""
    BHKdVAnsatz{T}(ϵ::T, v0::BHAnsatz{T})
    BHKdVAnsatz(ϵ::T, v0::BHAnsatz{T})

Construct a `BHKdVAnsatz` with the given `ϵ` and the tail term given
by that of `u0`.

The Clausen part of the tail is given by `u0.v0` and the Fourier part
by `u0.b`.

If `isnothing(u0.v0)` then an empty `v0` is created for the tail part
which doesn't contain any terms (technically one term with the
coefficient `0`).
"""
function BHKdVAnsatz{T}(ϵ::T, u0::BHAnsatz{T}) where {T}
    @assert ϵ > 0

    if isnothing(u0.v0)
        v0 = FractionalKdVAnsatz(-1 + ϵ, 0, 0)
        v0.a[0] = 0
    else
        v0 = deepcopy(u0.v0)
    end

    return BHKdVAnsatz{T}(ϵ, v0, copy(u0.b))
end

BHKdVAnsatz(ϵ::T, u0::BHAnsatz{T}) where {T} = BHKdVAnsatz{T}(ϵ, u0)

function Base.getproperty(u0::BHKdVAnsatz, name::Symbol)
    if name == :N
        return length(u0.b)
    elseif name == :w
        return x -> abs(x) * log(10 + inv(abs(x)))
    else
        return getfield(u0, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u0::BHKdVAnsatz)
    print(io, "$(typeof(u0)) N = $(u0.N), u0.v0.N0 = $(u0.v0.N0), ϵ = $(u0.ϵ)")
end
