export BHKdVAnsatz

"""
    BHKdVAnsatz{T}(ϵ::T, v0::BHAnsatz{T})
    BHKdVAnsatz(ϵ::T, v0::BHAnsatz{T})

Represents an ansatz used for the interval `α ∈ (-1, -1 + ϵ]`. The
ansatz is similar to that given by `v0` but uses a different main term
and a different weight.

The main term depends on `α` and is given by
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
by `v0.v0` and the Fourier terms by the vector `v0.b`. For more
details see `BHAnsatz`.

It assumes that `v0` uses the correct coefficient in front of the
leading Clausen function, i.e. `v0.a0 = 2 / π^2`.

The weight is fixed to be given by
```
u0.w(x) = abs(x) * log(10 + inv(abs(x)))
```
The reason for the `10` inside the `log` is to make sure that the
non-asymptotic value of the norm is sufficiently low.

If `isnothing(v0.v0)` then an empty `v0.v0` is created for the tail
part which doesn't contain any terms (technically one term with the
coefficient `0`). This is to avoid having to deal with this as a
special case.
"""
struct BHKdVAnsatz{T} <: AbstractAnsatz{T}
    ϵ::T
    v0::BHAnsatz{T}

    function BHKdVAnsatz{T}(ϵ::T, v0::BHAnsatz{T}) where {T}
        @assert ϵ > 0
        if T == arb
            @assert Arblib.overlaps(v0.a0, 2 / convert(T, π)^2)
        end

        if isnothing(v0.v0)
            v0v0 = FractionalKdVAnsatz(-1 + ϵ, 0, 0)
            v0v0.a[0] = 0
            v0 = BHAnsatz{T}(v0.a0, copy(v0.b), v0v0)
        else
            v0 = deepcopy(v0)
        end

        return new{T}(ϵ, v0)
    end
end

BHKdVAnsatz(ϵ::T, v0::BHAnsatz{T}) where {T} = BHKdVAnsatz{T}(ϵ, v0)

function Base.getproperty(u0::BHKdVAnsatz, name::Symbol)
    if name == :w
        return x -> abs(x) * log(10 + inv(abs(x)))
    else
        return getfield(u0, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u0::BHKdVAnsatz)
    print(io, "$(typeof(u0)) u0.v0.N = $(u0.v0.N), u0.v0.v0.N0 = $(u0.v0.v0.N0), ϵ = $(u0.ϵ)")
end
