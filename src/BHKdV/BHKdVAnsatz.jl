export BHKdVAnsatz

"""
    BHKdVAnsatz{T}(ϵ::T, v0::BHAnsatz{T})
    BHKdVAnsatz(ϵ::T, v0::BHAnsatz{T})

Represents an ansatz used for the interval `α ∈ (-1, -1 + ϵ]`. The
ansatz is similar to that given by `v0` but uses a different main term
and a different weight.

The main term depends on `α` and is given by
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
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

The weight is given by
```
u0.w(x) = abs(x)^(1 - u0.γ * (1 + α)) * log(u0.c + inv(abs(x)))
```
where `u0.c >= 1` and `0 <= u0.γ < 1`. By default `u0.c = 2ℯ` and
`u0.γ = 0`. The reason for the value of `u0.c` inside the `log` is
to make sure that the non-asymptotic value of the norm is sufficiently
low. The reason for the `u0.γ` is to make sure that the norm of the
operator is less than one at `x = 0`.

If `isnothing(v0.v0)` then an empty `v0.v0` is created for the tail
part which doesn't contain any terms (technically one term with the
coefficient `0`). This is to avoid having to deal with this as a
special case.
"""
struct BHKdVAnsatz{T} <: AbstractAnsatz{T}
    ϵ::T
    v0::BHAnsatz{T}
    γ::T
    c::T

    function BHKdVAnsatz{T}(
        ϵ::T,
        v0::BHAnsatz{T};
        γ::T = convert(T, 0),
        c::T = 2convert(T, ℯ),
    ) where {T}
        @assert 0 < ϵ
        @assert 0 <= γ < 1
        @assert 1 <= c
        T == Arb && @assert Arblib.overlaps(v0.a0, 2 / convert(T, π)^2)

        if isnothing(v0.v0)
            v0v0 = FractionalKdVAnsatz(-1 + ϵ, 0, 0)
            v0v0.a[0] = 0
            v0 = BHAnsatz{T}(v0.a0, copy(v0.b), v0v0)
        else
            # FIXME: This should be a deepcopy of v0 once Arblib
            # implements that.
            v0 = deepcopy(v0)
        end

        return new{T}(ϵ, v0, γ, c)
    end
end

BHKdVAnsatz(ϵ::T, v0::BHAnsatz{T}; γ::T = convert(T, 0), c::T = 2convert(T, ℯ)) where {T} =
    BHKdVAnsatz{T}(ϵ, v0; γ, c)

function Base.getproperty(u0::BHKdVAnsatz{T}, name::Symbol) where {T}
    if name == :w
        if iszero(u0.γ)
            return x -> abs(x) * log(u0.c + inv(abs(x)))
        else
            return x -> let αp1 = Arblib.nonnegative_part!(zero(u0.ϵ), Arb((0, u0.ϵ)))
                abspow(x, 1 - u0.γ * αp1) * log(u0.c + inv(abs(x)))
            end
        end
    elseif name == :wdivx
        if iszero(u0.γ)
            return x -> log(u0.c + inv(abs(x)))
        else
            return x -> let αp1 = Arblib.nonnegative_part!(zero(u0.ϵ), Arb((0, u0.ϵ)))
                abspow(x, -u0.γ * αp1) * log(u0.c + inv(abs(x)))
            end
        end
    else
        return getfield(u0, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u0::BHKdVAnsatz)
    println(
        io,
        "$(typeof(u0)) u0.v0.N = $(u0.v0.N), u0.v0.v0.N0 = $(u0.v0.v0.N0), ϵ = $(u0.ϵ)",
    )
    print(io, "γ = $(u0.γ), c = $(u0.c)")
end
