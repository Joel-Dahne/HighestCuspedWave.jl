export alpha0, alpha0_estimate

"""
    alpha0(u0::AbstractAnsatz{Arb})

Upper bound the value of `α₀` from the paper.

Specialized for different types of `AbstractAnsatz`.
"""
alpha0

"""
    alpha0_estimate(u0::AbstractAnsatz)

Estimate the value of `α₀` from the paper. Uses the observation that
the maximum is obtained at `x = π`. This always gives a lower bound of
`α₀`.

**TODO:** This is not necessarily true in all cases for all choices of
weights. But since it's only an estimate it doesn't really matter.
"""
alpha0_estimate(u0::AbstractAnsatz{T}) where {T} = u0.w(convert(T, π)) / (2u0(covert(T, π)))
