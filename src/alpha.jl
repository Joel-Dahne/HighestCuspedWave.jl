export alpha0, alpha0_estimate

"""
    alpha0(u0::AbstractAnsatz; M::Integer = 3, rtol::arb = 1e-5)
Upper bound the value of `α₀` from the paper.

Specialized for different types of ansatz.
"""
alpha0

"""
    alpha0_estimate(u0::AbstractAnsatz)
Estimate the value of `α₀` from the paper. Uses the observation that
the maximum is obtained at `x = π`. This always gives a lower bound of
`α₀`.
"""
function alpha0_estimate(u0::AbstractAnsatz{T}) where {T}
    if T == arb
        π = u0.parent(pi)
    else
        π = convert(T, pi)
    end
    return u0.w(π) / (2u0(π))
end
