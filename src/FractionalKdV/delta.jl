"""
    delta0(u0::FractionalKdVAnsatz; ϵ::arb = parent(u0.α)(0.1))
Upper bound the value of δ₀ from the paper. Uses an asymptotic
expansion on the interval [0, ϵ] and ball arithmetic on [ϵ, π].
"""
function delta0(u0::FractionalKdVAnsatz{arb};
                ϵ::arb = parent(u0.α)(0.1),
                )
    # Bound the value one [0, ϵ]
    # TODO: Implement this
    res1 = enclosemaximum(F0(u0, Asymptotic()), parent(u0.α)(0), ϵ,
                          absmax = true,
                          maxevals = 10^4,
                          )

    # Bound the value on [ϵ, π] by Ball evaluation
    # TODO: So far this naive approach only works for a very limited
    # range of parameters and needs to be improved.
    res2 = enclosemaximum(F0(u0), ϵ, parent(u0.α)(π),
                          absmax = true,
                          maxevals = 10^4,
                          )

    return max(res1, res2)
end

"""
    delta0_estimate(u0::FractionalKdVAnsatz; n::Integer = 100)
Estimate the value of δ₀ from the paper. Does this by evaluating F0 it
on n linearly spaced points. This always gives a lower bound of δ₀.
"""
function delta0_estimate(u0::FractionalKdVAnsatz{T};
                         n::Integer = 100,
                         ) where {T}
    xs = range(0, stop = π, length = n)[2:end]
    if T == arb
        xs = parent(u0.α).(xs)
    end
    return maximum(F0(u0).(xs))
end
