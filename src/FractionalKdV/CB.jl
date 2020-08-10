export CB, CB_estimate

"""
    CB(u0::FractionalKdVAnsatz{arb}; ϵ::arb = parent(u0.α)(0.1))
Upper bound the value of C_B from the paper. TODO: Uses an asymptotic
expansion on the interval [0, ϵ] and ball arithmetic on [ϵ, π].
"""
function CB(u0::FractionalKdVAnsatz{arb};
            ϵ::arb = parent(u0.α)(0.1),
            rtol::arb = parent(u0.α)(1e-2),
            show_trace = false,
            )
    # Bound the value one [0, ϵ]
    # TODO: Implement this
    n1 = enclosemaximum(T0(u0, Asymptotic()), parent(u0.α)(0), ϵ,
                        rtol = rtol,
                        absmax = true,
                        maxevals = 10^3,
                        show_trace = show_trace,
                        )

    # Bound the value one [ϵ, π]
    # TODO: This does not fully work yet
    n2 = enclosemaximum(T0(u0, rtol = 1e-2*Float64(rtol)), ϵ, parent(u0.α)(π),
                        rtol = rtol,
                        absmax = true,
                        maxevals = 10^3,
                        show_trace = show_trace,
                        )

    return max(n1, n2)
end

"""
    CB_estimate(u0::FractionalKdVAnsatz; n::Integer = 100)
Estimate the value of C_B from the paper. Does this by evaluating the
norm it on n linearly spaced points. This always gives a lower bound
of C_B. Currently it gives a very bad estimate.
"""
function CB_estimate(u0::FractionalKdVAnsatz{T};
                     n::Integer = 20,
                     show_trace = false,
                     ) where {T}
    xs = range(0, stop = π, length = n)[2:end]
    if T == arb
        xs = parent(u0.α).(xs)
    end

    res = abs.(T0(u0, rtol = 1e-3, show_trace = show_trace).(xs))
    m = zero(u0.α)
    for r in res
        m = max(m, r)
    end

    return m
end
