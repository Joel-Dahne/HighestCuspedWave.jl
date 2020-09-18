export CB, CB_bounded_by, CB_estimate

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
    tol = 1e-4*Float64(rtol)
    n2 = enclosemaximum(T0(u0, Ball(), rtol = tol, atol = tol), ϵ, parent(u0.α)(π),
                        rtol = rtol,
                        absmax = true,
                        maxevals = 10^3,
                        show_trace = show_trace,
                        )

    return max(n1, n2)
end

function CB_bounded_by(u0::FractionalKdVAnsatz{arb},
                       C::arb;
                       ϵ::arb = parent(u0.α)(0.1),
                       rtol = 1e-6,
                       atol = 1e-6,
                       show_trace = false,
                       )
    res1 = bounded_by(T0(u0, Asymptotic()), zero(u0.α), ϵ, C,
                      show_trace = show_trace,
                      )

    res1 || return false

    res2 = bounded_by(T0(u0, Ball(), rtol = rtol, atol = atol), ϵ, parent(u0.α)(π), C,
                      show_trace = show_trace, start_intervals = 256)
    return res2
end


"""
    CB_estimate(u0::FractionalKdVAnsatz; n::Integer = 100)
Estimate the value of C_B from the paper. Does this by evaluating the
norm it on n linearly spaced points. This always gives a lower bound
of C_B. Currently it gives a very bad estimate.

If `return_values = true` then also return the points and the computed
values.

TODO: Use the asymptotic version close to zero.
"""
function CB_estimate(u0::FractionalKdVAnsatz{T};
                     n::Integer = 20,
                     return_values = false,
                     show_trace = false,
                     ) where {T}
    xs = range(0, stop = π, length = n)[2:end]
    if T == arb
        xs = parent(u0.α).(xs)
    end

    res = similar(xs)
    f = T0(u0, Ball(), rtol = 1e-6, atol = 1e-6; show_trace)
    Threads.@threads for i in eachindex(xs)
        res[i] = f(xs[i])
    end

    m = zero(u0.α)
    for r in res
        m = max(m, abs(r))
    end

    if return_values
        return m, xs, res
    end
    return m
end
