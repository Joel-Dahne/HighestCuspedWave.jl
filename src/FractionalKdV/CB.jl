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
    # This is not implemented yet!
    # TODO: It might be worth it to try and find maximal value of ϵ
    # such that the bound can be shown with the asymptotic version.
    res1, enclosure1 = bounded_by(T0(u0, Asymptotic()), zero(u0.α), ϵ, C,
                              show_trace = false,
                              return_enclosure = true,
                      )
    res1 || return false, parent(C)(NaN)

    res2, enclosure2 = bounded_by(T0(u0, Ball(), rtol = rtol, atol = atol), ϵ, parent(u0.α)(π), C,
                              show_trace = show_trace, start_intervals = 256,
                              return_enclosure = true)
    return res2, max(enclosure1, enclosure2)
end
