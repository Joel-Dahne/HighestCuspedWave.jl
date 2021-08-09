export CB, CB_bounded_by, CB_estimate

"""
    CB(u0::AbstractAnsatz; ϵ = 0.1)
Upper bound the value of C_B from the paper. TODO: Uses an asymptotic
expansion on the interval [0, ϵ] and ball arithmetic on [ϵ, π].

Specialized for different types of ansatz.
"""
CB

"""
    CB_bounded_by(u0::AbstractAnsatz, C; ϵ = 0.1, ...)

TODO: Implement this

Specialized for different types of ansatz
"""
CB_bounded_by


"""
    CB_estimate(u0::FractionalKdVAnsatz; n::Integer = 100)
Estimate the value of C_B from the paper. Does this by evaluating the
norm of it on n linearly spaced points. This always gives a lower bound
of C_B. Currently it gives a very bad estimate.

If `return_values = true` then also return the points and the computed
values.

TODO: Use the asymptotic version close to zero.
"""
function CB_estimate(
    u0::AbstractAnsatz{T};
    n::Integer = 20,
    return_values = false,
    show_trace = false,
) where {T}
    xs = range(0, stop = π, length = n)[2:end]
    if T == arb
        xs = u0.parent.(xs)
    else
        xs = convert.(T, xs)
    end

    res = similar(xs)
    f = T0(u0, Ball(), rtol = 1e-6, atol = 1e-6; show_trace)
    Threads.@threads for i in eachindex(xs)
        res[i] = f(xs[i])
    end

    m = zero(xs[1])
    for r in res
        m = max(m, abs(r))
    end

    if return_values
        return m, xs, res
    end
    return m
end
