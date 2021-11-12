export CB, CB_bounded_by, CB_estimate

"""
    CB(u0::AbstractAnsatz; ϵ = 0.1)

Specialized for different types of ansatz.
"""
CB

"""
    CB_bounded_by(u0::AbstractAnsatz, C; ϵ = 0.1, ...)

Specialized for different types of ansatz
"""
CB_bounded_by


"""
    CB_estimate(u0::FractionalKdVAnsatz; n::Integer = 20, x_error = 0, include_zero = true)

Estimate the value of `C_B` from the paper. Does this by evaluating
the norm of it on `n` linearly spaced points. This always gives a
lower bound of `C_B`.

Notice that the default value for `n` is rather low since evaluation
of the integral is costly. Therefore it gives a rather poor estimate
of the norm.

If `return_values = true` then also return the points and the computed
values.

If `x_error` is non-zero and the type used is `Arb` then add the error
`x_error` to each value of `x`. This is useful if you want to assure
that the bounds holds even if the value for `x` is not a tight ball.

If `include_zero = true` then also evaluate the integral at `x = 0`
using the asymptotic expansion.

If `threaded = true` computes the evaluations in parallel.
"""
function CB_estimate(
    u0::AbstractAnsatz{T};
    n::Integer = 20,
    return_values = false,
    x_error = zero(T),
    include_zero = true,
    threaded = true,
) where {T}
    xs = collect(range(zero(T), π, length = n + 1)[2:end])

    if !iszero(x_error)
        @assert T == Arb
        Arblib.add_error!.(xs, x_error)
    end

    res = similar(xs)

    f = T0(u0, Ball())
    if threaded
        Threads.@threads for i in eachindex(xs)
            res[i] = f(xs[i])
        end
    else
        for i in eachindex(xs)
            res[i] = f(xs[i])
        end
    end

    if include_zero
        pushfirst!(xs, zero(T))
        pushfirst!(res, T0(u0, Asymptotic())(zero(T)))
    end

    m = maximum(abs.(res))

    return_values && return m, xs, res
    return m
end
