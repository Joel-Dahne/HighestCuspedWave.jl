export delta0, delta0_estimate

"""
    delta0(u0::AbstractAnsatz{Arb})

Upper bound the value of `δ₀` from the paper.

Specialized for different types of `AbstractAnsatz`.
"""
delta0

"""
    delta0_estimate(u0::AbstractAnsatz; ϵ = 0, n::Integer = 100, return_values = false)

Estimate the value of `δ₀` from the paper.

Does this by evaluating `F0` it on `n` linearly spaced points on the
interval. Uses an asymptotic expansion on the interval `[0, ϵ]` and
ball arithmetic on `[ϵ, π]`. This always gives a lower bound of `δ₀`.

If `return_values = true` then also return the points and the computed
values.
"""
function delta0_estimate(
    u0::AbstractAnsatz{T};
    ϵ = 0,
    M::Integer = 3,
    n::Integer = 100,
    return_values = false,
) where {T}
    xs = range(zero(T), π, length = n + 1)[2:end]
    res = similar(xs)

    # Asymptotic version might not be defined
    f1 = !iszero(ϵ) ? F0(u0, Asymptotic(); M) : missing
    f2 = F0(u0, Ball())
    Threads.@threads for i in eachindex(xs)
        x = xs[i]
        if x < ϵ
            res[i] = f1(x)
        else
            res[i] = f2(x)
        end
    end

    m = maximum(abs.(res))

    if return_values
        return m, xs, res
    end

    return m
end
