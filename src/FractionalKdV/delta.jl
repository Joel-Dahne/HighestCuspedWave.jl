export delta0, delta0_bounded_by, delta0_estimate

"""
    delta0(u0::FractionalKdVAnsatz; ϵ::arb = 0, M::Integer = 3, n::Integer = 6)
Upper bound the value of `δ₀` from the paper.

Uses an asymptotic expansion with `M` terms on the interval `[0, ϵ]`
and ball arithmetic together with Taylor expansions of degree `n` on
`[ϵ, π]`.

If `ϵ` is zero then it tries to determine an optimal value of it by
finding where the asymptotic expansion gives better bounds than ball
arithmetic. This is mostly relevant when `u0.α` is wide.

TODO: Look more at tuning for `n`, it seems to depend on the precise
values. It might be beneficial to reduce `ϵ` in smaller steps. It's
mainly the `[ϵ, π]` that takes time, could look more at what part of
it. Otherwise it's mainly important to look at the overestimation in
the evaluation due to the radius of `α`.
"""
function delta0(u0::FractionalKdVAnsatz{arb};
                ϵ::arb = zero(u0.α),
                M::Integer = 3,
                n::Integer = 6,
                rtol::arb = parent(u0.α)(1e-3),
                show_trace = false,
                )
    f = F0(u0, Asymptotic(), M = M)
    g = F0(u0)
    if iszero(ϵ)
        # Find a value of ϵ such that the asymptotic error is smaller
        # than the direct one
        ϵ = one(ϵ)
        while !(radius(f(ϵ)) < radius(g(ϵ))) && ϵ > 1e-3
            ϵ /= 2
        end
    end

    # Bound the value one [0, ϵ]
    # Estimate the value by evaluating it at ϵ
    estimate = abs(f(ϵ))
    res1 = enclosemaximum(f, parent(u0.α)(0), ϵ,
                          lower_bound = estimate,
                          rtol = rtol,
                          atol = 2radius(estimate),
                          absmax = true,
                          maxevals = 10^4,
                          show_trace = show_trace,
                          )

    # Bound the value on [ϵ, π] by Ball evaluation
    # Estimate the value by evaluating it at a number of points on the interval
    xs = ϵ .+ (parent(ϵ)(π) - ϵ).*range(0, stop = 1, length = 10)
    estimate = zero(ϵ)
    for x in xs
        estimate = max(estimate, abs(g(x)))
    end
    res2 = enclosemaximum(g, ϵ, parent(u0.α)(π),
                          evaltype = :taylor,
                          n = n,
                          lower_bound = estimate,
                          rtol = rtol,
                          atol = 4radius(estimate),
                          absmax = true,
                          maxevals = 10^4,
                          show_trace = show_trace,
                          )

    return max(res1, res2)
end

"""
    delta0_bounded_by(u0::FractionalKdVAnsatz, C::arb; ϵ::arb = parent(u0.α)(0.1))
Return true if `δ₀` is bounded by `C`. Uses an asymptotic expansion on
the interval [0, ϵ] and ball arithmetic on [ϵ, π].
"""
function delta0_bounded_by(u0::FractionalKdVAnsatz{arb},
                           C::arb;
                           ϵ::arb = parent(u0.α)(0.1),
                           M::Integer = 3,
                           show_trace = false,
                           )
    # TODO: Spawn these in separate threads?

    # Prove bound on [0, ϵ]
    res1 = bounded_by(F0(u0, Asymptotic(), M = M), parent(u0.α)(0), ϵ, C,
                      show_trace = show_trace,
                      )

    res1 || return false

    # Bound the value on [ϵ, π] by Ball evaluation
    # TODO: So far this naive approach only works for a very limited
    # range of parameters and needs to be improved.
    res2 = bounded_by(F0(u0), ϵ, parent(u0.α)(π), C,
                      show_trace = show_trace,
                      use_taylor = true,
                      )

    return res2
end

"""
    delta0_estimate(u0::FractionalKdVAnsatz; n::Integer = 100)
Estimate the value of δ₀ from the paper. Does this by evaluating F0 it
on n linearly spaced points on the interval. Uses an asymptotic
expansion on the interval [0, ϵ] and ball arithmetic on [ϵ, π]. This
always gives a lower bound of δ₀.

If `return_values = true` then also return the points and the computed
values.
"""
function delta0_estimate(u0::FractionalKdVAnsatz{T};
                         ϵ = 0,
                         M::Integer = 3,
                         n::Integer = 100,
                         return_values = false,
                         ) where {T}
    xs = range(0, stop = pi, length = n)[2:end]
    if T == arb
        xs = parent(u0.α).(xs)
    end

    res = similar(xs)
    Threads.@threads for i in eachindex(xs)
        x = xs[i]
        if x < ϵ
            res[i] = F0(u0, Asymptotic(), M = M)(x)
        else
            res[i] = F0(u0, Ball())(x)
        end
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
