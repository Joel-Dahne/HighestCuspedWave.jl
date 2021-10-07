"""
    delta0(u0::FractionalKdVAnsatz; ϵ::Arf = 0, M::Integer = 3, degree::Integer = 6)

Enclose the value of `δ₀` from the paper.

Uses an asymptotic expansion with `M` terms on the interval `[0, ϵ]`
and ball arithmetic together with Taylor expansions of degree `degree`
on `[ϵ, π]`.

If `ϵ` is zero then it tries to determine an optimal value for it by
finding where the asymptotic expansion gives better bounds than ball
arithmetic. This is mostly relevant when `u0.α` is wide.

If `threaded = true` it enables threading when calling
[`ArbExtras.maximum_enclosure`](@ref). If `verbose = true` print more
information about the process.

TODO: Look more at tuning for `n`, it seems to depend on the precise
values. It might be beneficial to reduce `ϵ` in smaller steps. It's
mainly the `[ϵ, π]` that takes time, could look more at what part of
it. Otherwise it's mainly important to look at the overestimation in
the evaluation due to the radius of `α`.
"""
function delta0(
    u0::FractionalKdVAnsatz{Arb};
    ϵ::Arf = zero(Arf),
    M::Integer = 3,
    degree::Integer = 6,
    rtol = 1e-3,
    threaded = true,
    verbose = false,
)
    # For asymptotic evaluation
    f = F0(u0, Asymptotic(); M)
    # For non-asymptotic evaluation
    g = F0(u0)

    if iszero(ϵ)
        # Find a value of ϵ such that the asymptotic error is smaller
        # than the direct one
        ϵ = one(ϵ)
        while !(Arblib.radius(f(Arb(ϵ))) < Arblib.radius(g(Arb(ϵ)))) && ϵ > 1e-3
            ϵ /= 2
        end

        verbose && @info "ϵ was determined to be" ϵ
    end

    # Bound the value one [0, ϵ]
    # Estimate the value by evaluating it at ϵ
    estimate = abs(f(Arb(ϵ)))
    res1 = ArbExtras.maximum_enclosure(
        f,
        zero(ϵ),
        ϵ,
        degree = -1,
        abs_value = true,
        point_value_max = estimate,
        atol = 4Arblib.radius(estimate); # We can expect to do better than this
        rtol,
        threaded,
        verbose,
    )

    # Bound the value on [ϵ, π] by Ball evaluation
    estimate = maximum(abs.(g.(range(Arb(ϵ), Arblib.ubound(Arb, Arb(π)), length = 10))))
    res2 = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        Arblib.ubound(Arb(π)),
        abs_value = true,
        point_value_max = estimate,
        atol = 4Arblib.radius(estimate); # We can expect to do better than this
        rtol,
        degree,
        threaded,
        verbose,
    )

    return max(res1, res2)
end

"""
    delta0_bounded_by(u0::FractionalKdVAnsatz, C::arb; ϵ::arb = parent(u0.α)(0.1))
Return true if `δ₀` is bounded by `C` together with an enclosure of
the maximum. Uses an asymptotic expansion on the interval [0, ϵ] and
ball arithmetic on [ϵ, π].
"""
function delta0_bounded_by(
    u0::FractionalKdVAnsatz{arb},
    C::arb;
    ϵ::arb = zero(u0.α),
    M::Integer = 3,
    n::Integer = 6,
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

    # Prove bound on [0, ϵ]
    res1, enclosure1 =
        bounded_by(f, parent(u0.α)(0), ϵ, C, return_enclosure = true; show_trace)

    res1 || return false, parent(C)(NaN)

    # Bound the value on [ϵ, π] by Ball evaluation
    res2, enclosure2 = bounded_by(
        F0(u0),
        ϵ,
        parent(u0.α)(π),
        C,
        use_taylor = true,
        return_enclosure = true;
        show_trace,
        n,
    )

    return res2, max(enclosure1, enclosure2)
end
