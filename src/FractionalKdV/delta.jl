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

**TODO:** Look more at tuning for `n`, it seems to depend on the
precise values. It might be beneficial to reduce `ϵ` in smaller steps.
It's mainly the `[ϵ, π]` that takes time, could look more at what part
of it. Otherwise it's mainly important to look at the overestimation
in the evaluation due to the radius of `α`.
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

    # Bound the value on [0, ϵ]
    # Estimate the value by evaluating it at ϵ
    estimate = abs(f(Arb(ϵ)))
    max_asymptotic = ArbExtras.maximum_enclosure(
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

    verbose && @info "Bound on [0, ϵ]" max_asymptotic

    # Bound the value on [ϵ, π] by Ball evaluation
    estimate = maximum(abs.(g.(range(Arb(ϵ), Arblib.ubound(Arb, Arb(π)), length = 10))))
    max_nonasymptotic = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        Arblib.ubound(Arb(π)),
        abs_value = true,
        point_value_max = estimate,
        atol = 3Arblib.radius(estimate); # We can expect to do better than this
        rtol,
        degree,
        threaded,
        verbose,
    )

    verbose && @info "Bound on [ϵ, π]" max_nonasymptotic

    res = max(max_asymptotic, max_nonasymptotic)
    # The result is always nonnegative
    return Arblib.nonnegative_part!(res, res)
end
