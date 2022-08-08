"""
    delta0_bound(u0::FractionalKdVAnsatz; ϵ, M, degree, rtol, threaded, verbose)

Enclose the value of `δ₀` from the paper.

Uses an asymptotic expansion with `M` terms on the interval `[0, ϵ]`
and ball arithmetic together with Taylor expansions of degree `degree`
on `[ϵ, π]`.

# Arguments
- `ϵ::Arf = zero(Arf)`: Determines the interval ``[0, ϵ]`` on which
  the asymptotic version is used. If `ϵ` is zero it automatically
  finds a suitable value for it.
- `M::Integer = 3`: Number of terms to use in the asymptotic
  expansions.
- `degree::Integer = 6`: Degree used for
  [`ArbExtras.maximum_enclosure`](@ref).
- `rtol = Arb(1e-3)`: Relative tolerance used when computing maximum.
- `threaded = true`: If true it enables threading when calling
  [`ArbExtras.maximum_enclosure`](@ref).
- `verbose = false`: Print information about the process.
"""
function delta0_bound(
    u0::FractionalKdVAnsatz{Arb};
    ϵ::Arf = zero(Arf),
    M::Integer = 3,
    degree::Integer = 6,
    rtol = Arb(1e-3),
    threaded = true,
    verbose = false,
)
    # For asymptotic evaluation
    f = F0(u0, Asymptotic(); M, ϵ = Arb(1.1))
    # For non-asymptotic evaluation
    g = F0(u0)

    if iszero(ϵ)
        ϵ = let ϵ = Arb(1)
            y = f(ϵ)
            z = g(ϵ)

            # Reduce ϵ until the value we get either satisfies the
            # required tolerance or is better than the non-asymptotic
            # version.
            while !ArbExtras.check_tolerance(y; rtol) && radius(y) > radius(z)
                ϵ /= 1.2

                y = f(ϵ)
                z = g(ϵ)
                if ϵ < 1e-3
                    @warn "could not determine working ϵ, last tried value" ϵ u0.α
                    break
                end
            end
            ubound(ϵ)
        end

        verbose && @info "ϵ was determined to be" ϵ
    end

    # Bound the value on [0, ϵ]
    # Estimate the value by evaluating it at ϵ.
    # Asymptotic evaluation is cheap and we are therefore not too
    # worried about extra evaluations, the absolute tolerance is
    # therefore set rather low.
    estimate = abs(f(Arb(ϵ)))
    max_asymptotic = ArbExtras.maximum_enclosure(
        f,
        zero(ϵ),
        ϵ,
        abs_value = true,
        point_value_max = estimate,
        atol = 1.1radius(estimate); # We can't expect to do much better than this
        degree,
        rtol,
        threaded,
        verbose,
    )

    verbose && @info "Bound on [0, ϵ]" max_asymptotic

    # Bound the value on [ϵ, π] by Ball evaluation
    # Ball evaluation is rather expensive and we are therefore worried
    # about extra evaluations, the absolute tolerance is therefore set
    # a bit higher.
    estimate = maximum(abs.(g.(range(Arb(ϵ), π, length = 10))))
    max_nonasymptotic = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        ubound(Arb(π)),
        abs_value = true,
        point_value_max = estimate,
        atol = max(max_asymptotic / 2, 2radius(estimate));
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
