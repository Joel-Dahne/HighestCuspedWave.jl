"""
    delta0_bound(u0::FractionalKdVAnsatz; ϵ, M, degree, rtol, threaded, verbose)

Enclose the value of `δ₀` from the paper. This is the supremum of
`F0(u0)` for `0 < x < π`.

Uses an asymptotic expansion with `M` terms close to zero and ball
arithmetic on the remaining.

The interval is split into three subintervals ``[0, ϵ1]``, ``[ϵ1,
ϵ2]`` and ``[ϵ2, π]``. On the first two the asymptotic expansion is
used and on the third ball arithmetic.

The values for `ϵ1` and `ϵ2` are determined automatically. The value
of `ϵ2` is taken such that evaluation with the asymptotic expansion
either satisfies the required tolerance or gives a smaller error than
evaluation with ball arithmetic. The value of `ϵ1` is taken such that
the interval ``[0, ϵ1]`` can be handled in one evaluation.

# Arguments
- `M::Integer = 10`: Number of terms to use in the asymptotic
  expansions.
- `degree::Integer = ifelse(u0.N0 > 100, 4, 6)`: Degree used for
  [`ArbExtras.maximum_enclosure`](@ref).
- `rtol = Arb(1e-3)`: Relative tolerance used when computing maximum.
- `ubound_tol = nothing: Any number less than this is determined to
  satisfy the tolerance. Useful if you only need to determine if the
  maximum is less than some given number. If set to `nothing` a
  default value is used depending on `u0`.
- `threaded = true`: If true it enables threading when calling
  [`ArbExtras.maximum_enclosure`](@ref).
- `verbose = false`: Print information about the process.
"""
function delta0_bound(
    u0::FractionalKdVAnsatz{Arb};
    M::Integer = 10,
    degree::Integer = ifelse(u0.N0 > 100, 4, 6),
    rtol = Arb(1e-3),
    ubound_tol = nothing,
    threaded = true,
    verbose = false,
)
    if isnothing(ubound_tol)
        # These values have been determined by testing what seems to
        # work well in practice.
        if u0.use_bhkdv
            if Float64(u0.α) < -0.99985
                ubound_tol = Arb(0.00025)
            elseif Float64(u0.α) < -0.99
                ubound_tol = Arb(0.0002)
            elseif Float64(u0.α) < -0.9
                ubound_tol = Arb(0.0004)
            else
                ubound_tol = Arb(-Inf)
            end
        else
            ubound_tol = Arb(-Inf)
        end
    end

    # For asymptotic evaluation
    f = F0(u0, Asymptotic(); M, ϵ = Arb(2.1))
    # For non-asymptotic evaluation
    g = F0(u0)

    # Find ϵ2 such that f(ϵ2) satisfies the required tolerance or has
    # a smaller radius than g(ϵ2)
    ϵ2 = Arb(2)
    fϵ2 = abs(f(ϵ2))
    gϵ2 = abs(g(ϵ2))
    while !ArbExtras.check_tolerance(fϵ2; rtol) &&
              !(fϵ2 < ubound_tol) &&
              radius(fϵ2) > radius(gϵ2)

        ϵ2 *= 0.8
        fϵ2 = abs(f(ϵ2))
        gϵ2 = abs(g(ϵ2))

        if ϵ2 < 1e-3
            verbose && @warn "Could not find working ϵ2" ϵ2 fϵ2 gϵ2
            break
        end
    end

    verbose && @info "Determined ϵ2" ϵ2

    # Set tolerances according to currently evaluated values
    atol = 1.1radius(fϵ2)
    ubound_tol = max(ubound_tol, ubound(Arb, fϵ2))

    # Find ϵ1 such that f(Arb(0, ϵ1)) satisfies the required tolerance
    ϵ1 = ϵ2
    fϵ1 = abs(f(Arb((0, ϵ1))))
    while !ArbExtras.check_tolerance(fϵ1; rtol, atol) && !(fϵ1 < ubound_tol)

        ϵ1 = ϵ1 < 0.5 ? ϵ1^2 : ϵ1 / 2
        fϵ1 = abs(f(Arb((0, ϵ1))))

        if ϵ1 < Arb("1e-100000")
            verbose && @error "Could not find working ϵ1" ϵ1 fϵ1
            return indeterminate(fϵ1)
        end
    end

    verbose && @info "Determined ϵ1" ϵ1

    # Bound the value on [0, ϵ]
    max_asymptotic = ArbExtras.maximum_enclosure(
        f,
        lbound(ϵ1),
        ubound(ϵ2),
        abs_value = true,
        degree = ifelse(u0.use_bhkdv, -1, degree),
        log_bisection = true,
        depth = ifelse(u0.use_bhkdv, 40, 20),
        maxevals = ifelse(u0.use_bhkdv, 50000, 5000);
        rtol,
        atol,
        ubound_tol,
        threaded,
        verbose,
    )

    # Take into account enclosure on [0, ϵ1]
    max_asymptotic = max(max_asymptotic, fϵ1)

    verbose && @info "Bound on [0, ϵ2]" max_asymptotic

    # Ball evaluation is expensive and we are worried about extra
    # evaluations, the tolerances are therefore set rather high.
    estimate = max(gϵ2, abs(g(Arb(π))))
    atol = max(max_asymptotic / 2, 1.5radius(gϵ2))
    ubound_tol = max(ubound_tol, ubound(Arb, max_asymptotic), ubound(Arb, estimate))

    # Bound the value on [ϵ, π] by Ball evaluation
    max_nonasymptotic = ArbExtras.maximum_enclosure(
        g,
        lbound(ϵ2),
        ubound(Arb(π)),
        abs_value = true,
        depth_start = 2;
        rtol,
        atol,
        ubound_tol,
        degree,
        threaded,
        verbose,
    )

    verbose && @info "Bound on [ϵ2, π]" max_nonasymptotic

    res = max(max_asymptotic, max_nonasymptotic)
    # The result is always nonnegative
    return Arblib.nonnegative_part!(res, res)
end
