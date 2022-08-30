"""
    prove(u0::KdVZeroansatz{Arb}; only_estimate_D0, threaded, verbose, extra_verbose)
"""
function prove(
    u0::KdVZeroAnsatz;
    M = 5, # Not used
    only_estimate_D0 = false,
    maxevals = 1000,
    threaded = true,
    verbose = false,
    extra_verbose = false,
)
    iszero(u0.α0) && @warn "this method is not intended for u0.α0 = 0"

    # Construct corresponding Fractionalkdvansatz
    u02 = FractionalKdVAnsatz{Arb}(
        u0.α,
        u0.p0(u0.α),
        map(a -> a(u0.α), u0.a),
        Arb[],
        one(Arb),
        Set{NTuple{3,Int}}([(2, 0, 0)]),
    )

    n₀_time = @elapsed n₀ = n0_bound(u0, verbose = extra_verbose)
    verbose && @info "Computed n₀" n₀ n₀_time

    # The maximum of T0 is in practice attained at x = 0
    D₀_estimate_time = @elapsed D₀_estimate = T0(u02, Asymptotic())(Arb(0))
    verbose && @info "Computed D₀ estimate" D₀_estimate D₀_estimate_time

    # Compute estimate of required bound of defect
    β_estimate = 1 / (1 - D₀_estimate)
    δ₀_goal = 1 / (4n₀ * β_estimate^2)

    verbose && @info "Required bound for defect $δ₀_goal"

    if !isfinite(δ₀_goal)
        proved = false
        proved_estimate = false
        δ₀ = indeterminate(Arb)
        D₀ = indeterminate(Arb)
        δ₀_time = NaN
        D₀_time = NaN
    else
        δ₀_time = @elapsed δ₀ = delta0_bound(
            u0,
            atol = 0.5δ₀_goal,
            maxevals = 1000,
            verbose = extra_verbose;
            threaded,
        )
        verbose && @info "Computed δ₀" δ₀ δ₀_time

        # Compute required bound for D₀, adding a little bit of head room
        D₀_goal = lbound(Arb, 1 - 2Arblib.sqrtpos!(zero(Arb), n₀ * δ₀)) - sqrt(eps())

        verbose && @info "Required bound for D₀" D₀_goal

        if !(D₀_estimate < D₀_goal)
            proved = false
            proved_estimate = false
            D₀ = indeterminate(Arb)
            D₀_time = NaN

            verbose && @warn "Required bound for D₀ not satisfied at x = 0"
        elseif only_estimate_D0
            proved = false
            proved_estimate = true
            D₀ = indeterminate(Arb)
            D₀_time = NaN
        else
            proved_estimate = true
            D₀ = D₀_goal
            D₀_time = @elapsed proved =
                D0_bounded_by(u02, lbound(D₀); maxevals, threaded, verbose)
        end
    end

    return (;
        proved,
        proved_estimate,
        n₀,
        δ₀,
        D₀_estimate,
        D₀,
        n₀_time,
        δ₀_time,
        D₀_estimate_time,
        D₀_time,
    )
end
