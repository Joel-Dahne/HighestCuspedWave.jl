"""
    prove(u0::KdVZeroansatz{Arb}; only_estimate_CB, threaded, verbose, extra_verbose)
"""
function prove(
    u0::KdVZeroAnsatz;
    M = 5, # Not used
    only_estimate_CB = false,
    maxevals = 1000,
    threaded = true,
    verbose = false,
    extra_verbose = false,
)
    iszero(u0.α0) && @warn "this method is not intended for u0.α0 = 0"

    # Construct corresponding Fractionalkdvansatz
    u02 = FractionalKdVAnsatz{Arb}(
        u0.α,
        u0.p0(u0.α - u0.α0),
        map(a -> a(u0.α - u0.α0), u0.a),
        Arb[],
        one(Arb),
        Set{NTuple{3,Int}}([(2, 0, 0)]),
    )

    α₀_time = @elapsed α₀ = alpha0(u0, verbose = extra_verbose)
    verbose && @info "Computed α₀" α₀ α₀_time

    # The maximum of T0 is in practice attained at x = 0
    C_B_estimate_time = @elapsed C_B_estimate = T0(u02, Asymptotic())(Arb(0))
    verbose && @info "Computed C_B estimate" C_B_estimate C_B_estimate_time

    # Compute estimate of required bound of defect
    β_estimate = 1 / (1 - C_B_estimate)
    δ₀_goal = 1 / (4α₀ * β_estimate^2)

    verbose && @info "Required bound for defect $δ₀_goal"

    if !isfinite(δ₀_goal)
        proved = false
        δ₀ = Arblib.indeterminate!(zero(Arb))
        C_B = Arblib.indeterminate!(zero(Arb))
        δ₀_time = NaN
        C_B_time = NaN
    else
        δ₀_time = @elapsed δ₀ = delta0(
            u0,
            atol = 0.5δ₀_goal,
            maxevals = 1000,
            verbose = extra_verbose;
            threaded,
        )
        verbose && @info "Computed δ₀" δ₀ δ₀_time

        # Compute required bound for C_B, adding a little bit of head room
        C_B_goal = lbound(Arb, 1 - 2Arblib.sqrtpos!(zero(Arb), α₀ * δ₀)) - sqrt(eps())

        verbose && @info "Required bound for C_B" C_B_goal

        if !(C_B_estimate < C_B_goal)
            proved = false
            C_B = Arblib.indeterminate!(zero(Arb))
            C_B_time = NaN

            verbose && @warn "Required bound for C_B not satisfied at x = 0"
        elseif only_estimate_CB
            proved = false
            C_B = Arblib.indeterminate!(zero(Arb))
            C_B_time = NaN
        else
            C_B = C_B_goal
            C_B_time = @elapsed proved =
                CB_bounded_by(u02, lbound(C_B); maxevals, threaded, verbose)
        end
    end

    return (;
        proved,
        α₀,
        δ₀,
        C_B_estimate,
        C_B,
        α₀_time,
        δ₀_time,
        C_B_estimate_time,
        C_B_time,
    )
end
