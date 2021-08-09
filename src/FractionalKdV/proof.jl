function prove(u0::FractionalKdVAnsatz{arb}; verbose = false)

    α₀ = alpha0(u0)
    δ₀ = delta0(u0)

    if verbose
        @show α₀ δ₀
    end

    C_B = CB_estimate(u0)

    β = 1 / (1 - C_B)
    C = 1 / (4α₀ * β^2)
    D = 1 - 2sqrt(α₀ * δ₀)

    if verbose
        @show C_B C
        println("Must have C_B ≤ 1 - 2√(α₀δ₀) = $D")
    end

    if !(δ₀ < C)
        return false
    end

    res, CB_bound = CB_bounded_by(u0, D, show_trace = true)
    CB_bound = ArbTools.ubound(CB_bound)
    β_bound = 1 / (1 - CB_bound)
    C_bound = 1 / (4α₀ * β_bound^2)
    @show C_bound

    return res, CB_bound
end

function prove2(α::arb; verbose = false)
    @time u0 = FractionalKdVAnsatz(α)

    # Enclose α₀
    @time α₀ = alpha0(u0)

    # Approximate δ₀
    @time δ₀_low = delta0(update_alpha(u0, ArbTools.lbound(α)))
    @time δ₀_mid = delta0(update_alpha(u0, midpoint(α))) # We don't need this in the end
    @time δ₀_upp = delta0(update_alpha(u0, ArbTools.ubound(α)))

    # See if it looks like we can prove it
    CB_est = CB_estimate(u0)
    D_est = 1 - 2sqrt(α₀ * max(δ₀_low, δ₀_upp))

    if !(CB_est < D_est)
        @error "bound doesn't hold, we don't have $CB_est < $D_est"
        return false, ArbTools.ubound(α₀), parent(α)(NaN), parent(α)(NaN)
    end

    # Bound and enclose CB
    p = 0.8 # This is an optimization parameter, p ∈ [0, 1] with p = 0
    # best for δ₀ and p = 1 best for CB
    D = ArbTools.ubound(CB_est) + p * (D_est - ArbTools.ubound(CB_est))
    @show ArbTools.lbound(CB_est) ArbTools.lbound(D_est) ArbTools.lbound(D)
    @time res, CB = CB_bounded_by(u0, D, show_trace = verbose)
    @show CB

    if !res
        @error "failed to prove bound for CB"
        return false, ArbTools.ubound(α₀), parent(α)(NaN), parent(α)(NaN)
    end

    # Bound and enclose δ₀
    # We need to beat whatever we get as the upper bound for CB
    β = 1 / (1 - ArbTools.ubound(CB))
    @show β (4α₀ * β^2)
    C = 1 / (4α₀ * β^2)

    @show max(δ₀_low, δ₀_upp) ArbTools.lbound(C)
    @assert max(δ₀_low, δ₀_upp) < C
    # TODO: We will have to bisect α to be able to bound this
    αs = mince(α, 10)
    δ₀s = similar(αs)
    δ₀ = parent(α)(-Inf)
    for i in eachindex(αs)
        @time res, δ₀s[i] =
            delta0_bounded_by(update_alpha(u0, αs[i]), C, show_trace = verbose)
        δ₀ = max(δ₀, δ₀s[i])
        if !res
            @error "failed to prove bound for δ₀ on α = $(αs[i])"
            return false, ArbTools.ubound(α₀), parent(α)(NaN), ArbTools.ubound(CB)
        end
    end

    α₀_ubound = ArbTools.ubound(α₀)
    δ₀_ubound = ArbTools.ubound(δ₀)
    CB_ubound = ArbTools.ubound(CB)

    @show δ₀_ubound (1 - CB_ubound)^2 / (4α₀_ubound)
    @assert δ₀_ubound < (1 - CB_ubound)^2 / (4α₀_ubound)

    return true, α₀_ubound, δ₀_ubound, CB_ubound
end
