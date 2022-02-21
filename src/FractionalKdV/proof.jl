"""
    prove(u0::FractionalKdVAnsatz{Arb}; M = 5; only_estimate_CB, threaded, verbose, extra_verbose)

Attempts to prove that the ansatz `u0` satisfies the requirements, that is
```
δ₀ <= 1 / (4α₀ * β^2)
```
with `β = inv(1 - C_B)` and `α₀` is given by [`alpha0`](@ref), `δ₀` by
[`delta0`](@ref) and `C_B` by [`CB`](@ref).

The most expensive part is the computation of `C_B`. Therefore it
first computes enclosures of `α₀` and `δ₀` but only an estimate of
`C_B` and checks if the condition holds. If the condition holds it
tries to prove that `C_B` is bounded by
```
1 - 2sqrt(α₀ * δ₀)
```
which is equivalent to the former inequality.

If `only_estimate_CB = true` it doesn't attempt to prove the bound on
`C_B` but only uses the estimate. This doesn't give a rigorous proof
but is useful if you only want to determine of the bound seems to
hold.

If `threaded = true` it uses threading to speed up the computations.

If `verbose = true` it prints information about the progress. For even
more information `extra_verbose` can also be set to be true.
"""
function prove(
    u0::FractionalKdVAnsatz{Arb};
    M = 5,
    only_estimate_CB = false,
    threaded = true,
    verbose = false,
    extra_verbose = false,
)
    # We don't wont to have to go further than a too large depth when
    # bounding C_B. Therefore we want to check the C_B values using x
    # values with a radius similar to that we would have at a depth we
    # are comfortable going to. Typically we have to go to around
    # depth = 12 but we could go a bit longer than that without
    # issues. So we take it to be 15.
    x_error = Arblib.mul_2exp!(zero(Arb), Arb(π), -15)

    α₀_time = @elapsed α₀ = alpha0(u0, verbose = extra_verbose; M, threaded)
    verbose && @info "Computed α₀" α₀ α₀_time

    δ₀_time = @elapsed δ₀ = delta0(u0, verbose = extra_verbose; M, threaded)
    verbose && @info "Computed δ₀" δ₀ δ₀_time

    C_B_estimate_time = @elapsed C_B_estimate = CB_estimate(u0; x_error, threaded)
    verbose && @info "Computed C_B estimate" C_B_estimate C_B_estimate_time

    β_estimate = 1 / (1 - C_B_estimate)

    # Bound that δ₀ needs to satisfy
    C_estimate = 1 / (4α₀ * β_estimate^2)

    # Bound that C_B needs to satisfy
    D = 1 - 2Arblib.sqrtpos!(zero(α₀), α₀ * δ₀)

    if verbose
        @info "Must have δ₀ ≤ 1 / (4α₀ * β^2) = C or equivalently C_B ≤ 1 - 2√(α₀δ₀) = D" lbound(
            C_estimate,
        ) δ₀ < C_estimate lbound(D) C_B_estimate < D
    end

    if !(δ₀ < C_estimate)
        proved = false
        proved_estimate = false
        C_B = Arblib.indeterminate!(zero(Arb))
        C_B_time = NaN
    elseif only_estimate_CB
        proved = false
        proved_estimate = true
        C_B = Arblib.indeterminate!(zero(Arb))
        C_B_time = NaN
    else
        proved_estimate = true
        C_B = lbound(Arb, D) - sqrt(eps()) # Add a little bit of head room
        C_B_time = @elapsed proved = CB_bounded_by(u0, lbound(C_B); threaded, verbose)

        verbose && @info "Bounded C_B" C_B C_B_time
    end

    return (;
        proved,
        proved_estimate,
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

function format_for_publishing(α₀, δ₀, C_B)
    α₀_float = Arblib.get_d(ubound(α₀), RoundUp)
    δ₀_float = Arblib.get_d(ubound(δ₀), RoundUp)
    C_B_float = Arblib.get_d(ubound(C_B), RoundUp)

    # Check that the inequality holds
    inequality_holds = Arb(δ₀_float) <= (1 - Arb(C_B_float))^2 / 4Arb(α₀)

    return inequality_holds, α₀_float, δ₀_float, C_B_float
end
