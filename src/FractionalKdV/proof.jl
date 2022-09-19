"""
    prove(u0::FractionalKdVAnsatz{Arb}; M = 5; only_estimate_D0, D0_maxevals, threaded, verbose, extra_verbose)

Attempts to prove that the ansatz `u0` satisfies the requirements, that is
```
δ₀ < 1 / (4n₀ * β^2)
```
with `β = inv(1 - D₀)` and `n₀` is given by [`n0_bound`](@ref), `δ₀`
by [`delta0_bound`](@ref) and `D₀` by [`D0_bound`](@ref).

The most expensive part is the computation of `D₀`. Therefore it
first computes enclosures of `n₀` and `δ₀` but only an estimate of
`D₀` and checks if the condition holds. If the condition holds it
tries to prove that `D₀` is bounded by
```
1 - 2sqrt(n₀ * δ₀)
```
which is equivalent to the former inequality.

If `only_estimate_D0 = true` it doesn't attempt to prove the bound on
`D₀` but only uses the estimate. This doesn't give a rigorous proof
but is useful if you only want to determine of the bound seems to
hold.

If `threaded = true` it uses threading to speed up the computations.

If `verbose = true` it prints information about the progress. For even
more information `extra_verbose` can also be set to be true.
"""
function prove(
    u0::FractionalKdVAnsatz{Arb};
    M = 10,
    only_estimate_D0 = false,
    D0_maxevals = 4000,
    threaded = true,
    verbose = false,
    extra_verbose = false,
)
    # We don't wont to have to go further than a too large depth when
    # bounding D₀. Therefore we want to check the D₀ values using x
    # values with a radius similar to that we would have at a depth we
    # are comfortable going to. Typically we have to go to around
    # depth = 12 but we could go a bit longer than that without
    # issues. So we take it to be 15.
    x_error = Arblib.mul_2exp!(zero(Arb), Arb(π), -15)

    n₀_time = @elapsed n₀ = n0_bound(u0, verbose = extra_verbose; M, threaded)
    verbose && @info "Computed n₀" n₀ n₀_time

    δ₀_time = @elapsed δ₀ = delta0_bound(u0, verbose = extra_verbose; M, threaded)
    verbose && @info "Computed δ₀" δ₀ δ₀_time

    D₀_estimate_time = @elapsed D₀_estimate = D0_estimate(u0; M, x_error, threaded)
    verbose && @info "Computed D₀ estimate" D₀_estimate D₀_estimate_time

    β_estimate = 1 / (1 - D₀_estimate)

    # Bound that δ₀ needs to satisfy
    C_estimate = 1 / (4n₀ * β_estimate^2)

    # Bound that D₀ needs to satisfy
    D = 1 - 2Arblib.sqrtpos!(zero(n₀), n₀ * δ₀)

    if verbose
        @info "Must have δ₀ < 1 / (4n₀ * β^2) = C or equivalently D₀ < 1 - 2√(n₀δ₀) = D" lbound(
            C_estimate,
        ) δ₀ < C_estimate lbound(D) D₀_estimate < D
    end

    if !(δ₀ < C_estimate)
        proved = false
        proved_estimate = false
        D₀ = indeterminate(Arb)
        D₀_time = NaN
    elseif only_estimate_D0
        proved = false
        proved_estimate = true
        D₀ = indeterminate(Arb)
        D₀_time = NaN
    else
        proved_estimate = true
        D₀ = lbound(Arb, D) - sqrt(eps()) # Add a little bit of head room
        D₀_time = @elapsed proved =
            D0_bounded_by(u0, lbound(D₀), maxevals = D0_maxevals; M, threaded, verbose)

        verbose && @info "Bounded D₀" D₀ D₀_time
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
        u0_N0 = u0.N0,
        u0_N1 = u0.N1,
    )
end
