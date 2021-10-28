"""
    prove(u0::FractionalKdVAnsatz{Arb}; M = 3, , only_estimate_CB = false, return_values = false, threaded = true, verbose = false, extra_verbose = false)

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

If `return_values = true` it returns the value of `δ₀` and the lower
bound of `1 / (4α₀ * β^2)` that is proved. If `only_estimate_CB =
true` it instead returns the estimate of it.

If `threaded = true` it uses threading to speed up the computations.

If `verbose = true` it prints information about the progress. For even
more information `extra_verbose` can also be set to be true.
"""
function prove(
    u0::FractionalKdVAnsatz{Arb};
    M = 5,
    only_estimate_CB = false,
    return_values = false,
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

    if verbose
        @time α₀ = alpha0(u0, verbose = extra_verbose; threaded)
        @show α₀

        @time δ₀ = delta0(u0, verbose = extra_verbose; M, threaded)
        @show δ₀

        @time C_B = CB_estimate(u0; x_error, threaded)
        @show C_B
    else
        α₀ = alpha0(u0; threaded)
        δ₀ = delta0(u0; M, threaded)
        C_B = CB_estimate(u0; x_error, threaded)
    end

    β = 1 / (1 - C_B)
    C = 1 / (4α₀ * β^2)
    D = 1 - 2Arblib.sqrtpos!(zero(C), α₀ * δ₀)

    if verbose
        @info "Must have δ₀ ≤ 1 / (4α₀ * β^2) = C" δ₀ lbound(C) δ₀ < C
        @info "Alternatively C_B ≤ 1 - 2√(α₀δ₀) = D" C_B lbound(D) C_B < D
    end

    if !(δ₀ < C)
        if return_values
            return false, δ₀, C
        else
            return false
        end
    end

    if only_estimate_CB
        return true, δ₀, C
    end

    if verbose
        @time res = CB_bounded_by(u0, lbound(D); threaded, verbose)
    else
        res = CB_bounded_by(u0, lbound(D); threaded, verbose)
    end

    C_B_bound = lbound(D)
    β_bound = 1 / (1 - C_B_bound)
    C_bound = 1 / (4α₀ * β_bound^2)

    return res, δ₀, C_bound
end
