"""
    prove(u0::FractionalKdVAnsatz{Arb}; M = 5; only_estimate_D0, D0_maxevals, threaded, verbose, extra_verbose)

Attempts to prove that the ansatz `u0` satisfies the requirements for
the fixed point proposition, that is
```
δ₀ < (1 - D₀)^2 / 4n₀
```
Where `n₀`, `δ₀` and `D₀` are given by the supremum for `0 < x < π` of
the following three functions
```
N(x) = u0.w(x) / 2u0(x)
F0(u0)
T0(u0)
```

The most expensive part is the computation of `D₀`. Therefore it
first computes enclosures of `n₀` and `δ₀`, using [`n0_bound`](@ref)
and [`delta0_bound`](@ref) respectively, but only an estimate of `D₀`,
using [`D0_estimate`](@ref), and checks if the condition holds. If the
condition holds it tries to prove that `D₀` is bounded by
```
1 - 2sqrt(n₀ * δ₀)
```
which is equivalent to the former inequality, using
[`D0_bounded_by`](@ref).


It returns the result together with a lot of metadata. More precisely
it returns a named tuple of the form
```
(;
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
    u0_N0,
    u0_N1,
)
```
where the fields have the following meaning
- `proved`: True if the inequality is true
- `proved_estimate`: True if the inequality is true with the estimate of `D₀`
- `n₀`: Upper bound of `n₀`
- `δ₀`: Upper bound of `δ₀`
- `D₀_estimate`: Estimate of `D₀`
- `D₀`: Upper bound of `D₀`
- `n₀_time`: Time for computing `n₀`
- `δ₀_time`: Time for computing `δ₀`
- `D₀_estimate_time`: Time for computing `D₀_estimate`
- `D₀_time`: Time for computing `D₀`
- `u0_N0`: The value of `u0.N0`
- `u0_N1`: The value of `u0.N1`

If `only_estimate_D0 = true` it doesn't attempt to prove the bound on
`D₀` but only uses the estimate. This doesn't give a rigorous proof
but is useful to determine if the bound seems to hold.

The maximum number of evaluations used by [`D0_bounded_by`](@ref) can
be set with `D0_maxevals`.

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
    n₀_time = @elapsed n₀ = n0_bound(u0, verbose = extra_verbose; M, threaded)
    verbose && @info "Computed n₀" n₀ n₀_time

    δ₀_time = @elapsed δ₀ = delta0_bound(u0, verbose = extra_verbose; M, threaded)
    verbose && @info "Computed δ₀" δ₀ δ₀_time

    # We don't want to have to go further than a too large depth when
    # bounding D₀. Therefore we want to estimate the D₀ value using x
    # values with a radius similar to that we would have at a depth we
    # are comfortable going to. Typically we have to go to around
    # depth = 12 but we could go a bit longer than that without
    # issues. So we take it to be 15.
    x_error = Arblib.mul_2exp!(zero(Arb), Arb(π), -15)

    D₀_estimate_time = @elapsed D₀_estimate = D0_estimate(u0; M, x_error, threaded)
    verbose && @info "Computed D₀ estimate" D₀_estimate D₀_estimate_time

    # Bound that δ₀ needs to satisfy
    δ₀_goal_estimate = (1 - D₀_estimate)^2 / 4n₀

    # Bound that D₀ needs to satisfy
    D₀_goal = 1 - 2Arblib.sqrtpos!(zero(n₀), n₀ * δ₀)

    if verbose
        @info "Must have δ₀ < (1 - D₀)^2 / 4n₀ or equivalently D₀ < 1 - 2√(n₀δ₀)" lbound(
            δ₀_goal_estimate,
        ) δ₀ < δ₀_goal_estimate lbound(D₀_goal) D₀_estimate < D₀_goal
    end

    if !(δ₀ < δ₀_goal_estimate)
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
        D₀ = lbound(Arb, D₀_goal) - sqrt(eps()) # Add a little bit of head room
        D₀_time = @elapsed proved =
            D0_bounded_by(u0, lbound(D₀), maxevals = D0_maxevals; M, threaded, verbose)

        if proved
            verbose && @info "Bounded D₀" D₀ D₀_time
        else
            verbose && @info "Failed to bound D₀" D₀ D₀_time
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
        u0_N0 = u0.N0,
        u0_N1 = u0.N1,
    )
end
