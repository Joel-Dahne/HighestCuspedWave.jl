"""
    D0_bound(u0::BHAnsatz; atol, verbose = false)

Compute an upper bound of `D₀` from the paper. This is the supremum of
```
abs(T0(u0)(x))
```
for `0 < x < π`.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use an asymptotic expansion of `T0(u0)` whereas on ´[ϵ,
π]` we use a non-asymptotic version.

In practice the maximum is attained on `[ϵ, π]` and for that reason we
compute the maximum on this part first and then only prove that the
value on `[0, ϵ]` is bounded by this value.

For the interval ``[ϵ, π]`` we compute `T0(u0)` using the argument
`skip_div_u0 = true` We then use the fact that `u0.v0(x)` gives a
**lower** bound for `u0(x)`. This follows from
[`lemma_bhkdv_main_term_limit`](@ref). This means that as long as
`u0.v0(x) > 0` division by `u0.v0(x)` instead of `u0(x)` gives an
upper bound of the value. To prove that `u0.v0(x)` is positive it is
enough to ensure that `u0.v0(π)` is positive and that the final bound
is finite, in this case `u0.v0(x)` can never be zero.

To efficiently compute a tight enclosure of `u0.v0(x)` we use that it
in practice is increasing in `x`. We therefore try to prove that it is
increasing on the interval ``[ϵ, b]``, with `b` slightly smaller than
`π`, by proving that the derivative is non-negative. This allows us to
only evaluate it on the endpoints if `x` lies in this interval.

In general we can get much better enclosures of `u0.v0(x)` than of
`T0(u0)(x)`. This means that we have to split a lot to get good bounds
for the latter but not as much to get good bounds for the former. We
therefore [`mince`](@ref) the interval into several smaller pieces on
which we evaluate `T0(u0)` and then put the result together. While
this doesn't decrease the number of times we have to evaluate `T0(u0)`
(in fact it probably increases it), it does decrease the number of
times we have to evaluate `u0.v0`, which turns out to be the more
costly part of the procedure.
"""
function D0_bound(u0::BHKdVAnsatz{Arb}; atol = Arb(1.5e-2), verbose = false)
    verbose && @info "Computing bound of D0"

    ϵ = midpoint(Arb("1e-1"))

    # Bound the value on [ϵ, π]

    # Assert that the lemma holds
    @assert lemma_bhkdv_main_term_limit(u0)

    # Check that u0.v0 is positive at x = π
    Arblib.ispositive(u0.v0(Arb(π))) || error("u0.v0(π) not positive")

    # In practice u0.v0 is increasing on [0, π]. Try to prove that
    # this is the case on the interval [ϵ, b] with b < π.
    b = Arf(3.14)

    verbose && @info "Trying to prove that u0 is increasing on [ϵ, b]" b

    # Try to prove that minus the derivative of u0 is non-positive,
    # meaning the derivative is non-negative, on [ϵ, b].
    u0v0_is_increasing = ArbExtras.bounded_by(
        ArbExtras.derivative_function(x -> -u0.v0(x)),
        ϵ,
        b,
        Arf(0),
        degree = 0,
        depth_start = 2,
        depth = 5,
        threaded = true;
        verbose,
    )

    if verbose
        if u0v0_is_increasing
            @info "Succeeded with proving that it is increasing"
        else
            # This means we fall back to using enclosure_series
            @info "Failed with proving that it is increasing"
        end
    end

    f = T0(u0, skip_div_u0 = true)

    g(x; minces = 5) =
        let
            if iswide(x)
                # Split x into several smaller intervals and evaluate
                # f on each, then put them together.
                x_minced = HighestCuspedWave.mince(x, minces)
                res = f(x_minced[1])
                for i = 2:length(x_minced)
                    isfinite(res) || break
                    res = union(res, f(x_minced[i]))
                end
            else
                res = f(x)
            end

            isfinite(res) || return res

            u0v0x = if iswide(x)
                if u0v0_is_increasing && ϵ < x < b
                    xₗ, xᵤ = ArbExtras.enclosure_getinterval(x)
                    u0v0xₗ = u0.v0(xₗ)
                    u0v0xᵤ = u0.v0(xᵤ)
                    Arb((u0v0xₗ, u0v0xᵤ))
                else
                    ArbExtras.enclosure_series(u0.v0, x)
                end
            else
                u0.v0(x)
            end

            return res / u0v0x
        end

    # Bound it on [ϵ, π]
    D0 = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        ubound(Arb(π)),
        degree = -1,
        abs_value = true,
        point_value_max = g(Arb(1.4)), # Maximum is near x = 1.4
        depth_start = 8,
        threaded = true;
        atol,
        verbose,
    )

    verbose && @info "Maximum on [$(Float64(ϵ)), π]" D0

    # Show that it is bounded by m on [0, ϵ]

    h = T0(u0, Asymptotic(), ϵ = Arb(1.1ϵ))
    ϵ2 = Arf(1e-10)

    # Handle the interval [0, ϵ2] with one evalution
    if !(h(Arb((0, ϵ2))) <= D0)
        @error "Bound doesn't hold on $((0, Float64(ϵ2)))"
        return indeterminate(Arb)
    end

    verbose && @info "Bound satisfied on $((0, Float64(ϵ2)))"

    # Bound it on [ϵ2, ϵ]
    check_bound = ArbExtras.bounded_by(
        h,
        ϵ2,
        ϵ,
        lbound(D0),
        degree = -1,
        abs_value = true,
        log_bisection = true,
        threaded = true;
        verbose,
    )

    if !check_bound
        @error "Bound doesn't hold on $((Float64(ϵ2), Float64(ϵ)))"
        return indeterminate(Arb)
    end

    verbose && @info "Computed bound" D0

    return D0
end
