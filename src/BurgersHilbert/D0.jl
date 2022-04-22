"""
    D0_bound(u0::BHAnsatz; atol = 1e-3, verbose = false)

Upper bound the value of `D₀` from the paper.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use an asymptotic expansion of `T0(u0)` whereas on ´[ϵ,
π]` we use a non-asymptotic version.

In practice the maximum is attained on `[ϵ, π]` and for that reason we
compute the maximum on this part first and then only prove that the
value on `[0, ϵ]` is bounded by this value.

For the non-asymptotic version we do a number of optimizations for
performance reasons.

The first optimization is to notice that even though `T0(u0)` doesn't
support evaluation on `ArbSeries`, and hence we can't use higher order
methods, `u0` does. We can thus get a much tighter enclosure of `u0`
by enclosing it with [`ArbExtras.enclosure_series`](@ref).

The second optimization comes from noticing that the enclosure for
`u0` that we get using the above strategy means that we get a much
tighter enclosure for `u0` than for `T0(u0)`. This means that we have
to split a lot to get good bounds for `T0(u0)` but not nearly as much
to get good bounds for `u0`. We therefore [`mince`](@ref) the interval
into several smaller pieces on which we evaluate `T0(u0)` and then put
the result together. While this doesn't decrease the number of times
we have to evaluate `T0(u0)` (in fact it probably increases it), it
does decrease the number of times we have to evaluate `u0`, which
turns out to be the more costly part of the procedure.
"""
function D0_bound(u0::BHAnsatz; atol = 1e-3, verbose = false)
    ϵ = Arf(1e-2)

    # Bound the value on [ϵ, π]

    # The evaluation of T0 doesn't support ArbSeries, but the
    # evaluation of u0 does. We can use ArbSeries to get a tighter
    # enclosure for u0.
    f = T0(u0, skip_div_u0 = true)

    g(x) = begin
        # Split x into several smaller intervals and evaluate f on
        # each, the put them together.
        xs = mince(x, 25)
        res = f(xs[1])
        for i = 2:length(xs)
            isfinite(res) || break
            Arblib.union!(res, res, f(xs[i]))
        end

        isfinite(res) || return res

        invu0 = inv(ArbExtras.enclosure_series(u0, x))

        return res * invu0#/ ArbExtras.enclosure_series(u0, x)
    end

    # Bound it on [ϵ, π]
    estimate = maximum(abs.(T0(u0).(range(Arb(ϵ), π, length = 10))))

    m = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        ubound(Arb(π)),
        degree = -1,
        abs_value = true,
        point_value_max = estimate,
        depth_start = 3,
        threaded = true;
        atol,
        verbose,
    )
    return m
    # Show that it is bounded by m on [0, ϵ]
    h = T0(u0, Asymptotic(), ϵ = Arb(1.1ϵ))
    ϵ2 = Arf(1e-100)

    # Handle the interval [0, ϵ2] with one evalution
    if !(h(Arb((0, ϵ2))) <= m)
        verbose && @error "Bounds not satisfied on [0, $(Float64(ϵ2))]"
        return Arblib.indeterminate!(zero(Arb))
    end

    verbose && @info "Bound satisfied on [0, $(Float64(ϵ2))]"

    # Bound it on [ϵ2, ϵ]
    bounded = ArbExtras.bounded_by(
        h,
        ϵ2,
        ϵ,
        lbound(m),
        degree = -1,
        log_bisection = true,
        threaded = true;
        verbose,
    )

    if !bounded
        verbose && @error "Could not prove bound on [$(Float64(ϵ2)), $(Float64(ϵ))]"
        if return_details
            return Arblib.indeterminate!(zero(Arb)), ϵ
        else
            return Arblib.indeterminate!(zero(Arb))
        end
    end

    return m
end
