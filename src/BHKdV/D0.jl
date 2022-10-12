"""
    D0_bound(u0::BHAnsatz; atol, verbose = false)

Compute an upper bound of `D₀` from the paper.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use an asymptotic expansion of `T0(u0)` whereas on ´[ϵ,
π]` we use a non-asymptotic version.

In practice the maximum is attained on `[ϵ, π]` and for that reason we
compute the maximum on this part first and then only prove that the
value on `[0, ϵ]` is bounded by this value.

This method uses the fact that `u0.v0(x)` gives a **lower** bound for
`u0(x)` for `0 <= x <= π` and any `α` in the interval `(-1, -1 +
u0.ϵ]`. This means that `inv(u0(x))` is **upper** bounded by
`inv(u0.v0(x))`.
- **PROVE:** That `u0(x)` indeed is lower bounded by `u0.v0(x)`. See
  [`n0_bound`](@ref).

For the non-asymptotic version we do a number of optimizations for
performance reasons.

The first optimization is to notice that even though `T0(u0)` doesn't
support evaluation on `ArbSeries`, and hence we can't use higher order
methods, `u0.v0` does. We can thus get a much tighter enclosure of
`u0.v0` by enclosing it with [`ArbExtras.enclosure_series`](@ref).

The second optimization comes from noticing that the enclosure for
`u0.v0` that we get using the above strategy means that we get a much
tighter enclosure for `u0` than for `T0(u0)`. This means that we have
to split a lot to get good bounds for `T0(u0)` but not nearly as much
to get good bounds for `u0`. We therefore [`mince`](@ref) the interval
into several smaller pieces on which we evaluate `T0(u0)` and then put
the result together. While this doesn't decrease the number of times
we have to evaluate `T0(u0)` (in fact it probably increases it), it
does decrease the number of times we have to evaluate `u0`, which
turns out to be the more costly part of the procedure.

One problem we have to deal with is that `T0(u0)` currently doesn't
fully allow evaluation on `x` values strictly larger than `π`. To
handle this we use [`ArbExtras.extrema_enclosure`](@ref) on the
interval `[ϵ, lbound(Arb(π))]`, which is strictly less than `π`, and
then evaluate on the remaining endpoint `[lbound(Arb(π)), π]`
separately.
"""
function D0_bound(u0::BHKdVAnsatz{Arb}; atol = Arb(1.5e-2), verbose = false)
    ϵ = midpoint(Arb("1e-1"))

    # Bound the value on [ϵ, π]

    # The evaluation of T0 doesn't support ArbSeries, but the
    # evaluation of u0 does. We can use ArbSeries to get a tighter
    # enclosure for u0.
    f = T0(u0, skip_div_u0 = true)

    minces = 5

    g(x) =
        let
            if iswide(x)
                # Split x into several smaller intervals and evaluate f on
                # each, then put them together.
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

            return res / ArbExtras.enclosure_series(u0.v0, x)
        end

    # Bound it on [ϵ, π]
    a = ϵ
    b = ubound(Arb(π))

    xs = range(Arb(a), Arb(b), length = 20)
    ys = similar(xs)
    Threads.@threads for i in eachindex(xs)
        ys[i] = g(xs[i])
    end

    m_approx = maximum(abs.(ys))

    verbose && @info "Approximate maximum" m_approx

    m = ArbExtras.maximum_enclosure(
        g,
        a,
        b,
        degree = -1,
        abs_value = true,
        point_value_max = m_approx,
        depth_start = 8,
        threaded = true;
        atol,
        verbose,
    )

    verbose && @info "Maximum on [$(Float64(ϵ)), π]" m

    # Show that it is bounded by m on [0, ϵ]

    h = T0(u0, Asymptotic(), ϵ = Arb(1.1ϵ))
    ϵ2 = Arf(1e-10)

    # Handle the interval [0, ϵ2] with one evalution
    if !(h(Arb((0, ϵ2))) <= m)
        @error "Bound doesn't hold on $((0, Float64(ϵ2)))"
        return indeterminate(Arb)
    end

    verbose && @info "Bound satisfied on $((0, Float64(ϵ2)))"

    # Bound it on [ϵ2, ϵ]
    check_bound = ArbExtras.bounded_by(
        h,
        ϵ2,
        ϵ,
        lbound(m),
        degree = -1,
        log_bisection = true,
        threaded = true;
        verbose,
    )

    if !check_bound
        @error "Bound doesn't hold on $((Float64(ϵ2), Float64(ϵ)))"
        return indeterminate(Arb)
    end

    return m
end
