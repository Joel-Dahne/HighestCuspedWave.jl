"""
    CB(u0::BHAnsatz; atol = 1e-3, verbose = false)

Upper bound the value of `C_B` from the paper.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use an asymptotic expansion of `T0(u0)` whereas on ´[ϵ,
π]` we use a non-asymptotic version.

For the non-asymptotic version we do a number of optimizations for
performance reasons.

The first optimization is to notice that even though `T0(u0)` doesn't
support evaluation on `ArbSeries`, and hence we can't use higher order
methods, `u0` does. We can thus get a much tighter enclosure of `u0`
but bounding it with [`ArbExtras.extrema_series`](@ref). We can also
notice that it's enough to do this using `degree = 1` since most of
the time it will pick up that `u0` is monotone on the interval and
only evaluate it on the endpoints.

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

One problem we have to deal with is that `T0(u0)` currently doesn't
fully allow evaluation on `x` values strictly larger than `π`. To
handle this we use [`ArbExtras.extrema_enclosure`](@ref) on the
interval `[ϵ, lbound(Arb(π))]`, which is strictly less than `π`, and
then evaluate on the remaining endpoint `[lbound(Arb(π)), π]`
separately.
"""
function CB(u0::BHAnsatz; atol = 1e-3, verbose = false)
    ϵ = Arb(1e-1)

    # TODO: Bound the value on [0, ϵ]
    @warn "evaluation on [0, ϵ] not implemented"
    m1 = zero(Arb)

    # Bound the value on [ϵ, π]

    # The evaluation of T0 doesn't support ArbSeries, but the
    # evaluation of u0 does. We can use ArbSeries to get a tighter
    # enclosure for u0.
    f = T0(u0, skip_div_u0 = true)

    g(x) = begin
        # Split x into several smaller intervals and evaluate f on
        # each, the put them together.
        xs = HighestCuspedWave.mince(x, 400)
        res = f(xs[1])
        for i = 2:length(xs)
            isfinite(res) || break
            res = union(res, f(xs[i]))
        end

        isfinite(res) || return res

        # Compute a tighter enclosure of u0 by using evaluation with
        # ArbSeries. In practice this will pick up that u0 is monotone
        # on the interval and therefore very efficient.
        invu0 = inv(
            Arb(
                ArbExtras.extrema_series(u0, Arblib.getinterval(Arf, x)..., degree = 1)[1:2],
            ),
        )

        return res * invu0
    end

    # Bound it on [ϵ, lbound(Arb(π))]
    a = Arblib.lbound(ϵ)
    b = Arblib.lbound(Arb(π))
    estimate = maximum(abs.(T0(u0).(range(Arb(a), Arb(b), length = 10))))

    m21 = ArbExtras.maximum_enclosure(
        g,
        a,
        b,
        degree = -1,
        abs_value = true,
        point_value_max = estimate,
        threaded = true;
        atol,
        verbose,
    )

    # Bound it on [lbound(Arb(π)), π]
    m22 = T0(u0)(union(Arb(b), Arb(π)))

    m2 = max(m21, m22)

    return max(m1, m2)
end
