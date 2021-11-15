"""
    alpha0(u0::BHKdVAnsatz{Arb})

Compute an upper bound of `α₀` from the paper.

This method relies on the fact that `u0.v0(x)` gives a **lower** bound
for `u0(x)` for `0 <= x <= π` and any `α` in the interval `(-1, -1 +
u0.ϵ]`.

This means that `u0.w(x) / 2u0(x)` is **upper** bounded by `u0.w(x) /
2u0.v0(x)`. And to get an upper bound for `α₀` it is enough to bound
this value. The approach is therefore the same as in the version for
`BHAnsatz` which the only difference being the different weight and
that we only get an upper bound instead of an enclosure.

**PROVE:** That `u0(x)` indeed is lower bounded by `u0.v0(x)`. The
tails are the same so it trivially holds for them. The only thing
remaining is thus to prove that it also holds for the main term. That
is, we have to prove,
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) >=
    u0.v0.a0 * clausencmzeta(x, 2, 1)
```
for all `0 <= x <= π` and all `α` in our interval.
"""
function alpha0(u0::BHKdVAnsatz{Arb}; rtol = Arb("1e-2"), verbose = false)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb("0.5")
        @assert isequal(u0.w(x), x * log(u0.c + inv(x)))
    end

    # This uses evaluation of u0.v0. It should always overlap with the
    # value from u0. As an extra precaution we check that this seems
    # to be the case.
    let x = Arb("0.5")
        @assert Arblib.overlaps(u0(x), u0.v0(x))
    end
    let x = Arb("1e-2")
        @assert Arblib.overlaps(u0(x, Asymptotic()), u0.v0(x, Asymptotic()))
    end

    ϵ = Arb(1e-1)

    # Bound the value on [0, ϵ] using an asymptotic expansion

    # The first step is to create a function which can be evaluated at
    # around zero. We do this by computing the asymptotic expansion of
    # u0 at zero and then manually cancelling factors between the
    # weight and the expansion. Finally care has to be take to make
    # sure that the weight with the cancelled factors can be evaluated
    # around zero, this is done using monotonicity.
    f = let u0v0_expansion = u0.v0(ϵ, AsymptoticExpansion())

        # Divide all terms in the expansion by x * log(x)
        u0v0_expansion_div_xlogx = empty(u0v0_expansion)
        for ((i, m, k, l), value) in u0v0_expansion
            u0v0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
        end

        numerator(x) = begin
            # PROVE: That this is monotonically decreasing on [0, 1]
            # and -1 at x = 0
            h(x) = log(u0.c + inv(x)) / log(x)

            iszero(x) && return -one(x)

            if Arblib.contains_zero(x) && abs(x) < 1
                return Arb((h(abs_ubound(Arb, x)), -one(x)))
            end

            return h(x)
        end
        denominator(x) = 2eval_expansion(u0.v0, u0v0_expansion_div_xlogx, x)

        x -> numerator(x) / denominator(x)
    end

    # While we could treat the whole interval [0, ϵ] directly with
    # ArbExtras.maximum_enclosure it goes a low quicker if we first
    # evaluate on an interval containing zero directly and then do the
    # rest with ArbExtras.maximum_enclosure. The reason for this is
    # the behaviour of the log-bisection for intervals containing
    # zero, see ArbExtras.bisect_interval.

    ϵ2 = midpoint(Arb("1e-1000"))

    # Evaluate on [0, ϵ2]
    m11 = abs(f(Arb((0, ϵ2))))

    verbose && @show m11

    m12 = ArbExtras.maximum_enclosure(
        f,
        ϵ2,
        ubound(ϵ),
        degree = -1,
        abs_value = true,
        log_bisection = true,
        point_value_max = m11, # Maximum is in practice attained at zero
        threaded = true;
        rtol,
        verbose,
    )

    verbose && @show m12

    m1 = max(m11, m12)

    # Bound the value on [ϵ, π]
    g(x) = u0.w(x) / (2u0.v0(x))

    m2 = ArbExtras.maximum_enclosure(
        g,
        lbound(ϵ),
        ubound(Arb(π)),
        degree = 1, # Enough to make use of monotonicity
        abs_value = true,
        point_value_max = abs(g(Arb(π))), # Maximum is in practice attained here
        depth_start = 4,
        threaded = true;
        rtol,
        verbose,
    )

    verbose && @show m2

    return max(m1, m2)
end
