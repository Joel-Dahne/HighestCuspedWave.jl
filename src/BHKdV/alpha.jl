"""
    alpha0(u0::BHKdVAnsatz{Arb})

Compute an upper bound of `α₀` from the paper.

This method relies on the fact that `u0.v0(x)` gives a **lower** bound
for `u0(x)` for `0 <= x <= π` and any `α` in the interval `(-1, -1 +
u0.ϵ]`.

This means that `u0.w(x) / 2u0(x)` is **upper** bounded by `u0.w(x) /
2u0.v0(x)`. and to get an upper bound for `α₀` it is enough to bound
this value. The approach is therefore the same as in the version for
`BHAnsatz` which the only difference being the different weight and
that we only get an upper bound instead of an enclosure.

**PROVE:** That `u0(x)` indeed is lower bounded by `u0.v0(x)`. The
tails are the same so it trivially holds for them. The only thing
remaining is thus to prove that it also holds for the main term. That
is, we have to prove,
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) >
    u0.v0.a0 * clausencmzeta(x, 2, 1)
```
for all `0 <= x <= π` and all `α` in our interval.
"""
function alpha0(
    u0::BHKdVAnsatz{Arb};
    rtol = 1e-5,
    degree = 1,
    maxevals = 2000,
    verbose = false,
)
    ϵ = Arb(1e-1)

    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * log(u0.c + inv(abs(x))))
    end

    # Bound the value on [0, ϵ] using an asymptotic expansion

    # The first step is to create a function which can be evaluated at
    # around zero. We do this by computing the asymptotic expansion of
    # u0 at zero and then manually cancelling factors between the
    # weight and the expansion. Finally care has to be take to make
    # sure that the weight with the cancelled factors can be evaluated
    # around zero, this is done using monotonicity.
    f = let u0v0_expansion = u0.v0(ϵ, AsymptoticExpansion())

        # Divide all terms in the expansion by abs(x) * log(abs(x))
        u0v0_expansion_div_xlogx = empty(u0v0_expansion)
        for ((i, m, k, l), value) in u0v0_expansion
            u0v0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
        end

        numerator(x) = begin
            # PROVE: That this is monotonically decreasing on [0, 1]
            # and -1 at x = 0
            h(x) = log(u0.c + inv(abs(x))) / log(abs(x))

            iszero(x) && return -one(x)

            if Arblib.contains_zero(x) && abs(x) < 1
                return Arb((h(Arblib.abs_ubound(Arb, x)), -one(x)))
            end

            return h(x)
        end
        denominator(x) = 2eval_expansion(u0.v0, u0v0_expansion_div_xlogx, x)

        x -> numerator(x) / denominator(x)
    end

    m1 = ArbExtras.maximum_enclosure(
        f,
        Arf(0),
        Arblib.ubound(ϵ),
        abs_value = true,
        point_value_max = abs(f(ϵ)), # Maximum is in practice attained here
        degree = -1,
        threaded = true,
        depth = 30,
        atol = 0.1; # This value is much smaller than m2 so we don't need good bounds
        maxevals,
        verbose,
    )

    # Bound the value on [ϵ, π]
    g(x) = u0.w(x) / (2u0.v0(x))

    m2 = ArbExtras.maximum_enclosure(
        g,
        Arblib.lbound(ϵ),
        Arblib.ubound(Arb(π)),
        abs_value = true,
        point_value_max = abs(g(Arb(π))), # Maximum is in practice attained here
        depth_start = 6,
        threaded = true;
        degree,
        maxevals,
        rtol,
        verbose,
    )

    @assert m1 < m2 # This should always hold in practice
    return max(m1, m2)
end
