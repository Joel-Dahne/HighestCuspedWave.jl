function alpha0(
    u0::BHAnsatz{Arb};
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
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    # Bound the value on [0, ϵ] using an asymptotic expansion

    # The first step is to create a function which can be evaluated at
    # around zero. We do this by computing the asymptotic expansion of
    # u0 at zero and then manually cancelling factors between the
    # weight and the expansion. Finally care has to be take to make
    # sure that the weight with the cancelled factors can be evaluated
    # around zero, this is done using monotonicity.
    f = let u0_expansion = u0(ϵ, AsymptoticExpansion())

        # Divide all terms in the expansion by abs(x) * log(abs(x))
        u0_expansion_div_xlogx = empty(u0_expansion)
        for ((i, m, k, l), value) in u0_expansion
            u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
        end

        numerator(x) = begin
            # PROVE: That this is monotonically decreasing on [0, 1]
            h(x) = sqrt(log((abs(x) + 1) / abs(x))) / log(abs(x))

            iszero(x) && return zero(x)

            if Arblib.contains_zero(x) && abs(x) < 1
                return union(zero(x), h(Arblib.abs_ubound(Arb, x)))
            end

            return h(x)
        end
        denominator(x) = 2eval_expansion(u0, u0_expansion_div_xlogx, x)

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
    g(x) = u0.w(x) / (2u0(x))

    m2 = ArbExtras.maximum_enclosure(
        g,
        Arblib.lbound(ϵ),
        Arblib.ubound(Arb(π)),
        abs_value = true,
        point_value_max = abs(g(Arb(π))), # Maximum is in practice attained here
        threaded = true;
        degree,
        maxevals,
        rtol,
        verbose,
    )

    @assert m1 < m2 # This should always hold in practice
    return max(m1, m2)
end
