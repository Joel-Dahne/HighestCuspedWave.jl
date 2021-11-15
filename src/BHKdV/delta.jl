"""
    delta0(u0::BHKdVAnsatz{Arb})

Enclose the value of `δ₀` from the paper.

The interval `[0, π]` is first split into three intervals, `[0, ϵ1]`,
`[ϵ1, ϵ2]` and `[ϵ2, π]`.
1. On the first interval we evaluate `F0` directly using the
   asymptotic expansion.
2. On the second interval we use the asymptotic expansion together
   with [`ArbExtras.maximum_enclosure`](@ref).
3. On the third interval we uses straight forward ball enclosures
   together with [`ArbExtras.maximum_enclosure`](@ref).

**IMPROVE:** We could make the choice of the tolerances used dynamic
and not hard coded like now.
"""
function delta0(u0::BHKdVAnsatz{Arb}; rtol = Arb(1e-1), degree = 4, verbose = false)
    ϵ1 = midpoint(Arb("1e-10000"))
    ϵ2 = midpoint(Arb("1e-5"))
    @show ϵ1 ϵ2
    f = F0(u0, Asymptotic())

    # Bound the value on [0, ϵ1] with one evaluation
    m1 = f(Arb((0, ϵ1)))

    verbose && @show m1

    # Bound on [ϵ1, ϵ2]

    # Evaluate f at ϵ2 to get a lower bound for the maximum
    fϵ2 = f(Arb(ϵ2))

    m2 = ArbExtras.maximum_enclosure(
        F0(u0, Asymptotic(), ϵ = 2Arb(ϵ2)),
        ϵ1,
        ϵ2,
        abs_value = true,
        log_bisection = true,
        point_value_max = fϵ2,
        threaded = true,
        degree = -1,
        depth = 30,
        maxevals = 20000,
        atol = Arb("2e-4");
        rtol,
        verbose,
    )

    verbose && @show m2

    return max(m1, m2)

    # Bound the value on [ϵ2, π]
    m3 = ArbExtras.maximum_enclosure(
        F0(u0),
        ϵ2,
        Arblib.ubound(Arb(π)),
        abs_value = true,
        threaded = true;
        degree = -1,
        rtol,
        verbose,
    )

    verbose && @show m3

    return max(m1, m2, m3)
end
