"""
    delta0(u0::BHKdVAnsatz{Arb})

Enclose the value of `δ₀` from the paper.

The interval `[0, π]` is first split into three intervals, `[0, ϵ1]`,
`[ϵ1, ϵ2]` and `[ϵ2, π]`.
1. On the first interval we evaluate `F0` directly using the
   asymptotic expansion.
2. On the second interval we use the asymptotic expansion together
   with [`ArbExtras.maximum_enclosure`](@ref).
3. On the third interval we use [`F0_upper_bound`](@ref) together with
   [`ArbExtras.maximum_enclosure`](@ref).

**IMPROVE:** We could make the choice of the tolerances used dynamic
and not hard coded like now.
"""
function delta0(u0::BHKdVAnsatz{Arb}; rtol = Arb("1e-1"), verbose = false)
    ϵ1 = midpoint(Arb("1e-10000"))
    ϵ2 = midpoint(Arb("1e-5"))

    f = F0(u0, Asymptotic())

    # Bound the value on [0, ϵ1] with one evaluation
    m1 = f(Arb((0, ϵ1)))

    verbose && @show m1

    # Bound on [ϵ1, ϵ2]
    m2 = ArbExtras.maximum_enclosure(
        F0(u0, Asymptotic(), ϵ = 2Arb(ϵ2)),
        ϵ1,
        ϵ2,
        abs_value = true,
        log_bisection = true,
        point_value_max = f(Arb(ϵ2)),
        threaded = true,
        degree = -1,
        depth = 30,
        maxevals = 20000,
        atol = Arb("2e-4");
        rtol,
        verbose,
    )

    verbose && @show m2

    # Bound the value on [ϵ2, π]
    g = F0_upper_bound(u0)

    # Compute an approximation of the maximum. Take log-spaced points
    # on [ϵ2, 1e-1] and linearly space on [1e-1, π].
    xs1 = exp.(range(log(Arb(ϵ2)), log(Arb("1e-1")), length = 20))
    xs2 = range(Arb("1e-1"), Arb(π), length = 20)
    xs = [xs1; xs2]
    ys = similar(xs)
    Threads.@threads for i in eachindex(xs)
        ys[i] = g(xs[i])
    end

    m3_approx = maximum(abs.(ys))

    verbose && @show m3_approx

    m3 = ArbExtras.maximum_enclosure(
        g,
        ϵ2,
        ubound(Arb(π)),
        degree = 4,
        abs_value = true,
        point_value_max = m3_approx,
        threaded = true,
        depth_start = 4,
        maxevals = 100000,
        depth = 30;
        rtol,
        verbose,
    )

    verbose && @show m3

    return max(m1, m2, m3)
end
