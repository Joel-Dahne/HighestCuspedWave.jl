"""
    delta0_bound(u0::BHKdVAnsatz{Arb}; atol = 0, rtol = 1e-1, ubound_tol = 5e-4, verbose = false)

Compute an upper bound of the value of `δ₀` from the paper. This is
the supremum of
```
abs(F0(u0))
```
for `0 < x < π`.

The interval `[0, π]` is split into three intervals, `[0, ϵ1]`, `[ϵ1,
ϵ]` and `[ϵ, π]`.
1. On the first interval we evaluate `F0` directly using the
   asymptotic expansion.
2. On the second interval we use the asymptotic expansion together
   with [`ArbExtras.maximum_enclosure`](@ref).
3. On the third interval we use [`F0_bound`](@ref) together with
   [`ArbExtras.maximum_enclosure`](@ref).

The value for `ϵ` is taken as large as possible but so that the
asymptotic evaluation still satisfies the required tolerance.

In practice the maximum is attained around `1e-2` so it computes an
approximation of the maximum on `[1e-10, ϵ]`. It then looks for `ϵ1`
such that the enclosure for the interval `[0, ϵ1]` is less than this
approximate maximum.

The tolerance is determined by `atol`, `rtol` and `ubound_tol`. For
`ubound_tol` an enclosure is determined to satisfy the tolerance if it
is smaller than `ubound_tol`. The default value for `ubound_tol` is
set so that the computed bound should satisfy the required inequality
with the expected values for `n₀` and `D₀`.
"""
function delta0_bound(
    u0::BHKdVAnsatz{Arb};
    atol = Arb(0),
    rtol = Arb(1e-1),
    ubound_tol = Arb(5e-4),
    verbose = false,
)
    verbose && @info "Computing bound of δ₀"

    # Determine a good choice of ϵ. Take it as large as possible so
    # that the asymptotic version still satisfies the specified
    # tolerance or it is better than the non-asymptotic version.
    ϵ = let ϵ = Arb(0.6), f = F0(u0, Asymptotic(), ϵ = 1.01ϵ, M = 10), g = F0_bound(u0)
        y = f(ϵ)
        z = g(ϵ)

        # Reduce ϵ until the value we get either satisfies the
        # required tolerance or is better than the non-asymptotic
        # version.
        while !ArbExtras.check_tolerance(y; atol, rtol) &&
                  !(y < ubound_tol) &&
                  radius(y) > radius(z)
            ϵ *= 0.8
            y = f(ϵ)
            z = g(ϵ)
            ϵ > 1e-2 || error("could not determine working ϵ, last tried value was $ϵ")
        end
        midpoint(ϵ)
    end

    verbose && @info "Determined ϵ" ϵ

    # We use two versions of the function for asymptotic evaluation,
    # with M = 3 and M = 10. We use the former for x < δ and the
    # latter otherwise
    δ = Arb(1e-2)
    f1 = F0(u0, Asymptotic(), ϵ = 1.01Arb(ϵ), M = 3)
    f2 = F0(u0, Asymptotic(), ϵ = 1.01Arb(ϵ), M = 10)
    f = x -> x < δ ? f1(x) : f2(x)

    verbose &&
        isfinite(ubound_tol) &&
        !(f(Arb(ϵ)) < ubound_tol) &&
        @warn "ubound_tol not satisfied at ϵ" f(Arb(ϵ)) ubound_tol

    # Compute an approximation of the maximum on [1e-10, ϵ]
    xs_asym = exp.(range(log(Arb(1e-10)), log(Arb(ϵ)), length = 24))
    ys_asym = similar(xs_asym)
    Threads.@threads for i in eachindex(xs_asym)
        ys_asym[i] = f(xs_asym[i])
    end
    m_asym_approx = maximum(abs.(ys_asym))

    verbose && @info "Approximate maximum on [1e-10, ϵ]" m_asym_approx

    verbose &&
        isfinite(ubound_tol) &&
        !(m_asym_approx < ubound_tol) &&
        @warn "ubound_tol not satisfied for approximate maximum" ubound_tol

    # Find ϵ1 such that the interval [0, ϵ1] can be handled in one
    # evaluation
    N = -100 # Use the interval [0, 10^-N]
    m1 = f(Arb((0, Arb(10)^N)))
    while !(m1 < ubound(Arb, m_asym_approx)) && !(m1 < ubound_tol)
        N -= 100
        m1 = f(Arb((0, Arb(10)^N)))
        N > -20000 || error("could not determine working N, last tried value was $N")
    end
    ϵ1 = lbound(Arb(10)^N)

    verbose && @info "Determined ϵ1" ϵ1

    # Bound on [0, ϵ]
    m2 = ArbExtras.maximum_enclosure(
        f,
        ϵ1,
        ϵ,
        abs_value = true,
        log_bisection = true,
        point_value_max = m_asym_approx,
        threaded = true,
        degree = -1,
        depth = 35,
        maxevals = 200000;
        atol,
        rtol,
        ubound_tol,
        verbose,
    )

    verbose && @info "Maximum on [ϵ1, ϵ]" m2

    # Bound the value on [ϵ, π]
    g = F0_bound(u0)

    m3 = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        ubound(Arb(π)),
        degree = 4,
        abs_value = true,
        threaded = true,
        depth_start = 4;
        atol,
        rtol,
        ubound_tol,
        verbose,
    )

    verbose && @info "Maximum on [ϵ, π]" m3

    δ₀ = max(m1, m2, m3)

    verbose && @info "Computed bound" δ₀

    return δ₀
end
