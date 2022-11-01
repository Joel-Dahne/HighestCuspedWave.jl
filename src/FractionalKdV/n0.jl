"""
    n0_bound(u0::FractionalKdVAnsatz; M::Integer, rtol, threaded, verbose)

Enclose the value of `n₀`. This is the supremum of
```
N(x) = u0.w(x) / 2u0(x)
```
for `0 < x < π`.

It uses an asymptotic expansion with `M` terms close to zero and ball
arithmetic on the remaining.

In practice the maximum is attained at `x = π` so the function is
evaluated there first. The interval is split into three subintervals
``[0, ϵ1]``, ``[ϵ1, ϵ2]`` and ``[ϵ2, π]``. On the first two
subintervals we use the asymptotic expansion and only prove that the
value is bounded by that at `x = π`, on the third subinterval we use
ball arithmetic compute an enclosure of the maximum.

The values for `ϵ1` and `ϵ2` are determined automatically. The value
of `ϵ1` is taken such that the interval ``[0, ϵ1]`` can be handled in
one evaluation. The value for `ϵ2` is taken such that the value at `x
= ϵ2` is less than the value at `x = π`, the interval ``[ϵ1, ϵ2]`` is
the handled with bisection.

On the interval `[ϵ, π]` we bound it using ball arithmetic with
[`ArbExtras.maximum_enclosure`](@ref). Notice that we do not have to
prove that the maximum is attained at `x = π`, it's only used to make
the procedure more efficient.

If `threaded = true` it enables threading when calling
[`ArbExtras.maximum_enclosure`](@ref). If `verbose = true` print more
information about the process.
"""
function n0_bound(
    u0::FractionalKdVAnsatz{Arb};
    M::Integer = 3,
    rtol = 1e-5,
    threaded = true,
    verbose = false,
)
    # For non-asymptotic evaluation
    f = x -> u0.w(x) / (2u0(x))
    # For asymptotic evaluation
    inv_u0 = inv_u0_normalised(u0; M, ϵ = Arb(π))
    g = x -> u0.wmulpow(x, u0.α) * inv_u0(x) / 2

    # The maximum is in practice attained at x = π
    m1 = f(Arb(π))

    verbose && @info "f(π) = $m1"

    # Find ϵ2 such that g(ϵ2) < m1
    ϵ2 = Arb(3)
    while !(g(ϵ2) < m1)
        ϵ2 *= 0.8
    end

    verbose && @info "Determined ϵ2" ϵ2

    # Find ϵ1 such that g(Arb(0, ϵ1)) < m1
    ϵ1 = ϵ2
    while !(g(Arb((0, ϵ1))) < m1)
        ϵ1 = ϵ1 < 0.5 ? ϵ1^2 : ϵ1 / 2
        if (!u0.use_bhkdv && ϵ1 < 1e-100) || (u0.use_bhkdv && ϵ1 < Arb("1e-10000"))
            verbose && @error "Could not prove bound on [0, ϵ1]" g(Arb((0, ϵ1)))
            return indeterminate(Arb)
        end
    end

    verbose && @info "Determined ϵ1" ϵ1

    # Prove the bound on [0, ϵ]
    bounded = ArbExtras.bounded_by(
        g,
        lbound(ϵ1),
        ubound(ϵ2),
        lbound(m1),
        abs_value = true,
        log_bisection = true,
        maxevals = ifelse(u0.use_bhkdv, 100000, 1000);
        threaded,
        verbose,
    )

    if !bounded
        verbose && @error "Could not prove bound on [ϵ1, ϵ2]"
        return indeterminate(Arb)
    end

    verbose && @info "Proved bound on [ϵ1, ϵ2]"

    # Bound the value on [ϵ, π]
    m = ArbExtras.maximum_enclosure(
        f,
        lbound(ϵ2),
        ubound(Arb(π)),
        abs_value = true,
        point_value_max = m1, # m1 is a lower bound of the maximum
        atol = 4radius(m1), # We cannot expect to do better than m1
        depth_start = 2;
        rtol,
        threaded,
        verbose,
    )

    return m
end
