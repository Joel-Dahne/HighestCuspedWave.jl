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
evaluated there first. Then starting at `ϵ = π/2` we halve `ϵ` until
the bound on `x = ϵ` is smaller than the value at `π`. We the prove
that the value on ``[0, ϵ]`` is bounded by the value at `x = π`.

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
    g = x -> abspow(x, u0.p + u0.α) * inv_u0(x) / 2

    # The maximum is in practice attained at x = π
    m1 = f(Arb(π))

    verbose && @info "f(π) = $m1"

    # Find ϵ such that g(ϵ) < m1
    ϵ = Arb(3)
    while !(g(ϵ) < m1)
        ϵ *= 0.8
    end

    verbose && @info "Determined ϵ" ϵ

    # Prove the bound on [0, ϵ]
    bounded = ArbExtras.bounded_by(
        g,
        Arf(0),
        ubound(ϵ),
        lbound(m1),
        abs_value = true;
        threaded,
        verbose,
    )

    if !bounded
        verbose && @error "Could not prove bound on [0, ϵ]"
        return indeterminate(Arb)
    end

    verbose && @info "Proved bound on [0, ϵ]"

    # Bound the value on [ϵ, π]
    m = ArbExtras.maximum_enclosure(
        f,
        lbound(ϵ),
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
