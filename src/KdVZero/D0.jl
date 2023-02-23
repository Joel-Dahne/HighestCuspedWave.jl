"""
    D0_bound(u0::KdVZeroAnsatz; verbose = false)

Compute a Taylor model in `α` of `D₀`, which is the supremum of
```
T0(u0)(x)
```
for `0 < x < π. The constant term in the expansion is `1`.

For a given value of `x` [`T0`](@ref) gives us a Taylor model in `α`
of the form
```
1 + ΔD(x) * α
```
, where `ΔD(x)` is the remainder term of the Taylor model, that gives
an enclosure of `T0(u0)(x)` for every `α ∈ u0.α`. We are then
interested in computing the minimum value of `ΔD(x)` for `x ∈ [0, π]`,
we take the minimum value since `α` is negative.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use the asymptotic version of `T0(u0)` whereas on ´[ϵ,
π]` we use the non-asymptotic version.

In practice the the minimum is attained at `x = 0`, we therefore
compute the minimum on `[0, ϵ]` first and then only prove that this is
a lower bound of the value on `[ϵ, π]`.

We take `ϵ = 1`, it might be possible to get better results by tuning
this choice but in practice it is good enough.
"""
function D0_bound(u0::KdVZeroAnsatz; rtol = Arb(1e-2), threaded = true, verbose = false)
    verbose && @info "Computing enclosure of ΔD"

    ϵ = Arf(1)

    # Compute the value at x = 0, this is the minimum in practice
    ΔD_zero = T0(u0, Asymptotic())(Arb(0)).p[1]

    verbose && @info "ΔD(0) = $ΔD_zero"

    # Function for computing ΔD(x) for x ∈ [0, ϵ]
    T0_asymptotic = T0(u0, Asymptotic(), ϵ = Arb(ϵ + 0.1))
    f(x) = T0_asymptotic(x).p[1]

    # Compute an enclosure on [0, ϵ]
    ΔD = ArbExtras.minimum_enclosure(
        f,
        zero(ϵ),
        ϵ,
        degree = -1,
        point_value_min = ΔD_zero; # Minimum in practice
        rtol,
        threaded,
        verbose,
    )

    verbose && @info "Enclosure of ΔD on [0, ϵ]" ΔD

    # To prove that the value on [ϵ, π] is lower bounded by Δ we use
    # ArbExtras.bounded_by. This only works for upper bounds so we
    # instead prove that -Δ(x) <= -Δ

    # Function for computing -Δ(x) for x ∈ [ϵ, π]
    T0_nonasymptotic = T0(u0, Ball())
    g(x) = -T0_nonasymptotic(x).p[1]

    # Prove that on the interval [ϵ, π] it is bounded by ΔD
    check_bound = ArbExtras.bounded_by(
        g,
        ϵ,
        ubound(Arb(π)),
        lbound(-ΔD),
        degree = -1,
        maxevals = 50000;
        threaded,
        verbose,
    )

    if !check_bound
        @error "Could not prove ΔD(x) on [ϵ, π] bounded by $ΔD"
        return TaylorModel(ArbSeries((1, NaN), degree = 1), u0.α, u0.α0)
    end

    verbose && @info "Proved ΔD(x) on [ϵ, π] bounded by $ΔD"

    verbose && @info "Computed enclosure" ΔD

    return TaylorModel(ArbSeries((1, ΔD), degree = 1), u0.α, u0.α0)
end
