"""
    delta0(u0::KdVZeroAnsatz; verbose)

Compute an expansion in `α` such that it gives an upper bound of
`abs(F0(u0)(x))` for all values of `α` in `u0.α`.

More precisely it returns an `ArbSeries` `p` of the form
```
p = 0 + 0 * α + p[2] * α^2
```
satisfying that `p(α)` gives an enclosure of \$δ_0\$ for every `α ∈
u0.α`.

For a given value of `x` [`F0`](@ref) gives us an expansion in `α`.
From this expansion we compute `p₂(x)` such that
```
0 + 0 * α + p₂(x) * α^2
```
gives an enclosure of `abs(F0(u0)(x))` for every `α ∈ u0.α`. We are then
interested in computing the maximum value of `p₂(x)` for `x ∈ [0, π]`.

To compute `p₂(x)` we note that `F0(u0)(x)` returns `q::ArbSeries`
such that `q(u0.α)` gives an enclosure for the specified `x`. We can
reduce this to a second order enclosure by rewriting `q` as
```
q = 0 + 0 * α + q[2] * α^2 + q[3] * α^3 + ...
  = (q[2] + q[3] * α + ...) * α^2
```
and we see that it is enough to compute an enclosure of `q[2] + q[3] *
α + ...`.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use the asymptotic version of `F0(u0)` whereas on ´[ϵ,
π]` we use the non-asymptotic version.
"""
function delta0(u0::KdVZeroAnsatz; threaded = true, verbose = false)
    ϵ = Arf(0.1)

    # In both the below methods we compute the expansion, divide by
    # α^2 and evaluate the resulting series at u0.α. This gives us an
    # enclosure of p₂(x). Note that the division by α checks so that
    # the first two coefficients are exactly zero.

    # Function for computing p₂(x) for x ∈ [0, ϵ]
    F0_asymptotic = F0(u0, Asymptotic())
    f(x) = (F0_asymptotic(x) << 2)(u0.α)

    # Function for computing p₂(x) for x ∈ [ϵ, π]
    F0_nonasymptotic = F0(u0, Ball())
    g(x) = (F0_nonasymptotic(x) << 2)(u0.α)


    # Bound the value on [0, ϵ]
    p2_asymptotic = ArbExtras.maximum_enclosure(
        f,
        zero(ϵ),
        ϵ,
        degree = -1,
        abs_value = true;
        threaded,
        verbose,
    )

    verbose && @info "Bound of p[2] from [0, ϵ]" p2_asymptotic

    # Bound the value on [ϵ, π]
    p2_nonasymptotic = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        ubound(Arb(π)),
        degree = -1,
        abs_value = true;
        threaded,
        verbose,
    )

    verbose && @info "Bound of p[2] on [ϵ, π]" p2_nonasymptotic

    res = max(p2_asymptotic, p2_nonasymptotic)

    return res
end
