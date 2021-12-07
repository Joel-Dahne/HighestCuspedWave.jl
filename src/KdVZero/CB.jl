"""
    CB(u0::KdVZeroAnsatz; verbose = false)

Compute an expansion in `α` such that it gives an upper bound of
`T0(u0)(x)` for all values of `α` in `u0.α`.

More precisely it returns an `ArbSeries` `p` of the form
```
p = 1 + p[1] * α
```
satisfying that `p(α)` gives an enclosure of \$C_B\$ for every `α ∈
u0.α`.

For a given value of `x` [`T0`](@ref) gives us an expansion in `α`.
From this expansion we compute `p₁(x)` such that
```
1 + p₁(x) * α
```
gives an enclosure of `T0(u0)(x)` for every `α ∈ u0.α`. We are then
interested in computing the minimum value of `p(x)` for `x ∈ [0, π]`,
we take the minimum value since `α` is negative.

To compute `p₁(x)` we note that `T0(u0)(x)` returns `q::ArbSeries`
such that `q(u0.α)` gives an enclosure for the specified `x`. We can
reduce this to a first order enclosure that rewriting `q` as
```
q = 1 + q[1] * α + q[2] * α^2 + ...
  = 1 + (q[1] + q[2] * α + ...) * α
```
and we see that it is enough to compute an enclosure of `(q[1] + q[2]
* α + ...)`.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use an asymptotic expansion of `T0(u0)` whereas on ´[ϵ,
π]` we use a non-asymptotic version.

In practice the maximum is attained at `x = 0` and for that reason we
compute the maximum on the interval `[0, ϵ]` first and then prove that
the value on `[ϵ, π]` is bounded by this value.
"""
function CB(u0::KdVZeroAnsatz; verbose = false)
    ϵ = Arf(0.1)

    # In both the below methods we compute the expansion, subtract 1
    # and divide by α, evaluate the resulting series at u0.α. This
    # gives us an enclosure of p₁(x).

    # Function for computing p₁(x) for x ∈ [0, ϵ]
    T0_asymptotic = T0(u0, Asymptotic())
    f(x) = ((T0_asymptotic(x) - 1) << 1)(u0.α)

    # Function for computing p₁(x) for x ∈ [ϵ, π]
    T0_nonasymptotic = T0(u0, Ball())
    g(x) = ((T0_nonasymptotic(x) - 1) << 1)(u0.α)

    # Compute an enclosure on [0, ϵ]
    # TODO: Implement this
    #p1 = ArbExtras.minimum_enclosure(f, Arf(0), ϵ, degree = -1; verbose)

    #verbose && @show p1

    # Check that the value on [ϵ, π] is bounded by that on [0, ϵ]
    # TODO: Implement this
    #ArbExtras..bounded_by(g, ϵ, ubound(Arb(π)), lbound(p1), degree = -1; verbose) ||
    #    error("could not show bound on [ϵ, π]")

    # FIXME: For now we only compute it on a discrete set of points
    xs = range(Arb(0), π, length = 100)[2:end-1]
    ys = g.(xs)

    p1 = minimum(ys)

    verbose && @show p1

    return ArbSeries((1, p1))
end
