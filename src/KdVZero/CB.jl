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
interested in computing the minimum value of `p₁(x)` for `x ∈ [0, π]`,
we take the minimum value since `α` is negative.

To compute `p₁(x)` we note that `T0(u0)(x)` returns `q::ArbSeries`
such that `q(u0.α)` gives an enclosure for the specified `x`. We can
reduce this to a first order enclosure by rewriting `q` as
```
q = 1 + q[1] * α + q[2] * α^2 + ...
  = 1 + (q[1] + q[2] * α + ...) * α
```
and we see that it is enough to compute an enclosure of `q[1] + q[2]
* α + ...`.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use the asymptotic version of `F0(u0)` whereas on ´[ϵ,
π]` we use the non-asymptotic version.
"""
function CB(
    u0::KdVZeroAnsatz;
    ϵ::Arf = Arf(0.1),
    rtol = Arb(1e-3), # TODO: Possibly tune this,
    maxevals = 1000, # TODO: This we can probably remove later
    threaded = true,
    verbose = false,
)
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
    p1_asymptotic = ArbExtras.minimum_enclosure(
        f,
        zero(ϵ),
        ϵ,
        degree = -1;
        rtol,
        maxevals,
        threaded,
        verbose,
    )

    verbose && @info "Bound of p[1] on [0, ϵ]" p1_asymptotic

    # Compute an enclosure on [0, ϵ]
    p1_nonasymptotic = ArbExtras.minimum_enclosure(
        g,
        ϵ,
        ubound(Arb(π)),
        degree = -1;
        rtol,
        maxevals,
        threaded,
        verbose,
    )

    verbose && @info "Bound of p[1] on [ϵ, π]" p1_nonasymptotic

    p1 = max(p1_asymptotic, p1_nonasymptotic)

    return ArbSeries((1, p1))
end

"""
    CB_estimate(u0::KdVZeroAnsatz)

Compute an estimate of `CB(u0)`. It returns an `ArbSeries` `p` of the
form
```
p = 1 + p[1] * α
```
such that `p(α)` gives an estimate of \$C_B\$ for the given `α ∈
u0.α`.
"""
function CB_estimate(
    u0::KdVZeroAnsatz;
    n::Integer = 100,
    return_values = false,
    include_zero = false,
    threaded = true,
)
    xs = collect(range(zero(Arb), π, length = n + 1)[2:end-1])

    ys = similar(xs)

    g = T0(u0, Ball())
    f(x) = ((g(x) - 1) << 1)(u0.α)
    if threaded
        Threads.@threads for i in eachindex(xs)
            ys[i] = f(xs[i])
        end
    else
        for i in eachindex(xs)
            ys[i] = f(xs[i])
        end
    end

    if include_zero
        # FIXME: Not properly implemented yet
        pushfirst!(xs, zero(Arb))
        pushfirst!(ys, ((T0(u0, Asymptotic())(zero(Arb)) - 1) << 1)(u0.α))
    end

    p1 = minimum(ys)

    res = ArbSeries((1, p1))

    return_values && return res, xs, ys
    return res
end
