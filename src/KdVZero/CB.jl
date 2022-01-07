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

For a given value of `x` [`T0`](@ref) gives us an expansion in `α` of
the form
```
1 + p₁(x) * α
```
that gives an enclosure of `T0(u0)(x)` for every `α ∈ u0.α`. We are
then interested in computing the minimum value of `p₁(x)` for `x ∈ [0,
π]`, we take the minimum value since `α` is negative.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use the asymptotic version of `T0(u0)` whereas on ´[ϵ,
π]` we use the non-asymptotic version.

In practice the the minimum is attained at `x = 0`, we therefore
compute the minimum on `[0, ϵ]` first and then only prove that this is
a lower bound of the value on `[ϵ, π]`.

The value for `ϵ` is chosen dynamically by first computing the value
at `x = 0` and then trying to pick `ϵ` as large but so that the
enclosure from the asymptotic version of `T0` is larger than the value
at `x = 0`
"""
function CB(u0::KdVZeroAnsatz; rtol = Arb(1e-2), threaded = true, verbose = false)
    # Compute the value at x = 0, this is the minimum in practice
    p1_zero = T0(u0, Asymptotic())(Arb(0))[1]

    # Determine a good choice of ϵ
    ϵ = let ϵ = Arb(π), f = T0(u0, Asymptotic(), ϵ = ubound(Arb, ϵ)), g = T0(u0)
        y = f(ϵ)[1]
        z = g(ϵ)[1]

        # Reduce ϵ until the value we get either has a lower bound
        # greater than the lower bound of the value at x = 0 or the
        # asymptotic version gives a better enclosure than the
        # non-asymptotic version.
        while !(isfinite(y) && lbound(y) > lbound(p1_zero)) && radius(y) > radius(z)
            if ϵ > 1
                ϵ -= 0.05
            else
                ϵ /= 1.2
            end
            y = f(ϵ)[1]
            z = g(ϵ)[1]
            ϵ > 0.1 || error("could not determine working ϵ, last tried value was $ϵ")
        end
        ubound(ϵ)
    end
    verbose && @info "Determined ϵ" ϵ

    # Function for computing p₁(x) for x ∈ [0, ϵ]
    T0_asymptotic = T0(u0, Asymptotic(), ϵ = Arb(ϵ + 0.1))
    f(x) = T0_asymptotic(x)[1]

    # Compute an enclosure on [0, ϵ]
    p1 = ArbExtras.minimum_enclosure(
        f,
        zero(ϵ),
        ϵ,
        degree = -1,
        point_value_min = p1_zero; # Minimum in practice
        rtol,
        threaded,
        verbose,
    )

    verbose && @info "Bound of p[1] on [0, ϵ]" p1

    # To prove that the value on [ϵ, π] is lower bounded by p1 we use
    # ArbExtras.bounded_by. This only works for upper bounds so we
    # instead prove that -p₁ < -p1

    # Function for computing -p₁(x) for x ∈ [ϵ, π]
    T0_nonasymptotic = T0(u0, Ball())
    g(x) = -T0_nonasymptotic(x)[1]

    # Prove that on the interval [ϵ, π] it is bounded by p1
    check_bound = ArbExtras.bounded_by(
        g,
        ϵ,
        ubound(Arb(π)),
        lbound(-p1),
        degree = -1,
        maxevals = 50000;
        threaded,
        verbose,
    )

    if !check_bound
        @error "Could not prove bound on [ϵ, π]"
        return ArbSeries((1, NaN))
    end

    verbose && @info "Proved bound on [ϵ, π]"

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
    xs = collect(range(zero(Arb), π, length = n + 1)[2:end])

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
        pushfirst!(xs, zero(Arb))
        pushfirst!(ys, ((T0(u0, Asymptotic())(zero(Arb)) - 1) << 1)(u0.α))
    end

    p1 = minimum(ys)

    res = ArbSeries((1, p1))

    return_values && return res, xs, ys
    return res
end
