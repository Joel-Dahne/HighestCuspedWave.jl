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
We truncate this expansion to degree `2` using
[`trunace_with_remainder`](@ref), giving us a polynomial of the form
```
0 + 0 * α + p₂(x) * α^2
```
that gives an enclosure of `abs(F0(u0)(x))` for every `α ∈ u0.α`. We
are then interested in computing the maximum value of `abs(p₂(x))` for
`x ∈ [0, π]`.

The interval `[0, π]` is split into two parts, `[0, ϵ]` and ´[ϵ, π]`.
On `[0 ϵ]` we use the asymptotic version of `F0(u0)` whereas on ´[ϵ,
π]` we use the non-asymptotic version.
"""
function delta0(
    u0::KdVZeroAnsatz;
    rtol = Arb(1e-2),
    atol = Arb(1e-3),
    threaded = true,
    verbose = false,
    maxevals = 10000, # A good value for this depends on u0.α
)
    # Determine a good choice of ϵ. Take it as large as possible so
    # that the asymptotic version still satisfies the required
    # tolerance when evaluated at ϵ
    ϵ = let ϵ = Arb(π), f = F0(u0, Asymptotic(), ϵ = ubound(Arb, ϵ))
        y = f(ϵ)[2]

        # Reduce ϵ until the value we get satisfies the required
        # tolerance
        while !ArbExtras.check_tolerance(y; rtol, atol)
            if ϵ > 1
                ϵ -= 0.05
            else
                ϵ /= 1.2
            end
            y = f(ϵ)[2]
            ϵ > 0.1 || error("could not determine working ϵ, last tried value was $ϵ")
        end
        ubound(ϵ)
    end

    verbose && @info "Determined ϵ" ϵ

    # Function for computing p₂(x) for x ∈ [0, ϵ]
    F0_asymptotic = F0(u0, Asymptotic(), ϵ = Arb(1.1ϵ))
    f(x) =
        let
            res = F0_asymptotic(x)
            @assert Arblib.valuation(res) == 2 # Check that p[0] and p[1] are zero
            truncate_with_remainder(res, u0.α, degree = 2)[2]
        end

    # Function for computing p₂(x) for x ∈ [ϵ, π]
    F0_nonasymptotic = F0(u0, Ball())
    g(x) =
        let
            res = F0_nonasymptotic(x)
            @assert Arblib.valuation(res) == 2 # Check that p[0] and p[1] are zero
            truncate_with_remainder(res, u0.α, degree = 2)[2]
        end

    # Compute an enclosure on [0, ϵ]
    p2_asymptotic = ArbExtras.maximum_enclosure(
        f,
        zero(ϵ),
        ϵ,
        degree = -1,
        abs_value = true;
        rtol,
        atol,
        threaded,
        verbose,
    )

    verbose && @info "Bound of p[2] from [0, ϵ]" p2_asymptotic

    # Compute an enclosure on [ϵ, π]
    p2_nonasymptotic = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        ubound(Arb(π)),
        degree = -1,
        abs_value = true;
        rtol,
        atol,
        maxevals,
        threaded,
        verbose,
    )

    verbose && @info "Bound of p[2] on [ϵ, π]" p2_nonasymptotic

    p2 = max(p2_asymptotic, p2_nonasymptotic)

    return ArbSeries((0, 0, p2))
end
