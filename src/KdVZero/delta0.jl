"""
    delta0_bound(u0::KdVZeroAnsatz; rtol, atol, maxevals, threaded, verbose)

Compute an upper bound of `abs(F0(u0)(x))`.

This method works differently depending on if `u0.α0 = 0` or not. In
the case that `u0.α0 = 0` any enclosure of the defect will contain
zero and hence not work. Instead we compute an expansion in `α` of the
defect in that case.

# `u0.α0 = 0`
In this case we compute a Taylor model in `α` such that it gives an
enclosure of `abs(F0(u0)(x))` for all `α ∈ u0.α`.

For a given value of `x` [`F0`](@ref) gives us a Taylor model in `α`.
We truncate this expansion to degree `2`, giving us a polynomial of
the form
```
0 + 0 * α + Δδ(x) * α^2
```
that gives an enclosure of `F0(u0)(x)` for every `α ∈ u0.α`. We are
then interested in computing the maximum value of `abs(Δδ(x))` for `x
∈ [0, π]`.

It uses the asymptotic version of `F0(u0)` on the entire interval `[0,
π]`.

# `u0.α0 < 0`
In this case we compute an enclosure of an upper bound of
`abs(F0(u0))` for `x ∈ [0, π]`.
"""
function delta0_bound(
    u0::KdVZeroAnsatz;
    rtol = Arb(1e-2),
    atol = ifelse(iszero(u0.α0), Arb(1e-3), Arb(0)),
    maxevals = 10000, # A good value for this depends on u0.α
    threaded = true,
    verbose = false,
)
    if iszero(u0.α0)
        verbose && @info "Computing enclosure of Δδ"

        # Function for computing Δ(x)
        F0_asymptotic = F0(u0, Asymptotic(), ϵ = Arb(3.2))
        f = x -> let
            res = F0_asymptotic(x)
            @assert Arblib.valuation(res.p) == 2 # Check that p[0] and p[1] are zero
            truncate(res, degree = 1).p[2]
        end

        # Compute an enclosure on [0, ϵ]
        Δδ = ArbExtras.maximum_enclosure(
            f,
            zero(Arf),
            ubound(Arb(π)),
            degree = -1,
            abs_value = true;
            maxevals,
            rtol,
            atol,
            threaded,
            verbose,
        )

        verbose && @info "Computed enclosure" Δδ

        return TaylorModel(ArbSeries((0, 0, Δδ), degree = 2), u0.α, u0.α0)
    else
        # Determine a good choice of ϵ. Take it as large as possible
        # so that the asymptotic version still satisfies the required
        # tolerance when evaluated at ϵ or gives a better bound than
        # the non-asymptotic version.
        ϵ = let ϵ = Arb(π), f = F0(u0, Asymptotic(), ϵ = ubound(Arb, ϵ)), g = F0(u0)
            y = f(ϵ)(u0.α)
            z = g(ϵ)(u0.α)

            # Reduce ϵ until the value we get either satisfies the
            # required tolerance or is better than the non-asymptotic
            # version.
            while !ArbExtras.check_tolerance(y; rtol, atol) && radius(y) > radius(z)
                if ϵ > 1
                    ϵ -= 0.05
                else
                    ϵ /= 1.2
                end

                y = f(ϵ)(u0.α)
                z = g(ϵ)(u0.α)
                ϵ > 0.1 ||
                    error("could not determine working ϵ, last tried value was $ϵ")
            end
            ubound(ϵ)
        end

        verbose && @info "Determined ϵ" ϵ

        # Function for computing F0(u0) for x ∈ [0, ϵ]
        f = F02(u0, Asymptotic(), ϵ = Arb(1.1ϵ))

        # Function for computing F0(u0) for x ∈ [ϵ, π]
        F0_nonasymptotic = F0(u0, Ball())
        g = x -> F0_nonasymptotic(x)(u0.α)

        # Compute an enclosure on [0, ϵ]
        res_asymptotic = ArbExtras.maximum_enclosure(
            f,
            zero(ϵ),
            ϵ,
            abs_value = true;
            rtol,
            atol,
            threaded,
            verbose,
        )

        verbose && @info "Bound on [0, ϵ]" res_asymptotic

        if ϵ >= Arb(π)
            # If ϵ = π then we don't need to compute any non-asymptotic
            # part
            verbose && @info "ϵ = π so [ϵ, π] is skipped"
            res = res_asymptotic
        else
            # Compute an enclosure on [ϵ, π]
            res_nonasymptotic = ArbExtras.maximum_enclosure(
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

            verbose && @info "Bound on [ϵ, π]" res_nonasymptotic

            res = max(res_asymptotic, res_nonasymptotic)
        end

        return res
    end
end
