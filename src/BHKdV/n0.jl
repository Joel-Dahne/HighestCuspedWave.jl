"""
    n0_bound(u0::BHKdVAnsatz{Arb})

Compute an upper bound of `n₀` from the paper. This is the supremum of
```
N(x) = u0.w(x) / 2u0(x)
```
for `0 < x < π`.

It splits the interval `[0, π]` into two parts, `[0, ϵ]` and `[ϵ, π]`
with `ϵ < 1`.

# The interval `[0, ϵ]`
As a first step we split this into three factors by multiplying and
dividing by `gamma(1 + α) * x^(-α) * (1 - x^p0)`, giving us.
```
u0.w(x) / 2u0(x) =  gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x) *
    u0.w(x) / (2gamma(1 + α) * x^(-α) * (1 - x^p0))
```
Let
```
F(x) = gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x)
G(x) = u0.w(x) / (2gamma(1 + α) * x^(-α) * (1 - x^p0))
```

## Handling `F(x)`
This is computed using [`inv_u0_bound`](@ref).

## Handling `G(x)`
Inserting
```
u0.w(x) = x^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x))
```
gives us
```
G(x) = x^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x)) / (2gamma(1 + α) * x^(-α) * (1 - x^p0))
     = x^((1 - u0.γ) * (1 + α)) * log(u0.c + inv(x)) / (2gamma(1 + α) * (1 - x^p0))
```
Multiplying and dividing by `log(inv(x))` and using that `gamma(1 + α)
= gamma(2 + α) / (1 + α)` we can split this into three factors as
```
G1 = inv(gamma(2 + α))
G2(x) = log(u0.c + inv(x)) / 2log(inv(x)) = -log(u0.c + inv(x)) / 2log(x)
G3(x) = x^((1 - u0.γ) * (1 + α)) * log(inv(x)) / (gamma(1 + α) * (1 - x^p0))
```
We can enclose `G1(x)` directly and `G2(x)` using that it is
increasing in `x` for `0 < x < 1` and that the limit as `x -> 0` is `1
/ 2`.

The function `G3(x)` occurs in
[`lemma_bhkdv_weight_div_asymptotic_enclosure`](@ref) and we can get
an enclosure from that.

# The interval `[ϵ, π]`
Here we use that `u0.v0(x)` gives a **lower** bound for `u0(x)`. This
is the statement of [`lemma_bhkdv_monotonicity_alpha`](@ref).

If `u0.v0(x)` is positive we then get that `u0.w(x) / 2u0(x)` is
**upper** bounded by `u0.w(x) / 2u0.v0(x)`. The approach is therefore
the same as in the version for `BHAnsatz` which the only difference
being the different weight and that we only get an upper bound instead
of an enclosure.

To prove that `u0.v0(x)` is positive it is enough to ensure that
`u0.v0(π)` is positive and that the bound for `u0.w(x) / 2u0.v0(x)` on
the interval is finite. In this case `u0.v0(x)` can never be zero.
"""
function n0_bound(u0::BHKdVAnsatz{Arb}; rtol = Arb("1e-2"), verbose = false)
    verbose && @info "Computing bound of n0"

    # Function for asymptotic evaluation

    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert Arblib.overlaps(
            u0.w(x),
            x^(1 - u0.γ * (1 + Arb((-1, -1 + u0.ϵ)))) * log(u0.c + inv(x)),
        )
    end
    F = inv_u0_bound(u0, ϵ = Arb(1 // 2))

    # inv(gamma(2 + α)) = rgamma(1 + u0.ϵ)
    G1 = rgamma(1 + u0.ϵ)

    G2(x) =
        if Arblib.contains_zero(x)
            xᵤ = ubound(Arb, x)
            Arb((1 // 2, -log(u0.c + inv(xᵤ)) / 2log(xᵤ)))
        else
            -log(u0.c + inv(x)) / 2log(x)
        end

    G3(x) = lemma_bhkdv_weight_div_asymptotic_enclosure(u0)

    f(x) = F(x) * G1 * G2(x) * G3(x)

    # Function for non-asymptotic evaluation

    # Assert that the lemma holds
    @assert lemma_bhkdv_monotonicity_alpha(u0)

    g(x) = u0.w(x) / 2u0.v0(x)

    # We expect the maximum to be attained at x = π. We therefore
    # first evaluate g at this point. We then find ϵ such that
    # f(Arb((0, ϵ))) is bounded by this value.

    gπ = g(Arb(π))

    # Check that the value at π is positive, so that u0.w(x) / 2u0(x)
    # indeed is upper bounded by u0.w(x) / 2u0.v0(x). This is still
    # provided that the resulting bound is finite. But if it is not
    # finite it doesn't matter anyway.
    Arblib.ispositive(gπ) || error("value at π not positive")

    verbose && @info "Value at x = π" gπ

    ϵ = Arf(0.45)
    while !(f(Arb((0, ϵ))) < gπ)
        # We use the factor 0.8 instead of, say, 0.5 to get
        # slightly better (higher) values for ϵ.
        ϵ *= 0.8

        # We use the non-asymptotic version on [ϵ, π] so we don't want
        # ϵ to be too small
        ϵ < 1e-3 && error("ϵ too small")
    end

    verbose && @info "ϵ was determined to be" ϵ

    # Finally we compute the maximum on the interval [ϵ, π]

    n0 = ArbExtras.maximum_enclosure(
        g,
        ϵ,
        ubound(Arb(π)),
        degree = 0, # Enough to make use of monotonicity
        abs_value = true,
        point_value_max = abs(g(Arb(π))), # Maximum is in practice attained here
        depth_start = 4,
        threaded = true;
        rtol,
        verbose,
    )

    verbose && @info "Computed bound" n0

    return n0
end
