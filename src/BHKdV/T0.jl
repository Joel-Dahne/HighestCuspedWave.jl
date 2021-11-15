function T0(
    u0::BHKdVAnsatz,
    evaltype::Ball;
    δ0::Arb = Arb(1e-5),
    δ1::Arb = Arb(1e-5),
    δ2::Arb = Arb(1e-5),
    skip_div_u0 = false,
)
    f = T01(u0, evaltype, skip_div_u0 = true; δ0, δ1)
    g = T02(u0, evaltype, skip_div_u0 = true; δ2)

    return x -> begin
        ## Integral on [0, x]
        part1 = f(x)

        # Short circuit on non-finite value
        isfinite(part1) || return part1

        ## Integral on [x, π]
        part2 = g(x)

        isfinite(part2) || return part2 # Avoid computing u0(x)

        if skip_div_u0
            return part1 + part2
        else
            return (part1 + part2) / u0(x)
        end
    end
end

function T0(u0::BHKdVAnsatz, evaltype::Asymptotic; non_asymptotic_u0 = false, ϵ = Arb(2e-1))
    f = T01(u0, evaltype; non_asymptotic_u0, ϵ)
    g = T02(u0, evaltype; non_asymptotic_u0, ϵ)

    return x -> f(x) + g(x)
end

"""
    T013plusT021(u0::BHKdVAnsatz{Arb}, ::Ball; δ::Arb)

Compute the part of the integral around the singularity at `y = x` by
factoring out part of the weight and integrating explicitly.

With the change of coordinates `t = y / x` we get that the integral
for the norm is given by
```
x / (π * u0(x) * log(u0.c + inv(x))) *
    ∫ abs(clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t * log(u0.c + inv(x * t)) dt
```
We are interested in computing this integral around the point `t = 1`
where it has a singularity. More precisely we compute it on the
interval `[1 - δ, 1 + δ]`.

Since `log(u0.c + inv(x * t))` is well behaved as long as `x` is not
very close to zero we can compute an enclosure of it for all `t ∈ [1 -
δ, 1 + δ]` and factor this out. This leaves us with the integral
```
I = ∫ (clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt
```
taken from `1 - δ` to `1 + δ`. Where we have removed the absolute
value since the expression inside is positive on the interval.
- **PROVE:** That the expression inside is positive, this holds as
  long as `δ` is not too large.

We can integrate this explicitly, see [`T0_p_one`](@ref), giving us
```
∫ (clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt =
    (clausenc(x * (1 - t), 2 - α) / x^2 - t * clausens(x * (1 - t), 1 - α) / x) +
    (clausenc(x * (1 + t), 2 - α) / x^2 + t * clausens(x * (1 + t), 1 - α) / x) -
    2(clausenc(x * t, 2 - α) / x^2 + t * clausens(x * t, 1 - α) / x)
```
Call this function `primitive(t)`. We then have
```
I = primitive(1 + δ) - primitive(1 - δ)
```

In our case `α` is the interval `[-1, -1 + u0.ϵ]` and we let `s1 = 1 -
α` and `s2 = 2 - α`.

We can notice that `primitive(t)` contains a division by `x` for all
terms and that we in the end want to multiply the whole integral by
`x`. To cancel these two we therefore let
```
primitive_mul_x(t) = (clausenc(x * (1 - t), 2 - α) / x - t * clausens(x * (1 - t), 1 - α)) +
    (clausenc(x * (1 + t), 2 - α) / x + t * clausens(x * (1 + t), 1 - α)) -
    2(clausenc(x * t, 2 - α) / x + t * clausens(x * t, 1 - α))
```
which we can rewrite as
```
primitive_mul_x(t) = (
        clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) - 2clausenc(x * t, 2 - α)
    ) / x + t * (
        -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) + 2clausens(x * t, 1 - α)
    )
```

Now specialising for the case `t = 1 - δ` and `t = 1 + δ` we get
```
(
    clausenc(x * δ, 2 - α) + clausenc(x * (2 - δ), 2 - α) - 2clausenc(x * (1 - δ), 2 - α)
) / x + (1 - δ) * (
    -clausens(x * δ, 1 - α) + clausens(x * (2 - δ), 1 - α) + 2clausens(x * (1 - δ), 1 - α)
)
```
and
```
(
    clausenc(x * δ, 2 - α) + clausenc(x * (2 + δ), 2 - α) - 2clausenc(x * (1 + δ), 2 - α)
) / x + (1 + δ) * (
    clausens(x * δ, 1 - α) + clausens(x * (2 + δ), 1 - α) + 2clausens(x * (1 + δ), 1 - α)
)
```
where we have used that `clausenc` is even and `clausens` is odd.

**IMPROVE:** We can use `ArbSeries` to expand in `x` and compute a
tighter enclosure for most of it.

# Handling `x` close to `π`
If `x` is sufficiently close to `π` then `1 + δ > π / x` and our
approach needs to be modified. In that case we are only want to
integrate up to `π / x` and we are hence interested in computing
`primitive_mul_x(π / x)`, for which we have
```
primitive_mul_x(π / x) = 2(clausenc(x + π, 2 - α) - clausenc(Arb(π), 2 - α)) / x
```
where `clausenc(Arb(π), 2 - α)` can also be given as the negated
alternating zeta function, `-eta(2 - α)`. Evaluating `clausenc(x + π,
2 - α)` directly when `x` overlaps `π` gives an indeterminate result.
Instead we use that it is `2π` periodic and even to get `clausenc(π -
x, 2 - α)`. Since `π - x` is expected to be small, we only need this
for `x` close to `π` we want to use the expansion at `x = 0`. However
in our case `2 - α` overlaps with `3` and hence the expansion is
singular. We therefore use
[`clausenc_expansion_odd_s_singular`](@ref).
"""
function T013plusT021(
    u0::BHKdVAnsatz{Arb},
    evaltype::Ball = Ball();
    δ::Arb = Arb(1e-3),
    skip_div_u0 = false,
)
    ϵ = u0.ϵ

    # s1 = 1 + α
    s1 = Arb((2 - ϵ, 2))
    # s2 = 2 - α
    s2 = Arb((3 - ϵ, 3))

    return x -> begin
        primitive(t) =
            (clausenc(x * (1 - t), s2) / x^2 - t * clausens(x * (1 - t), s1) / x) +
            (clausenc(x * (1 + t), s2) / x^2 + t * clausens(x * (1 + t), s1) / x) -
            2(clausenc(x * t, s2) / x^2 + t * clausens(x * t, s1) / x)

        primitive_mul_x(t) =
            (clausenc(x * (1 - t), s2) / x - t * clausens(x * (1 - t), s1)) +
            (clausenc(x * (1 + t), s2) / x + t * clausens(x * (1 + t), s1)) -
            2(clausenc(x * t, s2) / x + t * clausens(x * t, s1))

        primitive_1mδ_mul_x =
            (clausenc(x * δ, s2) + clausenc(x * (2 - δ), s2) - 2clausenc(x * (1 - δ), s2)) / x +
            (1 - δ) * (
                -clausens(x * δ, s1) + clausens(x * (2 - δ), s1) -
                2clausens(x * (1 - δ), s1)
            )

        primitive_1pδ_mul_x =
            (clausenc(x * δ, s2) + clausenc(x * (2 + δ), s2) - 2clausenc(x * (1 + δ), s2)) / x +
            (1 + δ) * (
                clausens(x * δ, s1) + clausens(x * (2 + δ), s1) -
                2clausens(x * (1 + δ), s1)
            )

        # IMPROVE: We don't need to compute this in general
        primitive_pidivx_mul_x = let
            # Compute an enclosure of clausenc(π - x, s2)
            if Arblib.overlaps(x, Arb(π))
                # Enclosure of π - x
                pimx = Arblib.nonnegative_part!(zero(x), π - x)

                # Compute the analytic terms in the expansion
                C, _, p, E = clausenc_expansion(ubound(Arb, pimx), s2, 3)

                @assert !isfinite(C)

                # Set the singular term to zero, we handle it separately.
                @assert !isfinite(p[2])
                p[2] = 0

                # Factor out x^1.5
                exponent = Arb(3 // 2)
                D = clausenc_expansion_odd_s_singular(ubound(Arb, pimx), s2, exponent)
                clausenc_pimx = D * pimx^exponent + p(pimx)
            else
                clausenc_pimx = clausenc(π - x, s2)
            end

            2(clausenc_pimx + real(eta(Acb(s2)))) / x
        end

        @assert Arblib.overlaps(primitive(1 - δ) * x, primitive_1mδ_mul_x)
        @assert Arblib.overlaps(primitive(1 + δ) * x, primitive_1pδ_mul_x)
        @assert Arblib.overlaps(primitive(π / x) * x, primitive_pidivx_mul_x)

        if 1 + δ < π / x
            # Standard case, when the upper integration limit is 1 + δ
            I_mul_x = primitive_1pδ_mul_x - primitive_1mδ_mul_x
        elseif 1 + δ > π / x
            # In this case we only want to integrate up to π / x
            # instead of 1 + δ
            I_mul_x = primitive_pidivx_mul_x - primitive_1mδ_mul_x
        elseif contains(π / x, 1 + δ)
            # In this case the integration limit depends on the
            # precise choice of x. However we can get an enclosure by
            # just using π / x
            I_mul_x = primitive_pidivx_mul_x - primitive_1mδ_mul_x
        else
            # In this case the integration limit also depends on the
            # precise choice of x. We could take the union of
            # primitive_1pδ_mul_x and primitive_pidivx_mul_x but the
            # first one of these will always be non-finite, so we just
            # return an indeterminate result directly.
            # IMPROVE: We could compute an enclosure in this case as
            # well, however I'm not sure it ever happens?
            @warn "π / x and 1 + δ overlaps" π / x 1 + δ
            return Arblib.indeterminate!(zero(x))
        end

        # Compute enclosure of log(u0.c + inv(x * t))
        t_interval = Arblib.add_error!(one(x), δ)
        weight_enclosure = log(u0.c + inv(x * t_interval))

        # Enclosure of factor in front of integral divided by x, excluding u0
        factor_div_x = inv(π * log(u0.c + inv(x)))

        res = factor_div_x * weight_enclosure * I_mul_x

        # Other version for computing it
        # T013 uses the variable t as well so δ1 = δ
        f = T013(u0, evaltype, δ1 = δ, skip_div_u0 = true)
        # T021 uses the variable y directly, so 1 + δ corresponds to x
        # + x * δ, hence δ2 = x * δ.
        g = T021(u0, evaltype, δ2 = x * δ, skip_div_u0 = true)
        part1 = abs(f(x))
        part2 = abs(g(x))
        res2 = part1 + part2

        @assert Arblib.overlaps(res, res2)

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
