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

function T0_alternative(
    u0::BHKdVAnsatz,
    evaltype::Ball = Ball();
    δ0::Arb = Arb(1e-3),
    δ::Arb = Arb(1e-1),
    skip_div_u0 = false,
)
    # Integration on [0, δ0]
    f1 = T011(u0, evaltype, skip_div_u0 = true; δ0)
    # Integration on [δ0, 1 - δ]
    f2 = T012(u0, evaltype, skip_div_u0 = true, δ1 = δ; δ0)
    # Integration on [1 - δ, b / x] with b picked later
    f3 = T0_primitive(u0, evaltype, skip_div_u0 = true)
    # Integration on [b, π] with b picked later
    f4 = T022(u0, evaltype, skip_div_u0 = true)
    # Construction depends on x

    return x -> begin
        # We take the tolerance used for the integration to be twice
        # the diameter of x but no less than 1e-5
        tol = max(Arb(1e-5), 4radius(Arb, x))

        part1 = f1(x)

        isfinite(part1) || return part1

        part2 = f2(x; tol)

        isfinite(part2) || return part2

        # We want to take the upper integration independent of x
        # in the variables used by T022. We take it to be ubound(x + δ)
        b_y = ubound(Arb, x + δ) # Upper bound in the variable y
        if b_y < π
            # Compute upper bound in the variable t = y / x
            b_t = b_y / x
            part3 = f3(x, 1 - δ, b_t)

            isfinite(part3) || return part3

            part4 = f4(x, b_y; tol)
        else
            # The value 1 + δ is not important
            part3 = f3(x, 1 - δ, 1 + δ, to_endpoint = true)

            part4 = zero(x)
        end

        res = part1 + part2 + part3 + part4

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T0_primitive(u0::BHKdVAnsatz{Arb}, ::Ball)

With the change of coordinates `t = y / x` the integral for the norm
is given by
```
x / (π * u0(x) * log(u0.c + inv(x))) *
    ∫ abs(clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t * log(u0.c + inv(x * t)) dt
```
This method returns a function `f` such that `f(x, a, b)` computes the
above, integrating `a` to `b`. Optionally it accepts the keyword
argument `to_endpoint` which if true makes it compute the integral
from `a` to `π / x`, i.e. `f(x, a, b, to_endpoint = true)`.

The computation is done by factoring out part of the weight and
integrating the rest explicitly.

More precisely, since `log(u0.c + inv(x * t))` is well behaved as long
as `x` is not very close to zero we can compute an enclosure of it for
all `t ∈ [1 - δ, 1 + δ]` and factor this out. This leaves us with the
integral
```
I = ∫ (clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt
```
taken from `a` to `b` (or `π / x`). Where we have removed the absolute
value since the expression inside is positive on the interval.
- **TODO:** We need to check the positivity of the integrand. This
  could be done by checking the endpoint `a`. However for it to work
  well for wide values of `x` we might need to do some more work.

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

By reorganising the terms we can rewrite `primitive(t)` as
```
primitive(t) =
    (clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) - 2clausenc(x * t, 2 - α)) / x^2 -
    t * (clausens(x * (1 - t), 1 - α) - clausens(x * (1 + t), 1 - α) + 2clausens(x * t, 1 - α)) / x
```
We can also improve the computed enclosures by replacing `clausenc`
with `clausencmzeta`. The constants always cancel out and
`clausencmzeta` has better monotonicity properties in the parameter,
giving tighter enclosures. Giving us
```
primitive(t) =
    (clausencmzeta(x * (1 - t), 2 - α) + clausencmzeta(x * (1 + t), 2 - α) - 2clausencmzeta(x * t, 2 - α)) / x^2 -
    t * (clausens(x * (1 - t), 1 - α) - clausens(x * (1 + t), 1 - α) + 2clausens(x * t, 1 - α)) / x
```

We can notice that `primitive(t)` contains a division by `x` for all
terms and that we in the end want to multiply the whole integral by
`x`. To cancel these two we therefore let
```
primitive_mul_x(t) =
    (clausencmzeta(x * (1 - t), 2 - α) + clausencmzeta(x * (1 + t), 2 - α) - 2clausencmzeta(x * t, 2 - α)) / x -
    t * (clausens(x * (1 - t), 1 - α) - clausens(x * (1 + t), 1 - α) + 2clausens(x * t, 1 - α))
```

Now, we are interested in computing `primitive_mul_x(b) -
primitive_mul_x(a)`. Since `a` and `b` are expect to be fairly close
together we can expect cancellations between the two terms. To help in
handling these cancellations we write up the full expression and
reorder the terms to put those parts with cancellations together.
This gives us
```
primitive_mul_x(b) - primitive_mul_x(a) =
    (
        (clausencmzeta(x * (1 - b), 2 - α) + clausencmzeta(x * (1 + b), 2 - α) - 2clausencmzeta(x * b, 2 - α)) -
        (clausencmzeta(x * (1 - a), 2 - α) + clausencmzeta(x * (1 + a), 2 - α) - 2clausencmzeta(x * a, 2 - α))
    ) / x - (
        b * (clausens(x * (1 - b), 1 - α) - clausens(x * (1 + b), 1 - α) + 2clausens(x * b, 1 - α)) -
        a * (clausens(x * (1 - a), 1 - α) - clausens(x * (1 + a), 1 - α) + 2clausens(x * a, 1 - α))
    )
```
There are two main terms which both are positive (if we include the
minus sign in the second term) so there are no cancellations between
them and we therefore enclose them separately. For most values of `x`
the two terms are monotone in `x` and we can therefore compute an
enclosure using `ArbExtras.extrema_series` with `degree = 0` which
picks up the monotonicity.

# Value at `t = π / x`
In the special case that `t = π / x` this simplifies to
```
primitive_mul_x(π / x) = 2(clausencmzeta(x + π, 2 - α) - clausencmzeta(Arb(π), 2 - α)) / x
    = 2(clausenc(x + π, 2 - α) - clausenc(Arb(π), 2 - α)) / x
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

# Handling wide values of `a` and `b`
If either `a` or `b` is a wide ball we can compute a tighter enclosure
by using that `primitive_mul_x(t)` is increasing in `t`. That it is
increasing is an easy consequence of being the primitive function of a
positive integrand.

**IMPROVE:** We can use `ArbSeries` to expand in `x` and compute a
tighter enclosure for most of it.
"""
function T0_primitive(u0::BHKdVAnsatz{Arb}, evaltype::Ball = Ball(); skip_div_u0 = false)
    #Enclosure of Arb((2 - u0.ϵ, 3)) computed in a way so that the
    #upper endpoint is exactly 2
    s1 = 2 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))
    #Enclosure of Arb((3 - u0.ϵ, 3)) computed in a way so that the
    #upper endpoint is exactly 3
    s2 = 3 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    return (x::Arb, a::Arb, b::Arb; to_endpoint = false) -> begin
        # Check that the integrand is positive on the left endpoint
        # PROVE: We would need to prove that it is enough to check the
        # lower bound of x
        let t = lbound(Arb, a), x = lbound(Arb, x), α = Arb((-1, -1 + u0.ϵ))
            @assert Arblib.ispositive(
                clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) -
                2clausenc(x * t, -α),
            )
        end

        primitive(t) =
            (
                clausencmzeta(x * (1 - t), s2) + clausencmzeta(x * (1 + t), s2) -
                2clausencmzeta(x * t, s2)
            ) / x^2 -
            t * (
                clausens(x * (1 - t), s1) - clausens(x * (1 + t), s1) +
                2clausens(x * t, s1)
            ) / x

        primitive_mul_x(t) =
            let
                # Compute enclosure of all clausenc terms using ArbSeries
                term1_f(x) =
                    (
                        clausencmzeta(x * (1 - t), s2) + clausencmzeta(x * (1 + t), s2) - 2clausencmzeta(x * t, s2)
                    ) / x

                term1 = Arb(
                    ArbExtras.extrema_series(term1_f, getinterval(x)..., degree = 0)[1:2],
                )

                # Compute enclosure of all clausens terms using ArbSeries
                term2_f(x) =
                    t * (
                        clausens(x * (1 - t), s1) - clausens(x * (1 + t), s1) +
                        2clausens(x * t, s1)
                    )

                term2 = Arb(
                    ArbExtras.extrema_series(term2_f, getinterval(x)..., degree = 0)[1:2],
                )

                return term1 - term2
            end

        # primitive_mul_x(b) - primitive_mul_x(a)
        primitive_bma_mul_x(a, b) =
            let
                # Compute enclosure of all clausenc terms using ArbSeries
                term1_f(x) =
                    (
                        (
                            clausencmzeta(x * (1 - b), s2) +
                            clausencmzeta(x * (1 + b), s2) - 2clausencmzeta(x * b, s2)
                        ) - (
                            clausencmzeta(x * (1 - a), s2) +
                            clausencmzeta(x * (1 + a), s2) - 2clausencmzeta(x * a, s2)
                        )
                    ) / x
                term1 = Arb(
                    ArbExtras.extrema_series(term1_f, getinterval(x)..., degree = 0)[1:2],
                )

                # Compute enclosure of all clausens terms using ArbSeries
                term2_f(x) = (
                    b * (
                        clausens(x * (1 - b), s1) - clausens(x * (1 + b), s1) +
                        2clausens(x * b, s1)
                    ) -
                    a * (
                        clausens(x * (1 - a), s1) - clausens(x * (1 + a), s1) +
                        2clausens(x * a, s1)
                    )
                )
                term2 = Arb(
                    ArbExtras.extrema_series(term2_f, getinterval(x)..., degree = 0)[1:2],
                )

                return term1 - term2
            end

        if b < π && !to_endpoint
            # We integrate on [a, b]
            weight_enclosure = log(u0.c + inv(x * Arb((a, b))))

            I_mul_x = if iswide(b)
                union(
                    primitive_bma_mul_x(a, lbound(Arb, b)),
                    primitive_bma_mul_x(a, ubound(Arb, b)),
                )
            else
                primitive_bma_mul_x(a, b)
            end
        elseif to_endpoint || b * x > π || contains(π / x, b)
            # We integrate on [a, π / x]
            weight_enclosure = log(u0.c + inv(x * Arb((a, π / x))))

            # Compute primitive_mul_x(a)
            primitive_a_mul_x = primitive_mul_x(a)

            # Enclosure of primitive_mul_x(π / x)
            primitive_pidivx_mul_x = let
                # Enclosure of π - x using that x <= π
                pimx = Arblib.nonnegative_part!(zero(x), π - x)

                # Compute an enclosure of clausenc(π - x, s2)
                if Arblib.overlaps(x, Arb(π))
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
                    clausenc_pimx = clausenc(pimx, s2)
                end

                2(clausenc_pimx + real(eta(Acb(s2)))) / x
            end

            I_mul_x = primitive_pidivx_mul_x - primitive_mul_x(a)
        else
            # In this case it is unclear what integration limit we
            # should use. Make a warning and return NaN. This should
            # not happen in practice.
            @warn "Unclear integration limit" b π / x
            return Arblib.indeterminate!(zero(x))
        end

        # Enclosure of factor in front of integral divided by x, excluding u0
        factor_div_x = inv(π * log(u0.c + inv(x)))

        res = factor_div_x * weight_enclosure * I_mul_x

        # Other version for computing it, can be used for checking.
        # Though it doesn't support all cases.
        if false
            # T013 uses the variable t as well so we get δ1 = 1 - a
            f = T013(u0, evaltype, δ1 = 1 - a, skip_div_u0 = true)
            # T021 uses the variable y directly, so b corresponds to x *
            # b, hence δ2 = x * (b - 1).
            g = T021(u0, evaltype, δ2 = x * (b - 1), skip_div_u0 = true)
            part1 = abs(f(x))
            part2 = abs(g(x))
            res2 = part1 + part2

            @assert Arblib.overlaps(res, res2)
        end

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
