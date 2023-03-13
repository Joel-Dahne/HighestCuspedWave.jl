function T0(u0::BHKdVAnsatz, evaltype::Ball; δ::Arb = Arb(1e-1), skip_div_u0 = false)
    # Integration on [0, 1 - δ]
    f1 = T012(u0, evaltype, skip_div_u0 = true; δ)
    # Integration on [1 - δ, b / x] with b picked later
    f2 = T0_primitive(u0, evaltype, skip_div_u0 = true)
    # Integration on [b, π] with b picked later
    f3 = T022(u0, evaltype, skip_div_u0 = true)

    return x -> begin
        # We take the tolerance used for the integration to be twice
        # the diameter of x but no less than 1e-3
        tol = max(Arb(1e-3), 4radius(Arb, x))

        part1 = f1(x; tol)

        isfinite(part1) || return part1

        # We want to take b, lower integration limit in T022,
        # independent of x in the variables used by T022. We take it
        # to be ubound((1 + δ) * x)
        b = ubound(Arb, (1 + δ) * x) # Upper bound in the variable y
        # We check (b / x) * x < π instead of b < π to make sure that
        # the inequality b_t * x < π holds even when x is a ball.
        if (b / x) * x < π
            # Compute upper bound of b in the variable t = y / x
            b_t = b / x
            part2 = f2(x, 1 - δ, b_t)

            isfinite(part2) || return part2

            part3 = f3(x, b; tol)
        else
            # In this case use b = π. This means that f3 is zero and
            # we want to integrate f2 all the way to the endpoint.

            # The valueargument 1 + δ is not important since it is not
            # used
            part2 = f2(x, 1 - δ, 1 + δ, to_endpoint = true)

            part3 = zero(x)
        end

        res = part1 + part2 + part3

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T0(u0::BHKdVAnsatz{Arb}, ::Asymptotic; ϵ::Arb = 0.5)

Return a function such that `T0(u0, Asymptotic())(x)` computes the
integral `T0(u0)(x)` using an evaluating strategy that works
asymptotically as `x` goes to zero. In general it only computes an
**upper bound** of it.

More precisely it computes
```
1 / (π * u0.w(x) * u0(x)) *
    ∫ abs(clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * u0.w(y) dy
```
with the integration going from `0` to `π`.

Then change of variables `t = y / x` gives us
```
x / (π * u0.w(x) * u0(x)) *
    ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * u0.w(x * t) dt
```
integrated from `0` to `π / x`. Using that
```
u0.w(x) = x^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x))
```
we can simplify this to
```
x / (π * log(u0.c + inv(x)) * u0(x)) *
    ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
This is the expression we are interested in computing an **upper
bound** of.

# Split into three factors
Similarly to in the asymptotic version of [`F0`](@ref) we split the
expression into three factors which we bound separately. We write it
as
```
inv(π)
* (gamma(1 + α) * x^(-α) * (1 - x^p0) / (π * u0(x)))
* (log(inv(x)) / log(u0.c + inv(x)))
* (
    x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(inv(x))) *
    ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
)
```
The first two factors, after `inv(π)`, are the same as in [`F0`](@ref)
and we handle them in the same way. For the third factor we use the
notation
```
W(x) = x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(inv(x)))
I(x) = ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```

# Computing `W(x) * I(x)`
As a first step we expand the integrand. We have the following
expansions for the Clausen functions in the integrand
```
clausenc(x * (1 - t), -α) = sinpi(-α / 2) * gamma(1 + α) * x^(-α - 1) * abs(1 - t)^(-α - 1) + P(x * abs(1 - t))
clausenc(x * (1 + t), -α) = sinpi(-α / 2) * gamma(1 + α) * x^(-α - 1) * (1 + t)^(-α - 1) + P(x * (1 + t))
clausenc(x * t, -α) = sinpi(-α / 2) * gamma(1 + α) * x^(-α - 1) * t^(-α - 1) + P(x * t)
```
where the remainder term `P` contains one constant term and `P(x *
abs(1 - t)) + P(x * (1 + t)) - 2P(x * t)` behaves like `O(x^2)`.

Inserting this into the integral allows us to split it into one main
integral
```
I_M(x) = sinpi(-α / 2) * gamma(1 + α) * x^(-α - 1) *
    ∫ abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
, where we have used that `sinpi(-α / 2) * gamma(1 + α) * x^(-α - 1)`
is positive to allow us to move it outside of the absolute value, and
one remainder integral
```
I_R(x) = x^(1 + α) * ∫ abs(P(x * abs(1 - t)) + P(x * (1 + t)) - 2P(x * t)) *
    t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
They satisfy `I(x) <= I_M(x) + I_R(x)`.

Letting
```
G1(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫_0^1 abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt

G2(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫_1^(π / x) abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt

R(x) = ∫ abs(P(x * abs(1 - t)) + P(x * (1 + t)) - 2P(x * t)) *
            t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
We have
```
W(x) * I(x) <= sinpi(-α / 2) * (G1(x) + G2(x)) +
    x^(1 + α) / (gamma(1 + α) * log(inv(x)) * (1 - x^p0)) * R(x)
```
Bounds of the functions `G1(x), G2(x), R(x)` are implemented in
[`_T0_asymptotic_main_1`](@ref), [`_T0_asymptotic_main_2`](@ref) and
[`_T0_asymptotic_remainder`](@ref) respectively. The only non-trivial
part remaining is bounding
```
inv(gamma(1 + α) * (1 - x^p0))
```
We rewrite it as
```
inv(gamma(2 + α) / (1 + α) * (1 - x^p0)) = inv(gamma(2 + α)) * inv((1 - x^p0) / (1 + α))
```
and handle the removable singularity with [`fx_div_x`](@ref).
"""
function T0(u0::BHKdVAnsatz{Arb}, ::Asymptotic; ϵ::Arb = Arb(0.5))
    # Enclosure of α
    α = Arb((-1, -1 + u0.ϵ))
    # Enclosure of α + 1, avoiding spurious negative parts
    αp1 = Arblib.nonnegative_part!(zero(α), Arb((0, u0.ϵ)))

    # Function for bounding gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x)
    f1 = inv_u0_bound(u0, M = 3; ϵ)

    # Function for enclosing log(inv(x)) / log(u0.c + inv(x))
    f2 = x -> if iszero(x)
        one(x)
    elseif Arblib.contains_zero(x)
        lower = let xᵤ = ubound(Arb, x)
            log(inv(xᵤ)) / log(u0.c + inv(xᵤ))
        end
        upper = one(x)
        Arb((lower, upper))
    else
        log(inv(x)) / log(u0.c + inv(x))
    end

    # Enclosure of sinpi(-α / 2)
    G_factor = sinpi(-α / 2)

    # Function for computing an upper bound of G1 and G2
    # respectively
    G1 = _T0_asymptotic_main_1(α, u0.γ, u0.c)
    G2 = _T0_asymptotic_main_2(α, u0.γ, u0.c)

    # Function for computing enclosure of
    # x^(1 + α) / (gamma(1 + α) * log(inv(x)) * (1 - x^p0))
    R_factor =
        x -> begin
            # Enclosure of inv(log(inv(x))) = -inv(log(x))
            invloginvx = if iszero(x)
                zero(x)
            elseif Arblib.contains_zero(x)
                -Arb((inv(log(ubound(Arb, x))), 0))
            else
                -inv(log(x))
            end

            # inv(gamma(1 + α) * (1 - x^p0))
            inv_gamma_1mxp0 = if iszero(x)
                rgamma(αp1)
            elseif Arblib.contains_zero(x)
                lower = rgamma(αp1)
                upper = let xᵤ = ubound(Arb, x)
                    rgamma(2 + α) / fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), αp1, extra_degree = 2)
                end
                Arb((lower, upper))
            else
                rgamma(2 + α) / fx_div_x(s -> 1 - x^(s + s^2 / 2), αp1, extra_degree = 2)
            end

            return abspow(x, αp1) * invloginvx * inv_gamma_1mxp0
        end

    # Function for computing enclosure of R
    R = _T0_asymptotic_remainder(α, u0.γ, u0.c)

    return x::Arb -> begin
        I_M = G_factor * (G1(x) + G2(x))

        I_R = R_factor(x) * R(x)

        I = I_M + I_R

        f1(x) * f2(x) / π * I
    end
end

"""
    T0_primitive(u0::BHKdVAnsatz{Arb}, ::Ball)

With the change of coordinates `t = y / x` the integral for the norm
is given by
```
x / (π * u0(x) * log(u0.c + inv(x))) *
    ∫ abs(clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
where for the full norm the integral is taken from `0` to `π / x`.
This method returns a function `f` such that `f(x, a, b)` computes the
above, integrating from `a` to `b`.

It accepts the keyword argument `to_endpoint` which if true makes it
compute the integral from `a` to `π / x`, i.e. `f(x, a, b, to_endpoint
= true)`. In that case the value for `b` is not used.

If `to_endpoint = false` then we require that `b * x < π`.

The value of `a` should be larger than the unique root of the
integrand on ``[0, 1]`` (see [`lemma_I_hat_root`](@ref)). This allows
us to remove the absolute value in the integral. To ensure that this
is the case we check that the integrand is positive at `t = a`. Since
the location of the zero on ``[0, 1]`` is decreasing in `x` it is
enough to check this for a lower bound of `x`, this makes it easier to
handle wide values of `x`.

The computation is done by factoring out part of the weight and
integrating the rest explicitly.

More precisely, since `t^(-u0.γ * (1 + α)) * log(u0.c + inv(x * t))`
is well behaved as long as `x` is not very close to zero we can
compute an enclosure of it for all `t` in the interval and factor this
out. This gives us
```
x / (π * u0(x) * log(u0.c + inv(x))) *
    Arb((a, b))^(-u0.γ * (1 + α)) * log(u0.c + inv(x * Arb((a, b)))) * I
```
where `I` is the integral
```
I = ∫ (clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt
```
taken from `a` to `b` (or `π / x` if `to_endpoint` is true).

# Computing `I`
We can integrate `I` explicitly, see [`T0_p_one`](@ref), giving us
```
∫ (clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt =
    (clausenc(x * (1 - t), 2 - α) / x^2 - t * clausens(x * (1 - t), 1 - α) / x) +
    (clausenc(x * (1 + t), 2 - α) / x^2 + t * clausens(x * (1 + t), 1 - α) / x) -
    2(clausenc(x * t, 2 - α) / x^2 + t * clausens(x * t, 1 - α) / x)
```
Call this function `primitive(t)`. We then have
```
I = primitive(b) - primitive(a)
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
them and we therefore enclose them separately. We compute a tighter
enclosure using [`enclosure_series`](@ref) which typically picks up
the monotonicity.

## Value at `t = π / x`
In the special case that `t = π / x` things simplify, we get
```
primitive_mul_x(π / x) = 2(clausencmzeta(x + π, 2 - α) - clausencmzeta(Arb(π), 2 - α)) / x
    = 2(clausenc(x + π, 2 - α) - clausenc(Arb(π), 2 - α)) / x
```
where `clausenc(Arb(π), 2 - α)` can also be given as the negated
alternating zeta function, `-eta(2 - α)`. Evaluating `clausenc(x + π,
2 - α)` directly when `x` overlaps `π` gives an indeterminate result.
Instead we use that it is `2π` periodic and even to get `clausenc(π -
x, 2 - α)`.

# Handling wide values of `a` and `b`
If either `a` or `b` is a wide ball we can compute a tighter enclosure
by using that `primitive_mul_x(t)` is increasing in `t`. That it is
increasing is an easy consequence of being the primitive function of a
positive integrand.
"""
function T0_primitive(u0::BHKdVAnsatz{Arb}, evaltype::Ball = Ball(); skip_div_u0 = false)
    #Enclosure of Arb((2 - u0.ϵ, 3)) computed in a way so that the
    #upper endpoint is exactly 2
    s1 = 2 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))
    #Enclosure of Arb((3 - u0.ϵ, 3)) computed in a way so that the
    #upper endpoint is exactly 3
    s2 = 3 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    return (x::Arb, a::Arb, b::Arb; to_endpoint = false) -> begin
        a < 1 < b || throw(ArgumentError("should have a < 1 < b"))
        b * x < π ||
            to_endpoint ||
            throw(ArgumentError("must have either b * x < π or integrate to endpoint"))

        # Check that the integrand is positive on the left endpoint.
        # It is enough to check this for a lower bound of x.
        let t = lbound(Arb, a), x = lbound(Arb, x), α = Arb((-1, -1 + u0.ϵ))
            Arblib.ispositive(
                clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) -
                2clausenc(x * t, -α),
            ) || throw(ArgumentError("integrand not positive at t = a = $a"))
        end

        primitive_mul_x(t) =
            let
                # Compute enclosure of all clausenc terms
                term1 = ArbExtras.enclosure_series(x) do x
                    (
                        clausencmzeta(x * (1 - t), s2) + clausencmzeta(x * (1 + t), s2) - 2clausencmzeta(x * t, s2)
                    ) / x
                end

                # Compute enclosure of all clausens terms
                term2 = ArbExtras.enclosure_series(x) do x
                    t * (
                        clausens(x * (1 - t), s1) - clausens(x * (1 + t), s1) +
                        2clausens(x * t, s1)
                    )
                end

                return term1 - term2
            end

        # primitive_mul_x(b) - primitive_mul_x(a)
        primitive_bma_mul_x(a, b) =
            let
                # Compute enclosure of all clausenc terms
                term1 = ArbExtras.enclosure_series(x) do x
                    (
                        (
                            clausencmzeta(x * (1 - b), s2) +
                            clausencmzeta(x * (1 + b), s2) -
                            2clausencmzeta(x * b, s2)
                        ) - (
                            clausencmzeta(x * (1 - a), s2) +
                            clausencmzeta(x * (1 + a), s2) -
                            2clausencmzeta(x * a, s2)
                        )
                    ) / x
                end

                # Compute enclosure of all clausens terms
                term2 = ArbExtras.enclosure_series(x) do x
                    (
                        b * (
                            clausens(x * (1 - b), s1) - clausens(x * (1 + b), s1) +
                            2clausens(x * b, s1)
                        ) -
                        a * (
                            clausens(x * (1 - a), s1) - clausens(x * (1 + a), s1) +
                            2clausens(x * a, s1)
                        )
                    )
                end

                return term1 - term2
            end

        if !to_endpoint
            @assert b * x < π # We have already checked this above

            # We integrate on [a, b]
            weight_enclosure = let α = Arb((-1, -1 + u0.ϵ)), t = Arb((a, b))
                t^(-u0.γ * (1 + α)) * log(u0.c + inv(x * t))
            end

            I_mul_x = if iswide(b)
                # Use that the integral is monotone in b
                union(
                    primitive_bma_mul_x(a, lbound(Arb, b)),
                    primitive_bma_mul_x(a, ubound(Arb, b)),
                )
            else
                primitive_bma_mul_x(a, b)
            end
        else
            # We integrate on [a, π / x]
            weight_enclosure = let α = Arb((-1, -1 + u0.ϵ)), t = Arb((a, π / x))
                t^(-u0.γ * (1 + α)) * log(u0.c + inv(x * t))
            end

            # Compute primitive_mul_x(a)
            primitive_a_mul_x = primitive_mul_x(a)

            # Enclosure of primitive_mul_x(π / x)
            primitive_pidivx_mul_x = let
                # Enclosure of π - x using that x <= π
                pimx = Arblib.nonnegative_part!(zero(x), π - x)

                2(clausenc(pimx, s2) + real(eta(Acb(s2)))) / x
            end

            I_mul_x = primitive_pidivx_mul_x - primitive_mul_x(a)
        end

        # Enclosure of factor in front of integral divided by x, excluding u0
        factor_div_x = inv(π * log(u0.c + inv(x)))

        res = factor_div_x * weight_enclosure * I_mul_x

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
