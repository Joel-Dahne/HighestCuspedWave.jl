"""
    T01(u0::BHKdVAnsatz, ::Ball; δ1, δ2, skip_div_u0)

Returns a function such that `T01(u0, Ball(); δ1, δ2)(x)` computes the
integral \$T_{0,1}\$ from the paper.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T01(
    u0::BHKdVAnsatz,
    evaltype::Ball;
    δ0::Arb = Arb(1e-5),
    δ1::Arb = Arb(1e-5),
    skip_div_u0 = false,
)
    f = T011(u0, evaltype, skip_div_u0 = true; δ0)
    g = T012(u0, evaltype, skip_div_u0 = true; δ0, δ1)
    h = T013(u0, evaltype, skip_div_u0 = true; δ1)

    if skip_div_u0
        return x -> f(x) + g(x) + h(x)
    else
        return x -> (f(x) + g(x) + h(x)) / u0(x)
    end
end

"""
    T01(u0::BHKdVAnsatz, ::Asymptotic)

Returns a function such that `T01(u0, Asymptotic())(x)` computes an
**upper bound** of the integral \$T_{0,1}\$ from the paper using an
evaluation strategy that works asymptotically as `x` goes to `0`.

It precomputes the expansions of `u0` and for that reason a number `ϵ`
has to be given, the resulting expansion will be valid for all `x <
ϵ`. The value of `ϵ` has to be less than `1 // 2`.

First of all the change of coordinates `t = y / x` leaves us with
```
x / (π * u0(x) * log(u0.c + inv(x))) *
    ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t * log(u0.c + inv(x * t)) dt
```
with the integration going from `0` to `1`.

Next the factor
```
F(x) = x^-α * log(x) / (π * u0(x) * log(u0.c + inv(x)))
```
is factored out from the whole expression and multiplied back in the
end. Notice that this factor is bounded in `x`. We can bound `log(x) /
log(u0.c + inv(x))` by using that it is `-1` at `x = 0` and
increasing.
- **PROVE:** That `log(x) / log(u0.c + inv(x))` is `-1` at `x = 0` and
    increasing.
To bound `x^-α / u0(x)` we compute the asymptotic expansion of
`u0(x)`, multiply that by `x^α` and evaluate.
- **TODO:** Deal with the fact that `u0(x) * x^α` blows up as `x ->
  0`. It's not a problem since we divide by it, but it needs to be
  handled. For now we assume monotonicity in `x` in this case and
  compute it at the upper bound for `x`

What we are left with computing is
```
W(x) * I
```
where `W(x) = x^(1 + α) / log(x)` and `I` the same integral as above.

Now consider the expansions
```
clausenc(x * (1 - t), -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * (1 - t)^(-α - 1) + R(x * (1 - t))
clausenc(x * (1 + t), -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * (1 + t)^(-α - 1) + R(x * (1 + t))
clausenc(x * t, -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * t^(-α - 1) + R(x * t)
```
where the error term `R` contains one constant term and `R(x * (1 -
t)) + R(x * (1 + t)) - 2R(x * t)` behaves like `O(x^2)`. Since we are
only looking for an upper bound we can split the absolute value and we
get the two integrals
```
I₁ = abs(gamma(1 + α) * sinpi(α / 2)) * x^(-α - 1) *
    ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(u0.c + inv(x * t)) dt
I₂ = ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(u0.c + inv(x * t)) dt
```
satisfying `I <= I₁ + I₂`. Furthermore we can split `log(u0.c + inv(x *
t))` as
```
log(u0.c + inv(x * t)) = log((u0.c * x * t + 1) / (x * t)) = log(1 + u0.c * x * t) - log(x) - log(t)
```
Which allows us to split the two above integrals into three integrals
each.
```
I₁₁ = -log(x) * abs(gamma(1 + α) * sinpi(α / 2)) * x^(-α - 1) *
    ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
I₁₂ = -abs(gamma(1 + α) * sinpi(α / 2)) * x^(-α - 1) *
    ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
I₁₃ = abs(gamma(1 + α) * sinpi(α / 2)) * x^(-α - 1) *
    ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(1 + u0.c * x * t) dt
I₂₁ = -log(x) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t dt
I₂₂ = -∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(t) dt
I₂₃ = ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(1 + u0.c * x * t) dt
```
Since we are interested in an upper bound of the absolute value of
their sum it is enough to upper bound the absolute value of each
integral separately.

# Handling `I₁`
Focusing on the first three integrals the first step is to understand
the behaviour of
```
gamma(1 + α) * sinpi(α / 2) * ((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))
```
The factor `sinpi(α / 2)` converges to `-1` and can be enclosed
directly, we can hence focus our attention on
```
gamma(1 + α) * ((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))
```
As `α -> -1` this converges to
```
-log(1 - t) - log(1 + t) + 2log(t)
```
- **TODO:** Finish this.

For the first integral we have
```
abs(W(x) * I₁₁) = abs(
    sinpi(α / 2) *
    ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
)
```
As mentioned `sinpi(α / 2)` can be enclosed directly. The integral
doesn't depend on `x` and can be enclosed.
- **TODO:** Enclose the integral
  `∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt`.

For the second integral we have
```
abs(W(x) * I₁₂) = abs(
    inv(log(x)) * sinpi(α / 2) *
    ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
)
```
where again the integral doesn't depend on `x` and can be enclosed.
The factor `sinpi(α / 2)` can be enclosed as above and `-inv(log(x))`
can be enclosed using that it is zero at `x = 0` and increasing.
- **TODO:** Enclose the integral
  `∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt`.

For the third integral we have
```
abs(W(x) * I₁₃) = abs(
    inv(log(x)) * sinpi(α / 2) *
    ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(1 + u0.c * x * t) dt
)
```
The factor outside the integral can be dealt with in the same way as
for `abs(W(x) * I₁₂)` The factor `log(1 + u0.c * x * t)` inside the
integral can be enclosed using that it is zero for `x = 0` and
increasing in both `x` and `t`. We can factor out the enclosure from
the integral and the integral we are left with is the same as for
`W(x) * I₁₁`.

# Handling `I₂`
Now for the remaining three terms the first step is to handle the term
`R₁(x) + R₂(x) - 2R₃(x)`. From the expansion of the Clausen functions
we get
```
R₁(x) + R₂(x) - 2R₃(x) = sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
Since we are interested in an upper bound of absolute value we get
```
abs(R₁(x) + R₂(x) - 2R₃(x)) <= x^2 * sum(
    abs(zeta(-α - 2m)) / factorial(2m) * abs((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
) <= x^2 * 2sum(abs(zeta(-α - 2m)) * 2^2m / factorial(2m) for m = 1:Inf)
```
where we have used that
```
abs((1 - t)^2m + (1 + t)^2m - 2t^2m) <= 2^2m + 3 <= 2^(2m + 1)
```
for `m >= 1`. The final sum doesn't depend on `x` and can be bounded,
call the bound `C`, i.e.,
```
c = 2sum(abs(zeta(-α - 2m)) * 2^2m / factorial(2m) for m = 1:Inf)
```
- **TODO:** Bound the tail of the sum
  `2sum(abs(zeta(-α - 2m)) * 2^2m / factorial(2m) for m = 1:Inf)`.

Now for the integrals we get
```
abs(W(x) * I₂₁) = abs(x^(1 + α) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t dt) <=
    x^(1 + α) * x^2 * C * ∫ t dt =
    x^(3 + α) * C / 2  <=
    x^2 * C / 2
```
```
abs(W(x) * I₂₂) = abs(x^(1 + α) / log(x) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(t) dt) <=
    abs(x^(1 + α) / log(x) * x^2 * C * ∫ t * log(t) dt) =
    abs(x^(3 + α) / log(x) * C * ∫ t * log(t) dt) <=
    abs(x^2 / log(x) * C * ∫ t * log(t) dt) =
    abs(x^2 / log(x) * C / 4)
```
```
abs(W(x) * I₂₃) = abs(x^(1 + α) / log(x) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(1 + u0.c * x * t) dt) <=
    abs(x^(1 + α) / log(x) * x^2 * C * ∫ t * log(1 + u0.c * x * t) dt) <=
    abs(x^(3 + α) / log(x) * C * log(1 + u0.c * x) / 2) <=
    abs(x^3 / log(x) * C * log(1 + u0.c * x) / 2)
```
where in the last step we have used that `log(1 + u0.c * x * t) <= log(1 + u0.c * x)`.
- **TODO:** Take a look at the sign for these three integrals. We
  should probably just consider the absolute value instead.

"""
function T01(u0::BHKdVAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), x * log(u0.c + inv(x)))
    end

    ϵ = convert(Arb, ϵ)
    @assert ϵ < 0.5

    u0_expansion = u0(ϵ, AsymptoticExpansion())
    u0_expansion_mul_xα = empty(u0_expansion)
    for ((p, q, i, j, k, l, m), value) in u0_expansion
        u0_expansion_mul_xα[(p, q, i + 1, j, k, l, m)] = value
    end

    # Enclosure of the factor
    # F(x) = x^-α * log(x) / (π * u0(x) * log(u0.c + inv(x)))
    F(x) = begin
        # Enclosure of log(x) / log(u0.c + inv(x))
        log_factor = if iszero(x)
            -one(x)
        elseif Arblib.contains_zero(x)
            # Use monotonicity, it is -1 at x = 0 and increasing
            xᵤ = ubound(Arb, x)
            Arb((-1, log(xᵤ) / log(u0.c + inv(xᵤ))))
        else
            log(x) / log(u0.c + inv(x))
        end

        # Enclosure of inv(u0(x) * x^α)
        # TODO: This behaves like inv(log(x)) and is hence increasing
        # in x. In case x overlaps with zero we therefore compute an
        # enclosure by evaluation at the upper bound for x. This is
        # not proved to be correct though
        invu0 = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            Arb((0, inv(eval_expansion(u0, u0_expansion_mul_xα, ubound(Arb, x)))))
        else
            inv(eval_expansion(u0, u0_expansion_mul_xα, x))
        end

        log_factor * invu0 / π
    end

    return x::Arb -> begin
        @assert x <= ϵ

        # Enclosure of inv(log(x))
        invlogx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            Arb((inv(log(ubound(Arb, x))), 0))
        else
            inv(log(x))
        end

        # Handle I₁

        # Enclosure of sinpi(α / 2)
        factor_I₁ = let α = Arb((-1, -1 + u0.ϵ))
            sinpi(α / 2)
        end

        # abs(W(x) * I₁₁)
        WxI₁₁ = begin
            # Enclosure of
            # ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
            # FIXME: Rigorously enclose this
            integral_I₁₁ = Arb("0.6931477765876257")

            abs(factor_I₁ * integral_I₁₁)
        end

        WxI₁₂ = begin
            # Enclosure of
            # ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
            # FIXME: Rigorously enclose this
            integral_I₁₂ = Arb("-0.46668696508451124")

            abs(invlogx * factor_I₁ * integral_I₁₂)
        end

        WxI₁₃ = begin
            # Enclosure of log(1 + u0.cx * t) for t ∈ [0, 1]
            factor_I₁₂ = Arb((0, log1p(u0.c * x)))

            # The integral is the same as for WxI₁₁
            integral_I₁₃ = integral_I₁₁

            abs(invlogx * factor_I₁ * factor_I₁₂ * integral_I₁₃)
        end

        WxI₁ = WxI₁₁ + WxI₁₂ + WxI₁₃

        # Handle I₂

        # Bound of 2sum((abs(zeta(-α - 2m))) * 2^2m / factorial(2m) for m = 1:Inf)
        # FIXME: Bound the tail of this sum.
        C = let α = Arb((-1, -1 + u0.ϵ))
            2sum((abs(zeta(-α - 2m))) * 2^2m / factorial(2m) for m = 1:10)
        end

        # Enclosure of inv(abs(gamma(1 + α) * sinpi(α / 2)))

        WxI₂₁ = x^2 * C / 2

        WxI₂₂ = abs(x^2 * invlogx * C / 4)

        WxI₂₃ = abs(x^2 * invlogx * C * log1p(u0.c * x) / 2)

        WxI₂ = WxI₂₁ + WxI₂₂ + WxI₂₃

        # Combine the above

        WxI = WxI₁ + WxI₂

        res = abs(F(x)) * WxI

        return res
    end
end

"""
    T011(u0::BHKdVAnsatz; δ0)

Computes the integral \$T_{0,1,1}\$ from the paper.

It uses the fact that the integrand is strictly increasing on the
interval `[0, 0.05]` for every value of `x` and 0 at `x = 0`. This
allows us to enclose the integrand on the interval which then easily
gives an enclosure of the integral by multiplying with the size of the
interval.

- **PROVE**: That the integrand indeed is increasing on the said
  interval.
"""
function T011(u0::BHKdVAnsatz, ::Ball = Ball(); δ0::Arb = Arb(1e-5), skip_div_u0 = false)
    δ0 < 0.05 || Throw(ArgumentError("δ0 must be less than 0.05, got $δ0"))

    #Enclosure of Arb((-1, -1 + u0.ϵ)) computed in a way so that
    #the lower endpoint is exactly -1
    α = -1 + Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    return x::Arb -> begin
        integrand(t) =
            abs(
                clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) -
                2clausenc(x * t, -α),
            ) *
            t *
            u0.wdivx(x * t)

        integral = δ0 * Arb((0, integrand(δ0)))

        res = integral * x / (π * u0.wdivx(x))
        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T012(u0::BHKdVAnsatz; δ0, δ1)

Returns a function such that `T012(u0; δ0, δ1)(x; tol)` computes the
integral \$T_{0,1,2}\$ from the paper using the prescribed tolerance
in the integration.
"""
function T012(
    u0::BHKdVAnsatz,
    ::Ball = Ball();
    δ0::Arb = Arb(1e-5),
    δ1::Arb = Arb(1e-5),
    skip_div_u0 = false,
)
    # Enclosure of -α so that the upper bound is exactly 1
    mα = 1 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    # Integration limits
    a = δ0
    b = 1 - δ1

    return (x::Arb; tol = Arb(1e-5)) -> begin
        integrand(t) = begin
            xt = x * t

            term = abs(
                clausenc(x * (1 - t), mα) + clausenc(x * (1 + t), mα) - 2clausenc(xt, mα),
            )

            return term * t * u0.wdivx(xt)
        end

        res = ArbExtras.integrate(integrand, a, b, atol = tol, rtol = tol)

        res *= x / (π * u0.wdivx(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T013(u0::BHKdVAnsatz; δ1)

Computes the integral \$T_{0,1,3}\$ from the paper.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval. We can also notice that the value inside the absolute value
is negative so we can remove the absolute value by putting a minus
sign, which we can bake in to the weight factor.

We are left with integrating the three Clausen terms
1. `clausenc(x * (1 - t), -α)`
2. `clausenc(x * (1 + t), -α)`
3. `2clausenc(x * t, -α)`
We have that the primitive functions for the three terms are given by
1. `-clausens(x * (1 - t), 1 - α) / x`
2. `clausens(x * (1 + t), 1 - α) / x`
3. `2clausens(x * t, 1 - α) / x`
Hence the integral from `1 - δ1` to `1` is
```
inv(x) * (
    (-clausens(0, 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α)) -
    (-clausens(x * δ1, 1 - α) + clausens(x * (2 - δ1), 1 - α) - 2clausens(x * (1 - δ1), 1 - α))
)
```
The multiplication by `inv(x)` can be cancelled by the multiplication
by `x` that is outside of the integral. Since `1 - α > 1` we have
`clausens(0, 1 - α) = 0`. If we also reorder the terms to more clearly
see which ones gives cancellations we get
```
clausens(x * δ1, 1 - α) +
(clausens(2x, 1 - α) - clausens(x * (2 - δ1), 1 - α)) -
2(clausens(x, 1 - α) - clausens(x * (1 - δ1), 1 - α))
```

In the case that `x` overlaps with `π` we get issues when evaluating
`(clausens(2x, 1 - α)` since it doesn't have a good implementation in
that case. We could subtract `2π` from the argument since it is `2π`
periodic, but that doesn't solve the issue since `clausens` currently
doesn't support evaluation on intervals containing zero.
- **FIXME:** Currently we do this by assuming that `clausens` is
  monotonic in `s`. In practice this is true for small enough
  arguments but not in general. If `x` is sufficiently close to `π`
  this will thus give a correct result, but it is not rigorously
  proved.

- **TODO:** Could improve enclosures by better handling cancellations
  for `clausens(2x, 1 - α) - clausens(x * (2 - δ1), 1 - α)` and
  `clausens(x, 1 - α) - clausens(x * (1 - δ1), 1 - α)`. Though this
  might not be needed.
"""
function T013(u0::BHKdVAnsatz, ::Ball = Ball(); δ1::Arb = Arb(1e-5), skip_div_u0 = false)
    return x -> begin
        x = convert(Arb, x)

        weight_factor = let t = Arb((1 - δ1, 1))
            -t * u0.wdivx(x * t)
        end

        # s = 1 - α
        s = Arb((2 - u0.ϵ, 2))

        integral = clausens(x * δ1, s) - 2(clausens(x, s) - clausens(x * (1 - δ1), s))
        if Arblib.overlaps(x, Arb(π))
            # FIXME: This assumes that clausens is monotonic on the
            # interval. In practice this is true for small enough
            # argument. But it is not true in general.

            # Compute an enclosure of clausens(2(x - π), s) on the
            # symmetric interval [-abs(2(x - π)), abs(2(x - π))] using
            # the oddness and assuming that the maximum is attained at
            # the endpoint.
            term = Arblib.add_error!(zero(x), clausens(abs_ubound(Arb, 2(x - Arb(π))), s))
            integral += term

            # Compute an enclosure of clausens(x * (2 - δ1) - 2π, s)
            # on the symmetric interval [-abs(x * (2 - δ1) - 2π),
            # abs(x * (2 - δ1) - 2π)] using the oddness and assuming
            # the maximum is attained at the endpoint.
            term = Arblib.add_error!(
                zero(x),
                clausens(abs_ubound(Arb, x * (2 - δ1) - 2Arb(π)), s),
            )
            integral -= term
        else
            integral += clausens(2x, s) - clausens(x * (2 - δ1), s)
        end
        integral *= weight_factor

        res = integral / (π * u0.wdivx(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
