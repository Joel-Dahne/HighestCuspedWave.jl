"""
    T02(u0::BHKdVAnsatz; δ2, skip_div_u0)

Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral \$T_{0,2}\$ from the paper.

The interval of integration is `[x, π]`. Since the integrand is
singular at `y = x` we split the interval into two parts, `[x, a]` and
`[a, π]`. In principle we want to take `a = x + δ2`, however this
gives issues if `x` is a wide ball. Instead we take `a` to always be a
thin ball between `x` and `π`.

In general we take `a = ubound(x + δ2)`. However if `x` is wide then
it's beneficial to take a larger value than `δ2`, depending on the
radius of `x`. In practice we take `8radius(x)` in that case.

If `x` is very close to `π`, so that `a` would be larger than `π`,
then we use the asymptotic version for the whole interval `[x, π]`.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T02(u0::BHKdVAnsatz, evaltype::Ball; δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    f = T021(u0, evaltype, skip_div_u0 = true; δ2)
    g = T022(u0, evaltype, skip_div_u0 = true; δ2)

    return x -> begin
        x = convert(Arb, x)
        δ2 = max(δ2, 8Arblib.radius(Arb, x))
        a = Arblib.ubound(Arb, x + δ2)

        if !(a < π)
            return f(x, π)
        end

        res = f(x, a) + g(x, a)

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T02(u0::BHKdVAnsatz, ::Asymptotic; ϵ = Arb(2e-1))

Returns a function such that `T02(u0, Asymptotic())(x)` computes an
**upper bound** of the integral \$T_{0,2}\$ from the paper using an
evaluation strategy that works asymptotically as `x` goes to 0.

It precomputes the expansions of `u0` and for that reason a number `ϵ`
has to be given, the resulting expansion will be valid for all `x <
ϵ`. The value of `ϵ` has to be less than `1 // 2`.

First of all the change of coordinates `t = y / x` leaves us with
```
x / (π * u0(x) * log(u0.c + inv(x))) *
    ∫ abs(clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t * log(u0.c + inv(x * t)) dt
```
with the integration going from `1` to `π / x`. Further we can notice
that the expression inside the absolute value is positive on the whole
interval and the absolute value can hence be removed.

Next the factor
```
F = abs(inv(π) * log(x) / log(u0.c + inv(x)) * gamma(1 + α) * x^-α * (1 - x^p0) / u0(x))
```
is factored out. Except for the addition of `inv(π)` this is the same
as the factor `F1` in [`F0_nonzer`](@ref) and we compute an enclosure
following the same procedure as described there.

What we are left with computing is
```
abs(W(x) * I)
```
where `W(x) = x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(x))` and
`I` the same integral as above (but without the absolute value in the
integrand).

Now consider the expansions
```
clausenc(x * (t - 1), -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * (t - 1)^(-α - 1) + R(x * (t - 1))
clausenc(x * (1 + t), -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * (1 + t)^(-α - 1) + R(x * (1 + t))
clausenc(x * t, -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * t^(-α - 1) + R(x * t)
```
where the error term `R` contains one constant term and `R(x * (1 -
t)) + R(x * (1 + t)) - 2R(x * t)` behaves like `O(x^2)`. We can split
the integral into the two integrals
```
I₁ = gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(u0.c + inv(x * t)) dt
I₂ = ∫ (R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(u0.c + inv(x * t)) dt
```
satisfying `I = I₁ + I₂`. Furthermore we can split `log(u0.c + inv(x *
t))` as
```
log(u0.c + inv(x * t)) = log((u0.c * x * t + 1) / (x * t)) = log(1 + u0.c * x * t) - log(x) - log(t)
```
Which allows us to split the two above integrals into three integrals
each.
```
I₁₁ = -log(x) * gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
I₁₂ = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
I₁₃ = gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(1 + u0.c * x * t) dt
I₂₁ = -log(x) * ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t dt
I₂₂ = -∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t * log(t) dt
I₂₃ = ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t * log(1 + u0.c * x * t) dt
```

# Computation of `W(x) * I₁`

## `W(x) * I₁₁`
For the first integral we get
```
W(x) * I₁₁ = -inv(1 - x^p0) * sinpi(α / 2) * ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
```
where the integral is taken from `1` to `π / x`. The factor `sinpi(α /
2)` converges to `-1` and can be enclosed directly. A primitive
function for the integral is given by
```
((t - 1)^-α * (1 - α * t) - (t + 1)^-α * (1 + α * t) + 2α * t^(1 - α)) / (α * (α - 1))
```
and we hence get that the integral is
```
((π / x - 1)^-α * (1 - α * π / x) - (π / x + 1)^-α * (1 + α * π /x) + 2α * (π / x)^(1 - α)) / (α * (α - 1)) -
(2α - 2^-α * (1 + α)) / (α * (α - 1))
```
The part `(2α - 2^-α * (1 + α)) / (α * (α - 1))` doesn't depend on
`x`, for the other part we can rewrite it as
```
x^(α - 1) * ((π - x)^-α * (x - α * π) - (π + x)^-α * (x + α * π) + 2α * π^(1 - α)) / (α * (α - 1))
```
and furthermore as
```
x^(α - 1) * (
    x * ((π - x)^-α - (π + x)^-α) -
    α * π * ((π - x)^-α + (π + x)^-α) +
    2α * π^(1 - α)
) / (α * (α - 1))
```
Expanding at `x = 0` we have
```
x * ((π - x)^-α - (π + x)^-α) = 2α * π^(-1 - α) * x^2 + O(x^4)
α * π * ((π - x)^-α + (π + x)^-α) = 2α * π^(1 - α) + α^2 * (1 + α) * π^(-1 - α) * x^2 + O(x^4)
```
This gives us
```
x^(α - 1) * (
    2α * π^(-1 - α) * x^2 + O(x^4) -
    α^2 * (1 + α) * π^(-1 - α) * x^2 - O(x^4)
) / (α * (α - 1))
```
Which simplifies to
```
x^(1 + α) * π^(-1 - α) * (
    2 + O(x^2) -
    α * (1 + α) - O(x^2)
) / (α - 1)
```
Ignoring the `O(x^2)` terms for now we get the following expression
for the integral
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt =
x^(1 + α) * π^(-1 - α) * (2 - α * (1 + α)) / (α - 1) - (2α - 2^-α * (1 + α)) / (α * (α - 1))
```
- **TODO:** Handle the `O(x^2)` terms.
Putting everything together we get
```
W(x) * I₁₁ =
    - inv(1 - x^p0) * sinpi(α / 2) * (
        x^(1 + α) * π^(-1 - α) * (2 - α * (1 + α)) / (α - 1) -
        (2α - 2^-α * (1 + α)) / (α * (α - 1))
    ) =
    - sinpi(α / 2) * (
        x^(1 + α) * π^(-1 - α) * (2 - α * (1 + α)) / (α - 1) -
        (2α - 2^-α * (1 + α)) / (α * (α - 1))
    ) / (1 - x^p0)
```
As mentioned above the factor `sinpi(α / 2)` can be enclosed directly.
The remaining part converges to one, which is easily seen by plugging
in `x = 0` and `α = -1` directly, and is increasing in both `x` and
`α`
- **PROVE:** The it is increasing in `x` and `α`
- **TODO:** That the limit is 1 could be proved better.

## `W(x) * I₁₂`
For the second term we have
```
W(x) * I₁₂ = -inv(log(x)) * inv(1 - x^p0) * sinpi(α / 2) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
```
where the integral is taken from `1` to `π / x`. The factor `sinpi(α /
2)` converges to `-1` and can be enclosed directly. A primitive
function for the integral is given by
```
-((1 - t)^α * (t - 1)^-α * beta_inc(1 - α, 0, 1 - t)) / α +
(
    (α - 1) * (1 - t)^α * (t - 1)^-α * beta_inc(2 - α, 0, 1 - t) +
    t^-α * (
        2α * t +
        (α - 1)^2 * t^α * beta_inc(1 - α, 0, 1 + t) +
        (1 - α) * (
            α * t^α * beta_inc(2 - α, 0, 1 + t) +
            (
                -2α * t +
                (t - 1)^-α * t^α * (α * t - 1) +
                (t / (1 + t))^α * (1 + α * t)
            ) * log(t)
        )
    ) / α
) / (α - 1)^2
```
This comes from
```
FullSimplify[Integrate[((t - 1)^(-a - 1) + (t + 1)^(-a - 1) - 2 t^(-a - 1))*t*Log[t], t]]
```
This gives us that the integral from `1` to `π / x` is given by
```
LONG EXPRESSION
```
From which we get
```
W(x) * I₁₂ =
```
- **TODO:** Finish this. Some more notes are given in the code below,
  implementation of the primitive function and also other variants of
  the primitive function. Write down the expression for the integral
  and check it. Simplify it and expand in `x` if needed. Determine the
  limit and figure out how to bound it. The expression uses complex
  arithmetic, if possible we should probably avoid that.

For the third term we have
```
W(x) * I₁₃ = inv(log(x)) * inv(1 - x^p0) * sinpi(α / 2) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(1 + u0.c * x * t) dt
```
The factor `log(1 + u0.c * x * t)` inside the integral can be enclosed
using that it is zero for `x = 0` and increasing in both `x` and `t`.
We can factor out the enclosure from the integral and the integral we
are left with is the same as for `W(x) * I₁₁`. In the end we thus get
```
W(x) * I₁₃ = -inv(log(x)) * log(1 + u0.c * x * [0, 1]) * W(x) * I₁₁
```
where the minus sign is added to cancel the one from `W(x) * I₁₁`.

Now for the remaining three terms the first step is to handle the term
`R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)`. From the expansion of
the Clausen functions we get
```
R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t) = sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
For `I₂₁` we want to multiply this by `t` and integrate from `1` to `π
/ x`. Switching the integration and summation gives us
```
sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ∫ ((1 - t)^2m + (1 + t)^2m - 2t^2m) * t dt for m = 1:Inf
)
```
The integral taken from `1` to `π / x` can be determined to be
```
∫ ((1 - t)^2m + (1 + t)^2m - 2t^2m) * t dt =
    2(π / x)^2m * sum(
        binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for k = 0:m-1
    )
```
Inserting this into the main sum we get
```
sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * (2(π / x)^2m * sum(
        binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for k = 0:m-1
    )) for m = 1:Inf
) =
2sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * π^2m * sum(
        binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for k = 0:m-1
    ) for m = 1:Inf
)
```
Call this sum `Σ₂₁`.
- **TODO:** This sum seems to be mostly constant in both `α` and `x`.
  We could therefore attempt to compute a uniform enclosure of it. The
  interior sum converges to`(2m - 1) / 2` as `x -> 0`. But does so
  from above so we don't get an immediate bound for it.

For `I₂₂` we instead of want to multiply by `t * log(t)`, still
integrating from `1` to `π / x`. Switching the integration and
summation this gives us
```
sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m *
    ∫ ((1 - t)^2m + (1 + t)^2m - 2t^2m) * t *log(t) dt for m = 1:Inf
)
```
By rewriting `(1 - t)^2m + (1 + t)^2m - 2t^2m` as
```
2sum(binomial(2m, 2k) * t^2k for k = 0:m-1)
```
the primitive function of the integrand can be determined to be
```
2sum(binomial(2m, 2k) * t^(2k + 2) * ((2k + 2) * log(t) - 1) / 4(1 + k)^2 for k = 0:m-1)
```
This gives us that the integral from `1` to `π / x` is
```
2sum(binomial(2m, 2k) / 4(1 + k)^2 * ((π / x)^(2k + 2) * ((2k + 2) * log(π / x) - 1) + 1) for k = 0:m-1)
```
To make the sum bounded we factor out `(π / x)^2m * log(π / x)` to get
```
(π / x)^2m * log(π / x) * 2sum(
    binomial(2m, 2k) / 4(1 + k)^2 * ((x / π)^(2(m - k - 1)) * ((2k + 2) - inv(log(π / x))) + (x / π)^2m / log(π / x)) for k = 0:m-1
)
```
The sum is now bounded by `2m - 1`.
- **PROVE:** That the above sum is bounded by `2m - 1`.
We can get an upper bound of the main sum by taking the absolute value
of each term and using the `2m - 1` bound, giving us
```
sum(
    abs(zeta(-α - 2m)) / factorial(2m) * x^2m * (π / x)^2m * log(π / x) * (2m - 1) for m = 1:Inf
) =
log(π / x) * sum(abs(zeta(-α - 2m)) / factorial(2m) * π^2m * (2m - 1) for m = 1:Inf)
```
Call the sum (not including the log-factor) `Σ₂₂`.

With this we get for `I₂₁`
```
W(x) * I₂₁ = -x^(1 + α) / (gamma(1 + α) * (1 - x^p0)) *
    ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t dt =
    -x^(1 + α) / (gamma(1 + α) * (1 - x^p0)) * Σ₂₁
```
With an enclosure of `Σ₂₁` we are left computing an enclosure of
`-x^(1 + α) / (gamma(1 + α) * (1 - x^p0))`.
- **TODO:** Compute an enclosure of `-x^(1 + α) / (gamma(1 + α) * (1 -
  x^p0))`. It seems to be increasing in `x` and for fixed `x` it is
  increasing but bounded as `α -> -1`.

For `I₂₂` we have
```
W(x) * I₂₂ = x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(x)) *
    ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t * log(t) dt
```
The integral is bounded by `log(π / x) * Σ₂₂`, giving us
```
W(x) * I₂₂ = x^(1 + α) / (gamma(1 + α) * (1 - x^p0)) * log(π / x) / log(x) * Σ₂₁
```
We can compute an enclosure of `x^(1 + α) / (gamma(1 + α) * (1 - x^p0)
* log(x))` as above. An enclosure of `log(π / x) / log(x)` is also
easily computed by rewriting it as `log(π) / log(x) - 1`.

Finally for the last integral we have
```
W(x) * I₂₃ = x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(x)) *
    ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t * log(1 + u0.c * x * t) dt
```
We can compute an enclosure of `x^(1 + α) / (gamma(1 + α) * (1 - x^p0)
* log(x))` as above. Similarly to for `I₁₃` we can factor out an
enclosure of `log(1 + u0.c * x * t)` from the integral, leaving us with the
same integral as for `I₂₁`, that is
```
W(x) * I₂₃ = inv(log(x)) * x^(1 + α) / (gamma(1 + α) * (1 - x^p0)) *
    log(1 + u0.c * x * [0, 1]) * Σ₂₁
```
"""
function T02(u0::BHKdVAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), x * log(u0.c + inv(x)))
    end

    ϵ = convert(Arb, ϵ)
    @assert ϵ < 0.5

    # Setup for enclosing F

    # Compute the expansion of u0 and remove the leading term, which
    # is handled separately.
    u0_expansion = u0(ϵ, AsymptoticExpansion())
    delete!(u0_expansion, (1, 0, 0, 0, 0, 0, 0))

    # Ensure that the tail of the expansion of u0 is positive, so that
    # we can remove it from the denominator of F1 and still get an
    # upper bound.
    expansion_ispositive(u0, u0_expansion, ϵ) ||
        error("expansion of u0 not proven to be positive, this should not happen")

    F(x) =
        let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, xᵤ = ubound(Arb, x)
            # Note that inside this statement α refers to -1 + u0.ϵ
            # and p0 to the corresponding p0 value.

            # Enclosure of abs(log(x) / log(u0.c + inv(x))), either by
            # direct evaluation or using monotonicity.
            F11 = if iszero(x)
                one(Arb)
            elseif Arblib.contains_zero(x)
                abs(Arb((-1, log(xᵤ) / log(u0.c + inv(xᵤ)))))
            else
                abs(log(x) / log(u0.c + inv(x)))
            end

            # Enclose F12
            F121_lower = -Arb(π)^2 / 2
            F122_upper = gamma(1 + α) / finda0(α)
            F121 = Arb((F121_lower, F122_upper))

            c(a) = gamma(a) * sinpi((1 - a) / 2)

            # Upper and lower bound of
            # (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
            F122_lower = (1 - xᵤ^p0) / (1 - c(α - p0) / c(α) * xᵤ^p0)
            F122_upper = one(Arb)
            # Combine upper and lower bound and multiply with
            # enclosure of inv(c(α)) to get an upper bound for F122.
            F122 = inv(Arb((c(α), -Arb(π) / 2))) * Arb((F122_lower, F122_upper))

            F12 = F121 * F122

            Arblib.ispositive(F12) ||
                error("leading term of u0 is not positive, this should not happen")

            inv(π) * F11 * F12
        end

    return x::Arb -> begin
        @assert x <= ϵ

        # Enclosure of sinpi(α / 2)
        factor_I₁ = sinpi(Arb((-1, -1 + u0.ϵ)) / 2)

        # Enclosure of inv(log(x))
        invlogx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            xᵤ = ubound(Arb, x)
            Arb((inv(log(xᵤ)), 0))
        else
            inv(log(x))
        end

        WxI₁₁ = begin
            # Enclosure of inv(1 - x^p0) * ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
            # FIXME: This doesn't enclose the full integral yet, there
            # are some O(x^2) terms left to handle.
            integral_I₁₁ =
                let xᵤ = ubound(Arb, x), α = -1 + u0.ϵ, p0 = (1 + α) + (1 + α)^2 / 2

                    # We use that it is 1 for α = -1, x = 0 and
                    # increasing in both variables. Notice that we
                    # always enclose it for x on the interval [0, xᵤ]
                    # even if x doesn't touch zero. This could
                    # possibly be improved but might not be needed.
                    lower = one(x)
                    upper =
                        (
                            xᵤ^(1 + α) * π^(-1 - α) * (2 - α * (1 + α)) / (α - 1) -
                            (2α - 2^-α * (1 + α)) / (α * (α - 1))
                        ) / (1 - xᵤ^p0)

                    Arb((lower, upper))
                end

            -factor_I₁ * integral_I₁₁
        end


        WxI₁₂ = let α = Arb(-1 + u0.ϵ) # FIXME: For testing TODO:
            # Finish this, the primitive functions seems to be correct
            # at least.

            primitive_I₁₂(t) =
                -((1 - t)^α * (t - 1)^-α * beta_inc(1 - α, zero(t), 1 - t)) / α +
                (
                    (α - 1) * (1 - t)^α * (t - 1)^-α * beta_inc(2 - α, zero(t), 1 - t) +
                    t^-α * (
                        2α * t +
                        (α - 1)^2 * t^α * beta_inc(1 - α, zero(t), 1 + t) +
                        (1 - α) * (
                            α * t^α * beta_inc(2 - α, zero(t), 1 + t) +
                            (
                                -2α * t +
                                (t - 1)^-α * t^α * (α * t - 1) +
                                (t / (1 + t))^α * (1 + α * t)
                            ) * log(t)
                        )
                    ) / α
                ) / (α - 1)^2

            # Primitive function of (t - 1)^(-α - 1) * t * log(t)
            primitive_I12_1(t) =
                t^-α * hypgeom_2f1(α, α, 1 + α, inv(t)) / ((α - 1) * α^2) -
                (t - 1)^(1 - α) / (α - 1)^2 -
                (t - 1)^-α * (α * t - 1) * log(t) / ((α - 1) * α)

            primitive_I12_1_one = hypgeom_2f1(α, α, 1 + α, Arb(1)) / ((α - 1) * α^2)

            # Primitive function of (t + 1)^(-α - 1) * t * log(t)
            primitive_I12_2(t) =
                -(t + 1)^(1 - α) * hypgeom_2f1(Acb(1), Acb(1 - α), Acb(2 - α), Acb(1 + t)) /
                ((α - 1) * α) +
                (t + 1)^(2 - α) * hypgeom_2f1(Acb(1), Acb(2 - α), Acb(3 - α), Acb(1 + t)) /
                ((α - 2) * (α - 1)) -
                (t + 1)^-α * (1 + α * t) * log(t) / ((α - 1) * α)

            primitive_I12_3(t) = 2t^(1 - α) * (1 + (α - 1) * log(t)) / (α - 1)^2

            @show a = primitive_I12_1(Arb(π / x)) - primitive_I12_1_one
            @show b = primitive_I12_2(Arb(π / x)) - primitive_I12_2(Arb(1))
            @show c = primitive_I12_3(Arb(π / x)) - primitive_I12_3(Arb(1))
            integral_I₁₂ = a + b + c

            term = let p0 = (1 + α) + (1 + α)^2 / 2
                integral_I₁₂ / ((1 - x^p0) * log(x))
            end

            @show term

            factor_I₁ * term
        end

        WxI₁₃ = begin
            # Enclosure of log(1 + u0.c * x * t) for t ∈ [0, 1]
            factor_I₁₂ = Arb((0, log1p(u0.c * x)))

            invlogx * factor_I₁₂ * WxI₁₁
        end

        WxI₁ = WxI₁₁ + WxI₁₂ + WxI₁₃

        # Enclosure of -x^(1 + α) / (gamma(1 + α) * (1 - x^p0))
        # FIXME: This currently computes with `α = -1 + u0.ϵ`
        # which doesn't given an upper bound, though it gives a
        # good estimate.
        factor_I₂ = let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2
            if Arblib.contains_zero(x)
                xᵤ = ubound(Arb, x)
                Arb((-xᵤ^(1 + α) / (gamma(1 + α) * (1 - xᵤ^p0)), 0))
            else
                -x^(1 + α) / (gamma(1 + α) * (1 - x^p0))
            end
        end

        WxI₂₁ = begin
            # Enclosure of
            #2sum(
            #    zeta(-α - 2m) * (-1)^m / factorial(2m) * π^2m * sum(
            #        binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for k = 0:m-1
            #    ) for m = 1:Inf
            #)
            # FIXME: This only computes the first 10 terms in the sum.
            # We need to bound the tail.
            Σ₂₁ = zero(x)
            let π = Arb(π), α = Arb((-1, -1 + u0.ϵ))
                for m = 1:10
                    term =
                        zeta(-α - 2m) * (-1)^m / factorial(2m) *
                        π^2m *
                        sum(
                            binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for
                            k = 0:m-1
                        )
                    Σ₂₁ += term
                end
            end
            Σ₂₁ *= 2

            factor_I₂ * Σ₂₁
        end

        WxI₂₂ = begin
            # Enclosure of
            # sum(abs(zeta(-α - 2m)) / factorial(2m) * π^2m * (2m - 1) for m = 1:Inf)
            # FIXME: This only computes the first 10 terms in the sum,
            # we need to bound the tail.
            Σ₂₂ = zero(x)
            let π = Arb(π), α = Arb((-1, -1 + u0.ϵ))
                for m = 1:10
                    Σ₂₂ += abs(zeta(-α - 2m)) / factorial(2m) * π^2m * (2m - 1)
                end
            end

            # Enclosure of log(π / x) / log(x) = log(π) / log(x) - 1
            factor_I₂₂ = log(Arb(π)) * invlogx - 1

            factor_I₂ * factor_I₂₂ * Σ₂₂
        end

        WxI₂₃ = begin
            # Enclosure of log(1 + u0.c * x * t) for t ∈ [0, 1]
            factor_I₂₃ = Arb((0, log1p(u0.c * x)))

            invlogx * factor_I₂ * factor_I₂₃ * Σ₂₁
        end

        WxI₂ = WxI₂₁ + WxI₂₂ + WxI₂₃

        #@show WxI₁₁ WxI₁₂ WxI₁₃ WxI₂₁ WxI₂₂ WxI₂₃

        return F(x) * (WxI₁ + WxI₂)
    end
end

"""
    T021(u0::BHKdVAnsatz)

Computes the (not yet existing) integral \$T_{0,2,1}\$ from the paper.

The interval of integration is `[x, a]`. Both `x` and `a` are assumed
to be less than or equal to `π`, if they are balls which overlap `π`
anything above `π` will be ignored.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval. We can also notice that the value inside the absolute value
is positive so we can remove the absolute value.

We are left with integrating the three Clausen terms
1. `clausenc(x - y, -α)`
2. `clausenc(x + y, -α)`
3. `2clausenc(y, -α)`
We have that the primitive functions for the three terms are given by
1. `-clausens(x - y, 1 - α)`
2. `clausens(x + y, 1 - α)`
3. `2clausens(y, 1 - α)`
Hence the integral from `x` to `a` is
```
(-clausens(x - a, 1 - α) + clausens(x + a, 1 - α) - 2clausens(a, 1 - α)) -
(-clausens(0, 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α))
```
Since `1 - α > 1` we have
`clausens(0, 1 - α) = 0`. If we also reorder the terms to more clearly
see which ones gives cancellations we get
```
-clausens(x - a, 1 - α) +
(clausens(x + a, 1 - α) - clausens(2x, 1 - α)) -
2(clausens(a, 1 - α) - clausens(x, 1 - α))
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
  for `clausens(x + a, 1 - α) - clausens(2x, 1 - α)` and `clausens(a,
  1 - α) - clausens(x, 1 - α)`. Though this might not be needed.
"""
function T021(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)
        δ = a - x

        interval = union(x, a)

        weight_factor = u0.w(interval)

        # s = 1 - α
        s = Arb((2 - u0.ϵ, 2))

        integral = -clausens(x - a, s) - 2(clausens(a, s) - clausens(x, s))
        if Arblib.overlaps(x, Arb(π))
            # FIXME: This assumes that clausens is monotonic on the
            # interval. In practice this is true for small enough
            # argument. But it is not true in general.

            # Compute an enclosure of clausens(x + a, s) on the
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
                clausens(abs_ubound(Arb, x * (2 - δ2) - 2Arb(π)), s),
            )
            integral -= term
        else
            integral += clausens(x + a, s) - clausens(2x, s)
        end
        integral *= weight_factor

        res = integral / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T022(u0::BHKdVAnsatz)

Computes the (not yet existing) integral \$T_{0,2,2}\$ from the paper.

The interval of integration is given by `[a, π]`. In practice `a`
should be a thin ball to not give problems with the integration.

This is done by directly computing the integral with the integrator in
Arb.

Notice that the expression inside the absolute value is always
positive, so we can remove the absolute value.

**FIXME:** This currently assumes that
```
clausenc(y - x, s) + clausenc(y + x, s) - 2clausenc(y, s)
```
and its derivatives up to the fourth one are monotonic in `s`. This is
true for most of the interval but there are some points where it
doesn't hold. One solution would be to prove that this only happens at
some places, isolate them and handle them separately. This might be
tedious though since we would have to do it for all required
derivatives. The point where it happens does depend on `x`.
"""
function T022(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * log(u0.c + inv(x)))
    end

    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)

        # FIXME: Currently we assume monotonicity in s, including for
        # all derivatives.
        s_l = 1 - u0.ϵ
        s_u = one(Arb)
        integrand(y) = begin
            term_l = clausenc(y - x, s_l) + clausenc(y + x, s_l) - 2clausenc(y, s_l)
            term_u = clausenc(y - x, s_u) + clausenc(y + x, s_u) - 2clausenc(y, s_u)

            if y isa ArbSeries
                coefficients = union.(Arblib.coeffs(term_l), Arblib.coeffs(term_u))
                term_union = ArbSeries(coefficients)
            else
                term_union = union(term_l, term_u)
            end

            return term_union * y * log(u0.c + inv(y))
        end

        res = ArbExtras.integrate(integrand, a, Arb(π), atol = 1e-5, rtol = 1e-5)

        res = res / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
