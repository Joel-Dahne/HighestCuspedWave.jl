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
F(x) = abs(inv(π) * log(x) / log(u0.c + inv(x)) * gamma(1 + α) * x^-α * (1 - x^p0) / u0(x))
```
is factored out. Except for the addition of `inv(π)` this is the same
as the factor `F1` in [`F0`](@ref) and we compute an upper bound
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

# Handling `I₁`
For `I₁` there are cancellations between the two terms `I₁₁` and `I₁₂`
and we therefore have to keep their sign. For `I₁₃` it is enough to
bound the absolute value.

## Handling `I₁₁`
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
`α`.
- **PROVE:** The it is increasing in `x` and `α`
- **TODO:** That the limit is 1 could be proved better.

## Handling `I₁₂`
For the second integral we have
```
W(x) * I₁₂ = -inv(log(x)) * inv(1 - x^p0) * sinpi(α / 2) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
```
where the integral is taken from `1` to `π / x`. The factor `sinpi(α /
2)` converges to `-1` and can be enclosed directly.

We can expand
```
(t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
```
in a Taylor around `α = -1` by writing it as
```
(t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1) =
    exp(-log(t - 1) * (α + 1)) + exp(-log(t + 1) * (α + 1)) - 2exp(-log(t) * (α + 1)) =
    sum((-1)^n * (α + 1)^n * (log(t - 1)^n + log(t + 1)^n - 2log(t)^n) / factorial(n) for n = 1:Inf)
```
Putting this inside the integral and switching the integral and the
sum we get
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt =
    sum((-1)^n * (α + 1)^n * ∫ (log(t - 1)^n + log(t + 1)^n - 2log(t)) * t * log(t) dt for n = 1:Inf)
```
where again the integral is taken from `1` to `π / x`. Now that the
integral no longer depend on `α` we can expand it around `x = 0`.

The remaining part is not finished, we here give the idea of what
should happen. We would like to show that the integral behaves like
```
∫ (log(t - 1)^n + log(t + 1)^n - 2log(t)) * t * log(t) dt = (-1)^n * n / (n + 1) * log(x)^(n + 1) + O(log(x)^n)
```
If we ignore the `O(log(x)^n)` term we would then get the sum
- **TODO:** Prove the above asymptotics and handle the `O(log(x)^n)` term.
```
sum((-1)^n * (α + 1)^n * (-1)^n * n / (n + 1) * log(x)^(n + 1) for n = 1:Inf) =
    sum((α + 1)^n * n / (n + 1) * log(x)^(n + 1) for n = 1:Inf)
```
This would give us
```
W(x) * I₁₂ =
    -inv(log(x)) * inv(1 - x^p0) * sinpi(α / 2) *
        sum((α + 1)^n * n / (n + 1) * log(x)^(n + 1) for n = 1:Inf) =
    -inv(1 - x^p0) * sinpi(α / 2) *
        sum((α + 1)^n * n / (n + 1) * log(x)^n for n = 1:Inf) =
```
The sum can be computed explicitly to be
```
x^(α + 1) + (1 - x^(α + 1)) / ((α + 1) * log(x))
```
Inserting this gives us
```
-inv(1 - x^p0) * sinpi(α / 2) * sum((α + 1)^n * log(x)^n for n = 1:Inf) =
    -inv(1 - x^p0) * sinpi(α / 2) * (x^(α + 1) + (1 - x^(α + 1)) / ((α + 1) * log(x))) =
    = -sinpi(α / 2) * (x^(α + 1) + (1 - x^(α + 1)) / ((α + 1) * log(x))) / (1 - x^p0)
```
- **TODO:** Finish this. Unfortunately the above seems to go to zero
  as `x` goes to `0`, whereas we would need to go to `-1 / 2`. Could
  this be because we have neglected the `O(log(x)^n)` terms in the
  integrals inside the sum? Or is there a chance that it actually goes
  to zero? In that case our whole approach fails :/

## Handling `I₁₃`
For the third integral we only need to bound the absolute value and we
have
```
abs(W(x) * I₁₃) = abs(
    inv(log(x)) * inv(1 - x^p0) * sinpi(α / 2) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(1 + u0.c * x * t) dt
)
```
The factor `log(1 + u0.c * x * t)` inside the integral can be enclosed
using that it is zero for `x = 0` and increasing in both `x` and `t`.
We can factor out the enclosure from the integral and the integral we
are left with is the same as for `W(x) * I₁₁`. In the end we thus get
```
abs(W(x) * I₁₃) = abs(inv(log(x)) * log(1 + u0.c * x * [0, 1]) * W(x) * I₁₁)
```

# Handling `I₂`
For the three integrals in `I₂` there are no cancellations we have to
keep track of and it is therefore enough to bound the absolute value
of each term separately.

In this case the only important part of `W(x)` is the `log(x)` factor
in the denominator, the rest we can treat by itself.

We also need to handle treat the term `R(x * (t - 1)) + R(x * (1 + t))
- 2R(x * t)`. From the expansion of the Clausen functions we get
```
R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t) = sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
Call this sum `S`.

## Handling `W(x) * log(x)`
This is the factor `W(x)` with the `log(x)` in the denominator
removed, taking the absolute value this gives us
```
W(x) * log(x) = x^(1 + α) / (gamma(1 + α) * (1 - x^p0))
```
Notice that this term is positive so we don't have to include the
absolute value.
- **TODO:** Finish bounding this

## Handling `I₂₁`
For the first integral we get
```
abs(inv(log(x)) * I₂₁) = abs(
    ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t dt
)
```
Replacing `R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)` with the sum
`S` defined above and switching the integration and summation gives us
```
abs(inv(log(x)) * I₂₁) = abs(sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ∫ ((1 - t)^2m + (1 + t)^2m - 2t^2m) * t dt for m = 1:Inf
)) <=
sum(
    abs(zeta(-α - 2m)) / factorial(2m) * x^2m * abs(∫ ((1 - t)^2m + (1 + t)^2m - 2t^2m) * t dt) for m = 1:Inf
)
```
The integral taken from `1` to `π / x` can be determined to be
```
∫ ((1 - t)^2m + (1 + t)^2m - 2t^2m) * t dt =
    2(π / x)^2m * sum(
        binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for k = 0:m-1
    )
```
Inserting this back into the sum we get
```
sum(
    abs(zeta(-α - 2m)) / factorial(2m) * x^2m * (2(π / x)^2m * abs(sum(
        binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for k = 0:m-1
    ))) for m = 1:Inf
) =
2sum(
    abs(zeta(-α - 2m)) / factorial(2m) * π^2m * abs(sum(
        binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for k = 0:m-1
    )) for m = 1:Inf
)
```
Call this sum `Σ₂₁`.
- **TODO:** This sum seems to be mostly constant in both `α` and `x`.
  We could therefore attempt to compute a uniform enclosure of it. The
  interior sum converges to`(2m - 1) / 2` as `x -> 0`. But does so
  from above so we don't get an immediate bound for it.

In the end this gives us
```
abs(W(x) * I₂₁) = (W(x) * log(x)) * (inv(log(x)) * I₂₁) <=
    (W(x) * log(x)) * Σ₂₁
```

## Handling `I₂₁`
For the second integral we get
abs(inv(log(x)) * I₂₂) = abs(
    inv(log(x)) *
    ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t * log(t) dt
)
Replacing `R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)` with the sum
`S` defined above and switching the integration and summation gives us
```
abs(inv(log(x)) * I₂₂) = abs(
    inv(log(x)) *
    sum(
        zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m *
        ∫ ((1 - t)^2m + (1 + t)^2m - 2t^2m) * t *log(t) dt for m = 1:Inf
    )
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
abs(W(x) * I₂₂) = ((W(x) * log(x)) * abs(inv(log(x)) * I₂₂) =
    ((W(x) * log(x)) * abs(log(π / x) / log(x) * Σ₂₂) =
```
An enclosure of `log(π / x) / log(x)` is easily computed by rewriting
it as `log(π) / log(x) - 1`.

## Handling I₂₃
For the third integral we have
```
abs(inv(log(x)) * I₂₃) = abs(
    inv(log(x)) *
    ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t * log(1 + u0.c * x * t) dt
)
```
Similarly to for `I₁₃` we can factor out an enclosure of `log(1 + u0.c
* x * t)` from the integral, leaving us with the same integral as for
`I₂₁`, that is
```
abs(inv(log(x)) * I₂₃) = abs(inv(log(x) * log(1 + u0.c * x * [0, 1]) * Σ₂₁)
```
and hence
```
abs(W(x) * I₂₃) = (W(x) * log(x)) * abs(inv(log(x) * log(1 + u0.c * x * [0, 1]) * Σ₂₁)
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

    # Setup for bounding F

    # Compute the expansion of u0 and remove the leading term, which
    # is handled separately.
    u0_expansion = u0(ϵ, AsymptoticExpansion())
    delete!(u0_expansion, (1, 0, 0, 0, 0, 0, 0))

    # Ensure that the tail of the expansion of u0 is positive, so that
    # we can remove it from the denominator of F1 and still get an
    # upper bound.
    expansion_ispositive(u0, u0_expansion, ϵ) ||
        error("expansion of u0 not proven to be positive, this should not happen")

    # Function for computing an upper bound of F
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

            return inv(π) * F11 * F12
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
        factor_I₁ = sinpi(Arb((-1, -1 + u0.ϵ)) / 2)

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

        # TODO: This is not finished and currently contains several
        # different implementations of the primitive function.
        WxI₁₂ = let α = Arb(-1 + u0.ϵ)
            # Primitive function of (t - 1)^(-α - 1) * t * log(t)
            primitive_I12_1(t) =
                beta_inc(Acb(α), Acb(1 - α), Acb(1 / t)) / ((α - 1) * α) -
                (t - 1)^(1 - α) / (α - 1)^2 -
                (t - 1)^-α * (α * t - 1) * log(t) / ((α - 1) * α)

            primitive_I12_1_one = beta_inc(α, 1 - α, Arb(1)) / ((α - 1) * α)

            # Primitive function of (t + 1)^(-α - 1) * t * log(t)
            primitive_I12_2(t) =
                beta_inc(Acb(1 - α), Acb(0), Acb(1 + t)) / α -
                beta_inc(Acb(2 - α), Acb(0), Acb(1 + t)) / (α - 1) -
                (1 + t)^(-α) * (1 + α * t) * log(t) / (α * (α - 1))

            primitive_I12_2_one =
                beta_inc(Acb(1 - α), Acb(0), Acb(2)) / α -
                beta_inc(Acb(2 - α), Acb(0), Acb(2)) / (α - 1)

            # Primitive function of -2t^(-α - 1) * t * log(t)
            primitive_I12_3(t) = 2t^(1 - α) / (α - 1)^2 + 2t^(1 - α) * log(t) / (α - 1)

            primitive_I12_3_one = 2 / (α - 1)^2

            # Integral on [1, t], putting similar terms together
            integral_I12(t) =
                let
                    # Integral of (t - 1)^(-α - 1) * t * log(t)
                    part1 =
                        (
                            (beta_inc(α, 1 - α, 1 / t) - beta_inc(α, 1 - α, one(t))) / ((α - 1) * α)
                        ) - (t - 1)^(1 - α) / (α - 1)^2 -
                        (t - 1)^-α * (α * t - 1) * log(t) / ((α - 1) * α)

                    # Integral of (t + 1)^(-α - 1) * t * log(t)
                    part2 =
                        (
                            beta_inc(Acb(1 - α), Acb(0), Acb(1 + t)) -
                            beta_inc(Acb(1 - α), Acb(0), Acb(2))
                        ) / α -
                        (
                            beta_inc(Acb(2 - α), Acb(0), Acb(1 + t)) -
                            beta_inc(Acb(2 - α), Acb(0), Acb(2))
                        ) / (α - 1) -
                        (1 + t)^(-α) * (1 + α * t) * log(t) / (α * (α - 1))

                    # Integral of t^(-α - 1) * t * log(t)
                    part3 = 2(t^(1 - α) - 1) / (α - 1)^2 + 2t^(1 - α) * log(t) / (α - 1)

                    return part1 + part2 + part3
                end

            # Integral of (t - 1)^(-α - 1) * t * log(t)
            part1 = primitive_I12_1(Arb(π / x)) - primitive_I12_1_one
            # Integral of (t + 1)^(-α - 1) * t * log(t)
            part2 = real(primitive_I12_2(Arb(π / x))) - primitive_I12_2_one
            # Integral of t^(-α - 1) * t * log(t)
            part3 = primitive_I12_3(Arb(π / x)) - primitive_I12_3_one

            # Full integral from parts
            I12_1 = part1 + part2 + part3

            I12_2 = integral_I12(π / x)

            integral_I₁₂ = I12_2

            term = let p0 = (1 + α) + (1 + α)^2 / 2
                integral_I₁₂ / ((1 - x^p0) * log(x))
            end

            # Alternative approach for computing the term, by
            # expanding the integrand in α.

            # Keep only leading term in the integral and approximate 1 - x^p0 = -(1 + α) * log(x)
            # Keeping only the leading terms i
            term_asym_1 =
                (6 - Arb(π)^2 - 12log(Arb(π))^2) / 24log(x)^2 + log(Arb(π)) / log(x) - 1 / 2

            # Keep only leading term in the integral but don't approximate 1 - x^p0
            term_asym_2 = let p0 = (1 + α) + (1 + α)^2 / 2
                -(1 + α) *
                inv(1 - x^p0) *
                ((6 - Arb(π)^2 - 12log(Arb(π))^2) / 24log(x) + log(Arb(π)) - log(x) / 2)
            end

            # Keep all terms in the integral but approximate the
            # integrals by their leading term
            term_asym_3 = let p0 = (1 + α) + (1 + α)^2 / 2
                (x^(α + 1) + (1 - x^(α + 1)) / ((α + 1) * log(x))) / (1 - x^p0)
            end

            -factor_I₁ * term_asym_3
        end

        WxI₁₃ = begin
            # Enclosure of log(1 + u0.c * x * t) for t ∈ [0, 1]
            factor_I₁₃ = Arb((0, log1p(u0.c * x)))

            abs(invlogx * factor_I₁₃ * WxI₁₁)
        end

        WxI₁ = WxI₁₁ + WxI₁₂ + WxI₁₃

        # Handle I₂

        # Enclosure of W(x) * log(x) = x^(1 + α) / (gamma(1 + α) * (1 - x^p0))
        # FIXME: This currently computes with `α = -1 + u0.ϵ` which
        # doesn't given an upper bound, though it gives a good
        # estimate.
        Wlogx = let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2
            if Arblib.contains_zero(x)
                xᵤ = ubound(Arb, x)
                Arb((0, xᵤ^(1 + α) / (gamma(1 + α) * (1 - xᵤ^p0))))
            else
                x^(1 + α) / (gamma(1 + α) * (1 - x^p0))
            end
        end

        WxI₂₁ = begin
            # Enclosure of
            #2sum(
            #    abs(zeta(-α - 2m)) / factorial(2m) * π^2m * abs(sum(
            #        binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for k = 0:m-1
            #    )) for m = 1:Inf
            #)
            # FIXME: This only computes the first 10 terms in the sum.
            # We need to bound the tail.
            Σ₂₁ = zero(x)
            let π = Arb(π), α = Arb((-1, -1 + u0.ϵ))
                for m = 1:10
                    term =
                        abs(zeta(-α - 2m)) / factorial(2m) *
                        π^2m *
                        abs(
                            sum(
                                binomial(2m, 2k) * ((x / π)^(2(m - 1 - k)) - (x / π)^2m) / (2k + 2) for
                                k = 0:m-1
                            ),
                        )
                    Σ₂₁ += term
                end
            end
            Σ₂₁ *= 2

            Wlogx * abs(Σ₂₁)
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

            Wlogx * abs(factor_I₂₂ * Σ₂₂)
        end

        WxI₂₃ = begin
            # Enclosure of log(1 + u0.c * x * t) for t ∈ [0, 1]
            factor_I₂₃ = Arb((0, log1p(u0.c * x)))

            Wlogx * abs(invlogx * factor_I₂₃ * Σ₂₁)
        end

        WxI₂ = WxI₂₁ + WxI₂₂ + WxI₂₃

        res = F(x) * (WxI₁ + WxI₂)

        return res
    end
end

"""
    T02_alternative(u0::BHKdVAnsatz, ::Asymptotic; ϵ = Arb(2e-1))

Returns a function such that `T02_alternative(u0, Asymptotic())(x)`
computes an **upper bound** of the integral \$T_{0,2}\$ from the paper
using an evaluation strategy that works asymptotically as `x` goes to
0.

This version uses the weight
```
x * (x^(-u0.γ * (1 + α)) + log(u0.c + inv(x)))
```

It precomputes the expansions of `u0` and for that reason a number `ϵ`
has to be given, the resulting expansion will be valid for all `x <
ϵ`. The value of `ϵ` has to be less than `1 // 2`.

First of all the change of coordinates `t = y / x` leaves us with
```
x / (π * u0(x) * (x^(-u0.γ * (1 + α)) + log(u0.c + inv(x)))) *
    ∫ abs(clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t * ((x * t)^(-u0.γ * (1 + α)) + log(u0.c + inv(x * t))) dt
```
with the integration going from `1` to `π / x`. Further we can notice
that the expression inside the absolute value is positive on the whole
interval and the absolute value can hence be removed.

Next the factor
```
F(x) = inv(π) *
    (x^(-u0.γ * (1 + α)) - log(x)) / (x^(-u0.γ * (1 + α)) + log(u0.c + inv(x))) *
     gamma(1 + α) * x^-α * (1 - x^p0) / u0(x)
```
is factored out. Except for the addition of `inv(π)` this is the same
as the factor `F1` in [`F0`](@ref) and we compute an upper bound
following the same procedure as described there.

What we are left with computing is
```
abs(W(x) * I)
```
where
```
W(x) = x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * (x^(-u0.γ * (1 + α)) - log(x)))
```
and `I` the same integral as above (but without the absolute value in
the integrand).

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
I1 = gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t *
        ((x * t)^(-u0.γ * (1 + α)) + log(u0.c + inv(x * t))) dt
I2 = ∫ (R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t *
         ((x * t)^(-u0.γ * (1 + α)) + log(u0.c + inv(x * t))) dt
```
satisfying `I = I1 + I2`. Furthermore we can split
```
(x * t)^(-u0.γ * (1 + α)) + log(u0.c + inv(x * t))
```
as
```
(x * t)^(-u0.γ * (1 + α)) + log(u0.c + inv(x * t)) =
    (x * t)^(-u0.γ * (1 + α)) -
    log(x) -
    log(t) +
    log(1 + u0.c * x * t)
```
Which allows us to split the two above integrals into four integrals
each.
```
I11 = gamma(1 + α) * sinpi(α / 2) * x^(-(1 + u0.γ) * (1 + α)) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(1 - u0.γ * (1 + α)) dt
I12 = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * log(x)
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
I13 = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
I14 = gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(1 + u0.c * x * t) dt
```
and
```
I21 = x^(-u0.γ * (1 + α)) * ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t^(1 - u0.γ * (1 + α))dt
I22 = -log(x) * ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t dt
I23 = -∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t * log(t) dt
I24 = ∫ (R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t)) * t * log(1 + u0.c * x * t) dt
```

# Handling `I1`
For `I1` there are important cancellations between `I11`, `I12` and
`I13` so we therefore have to keep track of their signs. For `I14` it
is enough to bound the absolute value.

## Handling `I11`
For the first integral we get
```
W(x) * I11 = sinpi(α / 2) *
    x^(-u0.γ * (1 + α)) / ((x^(-u0.γ * (1 + α)) - log(x))) *
    inv((1 - x^p0)) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(1 - u0.γ * (1 + α)) dt
```
where the integral is taken from `1` to `π / x`. The factor `sinpi(α /
2)` converges to `-1` and can be enclosed directly.

For `u0.γ = 1 / 2` a primitive function for the integral is given by
```
2t^(1 / 2 - 3α / 2) * (
    + 2 / (3α - 1)
    - t^(1 + α) * hypgeom_2f1(3 / 2 - α / 2, 1 + α, 5 / 2 - α / 2, -t) / (α - 3)
    + (1 - t)^α * (1 - t)^-α * t^(1 + α) * hypgeom_2f1(3 / 2 - α / 2, 1 + α, 5 / 2 - α / 2, t) / (α - 3)
)
```

Alternatively we can expand the integrand using that
```
(t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
```
can be Taylor expanded around `α = -1` by writing it as
```
(t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1) =
    exp(-log(t - 1) * (α + 1)) + exp(-log(t + 1) * (α + 1)) - 2exp(-log(t) * (α + 1)) =
    sum((-1)^n * (α + 1)^n * (log(t - 1)^n + log(t + 1)^n - 2log(t)^n) / factorial(n) for n = 1:Inf)
```
Putting this inside the integral and switching the integral and the
sum we get
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt =
    sum((-1)^n * (α + 1)^n * ∫ (log(t - 1)^n + log(t + 1)^n - 2log(t)) * t * log(t) dt for n = 1:Inf)
```
where again the integral is taken from `1` to `π / x`. Now that the
integral no longer depend on `α` we can expand it around `x = 0`.
- **TODO:** Fix the above for our case
It looks like the integrals behave like
```
(-1)^n * (3 / 2)^(n - 1) * log(x)^n
```
This would give us
```
2 / 3 * sum((α + 1)^n * (1 + u0.γ)^n * log(x)^n / factorial(n) for n = 1:Inf) =
    inv(1 + u0.γ) * (exp((α + 1) * (1 + u0.γ) * log(x)) - 1) =
    inv(1 + u0.γ) * (x^((1 + u0.γ) * (α + 1)) - 1)
```
Putting it back we get
```
sinpi(α / 2) *
    x^(-u0.γ * (1 + α)) / ((x^(-u0.γ * (1 + α)) - log(x))) *
    inv((1 - x^p0)) *
    inv(1 + u0.γ) * (x^((1 + u0.γ) * (α + 1)) - 1) =
-sinpi(α / 2) / (1 + u0.γ) *
    x^(-u0.γ * (1 + α)) / ((x^(-u0.γ * (1 + α)) - log(x))) *
    (x^((1 + u0.γ) * (α + 1)) - 1) / (x^p0 - 1)
```

## Handling `I12`
For the second integral we get
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
`α`.
- **PROVE:** The it is increasing in `x` and `α`
- **TODO:** That the limit is 1 could be proved better.

## Handling `I13`
For the third integral we have
```
W(x) * I13 = -inv(log(x)) * inv(1 - x^p0) * sinpi(α / 2) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
```
where the integral is taken from `1` to `π / x`. The factor `sinpi(α /
2)` converges to `-1` and can be enclosed directly.

We can expand
```
(t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
```
in a Taylor around `α = -1` by writing it as
```
(t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1) =
    exp(-log(t - 1) * (α + 1)) + exp(-log(t + 1) * (α + 1)) - 2exp(-log(t) * (α + 1)) =
    sum((-1)^n * (α + 1)^n * (log(t - 1)^n + log(t + 1)^n - 2log(t)^n) / factorial(n) for n = 1:Inf)
```
Putting this inside the integral and switching the integral and the
sum we get
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt =
    sum((-1)^n * (α + 1)^n * ∫ (log(t - 1)^n + log(t + 1)^n - 2log(t)) * t * log(t) dt for n = 1:Inf)
```
where again the integral is taken from `1` to `π / x`. Now that the
integral no longer depend on `α` we can expand it around `x = 0`.

The remaining part is not finished, we here give the idea of what
should happen. We would like to show that the integral behaves like
```
∫ (log(t - 1)^n + log(t + 1)^n - 2log(t)) * t * log(t) dt = (-1)^n * n / (n + 1) * log(x)^(n + 1) + O(log(x)^n)
```
If we ignore the `O(log(x)^n)` term we would then get the sum
- **TODO:** Prove the above asymptotics and handle the `O(log(x)^n)` term.
```
sum((-1)^n * (α + 1)^n * (-1)^n * n / (n + 1) * log(x)^(n + 1) for n = 1:Inf) =
    sum((α + 1)^n * n / (n + 1) * log(x)^(n + 1) for n = 1:Inf)
```
This would give us
```
W(x) * I₁₂ =
    -inv(log(x)) * inv(1 - x^p0) * sinpi(α / 2) *
        sum((α + 1)^n * n / (n + 1) * log(x)^(n + 1) for n = 1:Inf) =
    -inv(1 - x^p0) * sinpi(α / 2) *
        sum((α + 1)^n * n / (n + 1) * log(x)^n for n = 1:Inf) =
```
The sum can be computed explicitly to be
```
x^(α + 1) + (1 - x^(α + 1)) / ((α + 1) * log(x))
```
Inserting this gives us
```
-inv(1 - x^p0) * sinpi(α / 2) * sum((α + 1)^n * log(x)^n for n = 1:Inf) =
    -inv(1 - x^p0) * sinpi(α / 2) * (x^(α + 1) + (1 - x^(α + 1)) / ((α + 1) * log(x))) =
    = -sinpi(α / 2) * (x^(α + 1) + (1 - x^(α + 1)) / ((α + 1) * log(x))) / (1 - x^p0)
```
- **TODO:** Finish this. Unfortunately the above seems to go to zero
  as `x` goes to `0`, whereas we would need to go to `-1 / 2`. Could
  this be because we have neglected the `O(log(x)^n)` terms in the
  integrals inside the sum? Or is there a chance that it actually goes
  to zero? In that case our whole approach fails :/

## Handling `I₁₃`
We ignore this for now since it goes to zero.

# Handling `I₂`
We ignore these for now since they all go to zero.

"""
function T02_alternative(
    u0::BHKdVAnsatz,
    ::Asymptotic;
    non_asymptotic_u0 = false,
    ϵ = Arb(2e-1),
)
    ϵ = convert(Arb, ϵ)
    @assert ϵ < 0.5

    # Setup for bounding F

    # Compute the expansion of u0 and remove the leading term, which
    # is handled separately.
    u0_expansion = u0(ϵ, AsymptoticExpansion())
    delete!(u0_expansion, (1, 0, 0, 0, 0, 0, 0))

    # Ensure that the tail of the expansion of u0 is positive, so that
    # we can remove it from the denominator of F1 and still get an
    # upper bound.
    expansion_ispositive(u0, u0_expansion, ϵ) ||
        error("expansion of u0 not proven to be positive, this should not happen")

    # Function for computing an upper bound of F
    F(x) =
        let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, xᵤ = ubound(Arb, x)
            # Note that inside this statement α refers to -1 + u0.ϵ
            # and p0 to the corresponding p0 value.

            # Enclosure of
            # abs((x^(-u0.γ * (1 + α)) - log(x)) / (x^(-u0.γ * (1 + α)) + log(u0.c + inv(x))))
            # either by direct evaluation or using monotonicity in x and α.
            F11 =
                let f =
                        α ->
                            x -> abs(
                                (abspow(x, -u0.γ * (1 + α)) - log(x)) /
                                (abspow(x, -u0.γ * (1 + α)) + log(u0.c + inv(x))),
                            )
                    if iszero(x)
                        one(Arb)
                    elseif Arblib.contains_zero(x)
                        Arb((f(Arb(-1))(xᵤ), 1))
                    else
                        f(Arb((-1, -1 + u0.ϵ)))(x)
                    end
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

            return inv(π) * F11 * F12
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

        # Handle I1

        # Enclosure of sinpi(α / 2)
        factor_I1 = sinpi(Arb((-1, -1 + u0.ϵ)) / 2)

        Wx11 = begin
            zero(x)
        end

        # TODO: Rewrite this with updated weight
        WxI12 = begin
            # Enclosure of inv(1 - x^p0) * ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
            # FIXME: This doesn't enclose the full integral yet, there
            # are some O(x^2) terms left to handle.
            integral_I12 =
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

            -factor_I1 * integral_I12
        end

        # TODO: Rewrite this with updated weight
        WxI13 = let α = Arb(-1 + u0.ϵ)
            # Keep all terms in the integral but approximate the
            # integrals by their leading term
            term_asym = let p0 = (1 + α) + (1 + α)^2 / 2
                (x^(α + 1) + (1 - x^(α + 1)) / ((α + 1) * log(x))) / (1 - x^p0)
            end

            -factor_I1 * term_asym
        end

        WxI14 = zero(x) # Neglect this

        WxI1 = WxI11 + WxI12 + WxI13 + WxI14

        # Handle I2

        I2 = zero(x) # Neglect this

        res = F(x) * (WxI1 + WxI2)

        return res
    end
end

"""
    T02_alternative2(u0::BHKdVAnsatz, ::Asymptotic; ϵ = Arb(2e-1))

Returns a function such that `T02_alternative(u0, Asymptotic())(x)`
computes an **upper bound** of the integral \$T_{0,2}\$ from the paper
using an evaluation strategy that works asymptotically as `x` goes to
0.

This version uses the weight
```
w(x) = x^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x))
```

It precomputes the expansions of `u0` and for that reason a number `ϵ`
has to be given, the resulting expansion will be valid for all `x <
ϵ`. The value of `ϵ` has to be less than `1 // 2`.

The starting point is the integral
```
inv(π * u0(x) * w(x)) *
    ∫ (clausenc(x - y, - α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * w(y) dy
```
from `x` to `π`. Where we have removed an absolute value in the
integrand since it is positive. The change of variables `t = y / x`
leaves us with
```
x / (π * u0(x) * w(x)) *
    ∫ (clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * w(t * x) dt
```
with the integration going from `1` to `π / x`. Now
```
w(t * x) = (t * x)^(1 - u0.γ * (1 + α)) * log(u0.c + inv(t * x))
         = x^(1 - u0.γ * (1 + α)) * t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(t * x))
```
so this can be simplified to
```
x / (π * u0(x) * log(u0.c + inv(t * x))) *
    ∫ (clausenc(x * (t - 1), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(t * x)) dt
```

Next the factor
```
F(x) = abs(inv(π) * log(x) / log(u0.c + inv(x)) * gamma(1 + α) * x^-α * (1 - x^p0) / u0(x))
```
is factored out. Except for the addition of `inv(π)` this is the same
as the factor `F1` in [`F0`](@ref) and we compute an upper bound
following the same procedure as described there.

What we are left with computing is
```
abs(W(x) * I)
```
where
```
W(x) = x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(x))
```
and `I` the same integral as above.

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
I1 = gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t))) dt

I2 = ∫ (R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t))) dt
```
satisfying `I = I1 + I2`. Furthermore we can split
```
log(u0.c + inv(x * t))
```
as
```
log(u0.c + inv(x * t)) = -log(x) - log(t) + log(1 + u0.c * x * t)
```
Which allows us to split the two above integrals into three integrals
each.
```
I11 = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * log(x)
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(1 - u0.γ * (1 + α)) dt

I12 = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(1 - u0.γ * (1 + α)) * log(t) dt

I13 = gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(1 - u0.γ * (1 + α)) * log(1 + u0.c * x * t) dt
```
and similarly for `I2`

For now we are only interested in getting the asymptotic values and we
don't care about getting an enclosure. We will therefore neglect `I13`
and `I2` which are both of lower order in `x`.

# Handling `I1`
In all of these integrals the factor `sinpi(α / 2)` comes up which
converges to `-1` as `α -> -1` and can be enclosed directly.

We will also make heavy use of the expansion
```
((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-u0.γ * (1 + α))
    = (exp(-log(t - 1) * (α + 1)) + exp(-log(t + 1) * (α + 1)) - 2exp(-log(t) * (α + 1))) *
        exp(-u0.γ * log(t) * (1 + α))
    = exp(-(log(t - 1) + u0.γ * log(t)) * (α + 1)) +
        exp(-(log(t + 1) + u0.γ * log(t)) * (α + 1)) -
        2exp(-(1 + u0.γ) * log(t) * (α + 1))
    =  sum(
        (-1)^n * (α + 1)^n / factorial(n) *
        ((log(t - 1) + u0.γ * log(t))^n + (log(t + 1) + u0.γ * log(t))^n - 2(1 + u0.γ)^n * log(t)^n)
        for n = 1:Inf
    )
```
To reduce the complicated integrals to a sum of easier integrals.

## Handling `I11`
We get
```
W(x) * I11 = -sinpi(α / 2) / (1 - x^p0) *
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(1 - u0.γ * (1 + α)) dt
```

Using the expansion of
```
((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-u0.γ * (1 + α))
```
and switching the integration and summation we get
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-u0.γ * (1 + α)) * log(t) dt =
    sum(
        (-1)^n * (α + 1)^n / factorial(n) *
        ∫ ((log(t - 1) + u0.γ * log(t))^n + (log(t + 1) + u0.γ * log(t))^n - 2(1 + u0.γ)^n * log(t)^n) *
            t dt
        for n = 1:Inf
    )
```
Now that the integrals no longer depend on `α` we want to expand them
around `x = 0`.

As a first step we attempt to compute the leading term in the
asymptotic. The leading order in the integrals comes from the
integration for large values of `t`. We therefore want to expand the
integrand at infinity.

We focus on expanding
```
((log(t - 1) + u0.γ * log(t))^n + (log(t + 1) + u0.γ * log(t))^n - 2(1 + u0.γ)^n * log(t)^n)
```
since the multiplication by `t` is easily handled. Asymptotically as
`t -> Inf` this behaves like
```
- n * (1 + u0.γ)^(n - 1) * log(t)^(n - 1) / t^2
```
- **TODO:** Write down the calculations. It is done by rewriting
  `log(t ± 1) = log(t * (1 ± inv(t))) = log(t) + log(1 ± inv(t))` and
  using the binomial theorem, noticing that the two leading terms
  cancel and that the third one is `- n / (n + 1) * (1 + u0.γ)^(n - 1)
  * log(t)^(n + 1)`
This gives us the integrals
```
- n * (1 + u0.γ)^(n - 1) * ∫ log(t)^(n - 1) / t dt
```
The primitive function is given by
```
∫ log(t)^(n - 1) / t dt = log(t)^n / n
```
The integration is taken from `1` to `π / x` but since we only care
about what happens for small `x` (large `t`) we just let `t = 1 / x`,
giving us that
```
- n * (1 + u0.γ)^(n - 1) * ∫ log(t)^(n - 1) / t dt
```
asymptotically behaves like
```
-(-1)^n * (1 + u0.γ)^(n - 1) * log(1 / x)^n =
```
Inserting this into the sum we get
```
-sum((α + 1)^n / factorial(n) * (1 + u0.γ)^(n - 1) * log(x)^n for n = 1:Inf) =
-inv(1 + u0.γ) * sum((α + 1)^n / factorial(n) * (1 + u0.γ)^n * log(x)^n for n = 1:Inf) =
-inv(1 + u0.γ) * (exp((1 + u0.γ) * (1 + α) * log(x)) - 1) =
-(x^((1 + u0.γ) * (1 + α)) - 1) / (1 + u0.γ) =
(1 - x^((1 + u0.γ) * (1 + α))) / (1 + u0.γ)
```

Now adding the factor `-sinpi(α / 2) / (1 - x^p0)` in front of the
integral we get
```
W(x) * I11 = sinpi(α / 2) / (1 + u0.γ) * (1 - x^((1 + u0.γ) * (1 + α))) / (1 - x^p0)
```

## Handling `I12`
The procedure is very similar to that for `I11`. We get
```
W(x) * I12 = -sinpi(α / 2) / ((1 - x^p0) * log(x))
    ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(1 - u0.γ * (1 + α)) * log(t) dt
```

Using the expansion of
```
((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-u0.γ * (1 + α))
```
and switching the integration and summation we get
```
∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-u0.γ * (1 + α)) * log(t) dt =
    sum(
        (-1)^n * (α + 1)^n / factorial(n) *
        ∫ ((log(t - 1) + u0.γ * log(t))^n + (log(t + 1) + u0.γ * log(t))^n - 2(1 + u0.γ)^n * log(t)^n) *
            t * log(t) dt
        for n = 1:Inf
    )
```
Now that the integral no longer depend on `α` we can expand it around
`x = 0`.

As a first step we attempt to compute the leading term in the
asymptotic. The leading order in the integrals comes from the
integration for large values of `t`. We therefore want to expand the
integrand at infinity.

Similarly to above we focus on expanding
```
((log(t - 1) + u0.γ * log(t))^n + (log(t + 1) + u0.γ * log(t))^n - 2(1 + u0.γ)^n * log(t)^n)
```
and then add the multiplication by `t * log(t)`. As above we get that
asymptotically as `t -> Inf` this behaves like
```
- n * (1 + u0.γ)^(n - 1) * log(t)^(n - 1) / t^2
```
This gives us the integrals
```
- n * (1 + u0.γ)^(n - 1) * ∫ log(t)^n / t dt
```
The primitive function is given by
```
∫ log(t)^n / t dt = log(t)^(n + 1) / (n + 1)
```
The integration is taken from `1` to `π / x` but since we only care
about what happens for small `x` (large `t`) we just let `t = 1 / x`,
giving us that
```
- n * (1 + u0.γ)^(n - 1) * ∫ log(t)^(n - 1) / t dt
```
asymptotically behaves like
```
-(-1)^n * n / (n + 1) * (1 + u0.γ)^(n - 1) * log(1 / x)^(n + 1) =
```
Inserting this into the sum we get
```
-sum((α + 1)^n / factorial(n) * n / (n + 1) * (1 + u0.γ)^(n - 1) * log(x)^(n + 1) for n = 1:Inf) =
-log(x) / (1 + u0.γ) * sum(n / (n + 1) * (α + 1)^n / factorial(n) * (1 + u0.γ)^n * log(x)^n for n = 1:Inf) =
-log(x) / (1 + u0.γ) * (
    (1 - x^((1 + u0.γ) * (1 + α))) / ((1 + u0.γ) * (1 + α) * log(x)) +
    x^((1 + u0.γ) * (1 + α))
)
```

Now adding the factor `-sinpi(α / 2) / ((1 - x^p0) * log(x))` in front of the
integral we get
```
W(x) * I11 = sinpi(α / 2) / ((1 - x^p0) * log(x)) * log(x) / (1 + u0.γ) * (
    (1 - x^((1 + u0.γ) * (1 + α))) / ((1 + u0.γ) * (1 + α) * log(x)) +
    x^((1 + u0.γ) * (1 + α))
) =
sinpi(α / 2) / (1 + u0.γ) * (
    (1 - x^((1 + u0.γ) * (1 + α))) / ((1 + u0.γ) * (1 + α) * log(x)) +
    x^((1 + u0.γ) * (1 + α))
) / (1 - x^p0)
```

## Handling `I13`
We neglect this for now
- **TODO:** Handle this.

# Handling `I₂`
We ignore these for now since they all go to zero.
- **TODO:** Handle this.

"""
function T02_alternative2(
    u0::BHKdVAnsatz,
    ::Asymptotic;
    non_asymptotic_u0 = false,
    ϵ = Arb(2e-1),
)
    ϵ = convert(Arb, ϵ)
    @assert ϵ < 0.5

    # Setup for bounding F

    # Compute the expansion of u0 and remove the leading term, which
    # is handled separately.
    u0_expansion = u0(ϵ, AsymptoticExpansion())
    delete!(u0_expansion, (1, 0, 0, 0, 0, 0, 0))

    # Ensure that the tail of the expansion of u0 is positive, so that
    # we can remove it from the denominator of F1 and still get an
    # upper bound.
    expansion_ispositive(u0, u0_expansion, ϵ) ||
        error("expansion of u0 not proven to be positive, this should not happen")

    # Function for computing an upper bound of F
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

            return inv(π) * F11 * F12
        end

    return x::Arb -> begin
        @assert x <= ϵ

        sin_factor = let α = Arb((-1, -1 + u0.ϵ))
            sinpi(α / 2)
        end

        WxI11 = let α = -1 + u0.ϵ, p0 = (1 + α) + (1 + α)^2 / 2
            # TODO: Multiply by sin_factor
            1 / (1 + u0.γ) * (1 - x^((1 + u0.γ) * (1 + α))) / (1 - x^p0)
        end

        WxI12 = let α = -1 + u0.ϵ, p0 = (1 + α) + (1 + α)^2 / 2
            # TODO: Multiply by sin_factor
            1 / (1 + u0.γ) * (
                (1 - x^((1 + u0.γ) * (1 + α))) / ((1 + u0.γ) * (1 + α) * log(x)) +
                x^((1 + u0.γ) * (1 + α))
            ) / (1 - x^p0)
        end

        #@show WxI11 WxI12

        WxI1 = WxI11 + WxI12

        return WxI11, WxI12, WxI1

        # Handle I2

        WxI2 = zero(x) # Neglect this

        res = F(x) * (WxI1 + WxI2)

        return res
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
- **TODO:** Could improve enclosures by better handling cancellations
  for `clausens(x + a, 1 - α) - clausens(2x, 1 - α)` and `clausens(a,
  1 - α) - clausens(x, 1 - α)`. Though this might not be needed.
"""
function T021(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    # s = 1 - α computed such that the upper bound is exactly 2
    s = 2 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    return (x::Arb, a::Arb = x + δ2) -> begin
        integral =
            -clausens(x - a, s) + (clausens(x + a, s) - clausens(2x, s)) -
            2(clausens(a, s) - clausens(x, s))

        # Multiply by weight inside the integral that was factored out
        integral *= u0.w(union(x, a))

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

It returns a function `f` such that `f(x, a; tol)` computes the
integral on `[a, π]` for the given value of `x` and it uses the
prescribed tolerance for the integration. In practice `a` should be a
thin ball to not give problems with the integration.

This is done by directly computing the integral with the integrator in
Arb.

Notice that due to lemma [`lemma_integrand_2`](@ref) the expression
inside the absolute value is always positive, so we can remove the
absolute value.
"""
function T022(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    # Lower and upper bounds of s = -α
    s_l = 1 - u0.ϵ
    s_u = one(Arb)

    # Enclosure of -α so that the upper bound is exactly 1
    mα = 1 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    # Upper integration limit
    b = Arb(π)

    return (x::Arb, a::Arb = x + δ2; tol = Arb(1e-5)) -> begin
        integrand(y) =
            (clausenc(y - x, mα) + clausenc(y + x, mα) - 2clausenc(y, mα)) * u0.w(y)

        res = ArbExtras.integrate(integrand, a, b, atol = tol, rtol = tol)

        res /= (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
