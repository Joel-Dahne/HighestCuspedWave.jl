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

**TODO:** The primitive function of the integrand can be computed
explicitly and from this and explicit expression for the integral can
be computed. We expect `W(x) * I₁₁` to go to a non-zero limit so the
integral should thus behave similarly to `inv(log(x)) * inv(1 -
x^p0)`. However the expression you get contains several special
functions and is very involved, making it difficult to compute the
limit. The code below implements the primitive function in several
different versions but doesn't handle the limit yet.

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
))
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
            # First version
            begin
                # Primitive function of full integral
                primitive_I12_v1(t) =
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
            end

            # Second version
            begin
                # Primitive function of (t - 1)^(-α - 1) * t * log(t)
                primitive_I12_v2_1(t) =
                    t^-α * hypgeom_2f1(α, α, 1 + α, inv(t)) / ((α - 1) * α^2) -
                    (t - 1)^(1 - α) / (α - 1)^2 -
                    (t - 1)^-α * (α * t - 1) * log(t) / ((α - 1) * α)

                primitive_I12_v2_1_one =
                    hypgeom_2f1(α, α, 1 + α, Arb(1)) / ((α - 1) * α^2)

                # Primitive function of (t + 1)^(-α - 1) * t * log(t)
                primitive_I12_v2_2(t) =
                    -(t + 1)^(1 - α) *
                    hypgeom_2f1(Acb(1), Acb(1 - α), Acb(2 - α), Acb(1 + t)) /
                    ((α - 1) * α) +
                    (t + 1)^(2 - α) *
                    hypgeom_2f1(Acb(1), Acb(2 - α), Acb(3 - α), Acb(1 + t)) /
                    ((α - 2) * (α - 1)) -
                    (t + 1)^-α * (1 + α * t) * log(t) / ((α - 1) * α)

                primitive_I12_v2_2_one =
                    -2^(1 - α) * hypgeom_2f1(Acb(1), Acb(1 - α), Acb(2 - α), Acb(2)) /
                    ((α - 1) * α) +
                    2^(2 - α) * hypgeom_2f1(Acb(1), Acb(2 - α), Acb(3 - α), Acb(2)) /
                    ((α - 2) * (α - 1))

                # Primitive function of -2t^(-α - 1) * t * log(t)
                primitive_I12_v2_3(t) =
                    2t^(1 - α) / (α - 1)^2 + 2t^(1 - α) * log(t) / (α - 1)

                primitive_I12_v2_3_one = 2 / (α - 1)^2

                # Primitive function of whole integral
                primitive_I12_v2(t) = begin
                    # All terms with log(t) above.
                    # TODO: Explicitly handle the cancellations
                    # (t - 1)^-α - (t + 1)^-α is around -2
                    # (t - 1)^-α + (t + 1)^-α - 2t^-α is -0.00015?
                    a =
                        log(t) * (
                            ((t - 1)^-α - (t + 1)^-α) / α -
                            t * ((t - 1)^-α + (t + 1)^-α - 2t^-α)
                        ) / (α - 1)

                    # TODO: Handle cancellations between b1, b2 and c
                    b1 =
                        t^-α * hypgeom_2f1(α, α, 1 + α, inv(t)) / ((α - 1) * α^2) - (t - 1)^(1 - α) / (α - 1)^2

                    b2 =
                        -(t + 1)^(1 - α) * hypgeom_2f1(
                            Acb(1),
                            Acb(1 - α),
                            Acb(2 - α),
                            Acb(1 + t),
                        ) / ((α - 1) * α) +
                        (t + 1)^(2 - α) * hypgeom_2f1(
                            Acb(1),
                            Acb(2 - α),
                            Acb(3 - α),
                            Acb(1 + t),
                        ) / ((α - 2) * (α - 1))

                    c = 2t^(1 - α) / (α - 1)^2

                    a + b1 + b2 + c
                end

                primitive_I12_v2_one =
                    primitive_I12_v2_1_one + primitive_I12_v2_2_one + primitive_I12_v2_3_one

            end

            # Third version
            begin
                # Primitive function of (t - 1)^(-α - 1) * t * log(t)
                primitive_I12_v3_1(t) =
                    beta_inc(Acb(α), Acb(1 - α), Acb(1 / t)) / ((α - 1) * α) -
                    (t - 1)^(1 - α) / (α - 1)^2 -
                    (t - 1)^-α * (α * t - 1) * log(t) / ((α - 1) * α)

                primitive_I12_v3_1_one = beta_inc(α, 1 - α, Arb(1)) / ((α - 1) * α)

                # Primitive function of (t + 1)^(-α - 1) * t * log(t)
                primitive_I12_v3_2(t) =
                    beta_inc(Acb(1 - α), Acb(0), Acb(1 + t)) / α -
                    beta_inc(Acb(2 - α), Acb(0), Acb(1 + t)) / (α - 1) -
                    (1 + t)^(-α) * (1 + α * t) * log(t) / (α * (α - 1))

                primitive_I12_v3_2_one =
                    beta_inc(Acb(1 - α), Acb(0), Acb(2)) / α -
                    beta_inc(Acb(2 - α), Acb(0), Acb(2)) / (α - 1)

                # Primitive function of -2t^(-α - 1) * t * log(t)
                primitive_I12_v3_3(t) =
                    2t^(1 - α) / (α - 1)^2 + 2t^(1 - α) * log(t) / (α - 1)

                primitive_I12_v3_3_one = 2 / (α - 1)^2

                # Integral on [1, t], putting similar terms together
                integral_I12_v3(t) = begin
                    # Integral of (t - 1)^(-α - 1) * t * log(t)
                    part1_v3 =
                        (
                            (beta_inc(α, 1 - α, 1 / t) - beta_inc(α, 1 - α, one(t))) / ((α - 1) * α)
                        ) - (t - 1)^(1 - α) / (α - 1)^2 -
                        (t - 1)^-α * (α * t - 1) * log(t) / ((α - 1) * α)

                    # Integral of (t + 1)^(-α - 1) * t * log(t)
                    part2_v3 =
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
                    part3_v3 =
                        2(t^(1 - α) - 1) / (α - 1)^2 + 2t^(1 - α) * log(t) / (α - 1)
                    @show part2_v3
                    return part1_v3 + part2_v3 + part3_v3
                end

                # Like integral_I12_v3 above but with a different order for the terms
                integral_I12_v3_reorder(t, α = -1 + u0.ϵ) = begin
                    part_beta =
                        (
                            (beta_inc(α, 1 - α, 1 / t) - beta_inc(α, 1 - α, one(t))) / ((α - 1) * α)
                        ) +
                        (
                            beta_inc(Acb(1 - α), Acb(0), Acb(1 + t)) -
                            beta_inc(Acb(1 - α), Acb(0), Acb(2))
                        ) / α

                    part_log =
                        ((t - 1)^-α - (t + 1)^-α) * log(t) / (α * (α - 1)) -
                        ((t - 1)^-α + (t + 1)^-α - 2t^-α) * t * log(t) / (α - 1)

                    part_rest =
                        -(t - 1)^(1 - α) / (α - 1)^2 -
                        (
                            beta_inc(Acb(2 - α), Acb(0), Acb(1 + t)) -
                            beta_inc(Acb(2 - α), Acb(0), Acb(2))
                        ) / (α - 1) + 2t^(1 - α) / (α - 1)^2

                    part_constant = -2 / (α - 1)^2

                    return part_beta + part_log + part_rest + part_constant
                end

                # Like integral_I12_v3 above but replacing the beta
                # functions using beta_inc(a, b, x) = x^a / a *
                # hypgeom_2f1(a, 1 - b, a + 1, x)
                integral_I12_v3_2F1(t) =
                    let
                        # Integral of (t - 1)^(-α - 1) * t * log(t)
                        part1_v3 =
                            (
                                (
                                    t^-α * hypgeom_2f1(α, α, α + 1, 1 / t) -
                                    hypgeom_2f1(α, α, α + 1, one(t))
                                ) / ((α - 1) * α^2)
                            ) - (t - 1)^(1 - α) / (α - 1)^2 -
                            (t - 1)^-α * (α * t - 1) * log(t) / ((α - 1) * α)

                        # Integral of (t + 1)^(-α - 1) * t * log(t)
                        part2_v3 =
                            (
                                (1 + t)^(1 - α) * hypgeom_2f1(
                                    Acb(1 - α),
                                    Acb(1),
                                    Acb(2 - α),
                                    Acb(1 + t),
                                ) -
                                2^(1 - α) *
                                hypgeom_2f1(Acb(1 - α), Acb(1), Acb(2 - α), Acb(2))
                            ) / ((1 - α) * α) -
                            (
                                (1 + t)^(2 - α) * hypgeom_2f1(
                                    Acb(2 - α),
                                    Acb(1),
                                    Acb(3 - α),
                                    Acb(1 + t),
                                ) -
                                2^(2 - α) *
                                hypgeom_2f1(Acb(2 - α), Acb(1), Acb(3 - α), Acb(2))
                            ) / ((2 - α) * (α - 1)) -
                            (1 + t)^(-α) * (1 + α * t) * log(t) / (α * (α - 1))

                        # Integral of t^(-α - 1) * t * log(t)
                        part3_v3 =
                            2(t^(1 - α) - 1) / (α - 1)^2 + 2t^(1 - α) * log(t) / (α - 1)
                        @show part2_v3
                        return part1_v3 + part2_v3 + part3_v3
                    end
            end



            # Fourth version, with α = -1
            begin
                # Like integral_I12_v3_reorder above but with α = -1
                integral_I12_v4(t) = begin
                    part_beta = let α = Arb(-1 + 1e-5 * u0.ϵ)
                        (beta_inc(α, Arb(2), 1 / t) - beta_inc(α, Arb(2), one(t))) / 2 - (
                            beta_inc(Acb(2), Acb(0), Acb(1 + t)) -
                            beta_inc(Acb(2), Acb(0), Acb(2))
                        )
                    end

                    part_log = -log(t)

                    part_rest =
                        (t^2 + 2t) / 4 + real(beta_inc(Acb(3), Acb(0), Acb(1 + t))) / 2

                    part_constant = Arb(5 // 4)

                    return part_beta + part_log + part_rest + part_constant
                end
            end

            # Integral of (t - 1)^(-α - 1) * t * log(t)
            part1_v2 = primitive_I12_v2_1(Arb(π / x)) - primitive_I12_v2_1_one
            # Integral of (t + 1)^(-α - 1) * t * log(t)
            part2_v2 = real(primitive_I12_v2_2(Arb(π / x))) - primitive_I12_v2_2_one
            # Integral of t^(-α - 1) * t * log(t)
            part3_v2 = primitive_I12_v2_3(Arb(π / x)) - primitive_I12_v2_3_one

            # Integral of (t - 1)^(-α - 1) * t * log(t)
            part1_v3 = primitive_I12_v3_1(Arb(π / x)) - primitive_I12_v3_1_one
            # Integral of (t + 1)^(-α - 1) * t * log(t)
            part2_v3 = real(primitive_I12_v3_2(Arb(π / x))) - primitive_I12_v3_2_one
            # Integral of t^(-α - 1) * t * log(t)
            part3_v3 = primitive_I12_v3_3(Arb(π / x)) - primitive_I12_v3_3_one

            # Full integral from parts
            I12_1 = part1_v2 + part2_v2 + part3_v2

            # Full integral from full primitive functions
            I12_2 = primitive_I12_v2(π / x) - primitive_I12_v2_one

            I12_3 = integral_I12_v3(π / x)

            I12_4 = integral_I12_v3_reorder(π / x)

            I12_5 = integral_I12_v4(π / x)

            I12_6 = integral_I12_v3_2F1(π / x)

            # Compute the actual value of W(x) * I₁₂. The above is
            # mostly for testing.
            integral_I₁₂ = I12_2

            term = let p0 = (1 + α) + (1 + α)^2 / 2
                integral_I₁₂ / ((1 - x^p0) * log(x))
            end

            -factor_I₁ * term
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

It returns a function `f` such that `f(x, a; tol)` computes the
integral on `[a, π]` for the given value of `x` and it uses the
prescribed tolerance for the integration. In practice `a` should be a
thin ball to not give problems with the integration.

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
    # Lower and upper bounds of s = -α
    s_l = 1 - u0.ϵ
    s_u = one(Arb)

    # Upper integration limit
    b = Arb(π)

    return (x::Arb, a::Arb = x + δ2; tol = Arb(1e-5)) -> begin
        integrand(y) = begin
            # FIXME: Currently we assume monotonicity in s, including for
            # all derivatives.
            ymx = y - x
            ypx = y + x

            term_l = clausenc(ymx, s_l) + clausenc(ypx, s_l) - 2clausenc(y, s_l)
            term_u = clausenc(ymx, s_u) + clausenc(ypx, s_u) - 2clausenc(y, s_u)

            if y isa ArbSeries
                coefficients = union.(Arblib.coeffs(term_l), Arblib.coeffs(term_u))
                term_union = ArbSeries(coefficients)
            else
                term_union = union(term_l, term_u)
            end

            return term_union * u0.w(y)
        end

        res = ArbExtras.integrate(integrand, a, b, atol = tol, rtol = tol)

        res /= (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
