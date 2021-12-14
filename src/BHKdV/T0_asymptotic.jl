"""
    _bound_remainder_sum(α::Arb)

Compute an upper bound of the absolute value of the sum
```
sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^(2(m - 1)) * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
For `0 <= x < 1` and `0 <= t <= π / x`.

Taking the absolute value we have
abs(
    sum(
        zeta(-α - 2m) * (-1)^m / factorial(2m) * x^(2(m - 1)) * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
    )
) <=
abs(
    sum(
        abs(zeta(-α - 2m) * (-1)^m / factorial(2m) * x^(2(m - 1)) * ((1 - t)^2m + (1 + t)^2m - 2t^2m)) for m = 1:Inf
    )
)
```
Next we can notice that
```
(1 - t)^2m + (1 + t)^2m - 2t^2m =
    sum(binomial(2m, k) * (1 + (-1)^k) * t^k for k = 0:2m) - 2t^2m =
    2sum(binomial(2m, 2k) * t^2k for k = 0:m-1)
```
which is positive and increasing in `t`. An upper bound on the
interval `[0, π / x]` can thus be given by letting `t = π / x`.
Inserting this into the sum gives us
```
2sum(
    abs(zeta(-α - 2m)) / factorial(2m) * x^(2(m - 1)) *
        sum(binomial(2m, 2k) * (π / x)^2k for k = 0:m-1) for m = 1:Inf
) =
2sum(
    abs(zeta(-α - 2m)) * (-1)^m / factorial(2m) *
        sum(binomial(2m, 2k) * π^2k * x^(2(m - 1 - k)) for k = 0:m-1) for m = 1:Inf
)
```
Where we have also moved all factors which are positive outside the
absolute value. Now since `m - 1 - k >= 0` in all cases we get that
this sum is increasing in `x` and since `x < 1` an upper bound can be
computed by letting `x = 1`. This gives us
```
2sum(
    abs(zeta(-α - 2m)) * (-1)^m / factorial(2m) *
        sum(binomial(2m, 2k) * π^2k for k = 0:m-1) for m = 1:Inf
)
```
We have that
```
2sum(binomial(2m, 2k) * π^2k for k = 0:m-1) = (1 - π)^2m + (1 + π)^2m - 2π^2m
```
giving us
```
sum(
    abs(zeta(-α - 2m)) * (-1)^m / factorial(2m) * x^(2(m - 1)) *
        ((1 - π)^2m + (1 + π)^2m - 2π^2m) for m = 1:Inf
)
```
- **TODO:** Bound this sum. For now we just compute a finite number of
  terms.
"""
function _bound_remainder_sum(α::Arb)
    # FIXME: For now we just compute a finite number of terms in
    # the sum. We need to bound the remainder.

    C = zero(Arb)

    let π = Arb(π)
        for m = 1:10
            term = abs(zeta(-α - 2m)) / factorial(2m) * ((1 - π)^2m + (1 + π)^2m - 2π^2m)
            C += term
        end
    end

    return C
end

"""
    _T0_asymptotic_main(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
W(x) * I_M(x) = sinpi(α / 2) / ((1 - x^p0) * log(x)) *
                ∫ abs(abs(t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
                    t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
defined in [`T0_asymptotic`](@ref).

The function `sinpi(α / 2)` converges to `-1` and can be enclosed
directly, by factoring out `-sinpi(α / 2)` we can therefore focus on
computing an upper bound of
```
-1 / ((1 - x^p0) * log(x)) *
    ∫ abs(abs(t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
        t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```

# Split into two intervals
The function
```
abs(t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)
```
in the integrand has a singularity at `t = 1` so it is natural to
split the interval `[0, π / x]` into two parts, `[0, 1]` and `[1, π /
x]`.

We handle these two intervals separately.

# The interval `[0, 1]`
On this interval `t - 1` is non-positive so the problem reduces to
computing an upper bound of
```
G1(x) = -1 / ((1 - x^p0) * log(x)) *
            ∫ abs((1 - t)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
where the integration is taken from `0` to `1`.

This term tends to zero as `α -> -1` and `x -> 0` and we therefore
don't need to compute a very accurate upper bound.

**TODO:** Finish this. Possibly by multiplying and dividing by
`gamma(1 + α)` to make the integrand no longer go to zero. And then
splitting the log-term.

# The interval `[1, π / x]`
On this interval `t - 1` is non-negative so the problem reduces to
computing an upper bound of
```
G2(x) = -1 / ((1 - x^p0) * log(x)) *
            ∫ abs((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
where the integration is taken from `1` to `π / x`.

In this case the expression inside the absolute value
```
(t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)
```
can be checked to be positive on `[1, Inf]` and the absolute
value can hence be removed.
- **PROVE:** That `(t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)`
  is positive on `[1, Inf]`.
This gives us
```
G2(x) = -1 / ((1 - x^p0) * log(x)) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```

## Expanding integrand
We want to get rid of the fractional exponents inside the integral. To
do this we focus on the part of the integrand given by
```
((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-γ * (1 + α))
```
We can rewrite this as
```
((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-γ * (1 + α)) =
    (exp(-log(t - 1) * (1 + α)) + exp(-log(t + 1) * (1 + α)) - 2exp(-log(t) * (1 + α))) *
        exp(-γ * log(t) * (1 + α)) =
    exp(-(log(t - 1) + γ * log(t)) * (1 + α)) +
        exp(-(log(t + 1) + γ * log(t)) * (1 + α)) -
        2exp(-(1 + γ) * log(t) * (1 + α))
```
Expanding the exponentials and joining the three sums we arrive at
```
sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n)
    for n = 1:Inf
)
```
Inserting this into the integral and switching the order of
integration and summation gives us
```
G2(x) = -1 / ((1 - x^p0) * log(x)) * sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
        t * log(c + inv(x * t)) dt
    for n = 1:Inf
)
```
Where the integration is done from `1` to `π / x`. We are hence
interested in computing the integrals
```
I(n) = ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
    t * log(c + inv(x * t)) dt
```
for `n = 1, 2, ...`.

## Computing `I(n)`

### Split into parts
As a first step we split `I(n)` into three parts by using that
```
log(c + inv(x * t)) = log((c * x * t + 1) / (x * t)) = log(1 + c * x * t) - log(x) - log(t)
```
Letting
```
I1(n) = ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
    t dt

I2(n) = ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
    t * log(t) dt

I3(n) = ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
    t * log(1 + c * x * t) dt
```
we then get
```
I(n) = -log(x) * I1(n) - I2(n) + I3(n)
```

### Simplifying the integrand
To begin with we want to study the part of the integrand given by
```
(log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n
```
Using the binomial theorem we can write this as
```
sum(binomial(n, k) γ^k * log(t)^k * log(t - 1)^(n - k) for k = 0:n) +
sum(binomial(n, k) γ^k * log(t)^k * log(t + 1)^(n - k) for k = 0:n) -
2(1 + γ)^n * log(t)^n
```
By joining the sums we can write this as
```
sum(binomial(n, k) * γ^k * log(t)^k * (log(t - 1)^(n - k) + log(t + 1)^(n - k)) for k = 0:n) -
    2(1 + γ)^n * log(t)^n
```

For large values of `t` the factor
```
log(t - 1)^(n - k) + log(t + 1)^(n - k)
```
behaves almost like `2log(t)^(n - k)`. We can therefore split it into
this term plus a remainder and handle them separately. If we let
```
R(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l
```
we can split the sum as
```
sum(binomial(n, k) * γ^k * log(t)^k * (log(t - 1)^(n - k) + log(t + 1)^(n - k)) for k = 0:n) =
    sum(binomial(n, k) * γ^k * log(t)^k * 2log(t)^(n - k) for k = 0:n) +
    sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n)
```

We can now notice that the first of these sums
```
sum(binomial(n, k) * γ^k * log(t)^k * 2log(t)^(n - k) for k = 0:n)
```
is exactly equal to
```
2(1 + γ)^n * log(t)^n
```
and they hence cancel each other out. This leaves us with only one sum
remaining
```
sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n)
```
Noticing that `R(0, t) = 0` we can simplify this to
```
sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1)
```

Inserting this back into `I1(n)`, `I2(n)` and `I3(n)` we get
```
I1(n) = ∫ sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1) * t dt

I2(n) = ∫ sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1) * t * log(t) dt

I2(n) = ∫ sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1) * t * log(1 + c * x * t) dt
```

### Computing `I1(n)` and `I2(n)`
We start by computing `I1(n)` and `I2(n)` since `I3(n)` is of lower
order. Switching the summation and integration in `I1(n)` and `I2(n)`
we arrive at
```
I1(n) = sum(binomial(n, k) * γ^k * ∫ log(t)^k * R(n - k, t) * t dt for k = 0:n-1)

I2(n) = sum(binomial(n, k) * γ^k * ∫ log(t)^(k + 1) * R(n - k, t) * t dt for k = 0:n-1)
```
Using that `R(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l` we can
write these integrals as
```
∫ log(t)^k * (log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)) * t dt
```
and
```
∫ log(t)^(k + 1) * (log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)) * t dt
```
- **TODO:** Possibly split in three and integrate explicitly? We need
  to take care of cancellations between the integrals in that case.

#### Attempt using asymptotic behaviour of integrand
**NOTE:** This is an attempt which might not work out. In particular
this expansion only holds for large `t` and we would have to do
something more to handle `t` close to `1`.

Asymptotically as `t` goes to infinity we have that
```
log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)
```
behaves like
```
(n - k) * log(t)^(n - k - 2) * ((n - k - 1) - log(t)) / t^2 + O(1 / t^4)
```
- **TODO:** Prove this and bound remainder.

Inserting the main term into the integrals for `I1(n)` and `I2(n)`
gives us
```
∫ log(t)^k * (n - k) * log(t)^(n - k - 2) * ((n - k - 1) - log(t)) / t^2 * t dt
= (n - k) * ∫ log(t)^(n - 2) * ((n - k - 1) - log(t)) / t dt
= (n - k) * (n - k - 1) * ∫ log(t)^(n - 2) / t dt -
    (n - k) * ∫ log(t)^(n - 1) / t dt

∫ log(t)^(k + 1) * (n - k) * log(t)^(n - k - 2) * ((n - k - 1) - log(t)) / t^2 * t dt
= (n - k) * ∫ log(t)^(n - 1) * ((n - k - 1) - log(t)) / t dt
= (n - k) * (n - k - 1) * ∫ log(t)^(n - 1) / t dt -
    (n - k) * ∫ log(t)^n / t dt
```
The integrals can be computed exactly to be
```
(n - k) * (n - k - 1) * ∫ log(t)^(n - 2) / t dt =
    (n - k) * (n - k - 1) / (n - 1) * log(π / x)^(n - 1)

(n - k) * ∫ log(t)^(n - 1) / t dt =
    (n - k) / n * log(π / x)^n

(n - k) * (n - k - 1) * ∫ log(t)^(n - 1) / t dt =
    (n - k) * (n - k - 1) / n * log(π / x)^n

(n - k) * ∫ log(t)^n / t dt =
    (n - k) / (n + 1) * log(π / x)^(n + 1)
```
Inserting this into `I1(n)` and `I2(n)` gives us the approximations
```
I1(n) ≈ sum(
    binomial(n, k) * γ^k * (n - k) * log(π / x)^(n - 1) * ((n - k - 1) / (n - 1) - log(π / x) / n)
    for k = 0:n-1
) =
log(π / x)^(n - 1) * (
    1 / (n - 1) * sum(binomial(n, k) * γ^k * (n - k) * (n - k - 1)  for k = 0:n-1) -
    log(π / x) / n * sum(binomial(n, k) * γ^k * (n - k) for k = 0:n-1)
)

I2(n) ≈ sum(
    binomial(n, k) * γ^k * (n - k) * log(π / x)^n * ((n - k - 1) / n - log(π / x) / (n + 1))
    for k = 0:n-1
) =
log(π / x)^n * sum(
    1 / n * sum(binomial(n, k) * γ^k * (n - k) * (n - k - 1)  for k = 0:n-1) -
    log(π / x) / (n + 1) * sum(binomial(n, k) * γ^k * (n - k) for k = 0:n-1)
)
```

The sums in `I1(n)` and `I2(n)` above are the same and they can be
explicitly computed to be
```
sum(binomial(n, k) * γ^k * (n - k) * (n - k - 1)  for k = 0:n-1) =
    (1 + γ)^(n - 2) * (n - 1) * n

sum(binomial(n, k) * γ^k * (n - k) for k = 0:n-1) =
    (1 + γ)^(n - 1) * n
```
Inserting this backs gives us
```
I1(n) ≈ log(π / x)^(n - 1) * (1 + γ)^(n - 2) * (n - log(π / x) * (1 + γ))

I2(n) ≈ log(π / x)^n * (1 + γ)^(n - 2) * ((n - 1) - n / (n + 1) * log(π / x) * (1 + γ))
```

### Computing `I3(n)`
**TODO:** Do this.

## Inserting `I(n)` back into the sum
Recall that we are interested in computing
```
G2(x) = -1 / ((1 - x^p0) * log(x)) * sum((-1)^n * (1 + α)^n / factorial(n) * I(n) for n = 1:Inf)
```
With
```
I(n) = -log(x) * I1(n) - I2(n) + I3(n)
```
We can split `G2` into two sums as
```
G2(x) = -1 / ((1 - x^p0) * log(x)) * (
    - log(x) * sum((-1)^n * (1 + α)^n / factorial(n) * I1(n) for n = 1:Inf)
    - sum((-1)^n * (1 + α)^n / factorial(n) * I2(n) for n = 1:Inf)
)
```
We are now interested in computing these two sums.
- **TODO:** Add the sum coming from `I3(n)`.

Using the above approximations of `I1(n)` and `I2(n)` we get for the
first sum
```
sum((-1)^n * (1 + α)^n / factorial(n) * I1(n) for n = 1:Inf)

= sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^(n - 1) * (1 + γ)^(n - 2) * (n - log(π / x) * (1 + γ))
    for n = 1:Inf
)

= sum(
    (-1)^n * (1 + α)^n / factorial(n - 1) * log(π / x)^(n - 1) * (1 + γ)^(n - 2)
    for n = 1:Inf
)
- sum(
    (-1)^n * (1 + α)^n / factorial(n) * log(π / x)^(n) * (1 + γ)^(n - 1)
    for n = 1:Inf
)

= -(1 + α) / (1 + γ) * sum(
    (-1)^(n - 1) * (1 + α)^(n - 1) / factorial(n - 1) * log(π / x)^(n - 1) * (1 + γ)^(n - 1)
    for n = 1:Inf
)
- 1 / (1 + γ) * sum(
    (-1)^n * (1 + α)^n / factorial(n) * log(π / x)^(n) * (1 + γ)^n
    for n = 1:Inf
)

= -(1 + α) / (1 + γ) * (π / x)^(-(1 + α) * (1 + γ))
- 1 / (1 + γ) * ((π / x)^(-(1 + α) * (1 + γ)) - 1)

= (1 - (2 + α) * (π / x)^(-(1 + α) * (1 + γ))) / (1 + γ)
```
And for the second sum we get
```
sum((-1)^n * (1 + α)^n / factorial(n) * I2(n) for n = 1:Inf)

= sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^n * (1 + γ)^(n - 2) * ((n - 1) - n / (n + 1) * log(π / x) * (1 + γ))
    for n = 1:Inf
)

= sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^n * (1 + γ)^(n - 2) * (n - 1)
    for n = 1:Inf
)
- sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^n * (1 + γ)^(n - 2) * n / (n + 1) * log(π / x) * (1 + γ)
    for n = 1:Inf
)

= 1 / (1 + γ)^2 * sum(
    (-1)^n * (1 + α)^n / factorial(n) * log(π / x)^n * (1 + γ)^n * (n - 1) for n = 1:Inf
)
+ 1 / ((1 + α) * (1 + γ)^2) * sum(
    (-1)^(n + 1) * (1 + α)^(n + 1) / factorial(n + 1) * log(π / x)^(n + 1) * (1 + γ)^(n + 1) * n
    for n = 1:Inf
)

= 1 / (1 + γ)^2 * (1 - (π / x)^(-(1 + α) * (1 + γ)) * ((1 + α) * (1 + γ) * log(π / x) + 1))
+ 1 / ((1 + α) * (1 + γ)^2) * (1 - (π / x)^(-(1 + α) * (1 + γ)) * ((1 + α) * (1 + γ) * log(π / x) + 1))

= (2 + α) / ((1 + α) * (1 + γ)^2) *
    (1 - (π / x)^(-(1 + α) * (1 + γ)) * ((1 + α) * (1 + γ) * log(π / x) + 1))

= ((2 + α) * (1 - (π / x)^(-(1 + α) * (1 + γ)))) / ((1 + α) * (1 + γ)^2)
- (2 + α) / (1 + γ) * (π / x)^(-(1 + α) * (1 + γ)) * log(π / x)
```

If we let `s = (π / x)^(-(1 + α) * (1 + γ))` we can write `G2` as
```
G2(x) = -1 / ((1 - x^p0) * log(x)) * (
    - log(x) * (1 - (2 + α) * s) / (1 + γ)
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ)^2)
    + (2 + α) / (1 + γ) * s * (1 + α) * (1 + γ) * log(π / x)
)
```
- **TODO:** Add part from `I3(n)`
Writing `log(π / x) = log(π) - log(x)` and putting the log-terms
together, as well as factoring out `1 / (1 + γ)`, we get
```
G2(x) = -1 / ((1 - x^p0) * log(x) * (1 + γ) * (
    - log(x) * (1 - (2 + α) * s)
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ))
    + (2 + α) / (1 + γ) * s * (log(π) - log(x))
)

= -1 / ((1 - x^p0) * log(x) * (1 + γ)) * (
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ))
    + (2 + α) * s * log(π)
    - log(x) * (1 - (2 + α) * s)
    - log(x) * (2 + α) * s
)

= -1 / ((1 - x^p0) * log(x) * (1 + γ)) * (
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ))
    + (2 + α) * s * log(π)
    - log(x)
)
```
- **TODO:** Add part from `I3(n)`
Next we collect all occurrences of `s`, cancel the `log(x)` explicitly
and reorder the signs a bit
```
G2(x) = -1 / ((1 - x^p0) * log(x) * (1 + γ)) * (
    - (2 + α) / ((1 + α) * (1 + γ))
    + (2 + α) * (log(π) + 1 / ((1 + α) * (1 + γ))) * s
    - log(x)

= 1 / ((1 - x^p0) * (1 + γ)) * (
    + (2 + α) / ((1 + α) * (1 + γ)) / log(x)
    - (2 + α) * (log(π) + 1 / ((1 + α) * (1 + γ))) * s / log(x)
    + 1
)
```
- **TODO:** Add part from `I3(n)`

## Asymptotics in `α`
We now focus on determining the asymptotic behaviour of `G2(x)` as `α`
goes to `-1` to allow us to compute an enclosure for the full interval
of `α`.

To begin with we slightly rewrite `G2` as
```
G2(x) = 1 / ((1 - x^p0) * (1 + γ)) * (
    + (2 + α) / log(x) * ((1 - s) / (1 + α) / (1 + γ) - log(π) * s)
    + 1
)
```
The important parts are
1. `1 - x^p0`
2. `(1 - s) / (1 + α)`
3. `s`
Recall that `s = (π / x)^(-(1 + α) * (1 + γ))`, this give us
```
s = exp((1 + α) * (-(1 + γ) * log(π / x)))
  = sum((1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 0:Inf)
```
Looking at the factor
```
(2 + α) / log(x) * ((1 - s) / (1 + α) / (1 + γ) - log(π) * s) + 1
```
of `G2` we can write this as
```
(2 + α) / log(x) * (
    - 1 / (1 + γ) / (1 + α) * sum((1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 1:Inf)
    - log(π) * sum((1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 0:Inf)
) + 1

= (2 + α) / log(x) * (
    + log(π / x) * sum((1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n + 1) for n = 0:Inf)
    - log(π) * sum((1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 0:Inf)
) + 1

= (2 + α) / log(x) * sum(
        (log(π / x) / (n + 1) - log(π)) *
        (1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 0:Inf
) + 1
```
Extracting the first two terms in this sum gives us
```
= (2 + α) / log(x) * sum(
        (log(π / x) / (n + 1) - log(π)) *
        (1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 2:Inf
)
- (2 + α) / log(x) * (log(π / x) / 2 - log(π)) * (1 + α) * (1 + γ) * log(π / x)
+ (2 + α) / log(x) * (log(π / x) - log(π))
+ 1

= (2 + α) / log(x) * sum(
        (log(π / x) / (n + 1) - log(π)) *
        (1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 2:Inf
)
- (1 + γ) * (2 + α) / 2 * (log(x) - log(π)^2 / log(x)) * (1 + α)
+ (2 + α) / log(x) * (log(π / x) - log(π))
+ 1

= (2 + α) / log(x) * sum(
        (log(π / x) / (n + 1) - log(π)) *
        (1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 2:Inf
)
- (1 + γ) * (2 + α) / 2 * (log(x) - log(π)^2 / log(x)) * (1 + α)
- (1 + α)
```
Call the sum `Σ1`, i.e.
```
Σ1 = sum(
    (log(π / x) / (n + 1) - log(π)) *
    (1 + α)^n * (-(1 + γ) * log(π / x))^n / factorial(n) for n = 2:Inf
)
```

For `1 - x^p0` we have
```
1 - x^p0 = 1 - exp(p0 * log(x))
= -sum(p0^n * log(x)^n / factorial(n) for n = 1:Inf)
= -(1 + α) * (1 + (1 + α) / 2) * log(x) * sum(p0^n * log(x)^n / factorial(n + 1) for n = 0:Inf)
```
For the other factor of `G2`, `1 / ((1 - x^p0) * (1 + γ))`, this gives us
```
1 / ((1 - x^p0) * (1 + γ))
= -1 / ((1 + γ) * (1 + (1 + α) / 2) * log(x)) *
    1 / sum(p0^n * log(x)^n / factorial(n + 1) for n = 0:Inf) *
    1 / (1 + α)
```
Call the sum `Σ2`, i.e.
```
Σ2 = sum(p0^n * log(x)^n / factorial(n + 1) for n = 0:Inf)
```

For `G2` we then get
```
G2(x) = -1 / ((1 + γ) * (1 + (1 + α) / 2) * log(x)) * 1 / Σ2 * 1 / (1 + α) * (
    (2 + α) / log(x) * Σ1
    - (1 + γ) * (2 + α) / 2 * (log(x) - log(π)^2 / log(x)) * (1 + α)
    - (1 + α)
)

= -1 / ((1 + γ) * (1 + (1 + α) / 2) * Σ2) * (
    (2 + α) / log(x)^2 * Σ1 / (1 + α)
    - (1 + γ) * (2 + α) / 2 * (1 - log(π)^2 / log(x)^2)
    - 1 / log(x)
)
```
The problem now reduces to computing bounds for `Σ1 / (1 + α)` and
`Σ2` , where we will also need that `Σ2` is non-zero.

Unfortunately both `Σ1 / (1 + α)` and `Σ2` blow up as `x -> 0` so we
need cannot take the limit directly. For any fixed `x` which is not
too small, so that these terms are still relatively small, we can
compute explicit enclosures and from that get enclosures for the norm.

**IDEA:** Prove that it is increasing in `α` so that we only have to
bound it for an upper bound of `α`. Then we can fix `α` and do the
asymptotics in `x` only.
"""
function _T0_asymptotic_main(α::Arb, γ::Arb, c::Arb)
    # Construct function for computation of the term on the interval
    # [0, 1]
    G1 = let
        x -> zero(x)
    end

    # Construct function for computation of the term on the interval
    # [1, π / x]
    # FIXME: We currently do this for an upper bound of α
    G2(x) =
        let
            invlogx = if iszero(x)
                zero(x)
            elseif Arblib.contains_zero(x)
                Arb((inv(log(ubound(Arb, x))), 0))
            else
                inv(log(x))
            end

            # Compute enclosure of Σ1 / (1 + α)
            # TODO: Bound tail
            Σ1_div_onepα = let π = Arb(π)
                sum(
                    (log(π / x) / (n + 1) - log(π)) *
                    (1 + α)^(n - 1) *
                    (-(1 + γ) * log(π / x))^n / factorial(big(n)) for n = 2:10
                )
            end

            # Compute enclosure of Σ2
            # TODO: Bound tail
            Σ2 = let p0 = (1 + α) + (1 + α)^2 / 2
                sum(p0^n * log(x)^n / factorial(big(n + 1)) for n = 0:10)
            end

            res = let π = Arb(π)
                -1 / ((1 + γ) * (1 + (1 + α) / 2)) * 1 / Σ2 * (
                    (2 + α) * invlogx^2 * Σ1_div_onepα -
                    (1 + γ) * (2 + α) / 2 * (1 - log(π)^2 * invlogx^2) - invlogx
                )
            end

            return res
        end

    return x::Arb -> -sinpi(α / 2) * (G1(x) + G2(x))
end

"""
    _T0_asymptotic_remainder(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
W(x) * I_R(x) = -x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(x)) *
                ∫ abs(R(x * abs(1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
                    t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
defined in [`T0_asymptotic`](@ref).

# Bounding `W(x) * log(x)`
As a first step we take out the factor
```
W(x) * log(x) = x^(1 + α) / (gamma(1 + α) * (1 - x^p0))
```
and bound that separately.

- **TODO:** Finish bounding this

# Bound remainder terms
Next we want to compute a bound of the remainder terms from the
Clausen functions
```
abs(R(x * abs(1 - t)) + R(x * (1 + t)) - 2R(x * t))
```
More precisely we compute `C` such that
```
abs(R(x * abs(1 - t)) + R(x * (1 + t)) - 2R(x * t)) <= x^2 * C
```

From the expansion of the Clausen functions we get
```
R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t) = sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
Hence we have
```
R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t) = x^2 * sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^(2(m - 1)) * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
where we have removed the absolute value around `1 - t` since all
powers are even. We can thus get `C` by computing an upper bound of
the absolute value of this sum. This is implemented in the method
[`_bound_remainder_sum`](@ref).

# Tying it together
Using the above bound for the remainder terms we have
```
-1 / log(x) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
    t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt <=
C * x^2 / log(x) * ∫ t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
The integrand
```
t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t))
```
is increasing in `t`
- **PROVE:** That `t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t))` is
  increasing in `t`.
We can hence get an upper bound of the integrand by evaluating it
at `t = π / x`, for which we have
```
(π / x)^(1 - u0.γ * (1 + α)) * log(u0.c + inv(π))
```
and upper bound of the integral is thus given by
```
(π / x)^(1 - u0.γ * (1 + α)) * log(u0.c + inv(π)) * (π / x) =
    (π / x)^(2 - u0.γ * (1 + α)) * log(u0.c + inv(π)) =
```
Together with the factor in front this gives us
```
-C * x^2 / log(x) * (π / x)^(2 - u0.γ * (1 + α)) * log(u0.c + inv(π)) =
    -C * x^(u0.γ * (1 + α)) / log(x) * π^(2 - u0.γ * (1 + α)) * log(u0.c + inv(π)) <=
    -C * π^(2 - u0.γ * (1 + α)) * log(u0.c + inv(π)) / log(x)
```
where we have used that `x < 1` and `u0.γ * (1 + α) >= 0`.

- **IMPROVE:** This gives a fairly pessimistic bound of the remainder
  term, in particular for `x` not that small. Depending on where we
  want to start using the asymptotic version we might want to improve
  this. But it might not be needed.
"""
function _T0_asymptotic_remainder(α::Arb, γ::Arb, c::Arb)
    C = _bound_remainder_sum(α)

    return x::Arb -> begin
        # Enclosure of inv(log(x))
        invlogx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            Arb((inv(log(ubound(Arb, x))), 0))
        else
            inv(log(x))
        end

        # Enclosure of W(x) * log(x) = x^(1 + α) / (gamma(1 + α) * (1
        # - x^p0))
        # FIXME: Compute a proper upper bound of this
        Wxlogx =
            let αᵤ = ubound(Arb, α), p0ᵤ = (1 + αᵤ) + (1 + αᵤ)^2 / 2, xᵤ = ubound(Arb, x)
                xᵤ^(1 + αᵤ) / (gamma(1 + αᵤ) * (1 - x^p0ᵤ))
            end

        # Enclosure of integral
        I_R = C * Arb(π)^(2 - γ * (1 + α)) * log(c + inv(Arb(π)))

        # Enclosure of W(x) * I_R(x)
        Wx_I_R = -Wxlogx * invlogx * I_R

        return Wx_I_R
    end
end

"""
    T0_asymptotic(u0::BHKdVAnsatz{Arb}, ::Asymptotic)

Returns a function `f` such that `f(x)` computes an **upper bound** of
the integral \$T_0\$ from the paper using an evaluation strategy
that works asymptotically as `x` goes to `0`.

The integral in question is given by
```
1 / (π * u0.w(x) * u0(x)) *
    ∫ abs(clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * u0.w(y) dy
```
from `0` to `π`. Then change of variables `t = y / x` gives us
```
x / (π * u0.w(x) * u0(x)) *
    ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * u0.w(x * t) dt
```
from `0` to `π / x`. Using that
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

# Split into two factors
As a first step we split the above expression into two factors which
we bound separately.

The first factor is given by
```
F(x) = -inv(π) * log(x) / log(u0.c + inv(x)) * gamma(1 + α) * x^-α * (1 - x^p0) / u0(x)
```
Except for the addition of `-inv(π)` this is the same as the factor
`F1` in [`F0`](@ref) and we compute an upper bound following the same
procedure as described there.

The second factor is given by
```
x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(x)) *
    ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
where we will use the notation
```
W(x) = -x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(x))
I(x) = ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```

# Expand the integrand
We have the following expansions for the Clausen functions in the
integrand
```
clausenc(x * (1 - t), -α) = -sinpi(α / 2) * gamma(1 + α) * x^(-α - 1) * abs(1 - t)^(-α - 1) + R(x * abs(1 - t))
clausenc(x * (1 + t), -α) = -sinpi(α / 2) * gamma(1 + α) * x^(-α - 1) * (1 + t)^(-α - 1) + R(x * (1 + t))
clausenc(x * t, -α) = -sinpi(α / 2) * gamma(1 + α) * x^(-α - 1) * t^(-α - 1) + R(x * t)
```
where the error term `R` contains one constant term and `R(x * abs(1 -
t)) + R(x * (1 + t)) - 2R(x * t)` behaves like `O(x^2)`.

Inserting this into the integral allows us to split it into one main
integral
```
I_M(x) = -sinpi(α / 2) * gamma(1 + α) * x^(-α - 1) *
    ∫ abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
, where we have used that `-sinpi(α / 2) * gamma(1 + α) * x^(-α - 1)`
is positive to allow us to move it outside of the absolute value, and
one remainder integral
```
I_R(x) = ∫ abs(R(x * abs(1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
    t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
They satisfy `I(x) <= I_M(x) + I_R(x)`.

## Bounding the main integral
We need to bound `W(x) * I_M(x)`, this is implemented in
[`_T0_asymptotic_main`](@ref).

## Bounding the remainder integral
We need to bound `W(x) * I_R(x)`, this is implemented in
[`_T0_asymptotic_remainder`](@ref).
"""
function T0_asymptotic(
    u0::BHKdVAnsatz{Arb},
    evaltype::Asymptotic;
    non_asymptotic_u0 = false,
    ϵ::Arb = Arb(2e-1),
)
    # Enclosure of α
    α = Arb((-1, -1 + u0.ϵ))

    # Construct function for computing an upper bound of F(x)
    F = let
        # Function for computing an enclosure of -log(x) / log(u0.c +
        # inv(x))
        F1(x) =
            if iszero(x)
                one(x)
            elseif Arblib.contains_zero(x)
                xᵤ = ubound(Arb, x)
                -log(xᵤ) / log(u0.c + inv(xᵤ))
            else
                -log(x) / log(u0.c + inv(x))
            end

        # Function for computing an enclosure of gamma(1 + α) * x^-α * (1
        # - x^p0) / u0(x)
        F2 = inv_u0_bound(u0; ϵ)

        x -> F1(x) * F2(x) / π
    end

    # Function for computing an upper bound of main term
    Wx_I_M = _T0_asymptotic_main(α, u0.γ, u0.c)

    # Function for computing an upper bound of remainder term
    Wx_I_R = _T0_asymptotic_remainder(α, u0.γ, u0.c)

    return x::Arb -> F(x) * (Wx_I_M(x) + Wx_I_R(x))
end
