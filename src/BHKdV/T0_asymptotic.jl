"""
    _T0_asymptotic_main_1(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G1(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫ abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
defined in [`T0_asymptotic`](@ref), where the integration is taken
from `0` to `1`. Using that `1 - t >= 0` the problem reduces to
computing an upper bound of
```
G1(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```

As `α -> -1` and `x -> 0` we have that `G1(x)` tends to zero and we
therefore don't need to compute a very accurate upper bound.

As `α` goes to `-1` the factor `(1 - t)^(-α - 1) + (t + 1)^(-α - 1) -
2t^(-α - 1)` of the integrand goes to zero, by factoring out `(1 + α)`
we make it instead go to a non-zero function. This means we want to
analyse
```
G1(x) = (1 + α) / ((1 - x^p0) * log(inv(x))) *
            ∫ abs((1 - t)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
The factor `(1 + α) / ((1 - x^p0) * log(x))` goes to zero and it is
therefore enough to compute an enclosure of the integral. Let
```
J = ∫ abs((1 - t)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) *
        t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```

## Simplifying the integrand
We have
```
(1 - t)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)
= (exp(-log(t - 1) * (1 + α)) + exp(-log(t + 1) * (1 + α)) - 2exp(-log(t) * (1 + α)))
= sum((-1)^n * (1 + α)^n * (log(t - 1)^n + log(t + 1)^n - 2log(t)^n) / factorial(n) for n = 1:Inf)
```
Taking out the first term we can write this as
```
-(1 + α) * (log(t - 1) + log(t + 1) - 2log(t)) +
sum((-1)^n * (1 + α)^n * (log(t - 1)^n + log(t + 1)^n - 2log(t)^n) / factorial(n) for n = 2:Inf)
```
Call the sum `Σ`, we can the split the integral as
```
J1 = ∫ abs((log(t - 1) + log(t + 1) - 2log(t))) * t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt

J2 = ∫ abs(Σ) / (1 + α) * t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
with `J <= J1 + J2`.

## Computing `J1`
Splitting the log-term in the weight as
```
log(c + inv(x * t)) = log((c * x * t + 1) / (x * t)) = log(1 + c * x * t) - log(x) - log(t)
```
we can split `J1` into three integrals
```
J11 = -log(x) * ∫ abs((log(t - 1) + log(t + 1) - 2log(t))) * t^(1 - γ * (1 + α)) dt

J12 = -∫ abs((log(t - 1) + log(t + 1) - 2log(t))) * t^(1 - γ * (1 + α)) * log(t) dt

J13 = ∫ abs((log(t - 1) + log(t + 1) - 2log(t))) * t^(1 - γ * (1 + α)) * log(1 + c * x * t) dt
```

## Computing `J2`


## Finishing
**TODO:** We would need to cancel the `log(x)` and then enclose `(1 +
  α) / (1 - x^p0)`.
"""
function _T0_asymptotic_main_1(α::Arb, γ::Arb, c::Arb)
    # FIXME: Implement this
    return x::Arb -> zero(x)
end

"""
    _T0_asymptotic_main_2(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G2(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫ abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
defined in [`T0_asymptotic`](@ref), where the integration is taken
from `1` to `π / x`. Using that `1 - t <= 0` and that
```
(t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)
```
is positive on the interval ``[1, ∞]`` (see
[`lemma_integrand_2`](@ref)), the problem reduces to computing an upper
bound of
```
G2(x) = 1 / ((1 - x^p0) * log(inv(x))) *
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
Next we cancel the `log(x)` explicitly and reorder the signs a bit
```
G2(x) = 1 / ((1 - x^p0) * (1 + γ)) * (
    + ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ) * log(x))
    - (2 + α) * s * log(π) / log(x)
    + 1
)
```
Factoring out `(2 + α) / log(x)` inside the last factor gives
```
G2(x) = 1 / ((1 - x^p0) * (1 + γ)) * (
    + (2 + α) / log(x) * (
        (1 - s) / ((1 + α) * (1 + γ))
        - log(π) * s
    ) + 1
)
```

## Fixed `x`
To begin with we consider the case when `x` is fixed and non-zero,
then the only asymptotics we have to do is in `α`.

By multiplying and dividing by `1 + α` we can write `G2(x)` as
```
G2(x) = (1 + α) / ((1 - x^p0) * (1 + γ)) * (
    + (2 + α) / log(x) * (
        (1 - s) / ((1 + α) * (1 + γ))
        - log(π) * s
    ) + 1
    ) / (1 + α)
```
and then split it into the two factors
```
G21(x) = (1 + α) / ((1 - x^p0) * (1 + γ))

G22(x) = ((2 + α) / log(x) * ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s) + 1) / (1 + α)
```
Our goal will be to enclose these two separately.

### Enclosing `G21(x)`
Expanding `1 - x^p0` as
```
1 - x^p0 = 1 - exp(p0 * log(x))
= -sum(p0^n * log(x)^n / factorial(n) for n = 1:Inf)
= -p0 * sum(p0^(n - 1) * log(x)^n / factorial(n) for n = 1:Inf)
```
we can write `G21(x)` as
```
G21(x) = -1 / (1 + γ) * (1 + α) / p0 * 1 / sum(p0^(n - 1) * log(x)^n / factorial(n) for n = 1:Inf)
```
Since `p0 = (1 + α) * (1 + (1 + α) / 2)` we get
```
G21(x) = -1 / ((1 + γ) * (1 + (1 + α) / 2)) * 1 / sum(p0^(n - 1) * log(x)^n / factorial(n) for n = 1:Inf)
```
and what remains is to enclose the sum. Since `log(x) < 0` the sum is
alternating. For large enough `n` it is also decreasing, more
precisely we have
```
abs(p0^(n - 1) * log(x)^n / factorial(n)) > abs(p0^n * log(x)^(n + 1) / factorial(n + 1))
⟺
1 > abs(p0 * log(x) / (n + 1))
⟺
n > p0 * abs(log(x)) - 1
```
So for every `n` greater than `p0 * abs(log(x)) - 1` the terms in the
sum are decreasing. We can thus explicitly sum the first `n` and get
an error which is bounded by the magnitude of the `n + 1`-th term. To
get a slightly better enclosure we can explicitly sum a few extra
terms after the `n`-th.

### Enclosing `G22(x)`
To be able to explicitly cancel the division by `1 + α` we rewrite
`G22(x)` as.
```
G22(x) = (
    (1 + α) / log(x) * ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s)
    + 1 / log(x) * ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s)
    + 1
) / (1 + α)

= ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s) / log(x)
+ (1 / log(x) * ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s) + 1) / (1 + α)
```
Now, using that `s = (x / π)^((1 + α) * (1 + γ))` and letting `q0 = (1
+ α) * (1 + γ)` we can write this as
```
G22(x) = ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) / log(x)
+ (1 / log(x) * ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) + 1) / (1 + α)
```
The second term we can rewrite as
```
(1 / log(x) * ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) + 1) / (1 + α)

= ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0 + log(x)) / ((1 + α) * log(x))

= ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0 + log(x / π) + log(π)) / ((1 + α) * log(x))

= ((1 + q0 * log(x / π) - (x / π)^q0) / q0 + log(π) * (1 - (x / π)^q0)) / ((1 + α) * log(x))
```
Giving us
```
G22(x) = ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) / log(x)
+ ((1 + q0 * log(x / π) - (x / π)^q0) / q0 + log(π) * (1 - (x / π)^q0)) / ((1 + α) * log(x))

= inv(log(x)) * (
    (1 - (x / π)^q0) / q0
    - log(π) * (x / π)^q0
    + (1 + q0 * log(x / π) - (x / π)^q0) / (q0 * (1 + α))
    + log(π) * (1 - (x / π)^q0) / (1 + α)
)
```
We will now bound each of these four terms separately. The term
`log(π) * (x / π)^q0` can be enclosed directly. For the other three
terms we make use of the expansion
```
(x / π)^q0 = sum(q0^n * log(x / π)^n / factorial(n) for n = 0:Inf)
```
For the three terms this gives us
```
(1 - (x / π)^q0) / q0 = -sum(q0^(n - 1) * log(x / π)^n / factorial(n) for n = 1:Inf)

(1 + q0 * log(x / π) - (x / π)^q0) / (q0 * (1 + α)) =
    -sum(q0^(n - 1) * log(x / π)^n / factorial(n) for n = 2:Inf) / (1 + α) =
    -(1 + γ) * sum(q0^(n - 2) * log(x / π)^n / factorial(n) for n = 2:Inf) =

log(π) * (1 - (x / π)^q0) / (1 + α) =
    -log(π) * sum(q0^n * log(x / π)^n / factorial(n) for n = 1:Inf) / (1 + α) =
    -(1 + γ) * log(π) * sum(q0^(n - 1) * log(x / π)^n / factorial(n) for n = 1:Inf)
```
Note that the first and last sum are identical. Similar to the sum in
`G21(x)` all of these sums are alternating and to compute enclosures
it is enough to find when the terms start to decrease. In this case we
get
```
abs(q0^(n - 1) * log(x / π)^n / factorial(n)) > abs(q0^n * log(x / π)^(n + 1) / factorial(n + 1))
⟺
1 > abs(q0 * log(x / π) / (n + 1))
⟺
n > q0 * abs(log(x / π)) - 1
```
for the first sum and a similar calculation gives the same for the
second sum. So for every `n` greater than `q0 * abs(log(x / π)) - 1`
the terms in both sums are decreasing. We can thus explicitly sum the
first `n` and get an error which is bounded by the magnitude of the `n
+ 1`-th term. To get a slightly better enclosure we can explicitly sum
a few extra terms after the `n`-th.

## Monotonicity in `x`
Above we showed how to compute an enclosure of `G2(x)` for a fixed
non-zero `x`. We will now show a certain type of monotonicity in `x`,
so that it is enough to work with an upper bound of `x`. We want to
show that
```
G2(x) < 1 / (1 + γ)
```
We can rewrite this as
```
(2 + α) * ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) + log(x) >= (1 - x^p0) * log(x)
```
where we have multiplied both sides by `(1 + γ) * (1 - x^p0) *
log(x)`, since this is negative we also reversed the inequality. If we
move everything to one side we get
```
(2 + α) * ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) + x^p0 * log(x) >= 0
```
We now want to show that the left hand side is decreasing in `x`, so
that it is enough to check the inequality for an upper bound of `x`.

Differentiating the left hand side we get
```
(2 + α) * (-x^(q0 - 1) / π^q0 - q0 * log(π) * x^(q0 - 1) / π^q0) + x^(p0 - 1) * (p0 * log(x) + 1)
```
Grouping terms we have
```
x^(p0 - 1) * (
    - (2 + α) * x^(q0 - p0) / π^q0
    - (2 + α) * q0 * log(π) * x^(q0 - p0) / π^q0
    + p0 * log(x)
    + 1
)
```
It is therefore enough to show that
```
(
    - (2 + α) * x^(q0 - p0) / π^q0
    - (2 + α) * q0 * log(π) * x^(q0 - p0) / π^q0
    + p0 * log(x)
    + 1
)
```
is negative. We will prove that this is increasing in `x` and
non-positive at `x = 1`. At `x = 1` we get
```
(
    - (2 + α) / π^q0
    - (2 + α) * q0 * log(π) / π^q0
    + 1
) =
1 - (2 + α) * (1 + q0 * log(π)) / π^q0
```
which we want to prove is upper bounded by `0`. To prove move one of
the terms to the other side and multiply both sides by `π^q0` and use
the expansion `π^q0 = sum(q0^n * log(π)^n / factorial(n) for n =
0:Inf)` to get
```
sum(q0^n * log(π)^n / factorial(n) for n = 0:Inf) <= (2 + α) * (1 + q0 * log(π))
⟺
sum(q0^n * log(π)^n / factorial(n) for n = 2:Inf) <= (1 + α) * (1 + q0 * log(π))
⟺
(1 + α) * (1 + γ) * log(π) * sum(q0^(n - 1) * log(π)^n / factorial(n) for n = 2:Inf) <=
    (1 + α) * (1 + q0 * log(π))
⟺
(1 + γ) * log(π) * sum(q0^(n - 1) * log(π)^n / factorial(n) for n = 2:Inf) <=
    1 + q0 * log(π)
```
It therefore suffices to show that
```
sum(q0^(n - 1) * log(π)^n / factorial(n) for n = 2:Inf) <= 1 / (1 + γ) * log(π)
```
Since the sum is increasing in `q0` we only have to show this for an
upper bound of `q0`. The sum can be computed explicitly to be
```
sum(q0^(n - 1) * log(π)^n / factorial(n) for n = 2:Inf) = (π^q0 - q0 * log(π) - 1) / q0
```
and since we now know that this is increasing in `q0` we check that
```
(π^q0 - q0 * log(π) - 1) / q0 <= 1 / (1 + γ) * log(π)
```
holds for an upper bound of `q0`.

For the derivative we get
```
(
    - (q0 - p0) * (2 + α) * x^(q0 - p0 - 1) / π^q0
    - (q0 - p0) * (2 + α) * q0 * log(π) * x^(q0 - p0 - 1) / π^q0
    + p0 / x
)
```
Multiplying by `x`, which doesn't change the sign, and reordering we
get
```
p0 - (q0 - p0) * (2 + α) * (1 + q0 * log(π)) / π^q0 * x^(q0 - p0)
```
Since `q0 - p0 > 0` and hence `0 <= x^(q0 - p0) <= 1`, it is enough to
show
```
p0 - (q0 - p0) * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
Factoring out `1 + α` and cancelling it we get
```
(1 + (1 + α) / 2) - ((1 + γ) - (1 + (1 + α) / 2)) * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
⟺
1 + (1 + α) / 2 - (γ - (1 + α) / 2) * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
Since `(1 + α) / 2` is positive we can remove it, giving us
```
1 - (γ - (1 + α) / 2) * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
Similarly we can replace `(γ - (1 + α) / 2)` with `γ` without losing
anything.
```
1 - γ * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
Finally `γ < 1` so it is enough to show
```
1 - (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
which is exactly what we already had above for the `x = 1` case.
"""
function _T0_asymptotic_main_2(α::Arb, γ::Arb, c::Arb)
    # We use 1 + α in many places and this is assumed to be
    # non-negative, we therefore compute such an enclosure
    onepα = Arblib.nonnegative_part!(zero(α), 1 + α)

    # TODO: Add the two remainder terms which we for now ignore.
    return x::Arb -> begin
        # This function assumes that x is less than or equal to 1
        @assert x <= 1

        if Arblib.contains_zero(x)
            # In this case we prove that G2 is bounded by 1 / (1 +
            # γ). By above we only have to check two things
            # 1. 1 - (2 + α) * (1 + q0 * log(π)) / π^q0 <= 0
            # 2. 0 <= p0 - (q0 - p0) * (2 + α) * (1 + q0 * log(π)) / π^q0

            # 1. This was reduced to checking that
            # (π^q0 - q0 * log(π) - 1) / q0 <= 1 / (1 + γ) * log(π)
            # holds for an upper bound of q0
            lemma1 = let q0 = ubound(Arb, (1 + α) * (1 + γ)), π = Arb(π)
                (π^q0 - q0 * log(π) - 1) / q0 <= 1 / (1 + γ) * log(π)
            end

            # 2. This was reduced to checking the same inequality
            # as for 1. and using that γ <= 1
            lemma2 = γ <= 1 && lemma1

            if lemma1 && lemma2
                return 1 / (1 + γ)
            else
                return Arblib.indeterminate!(zero(α))
            end
        else
            # Compute G21
            G21 = let
                # Enclosure of sum(p0^(n - 1) * log(x)^n / factorial(n) for n = 1:Inf)
                G21_sum = let p0 = onepα + onepα^2 / 2
                    # N such that the terms are decreasing
                    N = Int(Arblib.ceil!(zero(Arf), ubound(p0 * abs(log(x)) - 1)))
                    # IMPROVE: If we need to go to very small x we
                    # will have to take more extra terms
                    extra = 10 # Take 10 extra terms
                    term(n) = p0^(n - 1) * log(x)^n / factorial(big(n))
                    Σ = sum(term(n) for n = 1:N+extra)
                    # Error is bounded by magnitude of next term
                    Arblib.add_error!(Σ, term(N + extra + 1))
                    Σ
                end

                G21 = -1 / ((1 + γ) * (1 + onepα / 2) * G21_sum)
            end

            G22 = let q0 = onepα * (1 + γ), π = Arb(π)
                # Compute N such that the terms in the sum are
                # decreasing for n >= N. This is the same number
                # of all sums
                N = Int(Arblib.ceil!(zero(Arf), ubound(q0 * abs(log(x / π)) - 1)))

                # IMPROVE: If we need to go to very small x we
                # will have to take more extra terms
                extra = 10 # Take 10 extra terms

                # Compute term1 and term4 since they share the same sum
                term1, term4 = let
                    # Enclosure of sum(q0^(n - 1) * log(x / π)^n / factorial(n) for n = 1:Inf)
                    term14_sum = let
                        term(n) = q0^(n - 1) * log(x / π)^n / factorial(big(n))
                        Σ = sum(term(n) for n = 1:N+extra)
                        # Error is bounded by magnitude of next term
                        Arblib.add_error!(Σ, term(N + extra + 1))
                        Σ
                    end

                    -term14_sum, -(1 + γ) * log(π) * term14_sum
                end

                # This we enclose directly
                term2 = -log(π) * abspow(x / π, q0)

                term3 = let
                    term3_sum = let
                        term(n) = q0^(n - 2) * log(x / π)^n / factorial(big(n))
                        Σ = sum(term(n) for n = 2:N+extra)
                        # Error is bounded by magnitude of next term
                        Arblib.add_error!(Σ, term(N + extra + 1))
                        Σ
                    end

                    -(1 + γ) * term3_sum
                end

                (term1 + term2 + term3 + term4) / log(x)
            end

            return G21 * G22
        end
    end
end

"""
    _T0_asymptotic_remainder(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
R(x) = ∫ abs(R(x * abs(1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
            t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
defined in [`T0_asymptotic`](@ref).

From the expansion of the Clausen functions we have
```
R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t) = sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
where we have removed the absolute value around `1 - t` since all
powers are even. We have that
```
(1 - t)^2m + (1 + t)^2m - 2t^2m =
    sum(binomial(2m, k) * (1 + (-1)^k) * t^k for k = 0:2m) - 2t^2m =
    2sum(binomial(2m, 2k) * t^2k for k = 0:m-1)
```
Giving us
```
R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t) = 2sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * sum(binomial(2m, 2k) * t^2k for k = 0:m-1) for m = 1:Inf
)
```
The factor `zeta(-α - 2m) * (-1)^m` is positive since `0 < -α < 1` and
`m >= 1`. All terms are therefore positive and adding the integrating
we can write it as
```
R(x) = 2sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m *
        sum(binomial(2m, 2k) * ∫ t^(2k + 1 - γ * (1 + α)) * log(c + inv(x * t)) dt for k = 0:m-1) for m = 1:Inf
)
```
where the integration is taken from `0` to `π / x`. For the integral
we first split the log-factor as
```
log(1 + inv(x * t)) = log(1 + c * x * t) - log(x * t)
```
giving us
```
∫ t^(2k + 1 - γ * (1 + α)) * log(c + inv(x * t)) dt =
    ∫ t^(2k + 1 - γ * (1 + α)) * log(1 + c * x * t) dt
    - ∫ t^(2k + 1 - γ * (1 + α)) * log(x * t) dt
```
Using that `log(1 + c * x * t)` is bounded by `log(1 + c * π)` and
integrating explicitly gives us
```
∫ t^(2k + 1 - γ * (1 + α)) * log(1 + c * x * t) dt <=
    log(1 + c * π) * ∫ t^(2k + 1 - γ * (1 + α)) dt =
    log(1 + c * π) * (π / x)^(2(k + 1) - γ * (1 + α)) / (2(k + 1) - γ * (1 + α)) <=
    log(1 + c * π) * (π / x)^(2(k + 1)) / (2(k + 1) - γ * (1 + α))
```
```
∫ t^(2k + 1 - γ * (1 + α)) * log(x * t) dt =
    ((2(k + 1) - γ * (1 + α)) * log(π) - 1) * (π / x)^(2(k + 1) - γ * (1 + α)) / (2(k + 1) - γ * (1 + α))^2
    <= log(π) * (π / x)^(2(k + 1)) / (2(k + 1) - γ * (1 + α))
```
And hence
```
∫ t^(2k + 1 - γ * (1 + α)) * log(c + inv(x * t)) dt <=
    log(c + inv(π)) * (π / x)^(2(k + 1)) / (2(k + 1) - γ * (1 + α))
```
Inserting this back into `R(x)` and factoring out `(π / x)^2m` we get
```
R(x) <= 2log(c + inv(π)) * sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * π^2m *
        sum(binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α)) for k = 0:m-1) for m = 1:Inf
)
```
We split this sum into `m = 1:M-1` and `m = M:Inf`. For the finite
part we sum directly. For the infinite tail we note that
```
sum(binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α)) for k = 0:m-1) <=
    sum(binomial(2m, 2k) * (1 / π)^(2(m - 1 - k)) / (2k + 1) for k = 0:m-1) <=
    sum(binomial(2m, 2k) * (1 / π)^(2(m - 1 - k)) / (2k + 1) for k = 0:m) =
    π^(2 - 2m) * ((π - 1)^2m + (π + 1)^2m) / 2 <=
    π^(2 - 2m) * (π + 1)^2m <=
    = π^2 * (1 + inv(π))^2m
```
Giving us
```
sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * π^2m *
        sum(binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α)) for k = 0:m-1) for m = M:Inf
) <=
π^2 * sum(zeta(-α - 2m) * (-1)^m / factorial(2m) * (π * (1 + inv(π)))^2m for m = M:Inf)
```
This sum is exactly the same as in the remainder for `clausenc` in
[`clausenc_expansion_remainder`](@ref), substituting `x` for `π * (1 +
inv(π))`, and we can get a bound from that function.
"""
function _T0_asymptotic_remainder(α::Arb, γ::Arb, c::Arb; M = 20)
    # Precompute the factor
    # (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m))
    # in the finite sum
    factors = [(-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) for m = 1:M-1]

    # Compute tail of sum, which is the same as for clausenc_expansion
    # It is valid for x < 1
    tail = let y = π * (1 + inv(Arb(π)))
        Arb(π)^2 * clausenc_expansion_remainder(y, -α, M) * y^2M
    end

    return x::Arb -> begin
        x < 1 || throw(DomainError(x, "x must be less than 1"))

        # Sum directly for m = 1:M-1
        main = sum(1:M-1, init = zero(x)) do m
            factors[m] * sum(
                binomial(2m, 2k) * abspow(x / π, 2(m - 1 - k)) / (2(k + 1) - γ * (1 + α)) for k = 0:m-1
            )
        end

        return 2log(c + inv(Arb(π))) * (main + tail)
    end
end

"""
    T0_asymptotic(u0::BHKdVAnsatz{Arb}, ::Asymptotic)

Returns a function `f` such that `f(x)` computes an **upper bound** of
the integral ``T_0`` from the paper using an evaluation strategy
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

# Split into three factors
Similarly to in the asymptotic version of [`F0`](@ref) we split the
expression into three factors which we bound separately. We write it
as
```
(gamma(1 + α) * x^(-α) * (1 - x^p0) / (π * u0(x)))
* (log(inv(x)) / log(u0.c + inv(x)))
* (
    x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(inv(x))) *
    ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
)
```
Except for the addition of `1 / π` the first two factors are the same
as in [`F0`](@ref) and we handle them in the same way. For the third
factor we use the notation
```
W(x) = x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(inv(x)))
I(x) = ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```

# Expand the integrand
We have the following expansions for the Clausen functions in the
integrand
```
clausenc(x * (1 - t), -α) = sinpi(-α / 2) * gamma(1 + α) * x^(-α - 1) * abs(1 - t)^(-α - 1) + R(x * abs(1 - t))
clausenc(x * (1 + t), -α) = sinpi(-α / 2) * gamma(1 + α) * x^(-α - 1) * (1 + t)^(-α - 1) + R(x * (1 + t))
clausenc(x * t, -α) = sinpi(-α / 2) * gamma(1 + α) * x^(-α - 1) * t^(-α - 1) + R(x * t)
```
where the error term `R` contains one constant term and `R(x * abs(1 -
t)) + R(x * (1 + t)) - 2R(x * t)` behaves like `O(x^2)`.

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
I_R(x) = x^(1 + α) * ∫ abs(R(x * abs(1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
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

R(x) = ∫ abs(R(x * abs(1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
            t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
We have
```
W(x) * I(x) <= sinpi(-α / 2) * (G1(x) + G2(x)) +
    x^(1 + α) / (gamma(1 + α) * log(inv(x)) * (1 - x^p0)) * R(x)
```
Bounds of the functions `G1(x), G2(x), R(x)` are implemented in
[`_T0_asymptotic_main_2`](@ref), [`_T0_asymptotic_main_1`](@ref) and
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
function T0_asymptotic(
    u0::BHKdVAnsatz{Arb},
    evaltype::Asymptotic;
    non_asymptotic_u0 = false,
    ϵ::Arb = Arb(2e-1),
)
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

    # Enclosure of sin(-α / 2)
    G_factor = sinpi(-α / 2)

    # Function for computing an upper bound of G1, G2 and R
    # respectively
    G1 = _T0_asymptotic_main_1(α, u0.γ, u0.c)
    G2 = _T0_asymptotic_main_2(α, u0.γ, u0.c)

    # Function for computing enclosure of
    # x^(1 + α) / (gamma(1 + α) * log(inv(x)) * (1 - x^p0))
    R_factor =
        x -> begin
            # Enclosure of inv(log(inv(x))) = -inv(log(x))
            invlogx = if iszero(x)
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

            return abspow(x, αp1) * invlogx * inv_gamma_1mxp0
        end

    R = _T0_asymptotic_remainder(α, u0.γ, u0.c)

    return x::Arb -> f1(x) * f2(x) / π * (G_factor * (G1(x) + G2(x)) + R_factor(x) * R(x))
end
