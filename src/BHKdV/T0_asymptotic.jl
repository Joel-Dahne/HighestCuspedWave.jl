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

`G1(x)` tends to zero as `α -> -1` and `x -> 0` and we therefore don't
need to compute a very accurate upper bound.

As `α` goes to `-1` the factor `(1 - t)^(-α - 1) + (t + 1)^(-α - 1) -
2t^(-α - 1)` of the integrand goes to zero, by factoring out `(1 + α)`
we make it instead go to a non-zero function. This means we want to
analyse
```
G1(x) = -(1 + α) / ((1 - x^p0) * log(x)) *
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
and what remains is to enclose the sum. Since `log(x) < 1` the sum is
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
We will not bound each of these four terms separately. The term
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
function _T0_asymptotic_main(α::Arb, γ::Arb, c::Arb)
    # This function assumes that x is less than or equal to 1
    @assert x <= 1

    # We use 1 + α in many places and this is assumed to be
    # non-negative, we therefore compute such an enclosure
    onepα = Arblib.nonnegative_part!(zero(α), 1 + α)

    # Construct function for computation of the term on the interval
    # [0, 1]
    G1 = let
        x -> zero(x)
    end

    # Construct function for computation of the term on the interval
    # [1, π / x]
    # TODO: Add the two remainder terms which we for now ignore.
    G2(x) =
        let
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
