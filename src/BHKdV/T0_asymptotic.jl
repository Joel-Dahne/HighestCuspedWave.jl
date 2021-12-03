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


"""
function _T0_asymptotic_main(α::Arb, γ::Arb, c::Arb)
    return x::Arb -> begin
        Wx_I_M = zero(x)

        return Wx_I_M
    end
end

"""
    _T0_asymptotic_remainder(α::Arb, γ::Arb, c::Arb)

Compute an enclosure of
```
W(x) * I_R(x) = -x^(1 + α) / (gamma(1 + α) * (1 - x^p0) * log(x)) *
                ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
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
abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t))
```
More precisely we compute `C` such that
```
abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) <= x^2 * C
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
We can thus get `C` by computing an upper bound of the absolute value
of this sum. This is implemented in the method
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
the integral \$T_{0}\$ from the paper using an evaluation strategy
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
This is the expression we are interested in computing an *upper
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
clausenc(x * (t - 1), -α) = -sinpi(α / 2) * gamma(1 + α) * x^(-α - 1) * (t - 1)^(-α - 1) + R(x * (t - 1))
clausenc(x * (1 + t), -α) = -sinpi(α / 2) * gamma(1 + α) * x^(-α - 1) * (1 + t)^(-α - 1) + R(x * (1 + t))
clausenc(x * t, -α) = -sinpi(α / 2) * gamma(1 + α) * x^(-α - 1) * t^(-α - 1) + R(x * t)
```
where the error term `R` contains one constant term and `R(x * (1 -
t)) + R(x * (1 + t)) - 2R(x * t)` behaves like `O(x^2)`.

Inserting this into the integral allows us to split it into one main
integral
```
I_M(x) = -sinpi(α / 2) * gamma(1 + α) * x^(-α - 1) *
    ∫ abs((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
        t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
, where we have used that `-sinpi(α / 2) * gamma(1 + α) * x^(-α - 1)`
is positive to allow us to move it outside of the absolute value, and
one remainder integral
```
I_R(x) = ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
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
