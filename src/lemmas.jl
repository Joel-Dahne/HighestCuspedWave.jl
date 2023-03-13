"""
    lemma_kdv_asymptotic_expansion

Corresponds to Lemma 4.2 in the paper.

# Statement

TODO
"""
function lemma_kdv_asymptotic_expansion end

"""
    lemma_kdv_asymptotic_expansion_defect

Corresponds to Lemma 4.3 in the paper.

# Statement

TODO
"""
function lemma_kdv_asymptotic_expansion_defect end

"""
    lemma_bhkdv_main_term_limit

Corresponds to Lemma 4.5 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_main_term_limit end

"""
    lemma_bhkdv_asymptotic_expansion

Corresponds to Lemma 4.6 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_asymptotic_expansion end

"""
    lemma_kdvzero_as_asymptotics

Corresponds to Lemma 4.7 in the paper.

# Statement

TODO
"""
function lemma_kdvzero_as_asymptotics end

"""
    lemma_kdvzero_asymptotic_expansion

Corresponds to Lemma 4.8 in the paper.

# Statement

TODO
"""
function lemma_kdvzero_asymptotic_expansion end

"""
    lemma_inv_u0_normalised

Corresponds to Lemma 8.1 in the paper.

# Statement

TODO
"""
function lemma_inv_u0_normalised end

"""
    lemma_bhkdv_inv_u0_normalised

Corresponds to Lemma 8.2 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_inv_u0_normalised end

"""
    lemma_bhkdv_weight_factor_enclosure

Corresponds to Lemma 9.1 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_weight_factor_enclosure end

"""
    lemma_kdvzero_F0_expansion

Corresponds to Lemma 10.1 in the paper.

# Statement

TODO
"""
function lemma_kdvzero_F0_expansion end

"""
    lemma_I_hat_root_zero

Corresponds to Lemma 11.1 in the paper.

# Statement

TODO
"""
function lemma_I_hat_root_zero end

"""
    lemma_I_hat_root

Corresponds to Lemma 11.2 in the paper.

# Statement

TODO
"""
function lemma_I_hat_root end

"""
    lemma_I_hat_root_alpha

Corresponds to Lemma 11.3 in the paper.

# Statement

TODO
"""
function lemma_I_hat_root_alpha end

"""
    lemma_I_positive

Corresponds to Lemma 11.4 in the paper.

# Statement

TODO
"""
function lemma_I_positive end

"""
    lemma_U0_primitive_weight_x

Corresponds to Lemma 11.5 in the paper.

# Statement

TODO
"""
function lemma_U0_primitive_weight_x end

"""
    lemma_kdv_U0_asymptotic

Corresponds to Lemma 11.6 in the paper.

# Statement

TODO
"""
function lemma_kdv_U0_asymptotic end

"""
    lemma_bhkdv_U0_asymptotic_split

Corresponds to Lemma 11.7 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_T0_asymptotic_split end

"""
    lemma_kdvzero_U0_constant_term

Corresponds to Lemma 11.8 in the paper.

# Statement

TODO
"""
function lemma_kdvzero_U0_constant_term end

"""
    lemma_kdvzero_U0_primitive_K_expansion

Corresponds to Lemma 8.2 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_inve_u0_normalised end

"""
    lemma_removable_singularities

Corresponds to Lemma A.1 in the paper.

# Statement

TODO
"""
function lemma_removable_singularities end

"""
    lemma_clausenc_monotone

Corresponds to Lemma C.1 in the paper.

# Statement
For all real `s` the Clausen function `clausenc(x, s)` is monotone in
`x` for `0 < x < π`. For `s > 0` it is non-increasing. For `s <= 0`
the sign of the derivative is the same as that of `-sinpi(s / 2)`.
"""
function lemma_clausenc_monotone end

"""
    lemma_clausens_monotone

Corresponds to Lemma C.2 in the paper.

# Statement
For `s <= 1` the Clausen function `clausens(x, s)` is monotone in
`x` for `0 < x < 2π`. The sign of the derivative is the same as that
of `-cospi(s / 2)`.
"""
function lemma_clausens_monotone end

"""
    lemma_clausencmzeta_monotone

Corresponds to Lemma C.3 in the paper.

# Statement
For `s > 1` and `x` real the function `clausencmzeta(x, s)` is
non-decreasing in `s`.
"""
function lemma_clausencmzeta_monotone end

"""
    lemma_clausenc_expansion_odd_integer

Corresponds to Lemma C.4 in the paper.

# Statement

TODO
"""
function lemma_clausenc_expansion_odd_integer end

"""
    lemma_x_pow_t_m1_div_t

Corresponds to Lemma C.5 in the paper.

# Statement

TODO
"""
function lemma_absxt_m1_div_t end

"""
    lemma_clausen_remainder

Corresponds to Lemma C.6 in the paper.

# Statement
Let `s >= 0, 2M >= s + 1` and `abs(x) < 2π`, we then have the
following bounds for the tails in
[`equation_clausenc_expansion`](@ref) and
[`equation_clausens_expansion`](@ref)
```
abs(sum((-1)^m * zeta(s - 2m) * x^2m / factorial(2m) for m = M:Inf)) <=
    2(2π)^(1 + s - 2M)abs(sinpi(s / 2)) * zeta(2M + 1 - s) * x^2M / (4π^2 - x^2)

abs(sum((-1)^m * zeta(s - 2m - 1) * x^(2m + 1) / factorial(2m + 1) for m = M:Inf)) <=
    2(2π)^(s - 2M)abs(cospi(s / 2)) * zeta(2M + 2 - s) * x^(2M + 1) / (4π^2 - x^2)
```
"""
function lemma_clausen_remainder end

"""
    lemma_clausen_derivative_remainder

Corresponds to Lemma C.7 in the paper.

# Statement
**TODO:** Check that this corresponds to the version in the paper.

Let `d(f, k)(x)` denote the `k`th derivative of `f` evaluated at `x`.

Let `β >= 1, s >= 0, 2M >= s + 1` and `abs(x) < 2π`, we then have the
following bounds for the tails in ... and ...
```
abs(sum((-1)^m * d(zeta, β)(s - 2m) * x^2m / factorial(2m) for m = M:Inf)) <=
    sum(j1 + j2 + j3 = β) do j1, j2, j3
        binomial(β, (j1, j2, j3)) *
        2(log(2π) + π / 2)^j1 *
        (2π)^(s - 1) *
        abs(d(zeta, j3)(1 + 2M - s)) *
        sum(p[j2](1 + 2m - s) * (x / 2π)^2m for m = M:Inf)
    end

abs(sum((-1)^m * d(zeta, β)(s - 2m - 1) * x^(2m + 1) / factorial(2m + 1) for m = M:Inf)) <=
    sum(j1 + j2 + j3 = β) do j1, j2, j3
        binomial(β, (j1, j2, j3)) *
        2(log(2π) + π / 2)^j1 *
        (2π)^(s - 1) *
        abs(d(zeta, j3)(2 + 2M - s)) *
        sum(p[j2](2 + 2m - s) * (x / 2π)^(2m + 1) for m = M:Inf)
    end
```

Here `p[j2]` is recursively given by
```
p[k + 1](s) = polygamma(0, s) * p[k](s) + d(p[k], 1)(s)
p[0] = 1
```
It is given by a linear combination of terms of the form
```
polygamma(0, s)^q[0] * polygamma(1, s)^q[1] * ⋯ * polygamma(j2 - 1, s)^q[j2 - 1]
```
We have the following bounds
```
sum(
    polygamma(0, 1 + 2m - s)^q[0] * ⋯ * polygamma(j2 - 1, 1 + 2m - s)^q[j2 - 1] (x / 2π)^2m
    for m = M:Inf
) <= abs(polygamma(1, 1 + 2M - s)^q[1] * ⋯ * polygamma(j2 - 1, 1 + 2M - s)^q[j2 - 1]) *
    1 / 2^(q[0] / 2) * (2π)^(-2M) * lerch_phi(x^2 / 4π^2, -q[0] / 2, M + 1) * x^2M

sum(
    polygamma(0, 2 + 2m - s)^q[0] * ⋯ * polygamma(j2 - 1, 2 + 2m - s)^q[j2 - 1] (x / 2π)^(2m + 1)
    for m = M:Inf
) <= abs(polygamma(1, 2 + 2M - s)^q[1] * ⋯ * polygamma(j2 - 1, 2 + 2M - s)^q[j2 - 1]) *
    1 / 2^(q[0] / 2) * (2π)^(-2M - 1) * lerch_phi(x^2 / 4π^2, -q[0] / 2, M + 1) * x^(2M + 1)
```
"""
function lemma_clausen_derivative_remainder end

"""
    lemma_clausen_derivative_remainder_interval

Corresponds to Lemma C.9 in the paper.

# Statement

TODO
"""
function lemma_clausen_derivative_remainder_interval end

"""
    lemma_bhkdv_integrand_enclosure_zero

Corresponds to Lemma D.1 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_integrand_enclosure_zero end

"""
    lemma_bhkdv_U0_G1

Corresponds to Lemma F.1 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_U0_G1 end

"""
    lemma_bhkdv_U0_G2_split

Corresponds to Lemma F.2 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_U0_G2_split end

"""
    lemma_bhkdv_U0_G2_M

Corresponds to Lemma F.3 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_U0_G2_M end

"""
    lemma_bhkdv_U0_G2_R_n1

Corresponds to Lemma F.4 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_U0_G2_R_n1 end

"""
    lemma_bhkdv_U0_G2_R_1

Corresponds to Lemma F.5 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_U0_G2_R_1 end

"""
    lemma_bhkdv_U0_G2_R_2

Corresponds to Lemma F.6 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_U0_G2_R_2 end

"""
    lemma_bhkdv_U0_G21

Corresponds to Lemma F.7 in the paper.

# Statement

TODO
"""
function lemma_bhkdv_U0_G21 end

"""
    lemma_integrand_1()

# Statement
Let
```
I(α, x, y) = clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)
```
For every `α` in ``[-1, 0)`` and `x` in ``(0, π)`` the function `I(α,
x, y)` is increasing and continuous in `y` for `y` in ``(0, x)`` and
has the limit `-Inf` at `y = 0` and `Inf` at `y = x`.

Moreover, the unique root, `r(α, x)`, in `y` satisfies the inequality
```
x / 2 <= r(α, x) <= B(α) * x
```
and `hat_r(α, x) = r(α, x) / x` is decreasing in `x`.

This can also be stated for
```
hat_I(α, x, y) = clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
on the interval ``(0, 1)`` in `t`. In which case the root is given by
`hat_r(α, x)`.

Taking the limit as `x` goes to zero we can expand in `x` and get a
similar result. For `α` in the interval ``(-1, 0)`` we get
```
(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
```
which has the unique root `hat_r(α, 0)`. For `α = -1` we instead get
```
log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2) =
    log(1 / t^2 - 1)
```
with the unique root `hat_r(-1, 0) = 1 / sqrt(2)`.

# Proof
For a proof see the corresponding lemma in the paper.
"""
lemma_integrand_1() = true

"""
    lemma_integrand_2()

# Statement
Let
```
I(α, x, y) = clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)
```
For every `α` in ``[-1, 0)``, `x` in ``(0, π)`` and `y` in ``(x, π)``
the function `I(α, x, y)` is positive.

This can also be stated for
```
hat_I(α, x, y) = clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
on the interval ``(1, π / x)`` in `t`.

Taking the limit as `x` goes to zero we can expand in `x` and get a
similar result. Let `t` be in the interval ``(1, π / x)``, for `α` in
the interval ``(-1, 0)`` we get that
```
(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
```
is positive and for `α = -1`
```
log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2) =
    log(1 / t^2 - 1)
```
is positive.

# Proof
For a proof see the corresponding lemma in the paper.
"""
lemma_integrand_2() = true

"""
    lemma_bhkdv_monotonicity_alpha(u0::BHKdVAnsatz)

# Statement
The function `u0(x)` is **lower bounded** by `u0.v0(x)` for all `x ∈
[0, π]`. This holds for all `u0.ϵ ∈ (0, 1 / 4)`.

# Use
This is used when upper bounding terms where `u0` occurs in the
denominator since we can then replace it by `u0.v0` and still get an
upper bound. For example this is used in
- [`n0_bound`](@ref)

# Proof
The only difference between `u0` and `u0.v0` is the leading term. The
Clausen terms in the tail as well as the Fourier terms are exactly the
same. It hence suffices to prove the inequality for the main term.

The main term for `u0` is for `α ∈ [-1, -1 + u0.ϵ]` given by
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
where `p0 = 1 + α + (1 + α)^2 / 2` and
```
a0 = finda0(α) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))^2
```
The main term for `u0.v0` is
```
2 / π^2 * clausencmzeta(x, 2, 1)
```

That
```
2 / π^2 * clausencmzeta(x, 2, 1) <=
    a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
for all `α ∈ [-1, -1 + ϵ]` and `x ∈ [0, π]` follows easily from
[`lemma_bhkdv_monotonicity_main_term`](@ref) and
[`lemma_bhkdv_limit_main_term`](@ref).
"""
lemma_bhkdv_monotonicity_alpha(u0::BHKdVAnsatz) =
    lemma_bhkdv_main_term_monotonicity(u0) && lemma_bhkdv_main_term_limit()

"""
    lemma_bhkdv_main_term_monotonicity()

# Statement
The function
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
where `p0 = 1 + α + (1 + α)^2 / 2` and
```
a0 = finda0(α) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))^2
```
is increasing in `α` for all real `x` and all `α` in the interval
`(-1, -3 / 4)`.

# Use
This is used in the proof of [`lemma_bhkdv_monotonicity_alpha`](@ref).

# Proof
We split the function as
```
-(1 + α) * a0 *
    (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) / (1 + α)
```
To prove that it is increasing in `α` it is enough to prove that both
factors are positive and increasing in `α`.

For the proof that
```
(clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) / (1 + α)
```
is positive and increasing in `α` see the lemma in the paper.

For `-(1 + α) * a0` we explicitly compute a lower bound for the value
and the derivative and check so that they are both positive.
"""
function lemma_bhkdv_main_term_monotonicity(u0::BHKdVAnsatz)
    @assert u0.ϵ < 1 // 4

    a0_minimum =
        ArbExtras.minimum_enclosure(Arf(-1), Arf(-3 // 4), lbound_tol = Arf(0)) do α
            -finda0αp1(α)
        end

    Arblib.ispositive(a0_minimum) || error("failed proving that a0 is positive")

    a0αp1_deriv = ArbExtras.derivative_function(finda0αp1)

    a0_deriv_minimum =
        ArbExtras.minimum_enclosure(Arf(-1), Arf(-3 // 4), lbound_tol = Arf(0)) do α
            -a0αp1_deriv(α)
        end

    Arblib.ispositive(a0_deriv_minimum) || error("failed proving that a0 is increasing")

    return true
end

"""
    lemma_bhkdv_main_term_limit()

# Statement
The function
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
where `p0 = 1 + α + (1 + α)^2 / 2` and
```
a0 = finda0(α) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))^2
```
converges to
```
2 / π^2 * clausencmzeta(x, 2, 1)
```
from above as `α -> -1` for all real values of `x`.

# Use
This is used in the proof of [`lemma_bhkdv_monotonicity_alpha`](@ref).

# Proof
That the convergence has to be from above is an easy consequence of
[`lemma_bhkdv_monotonicity_main_term`](@ref). It is therefore enough
to prove that the limit is what is claimed.

We split the function as
```
-(1 + α) * a0 *
    (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) / (1 + α)
```
The first factor converges to `2 / π^2`. The second factor converges
to the derivative of `clausencmzeta` w.r.t. `α` at `α = -1`, which is
exactly `clausencmzeta(x, 2, 1)`
"""
lemma_bhkdv_main_term_limit() = true

"""
    lemma_bhkdv_weight_div_asymptotic_enclosure()

# Statement
Let
```
γ = 1 // 2
-1 < α < 0
p0 = (1 + α) + (1 + α)^2 / 2
```
then the function
```
x^((1 + α) / 2) * log(inv(x)) * (1 + α) / (1 - x^p0)
```
satisfies the inequality
```
0 < x^((1 + α) / 2) * log(inv(x)) * (1 + α) / (1 - x^p0) < 1
```

# Proof
See paper
"""
function lemma_bhkdv_weight_div_asymptotic_enclosure(u0::BHKdVAnsatz{Arb})
    @assert u0.γ == 1 // 2
    return Arblib.unit_interval!(zero(u0.ϵ))
end
