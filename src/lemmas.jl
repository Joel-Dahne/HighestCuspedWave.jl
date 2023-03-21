"""
    lemma_kdv_asymptotic_expansion

Corresponds to Lemma 4.2 in the paper.

# Statement

Let `u0` be as in [`equation_kdv_u0`](@ref). Then the following
asymptotic expansions hold near `x = 0`
```
u0(x) = sum(a⁰[j] * abspow(x, -α + j * p0) for j = 0:N0) +
    sum(1:Inf) do m
        (-1)^m / factorial(2m) * (
            sum(a[j] * zeta(1 - α + j * p0 - 2m) for j = 0:N0)
            + sum(b[n] * n^2m for n = 1:N1)
        ) * x^2m
    end
```
and
```
H(u0)(x) = -sum(A⁰[j] * abspow(x, -2α + j * p0) for j = 0:N0) -
    sum(1:Inf) do m
        (-1)^m / factorial(2m) * (
            sum(a[j] * zeta(1 - 2α + j * p0 - 2m) for j = 0:N0)
            + sum(b[n] * n^(2m + α) for n = 1:N1)
        ) * x^2m
    end
```
valid for `abs(x) < 2π`. The following notation is used
```
p0 = u0.p0
N0 = u0.N0
N1 = u0.N1
a[j] = u0.a[j]
b[n] = u0.b[n]

a⁰[j] = gamma(α - j * p0) * cospi((α - j * p0) / 2) * a[j]
A⁰[j] = gamma(α - 1 - j * p0) * cospi((α - 1 - j * p0) / 2) * a[j]
```
"""
function lemma_kdv_asymptotic_expansion end

"""
    lemma_kdv_asymptotic_expansion_defect

Corresponds to Lemma 4.3 in the paper.

# Statement

Let `u0` be as in [`equation_kdv_u0`](@ref). Then the following
asymptotic expansions hold near `x = 0`
```
H(u0)(x) + u0(x)^2 / 2 = ... # Not written out, see paper
```
The code doesn't use the expression of the expansion given in the
paper directly and we therefore don't write it out above. See the
asymptotic version of [`defect`](@ref) for [`BHAnsatz`](@ref) for how
the expansion is computed.
"""
function lemma_kdv_asymptotic_expansion_defect end

"""
    lemma_bhkdv_main_term_limit(u0::BHKdVAnsatz{Arb})

Corresponds to Lemma 4.5 in the paper.

# Statement

The function
```
a0(α) * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 2 + (1 + α)^2 / 2))
```
with
```
a0(α) = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```
converges to
```
2 / π^2 * clausencmzeta(x, 2, 1)
```
pointwise in `x` as `α` goes to `-1`. Moreover, for `-1 < α < -3 / 4`
it satisfies the inequality
```
a0(α) * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 2 + (1 + α)^2 / 2)) >=
    2 / π^2 * clausencmzeta(x, 2, 1)
```

# Computer assisted part of proof
Part of the proof is computer assisted and handled by this function.

The computer assisted part is proving that `-(1 + α) * a0(α)` is
positive and increasing for `-1 < α < -3 / 4`.

This is done by checking that the minimum values of `a0` and its
derivative on the interval are positive. If the proof fails it throws
an error.

It also checks that `u0.ϵ < 1 / 4` so that the lemma actually applies
to all values of `α` covered by `u0`.
"""
function lemma_bhkdv_main_term_limit(u0::BHKdVAnsatz{Arb})
    u0.ϵ < 1 // 4 || error("requires that u0.ϵ < 1 / 4")

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
    lemma_bhkdv_asymptotic_expansion

Corresponds to Lemma 4.6 in the paper.

# Statement

The approximation for `α` in `I₁`, with `u0` as in
[`equation_bhkdv_u0`](@ref) and `H(u0)` as in
[`equation_bhkdv_Hu0`](@ref) has the asymptotic expansion
```
u0(x) = a0 * (
        gamma(α) * cospi(α / 2)
        - gamma(-1 - (1 + α)^2 / 2) * cospi((1 + (1 + α)^2 / 2) / 2) * abspow(x, 1 + α + (1 + α)^2 / 2)
    ) * abspow(x, -α) +
    a0 * sum(1:Inf) do m
        (-1)^m / factorial(2m) *
        (zeta(1 - α - 2m) - zeta(2 + (1 + α)^2 / 2 - 2m)) * x^2m
    end +
    sum(a⁰[j] * abspow(x, -α_hat + j * p0) for j = 1:N0) +
    sum(1:Inf) do m
        (-1)^m / factorial(2m) * (
            sum(a[j] * zeta(1 - α_hat + j * p0 - 2m) for j = 1:N0)
            + sum(b[n] * n^2m for n = 1:N1)
        ) * x^2m
    end
```
and
```
H(u0)(x) = -a0 * (
        gamma(2α) * cospi(α)
        - gamma(-1 + α - (1 + α)^2 / 2) * cospi((1 - α + (1 + α)^2 / 2) / 2) * abspow(x, 1 + α + (1 + α)^2 / 2)
    ) * abspow(x, -2α) -
    a0 * sum(
        (-1)^m / factorial(2m) *
        (zeta(1 - 2α - 2m) - zeta(2 - α + (1 + α)^2 / 2 - 2m)) * x^2m
        for m = 1:Inf
    ) -
    sum(A⁰[j] * abspow(x, -α - α_hat + j * p0) for j = 0:N0) -
    sum(m = 1:Inf) do m
        (-1)^m / factorial(2m) * (
            sum(a[j] * zeta(2 - α - α_hat + j * p0 - 2m) for j = 0:N0)
            + sum(b[n] * n^(2m + α) for n = 1:N1)
        ) * x^2m
    end
```
valid for `abs(x) < 2π`. The following notation is used
```
a0 = finda0(α)
α_hat = u0.v0.α
p0 = u0.p0
N0 = u0.N0
N1 = u0.N1
a[j] = u0.a[j]
b[n] = u0.b[n]

a⁰[j] = gamma(α_hat - j * p0) * cospi((α_hat - j * p0) / 2) * a[j]
A⁰[j] = gamma(α_hat - 1 - j * p0) * cospi((α_hat - 1 - j * p0) / 2) * a[j]
```
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

For `-1 < α < 0` the function
```
abs(x)^-α / u0(x)
```
is non-zero and bounded at `x = 0`.
"""
function lemma_inv_u0_normalised end

"""
    lemma_bhkdv_inv_u0_normalised

Corresponds to Lemma 8.2 in the paper.

# Statement

For `α ∈ I1` the function
```
gamma(1 + α) * abspow(x, -α) * (1 - abspow(x, 1 + α + (1 + α)^2 / 2)) / u0(x)
```
is non-zero and uniformly bounded in `α` at `x = 0`.
"""
function lemma_bhkdv_inv_u0_normalised end

"""
    lemma_bhkdv_weight_factor_enclosure

Corresponds to Lemma 9.1 in the paper.

# Statement

Let `-1 < α < 0` and `0 < x < 1`, then
```
0 < x^(1 + α) / 2 * log(1 / x) * (1 + α) / (1 - x^(1 + α + (1 + α)^2 / 2)) < 1
```
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

For all `-1 < α < 0` the function
```
(1 - t)^(-α - 1) + (1 - t)^(-α - 1) - 2t^(-α - 1)
```
is increasing and continuous in `t` for `0 < t < 1` and has the unique
root `r0(α)`. For `0 < x < π` and `0 < t < r0(α)` the function
`I_hat(x, t, α)` (see [`equation_I_hat`](@ref)) is increasing in `x`.
"""
function lemma_I_hat_root_zero end

"""
    lemma_I_hat_root

Corresponds to Lemma 11.2 in the paper.

# Statement

For all `-1 < α < 0` and `0 < x < π` the function `I_hat(x, t, α)`
(see [`equation_I_hat`](@ref)) is increasing and continuous in `t` for
`0 < t < 1` and has the limits `-Inf` and `Inf` for `t` going to `0`
and `1` respectively.

Moreover, the unique root, `r(α, x)`, in `t` is decreasing in `x` and
satisfies the inequality
```
1 / 2 < r(α, x) < r0(α)
```
with `r0(α)` as in [lemma_I_hat_root_zero`](@ref).
"""
function lemma_I_hat_root end

"""
    lemma_I_hat_root_alpha

Corresponds to Lemma 11.3 in the paper.

# Statement

Let `α_interval = [αₗ, αᵤ] ⊆ (-1, 0)` and `0 < x < π`. Let `I_hat_dα`
denote the derivative of `I_hat(x, t, α)` (see
[`equation_I_hat`](@ref)) with respect to `α`. If
```
I_hat_dα(x, r(αₗ, x), α) > 0
```
, with `r(α, x)` as in [`lemma_I_hat_root`](@ref), for all `α ∈
α_interval` then `r(α, x) <= r(αₗ, x)` for all `α ∈ α_interval`.
Similarly, if `I_hat_dα(x, r(αᵤ, x), α) > 0` for all `α ∈ α_interval`,
then `r(α, x) >= r(αᵤ, x)` for all `α ∈ α_interval`.
"""
function lemma_I_hat_root_alpha end

"""
    lemma_I_positive

Corresponds to Lemma 11.4 in the paper.

# Statement

For all `-1 < α < 0` and `0 < x < π` we have `I(x, y, α) > 0` for `x <
y < π` and `I_hat(x, t, α) > 0` for `1 < t < π / x`. Here `I(x, y, α)`
and `I_hat(x, t, α)` are as in [`equation_I`](@ref) and
[`equation_I_hat`](@ref) respectively,
"""
function lemma_I_positive end

"""
    lemma_U0_primitive_weight_x

Corresponds to Lemma 11.5 in the paper.

# Statement

For `-1 < α < 0`, `0 < x < π` and `w(x) = abs(x)` we have
```
U0(x) = 2clausencmzeta(x, 2 - α) + 2(clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) -
    2(clausenc(x * (1 - r), 2 - α) + clausenc(x * (1 + r), 2 - α) - 2clausenc(x * r, 2 - α)) -
    2x * r * (-clausens(x * (1 - r), 1 - α) + clausens(x * (1 + r), 1 - α) - 2clausens(x * r, 1 - α))
```
Here `r` is the root occurring in [`lemma_I_hat_root`](@ref).
"""
function lemma_U0_primitive_weight_x end

"""
    lemma_kdv_U0_asymptotic

Corresponds to Lemma 11.6 in the paper.

# Statement

Let `0 < ϵ < π / 2`, for `-1 < α < 0`, `0 < x < ϵ` and `u0.w(x) =
abspow(x, p)` with `-α < p < 1` and `1 + α != p` we have
```
U01 / x^(-α + p) <= c1 + d1 * x^(3 + α)
```
and
```
U02 / x^(-α + p) <= c2 + d2 * x^(2 + α - p)
```
where
```
c = gamma(1 + α) * sinpi(-α / 2) * (
    2 / (α - p) +
    gamma(-α) * gamma(1 + p) / gamma(1 - α + p) +
    hypgeom_2f1(1 + α, 1 + p, 2 + p, -1) / (1 + p) -
    2r^p * (
        2r^-α / (α - p) +
        r * hypgeom_2f1(1 + α, 1 + p, 2 + p, -r) / (1 + p) +
        r * hypgeom_2f1(1 + α, 1 + p, 2 + p, r) / (1 + p)
    )
)

c2 = gamma(1 + α) * sinpi(-α / 2) * (
    gamma(-α) * gamma(α - p) / gamma(-p) +
    (hypgeom_2f1(1 + α, α - p, 1 + α - p, -1) - 2) / (α - p)
)

d1 = 2sum(1:M-1) do m
        (-1)^m * zeta(-α - 2m) * ϵ^(2m - 2) / factorial(2m) *
            sum(binomial(2m, 2k) / (2k + 1 + p) for k = 0:m-1)
    end +
    1 / ϵ^2 * sum((-1)^m * zeta(-α - 2m) * (2ϵ)^2m / factorial(2m) for m = M:Inf)

d2 = -gamma(1 + α) * sinpi(-α / 2) * (1 + α) * (2 + α) / (2 + α - p) / π^(2 + α - p) +
    2π^(p - 1) * sum(1:M-1) do m
        (-1)^m *
        zeta(-α - 2m) *
        Arb(π)^2m /
        factorial(2m) *
        sum(binomial(2m, 2k) * (ϵ / π)^(2(m - 1 - k)) / (2k + 1 + p) for k = 0:m-1)
    end +
    6π^(p - 1) * sum((-1)^m * zeta(-α - 2m) * (3π / 2)^2m / factorial(2m) for m = M:Inf)
```
The tails for `d1` and `d2` are of the same form as in those for the
Clausen functions and can be bounded using
[`lemma_clausen_remainder`](@ref).
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

Let `d(f, k)(x)` denote the `k`th derivative of `f` evaluated at `x`.

Let `m` be a non-negative integer and let `I` be an interval
containing zero. Consider a function `f(x)` with a zero of order `n`
at `x = 0` and such that `d(f, m + n)` is absolutely continuous on
`I`. Then for all `x ∈ I` we have
```
f(x) / x^n = sum(f[k + n](0) * x^k for k = 0:m) + f[m + n + 1](ξ) * x^(m + 1)
```
for some `ξ` between `0` and `x`. Here we use `f[k](x) = d(f, k)(x) /
factorial(k)`.

Furthermore, if `d(f, m + n + p)` is absolutely continuous for some
non-negative integer `p` we have
```
d(x -> f(x) / x^n, p)(x) =
    sum(factorial(k + p) / factorial(k) * f[k + n + p](0) * x^k for k = 0:m) +
    factorial(m + p + 1) / factorial(m + 1) * f[m + n + p + 1](ξ) * x^(m + 1)
```
for some `ξ` between `0` and `x`.
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

For `x != 0` the function
```
(abs(x)^t - 1) / t
```
is non-decreasing in `t`. In the limit as `t` goes to zero it is equal
to `log(abs(x))`.
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
Let `d(f, k)(x)` denote the `k`th derivative of `f` evaluated at `x`.

Let `β >= 1, s >= 0, 2M >= s + 1` and `abs(x) < 2π`, we then have the
following bounds:
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
    lemma_clausen_remainder_bhkdv

Corresponds to Lemma C.8 in the paper.

# Statement

For `-1 < α < 0` and `M1 >= 1` we have
```
abs(
    sum(M1:Inf) do m
        (-1)^m / factorial(2m) *
        (zeta(1 - α - 2m) - zeta(2 + (1 + α)^2 / 2 - 2m)) / (1 + α) * x^2m
    end
) <=
abs(sum((-1)^m / factorial(2m) * dzeta(s1[m] - 2m) * x^2m for m = M1:Inf)) +
(1 + α) / 2 * abs(sum((-1)^m / factorial(2m) * dzeta(s2[m] - 2m) * x^2m for m = M1:Inf))
```
with `1 - α < s1[m] < 2` and `2 < s2[m] < 2 + (1 + α)^2 / 2`.
Similarly, for `M2 >= 2`,
```
abs(
    sum(M2:Inf) do m
        (-1)^m / factorial(2m) *
        (zeta(1 - 2α - 2m) - zeta(2 - α +  (1 + α)^2 / 2 - 2m)) / (1 + α) * x^2m
    end
) <=
2abs(sum((-1)^m / factorial(2m) * dzeta(s3[m] - 2m) * x^2m for m = M2:Inf)) +
(1 - (1 + α) / 2) * abs(sum((-1)^m / factorial(2m) * dzeta(s4[m] - 2m) * x^2m for m = M2:Inf))
```
with `1 - 2α < s3[m] < 3` and `2 - α + (1 + α)^2 / 2 < s4[m] < 3`.
"""
function lemma_clausen_remainder_bhkdv end

"""
    lemma_clausen_derivative_remainder_interval

Corresponds to Lemma C.9 in the paper.

# Statement

Let `d(f, k)(x)` denote the `k`th derivative of `f` evaluated at `x`.

Let `β >= 1`, `s` be the interval ``[sₗ, sᵤ]`` with `sₗ >= 0`, 2M >=
sᵤ + 1`, `s[m] ∈ s` for `m >= M` and `abs(x) < 2π`, we then have the
following bounds
```
abs(sum((-1)^m * d(zeta, β)(s[m] - 2m) * x^2m / factorial(2m) for m = M:Inf)) <=
    sum(j1 + j2 + j3 = β) do j1, j2, j3
        binomial(β, (j1, j2, j3)) *
        2(log(2π) + π / 2)^j1 *
        (2π)^(sᵤ - 1) *
        abs(d(zeta, j3)(1 + 2M - sᵤ)) *
        sum(p[j2](1 + 2m - s[m]) * (x / 2π)^2m for m = M:Inf)
    end

abs(sum((-1)^m * d(zeta, β)(s[m] - 2m - 1) * x^(2m + 1) / factorial(2m + 1) for m = M:Inf)) <=
    sum(j1 + j2 + j3 = β) do j1, j2, j3
        binomial(β, (j1, j2, j3)) *
        2(log(2π) + π / 2)^j1 *
        (2π)^(sᵤ - 1) *
        abs(d(zeta, j3)(2 + 2M - sᵤ)) *
        sum(p[j2](2 + 2m - s[m]) * (x / 2π)^(2m + 1) for m = M:Inf)
    end
```

Here `p[j2]` is the same as in
[`lemma_clausen_derivative_remainder`](@ref) and given recursively by
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
    polygamma(0, 1 + 2m - s[m])^q[0] * ⋯ * polygamma(j2 - 1, 1 + 2m - s[m])^q[j2 - 1] (x / 2π)^2m
    for m = M:Inf
) <= abs(polygamma(1, 1 + 2M - sᵤ)^q[1] * ⋯ * polygamma(j2 - 1, 1 + 2M - sᵤ)^q[j2 - 1]) *
    1 / 2^(q[0] / 2) * (2π)^(-2M) * lerch_phi(x^2 / 4π^2, -q[0] / 2, M + 1) * x^2M

sum(
    polygamma(0, 2 + 2m - s[m])^q[0] * ⋯ * polygamma(j2 - 1, 2 + 2m - s[m])^q[j2 - 1] (x / 2π)^(2m + 1)
    for m = M:Inf
) <= abs(polygamma(1, 2 + 2M - sᵤ)^q[1] * ⋯ * polygamma(j2 - 1, 2 + 2M - sᵤ)^q[j2 - 1]) *
    1 / 2^(q[0] / 2) * (2π)^(-2M - 1) * lerch_phi(x^2 / 4π^2, -q[0] / 2, M + 1) * x^(2M + 1)
```
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
    lemma_bhkdv_F0_P_factor

Corresponds to Lemma E.1 in the paper.

# Statement

The factor
```
a0 * (c(α) - c(α - p0) * x^p0)
```
can be written as
```
C1 * x^p0 - C2 * (x^p0 - 1) / x^p0
```
with
```
C1 = a0 * (c(α) - c(α - p0))
C2 = a0 * c(α) * p0
```

Here
```
a0 = finda0(α)
c(s) = gamma(1 + s) * cospi(s  2)
p0 = 1 + α + (1 + α)^2 / 2
```
"""
function lemma_bhkdv_F0_P_factor end

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
    lemma_kdv_U0_hybrid_R

Corresponds to Lemma G.1 in the paper.

# Statement

Let `0 < ϵ < π / 2`, for `-1 < α < 0`, `0 <= x <= ϵ` and `-α < p < 1`
with `1 + α != p` we have
```
U03R(x) / (x^(-α + p) * log(1 / x)) <=
    2x^(2 + α - p) * π^(p - 1) / log(1 / x) * sum(1:M-1) do m
        binomial(2m, 2k) / (2k + 1 + p)^2 * (x / π)^(2(m - 1 - l))
    end -
    2x^(2 + α - p) * π^(p - 1) * log(π / x) / log(1 / x) * sum(1:M-1) do m
        binomial(2m, 2k) / (2k + 1 + p) * (x / π)^(2(m - 1 - l))
    end +
    6x^(2 + α - p) * π^(p - 1) / log(1 / x) * sum(M:Inf) do m
        (-1)^m * zeta(-α - 2m) * (3π / 2)^2m / factorial(2m)
    end
```

Here
```
U03R(x) = 2x^(1 + p) * sum(1:Inf) do m
    (-1)^m * zeta(-α - 2m) * x^2m / factorial(2m) *
        sum(binomial(2m, 2k) * ∫ t^(2k + p) * log(inv(t)) dt for k = 0:m-1)
end
```
where the integration is taken from `0` to `π / x`.
"""
function lemma_kdv_U0_hybrid_R end
