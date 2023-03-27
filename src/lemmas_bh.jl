"""
    lemma_bh_asymptotic_expansion

Corresponds to Lemma 4.3 in the Burgers-Hilbert paper.

# Statement

The approximation `u0` has the asymptotic expansions
```
u0(x) = -1 / π * abs(x) * log(abs(x)) - (γ - 1) / π * abs(x) +
    sum(a⁰[j] * abspow(x, -α + j * p0) for j = 1:N0) +
    sum(m = 1:Inf) do m
        (-1)^m / factorial(2m) * (
            2 / π^2 * dzeta(2 - 2m)
            + sum(a[j] * zeta(1 - α + j * p0 - 2m) for j = 1:N0)
            + sum(b[n] * n^2m for n = 1:N1)
        ) * x^2m
    end
```
and
```
H(u0)(x) = -1 / 2π^2 * x^2 * log(abs(x))^2 + (3 - 2γ) / 2π^2 * x^2 * log(abs(x)) -
    sum(A⁰[j] * abspow(x, 1 - α + j * p0) for j = 1:N0) +
    1 / 2 *(
        (36γ - 12γ^2 - 24γ₁ - 42 + π^2) / 12π^2
        + sum(a[j] * zeta(-α + j * p0) for j = 1:N0)
        + sum(b[n] * n for n = 1:N1)
    ) * x^2 -
    sum(m = 2:Inf) do m
        (-1)^m / factorial(2m) * (
            2 / π^2 * dzeta(3 - 2m)
            + sum(a[j] * zeta(2 - α + j * p0 - 2m) for j = 1:N0)
            + sum(b[n] * n^(2m - 1) for n = 1:N1)
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
`γ₁` is the Stieltjes constant and `γ = γ₀` is the Euler's constant.
"""
function lemma_bh_asymptotic_expansion end

"""
    lemma_bh_inv_u0_normalised

Corresponds to Lemma 5.1 in the Burgers-Hilbert paper.

# Statement
The function
```
-abs(x) * log(abs(x)) / u0(x)
```
is positive and bounded at `x = 0` and for `abs(x) < 1` it has the
expansion
```
abs(x) * log(abs(x)) / u0(x) = -inv(
    -1 / π
    -(γ - 1) / π / log(abs(x))
    + sum(a⁰[j] * abspow(x, -1 - α + j * p0) / log(abs(x)) for j = 1:N0)
    + sum(m = 1:Inf) do m
        (-1)^m / factorial(2m) * (
            2 / π^2 * dzeta(2 - 2m)
            + sum(a[j] * zeta(1 - α + j * p0 - 2m) for j = 1:N0)
            + sum(b[n] * n^2m for n = 1:N1)
        ) * abs(x)^(2m - 1) / log(abs(x))
    end
)
```
"""
function lemma_bh_inv_u0_normalised end

"""
    lemma_bh_defect_normalised

Corresponds to Lemma 5.2 in the Burgers-Hilbert paper.

# Statement
The function
```
(H(u0)(x) + 1 / 2 * u0(x)^2) / (x^2 * log(abs(x)))
```
is non-zero and bounded at `x = 0` and for `abs(x) < 1`.

The lemma also gives the expansion at `x = 0`. We do however not use
this expression of the expansion directly. See the asymptotic version
of [`defect`](@ref) for [`BHAnsatz`](@ref) for how the expansion is
computed.
"""
function lemma_bh_defect_normalised end

"""
    lemma_bh_I_hat

Corresponds to Lemma 6.1 in the Burgers-Hilbert paper.

# Statement
For all `0 < x < π` the function
```
I_hat(x, t) = log(sin(abs(x * (1 - t) / 2)) * sin(x * (1 + t) / 2) / sin(x * t / 2)^2)
```
is decreasing and continuous in `t` for `0 < t < 1` and has the limit
`Inf` for `t` tending to `0` and `-Inf` for `t` tending to `1`.
Moreover, the unique root, `r(x)`, is decreasing in `x` and satisfies
the inequality
```
1 / 2 < r(x) < 1 / sqrt(2)
```
"""
function lemma_bh_I_hat end

"""
    lemma_bh_I

Corresponds to Lemma 6.2 in the Burgers-Hilbert paper.

# Statement
For `0 < x < π` we have
```
I(x, y) = log(sin(abs(x - y) / 2) * sin((x + y) / 2) / sin(y / 2)^2) < 0
```
for all `x < y < π`.
"""
function lemma_bh_I end

"""
    lemma_bh_U1_asymptotic

Corresponds to Lemma 6.3 in the Burgers-Hilbert paper.

# Statement
Let
```
U1(x) = x^2 * ∫ abs(I_hat(x, t)) * t * sqrt(log(1 + inv(x * t))) dt
```
where the integration is from `0` to `1` and `I_hat(x, t)` is as in
[`lemma_bh_I_hat`](@ref).

For `0 <= x <= ϵ` with `ϵ < 1` we have
```
U1(x) / (-x^2 * log(x) * sqrt(log(1 + inv(x)))) <= 1 / sqrt(log(1 + inv(x))) * (
    log(2) / sqrt(log(inv(x)))
    + (c1 + log(2) * sqrt(log(1 + x))) / log(inv(x))
    + 3R1 / 8 * x^2 * (2 / sqrt(log(inv(x))) + (sqrt(π / 2) + 2sqrt(log(1 + x))) / log(inv(x)))
)
```
where
```
c1 = ∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt
```
integrated from `0` to `1` and `R1` is the supremum for `0 <= y <= ϵ`
of
```
1 / 2 * abs(d²/dy²(log(sin(y) / y)))
```
where we use `d²/dy²` to denote differentiating twice w.r.t. `y`.
"""
function lemma_bh_U1_asymptotic end

"""
    lemma_bh_U2_asymptotic

Corresponds to Lemma 6.4 in the Burgers-Hilbert paper.

# Statement
Let
```
U2(x) = ∫ abs(I(x, y)) * y * sqrt(log(1 + inv(y))) dy
```
where the integration is from `x` to `π` and `I(x, t)` is as in
[`lemma_bh_I`](@ref).

For `0 <= x <= ϵ` with `ϵ < 1 / 2` we have
```
U2(x) / (-x^2 * log(x) * sqrt(log(1 + inv(x)))) <= log(16 / 3sqrt(3)) / log(inv(x)) +
    (
        2 / 3 * log(inv(2x))^(3 / 2) / (log(inv(x)) * sqrt(log(1 + inv(x))))
        + R2 * sqrt(log(inv(2x))) / (8log(inv(x)) * sqrt(log(1 + inv(x))))
        + sqrt(log(2)) / sqrt(log(inv(1 + inv(x))))
        - log(2)^(3 / 2) / (log(inv(x)) * sqrt(log(1 + inv(x))))
        + R2 * sqrt(log(2)) * (1 - 4x^2) / (8log(inv(x)) * sqrt(log(1 + inv(x))))
    ) + sqrt(log(2)) / (log(inv(x)) * sqrt(log(1 + inv(x)))) * (
        1 / 2 * log((π^2 - x^2) / (1 - x^2))
        + log(1 - x^2) / 2x^2
        - π^2 * log(1 - x^2 / π^2) / 2x^2
    ) + D1 * c2 / (log(inv(x)) * sqrt(log(1 + inv(x))))
```
where
```
c2 = ∫ y * sqrt(log(1 + inv(y))) dt
```
integrated from `0` to `π`, `R2` is the supremum for `0 <= y <= 1 / 4`
of
```
1 / 2 * abs(d²/dy²(log(1 - y)))
```
where we use `d²/dy²` to denote differentiating twice w.r.t. `y` and
`D1` is the supremum for `0 < x < ϵ` of
```
-(log(sinc((1 - x / π) / 2)) + log(sinc((1 + x / π) / 2)) - 2log(sinc(1 / 2))) / x^2
```
Note that we are here using the Julia convention that `sinc(x) = sin(π
* x) / (π * x)`, which is different from the one used in the paper.
"""
function lemma_bh_U2_asymptotic end

"""
    lemma_bh_U2_term_bound

Corresponds to Lemma 6.5 in the Burgers-Hilbert paper.

# Statement
The function
```
f(x) = 2 / 3 * log(inv(2x))^(3 / 2) / (log(inv(x)) * sqrt(log(1 + inv(x))))
```
is bounded from above by `2 / 3` and is decreasing in `x` for `0 < x <
1 / 2`.
"""
function lemma_bh_U2_term_bound end

"""
    lemma_bh_bound_delta0

Corresponds to Lemma 7.2 in the Burgers-Hilbert paper.

We don't use the statement of the lemma in the code but some of the
details of the asymptotic version of [`F0`](@ref) are discussed in the
proof of the lemma.
"""
function lemma_bh_bound_delta0 end

"""
    lemma_bh_U_singular_integrals

Corresponds to Lemma C.1 in the Burgers-Hilbert paper.

# Statement
Let
```
U122(x) = -∫ I_hat(x, t) * t * sqrt(log(1 + inv(t))) dt
```
integrated from `1 - δ1` to `1` and
```
U21(x) = -∫ I(x, y) * y * sqrt(log(1 + inv(y))) dy
```
integrated from `x` to `x + δ2`.

We have
```
U122(x) = ξ1 / x * (
    (-clausens(0, 2) + clausens(2x, 2) - 2clausens(x, 2)) -
    (-clausens(x * δ1, s) + clausens(x * (2 - δ1), s) - 2clausens(x * (1 - δ1), s))
)
```
for some `ξ1` in the image of `t * sqrt(1 + inv(x * t))` for `1 - δ1
<= t <= 1` Furthermore
```
U21(x) = ξ2 * (
    (clausens(δ2, 2) + clausens(2x + δ2, 2) - 2clausens(x + δ2, 2)) -
    (-clausens(0, s) + clausens(2x, s) - 2clausens(x, s))
)
```
for some `ξ2` in the image of `y * sqrt(1 + inv(y))` for `x <= y <= x
+ δ2`.
"""
function lemma_bh_U_singular_integrals end

"""
    lemma_bh_U11_bound_zero

Corresponds to Lemma C.2 in the Burgers-Hilbert paper.

# Statement
For `0 < x < π` and `0 < t < 1 / 10x` the function
```
log(sin(x * t / 2)) * t * sqrt(log(1 + inv(x * t)))
```
is decreasing in `t`.
"""
function lemma_bh_U11_bound_zero end
