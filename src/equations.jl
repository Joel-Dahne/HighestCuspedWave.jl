"""
    equation_kdv_u0

Corresponds to Equation 12 in the paper.

# Equation
```
u0(x) = sum(a[j] * clausencmzeta(x, 1 - α + j * p0) for j = 0:N0) +
    sum(b[n] * (cos(n * x) - 1) for n = 1:N1)
```
The following notation is used
```
p0 = u0.p0
N0 = u0.N0
N1 = u0.N1
a[j] = u0.a[j]
b[n] = u0.b[n]
```
"""
function equation_kdv_u0 end

"""
    equation_a0

Corresponds to Equation 13 in the paper.

# Equation
```
a0 = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```
"""
function equation_a0 end

"""
    equation_kdv_Hu0

Corresponds to Equation 14 in the paper.

# Equation
```
H(u0)(x) = -sum(a[j] * clausencmzeta(x, 1 - 2α + j * p0) for j = 0:N0) -
    sum(b[n] * n^α * (cos(n * x) - 1) for n = 1:N1)
```
The following notation is used
```
p0 = u0.p0
N0 = u0.N0
N1 = u0.N1
a[j] = u0.a[j]
b[n] = u0.b[n]
```
"""
function equation_kdv_Hu0 end

"""
    equation_p0

Corresponds to Equation 17 in the paper.

# Equation
```
gamma(2α - p0) * cospi((2α - p0) / 2) / (gamma(α - p0) * cospi((α - p0) / 2)) =
    2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```
"""
function equation_p0 end

"""
    equation_bhkdv_u0

Corresponds to Equation 20 in the paper.

# Equation
```
u0(x) = a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 2 + (1 + α)^2 / 2)) +
    sum(a[j] * clausencmzeta(x, 1 - α_hat + j * p0) for j = 1:N0) +
    sum(b[n] * (cos(n * x) - 1) for n = 1:N1)
```
Note that the notation is slightly different than in the paper. In the
paper there are two terms of `clausecmzeta(x, 1 - α + p0)`, here we
have merged them into one term with the corresponding value for
`a[1]`. The notation used is
```
a0 = finda0(α)
α_hat = u0.v0.α
p0 = u0.v0.p0
N0 = u0.v0.N0
N1 = u0.v0.N1
a[j] = u0.v0.a[j]
b[n] = u0.v0.b[n]
```
"""
function equation_bhkdv_u0 end

"""
    equation_bhkdv_Hu0

Corresponds to Equation 22 in the paper.

# Equation
```
H(u0)(x) = -a0 * (clausencmzeta(x, 1 - 2α) - clausecmzeta(x, 2 - α + (1 + α)^2 / 2))
    -sum(a[j] * clausencmzeta(x, 1 - α - α_hat + j * p0) for j = 1:N0) -
    sum(b[n] * n^α * (cos(n * x) - 1) for n = 1:N1)
```
Note that as in [`equation_bhkdv_u0`](@ref) the notation is slightly
different than in the paper. In the paper there are two terms of
`clausecmzeta(x, 1 - α + p0)`, here we have merged them into one term
with the corresponding value for `a[1]`. The notation used is
```
a0 = finda0(α)
α_hat = u0.v0.α
p0 = u0.v0.p0
N0 = u0.v0.N0
N1 = u0.v0.N1
a[j] = u0.v0.a[j]
b[n] = u0.v0.b[n]
```
"""
function equation_bhkdv_Hu0 end

"""
    equation_bhkdv_F0_factor

Corresponds to Equation 28 in the paper.

# Equation
```
(H(u0)(x) + u0(x)^2 / 2) /
    (gamma(1 + α) * log(inv(x)) * (1 - x^(1 + α + (1 + α)^2 / 2)) * x^(1 - α))
```
"""
function equation_bhkdv_F0_factor end

"""
    equation_I

This equation is not numbered in the paper but given shortly after
Equation 29

# Equation
```
I(x, y, α) = clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)
```
"""
function equation_I end

"""
    equation_I_hat

This equation is not numbered in the paper but given right before
Lemma 11.1.

# Equation
```
I_hat(x, t, α) = clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
"""
function equation_I_hat end

"""
    equation_clausenc_zeta

Corresponds to Equation 43 in the paper.

# Equation
```
clausenc(x, s) = gamma(1 - s) / (2π)^(1 - s) * cospi((1 - s) / 2) *
    (zeta(1 - s, x / 2π) + zeta(1 - s, 1 - x / 2π))
```

# Notes
When `s` is a non-negative integer the expression has a removable
singularity. See the discussions in the paper for how these removable
singularities can be handled.
"""
function equation_clausenc_zeta end

"""
    equation_clausens_zeta

Corresponds to Equation 44 in the paper.

# Equation
```
clausens(x, s) = gamma(1 - s) / (2π)^(1 - s) * sinpi((1 - s) / 2) *
    (zeta(1 - s, x / 2π) - zeta(1 - s, 1 - x / 2π))
```

# Notes
When `s` is a non-negative integer the expression has a removable
singularity. See the discussions in the paper for how these removable
singularities can be handled.
"""
function equation_clausens_zeta end

"""
    equation_clausenc_expansion

Corresponds to Equation 46 in the paper.

# Equation
```
clausenc(x, s) = gamma(1 - s) * sinpi(s / 2) * abspow(x, s - 1) +
    sum((-1)^m * zeta(s - 2m) * x^2m / factorial(2m) for m = 0:Inf)
```
"""
function equation_clausenc_expansion end

"""
    equation_clausens_expansion

Corresponds to Equation 47 in the paper.

# Equation
```
clausens(x, s) = gamma(1 - s) * cospi(s / 2) * sign(x) * abspow(x, s - 1) +
    sum((-1)^m * zeta(s - 2m - 1) * x^(2m + 1) / factorial(2m + 1) for m = 0:Inf)
```
"""
function equation_clausens_expansion end

"""
    equation_clausenc_derivative_expansion

Corresponds to Equation 49 in the paper.

# Equation
```
clausenc(x, s, β) = d(s -> gamma(1 - s) * sinpi(s / 2) * abspow(x, s - 1), β)(s) +
    sum((-1)^m * d(zeta, β)(s - 2m) * x^2m / factorial(2m) for m = 0:Inf)
```
where `d(f, k)(x)` denote the `k`th derivative of `f` evaluated at `x`.
"""
function equation_clausenc_derivative_expansion end

"""
    equation_clausens_derivative_expansion

Corresponds to Equation 50 in the paper.

# Equation
```
clausens(x, s, β) = d(s -> gamma(1 - s) * cospi(s / 2) * abspow(x, s - 1), β)(s) +
    sum((-1)^m * d(zeta, β)(s - 2m - 1) * x^(2m + 1) / factorial(2m + 1) for m = 0:Inf)
```
where `d(f, k)(x)` denote the `k`th derivative of `f` evaluated at `x`.
"""
function equation_clausens_derivative_expansion end
