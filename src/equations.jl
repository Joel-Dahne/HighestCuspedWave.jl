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
