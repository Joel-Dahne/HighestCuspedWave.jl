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
