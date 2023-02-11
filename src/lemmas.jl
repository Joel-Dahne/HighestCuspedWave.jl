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
