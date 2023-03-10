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
