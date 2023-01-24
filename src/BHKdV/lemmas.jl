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

    a0_deriv_minimum =
        ArbExtras.minimum_enclosure(Arf(-1), Arf(-3 // 4), lbound_tol = Arf(0)) do α
            if α isa Arb
                -finda0αp1(ArbSeries((α, 1)))[1]
            else # α isa ArbSeries
                Arblib.derivative(-finda0αp1(ArbSeries(α, degree = Arblib.degree(α) + 1)))
            end
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
-x^((1 + α) / 2) * log(x) * (1 + α) / (1 - x^p0)
```
satisfies the inequality
```
0 < -x^((1 + α) / 2) * log(x) * (1 + α) / (1 - x^p0) < 1
```

# Proof
See paper
"""
function lemma_bhkdv_weight_div_asymptotic_enclosure(u0::BHKdVAnsatz{Arb})
    @assert u0.γ == 1 // 2
    return Arblib.unit_interval!(zero(u0.ϵ))
end
