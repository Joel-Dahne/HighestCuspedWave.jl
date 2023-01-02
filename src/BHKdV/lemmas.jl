"""
    lemma_bhkdv_monotonicity_alpha(u0::BHKdVAnsatz)

# Statement
The function `u0(x)` is **lower bounded** by `u0.v0(x)` for all `x ∈
[0, π]`. This holds for all `ϵ ∈ (0, 1)`.

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

- **PROVE:** We probably need to prove that `u0` is positive as well.
  So that the inequality holds after taking the absolute value. This
  could be done somewhere else though.
- **IMRPOVE:** We could add a check that this bound indeed holds for
  some points. To avoid potential bugs.
"""
lemma_bhkdv_monotonicity_alpha(u0::BHKdVAnsatz) =
    lemma_bhkdv_main_term_monotonicity() && lemma_bhkdv_main_term_limit()

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
`(-1, 0)`.

# Use
This is used in the proof of [`lemma_bhkdv_monotonicity_alpha`](@ref).

# Proof
Since `clausencmzeta` is `2π` periodic and even it is enough to prove
it for the `x` in the interval `[0, π]`.

**PROVE:** Finish this proof.
"""
lemma_bhkdv_main_term_monotonicity() = true

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
That the convergence has to be from below is an easy consequence of
[`lemma_bhkdv_monotonicity_main_term`](@ref). It is therefore enough
to prove that the limit is what is claimed.

Since `clausencmzeta` is `2π` periodic and even it is enough to prove
it for the `x` in the interval `[0, π]`.

**PROVE:** Finish this proof.
"""
lemma_bhkdv_main_term_limit() = true

"""
    lemma_bhkdv_main_term_H_monotonicity()

# Statement
The function
```
-a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
```
, which is `H` to the term in
[`lemma_bhkdv_main_term_H_monotonicity`](@ref), where `p0 = 1 + α + (1
+ α)^2 / 2` and
```
a0 = finda0(α) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))^2
```
is decreasing in `α` for all real `x` and all `α` in the interval
`(-1, 0)`.

# Proof
Since `clausencmzeta` is `2π` periodic and even it is enough to prove
it for the `x` in the interval `[0, π]`.

**PROVE:** Finish this proof. We might be able to use
  [`lemma_bhkdv_main_term_H_monotonicity`](@ref), or at the very least
  the proof should be similar.
"""
lemma_bhkdv_main_term_H_monotonicity() = true

"""
    lemma_bhkdv_main_term_H_limit()

# Statement
The function
```
-a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
```
, which is `H` to the term in
[`lemma_bhkdv_main_term_H_monotonicity`](@ref), where `p0 = 1 + α + (1
+ α)^2 / 2` and
```
a0 = finda0(α) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))^2
```
converges to
```
-2 / π^2 * clausencmzeta(x, 3, 1)
```
from above as `α -> -1` for all real values of `x`.

# Proof
That the convergence has to be from above is an easy consequence of
[`lemma_bhkdv_monotonicity_H_main_term`](@ref). It is therefore enough
to prove that the limit is what is claimed.

Since `clausencmzeta` is `2π` periodic and even it is enough to prove
it for the `x` in the interval `[0, π]`.

**PROVE:** Finish this proof.
"""
lemma_bhkdv_main_term_H_limit() = true

"""
    lemma_bhkdv_weight_div_asymptotic_bound()

# Statement
Let
```
0 <= γ < 1
-1 < α < 0
p0 = (1 + α) + (1 + α)^2 / 2
```
then the function
```
-x^((1 - γ) * (1 + α)) * log(x) / (gamma(1 + α) * (1 - x^p0))
```
is bounded by
```
1.1
```
- **PROVE:** This bound can be checked numerically to be valid but
  should be updated with a better bound once we have that. Annoyingly
  enough the upper bound is not quite one, but some number slightly
  larger than one.

# Proof
**PROVE:** Prove this.
"""
lemma_bhkdv_weight_div_asymptotic_bound(u0::BHKdVAnsatz{T}) where {T} = convert(T, 1.1)
