"""
    lemma_bhkdv_monotonicity_alpha(u0::BHKdV)

# Statement
The function `u0(x)` is **lower bounded** by `u0.v0(x)` for all `x ∈
[0, π]`. This holds for all `ϵ ∈ (0, 1)`.

# Use
This is used when upper bounding terms where `u0` occurs in the
denominator since we can then replace it by `u0.v0` and still get an
upper bound. For example this is used in
- [`alpha0`](@ref)

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
for all `α ∈ [-1, -1 + ϵ]` and `x ∈ [0, π]` is the result of
[`lemma_bhkdv_monotonicity_main_term`](@ref).

- **TODO:** We probably need to prove that `u0` is positive as well.
  So that the inequality holds after taking the absolute value. This
  could be done somewhere else though.
- **IMRPOVE:** We could add a check that this bound indeed holds for
  some points. To avoid potential bugs.
"""
lemma_bhkdv_monotonicity_alpha(u0::BHKdV) = lemma_bhkdv_monotonicity_main_term(u0)

"""
    lemma_bhkdv_monotonicity_main_term()

# Statement
The inequality
```
2 / π^2 * clausencmzeta(x, 2, 1) <=
    a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
where `p0 = 1 + α + (1 + α)^2 / 2` and
```
a0 = finda0(α) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))^2
```
holds for real `x` and all `α` in the interval `(-1, 0)`.

# Use
This is used in the proof of [`lemma_bhkdv_monotonicity_alpha`](@ref).

# Proof
Since `clausencmzeta` is `2π` periodic and even it is enough to prove
it for the `x` in the interval `[0, π]`.

**PROVE:** Finish this proof.
"""
lemma_bhkdv_monotonicity_main_term() = true
