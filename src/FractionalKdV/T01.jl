"""
    _integrand_compute_root(u0::Fractionalkdvansatz, x::Arb)

Compute the unique root of
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
in `t` on the interval `[0, 1]`. It assumes that `0 <= x <= π`.

The existence and uniqueness of the root is based on lemma
[`lemma_integrand_1`](@ref).

For wide values of `x` it uses that the root is decreasing in `x` to
only have to evaluate at the endpoints.

If the lower bound of `x` is zero or close to zero (smaller than
`eps(x)`) it computes an upper bound of the root by considering the
limiting case as `x` goes to zero. The limiting root can then be bound
by computing the root of
```
(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
```

If the upper bound of `x` is close to zero, smaller than `eps(x)`, we
compute the root at `eps(x)` and use that as a lower bound. This
avoids computing with very small values of `x`.
"""
function _integrand_compute_root(u0::FractionalKdVAnsatz, x::Arb)
    compute_root(x) =
        let
            # Function we are computing the root of
            f =
                t ->
                    clausenc(x * (1 - t), -u0.α) + clausenc(x * (1 + t), -u0.α) -
                    2clausenc(x * t, -u0.α)

            # The root is lower bounded by 1 / 2, take a value
            # slightly larger so that we can still isolate it even if
            # it touches 1 / 2.
            root_lower = Arf(0.5) - sqrt(eps(Arf))

            # Find a crude upper bound for the root
            δ = Arb(0.4)
            # IMPROVE: We can remove this check if we can prove that
            # the root is less than x + δ.
            Arblib.ispositive(f(root_lower + δ)) || return Arblib.indeterminate!(zero(x))
            while Arblib.ispositive(f(root_lower + δ / 2)) && δ > 1e-5
                Arblib.mul_2exp!(δ, δ, -1)
            end
            root_upper = ubound(root_lower + δ)

            # Improve the enclosure of the root
            roots, flags = ArbExtras.isolate_roots(f, root_lower, root_upper)
            if length(flags) == 1 && flags[1]
                # Refine the unique root
                root = ArbExtras.refine_root(f, Arb(only(roots)))
            else
                root = Arb((roots[1][1], roots[end][2]))
            end

            return root
        end

    compute_root_zero() =
        let
            # Function we are computing root of
            f = t -> (1 - t)^(-u0.α - 1) + (1 + t)^(-u0.α - 1) - 2t^(-u0.α - 1)

            roots, flags = ArbExtras.isolate_roots(f, Arf(0.5), Arf(0.9))
            length(flags) == 1 && flags[1] || error("could not isolate root for x = 0")

            ArbExtras.refine_root(f, Arb(only(roots)))
        end

    xₗ, xᵤ = getinterval(Arb, x)
    xᵤ = min(Arb(π), xᵤ) # We assume that xᵤ <= π
    ϵ = eps(x)

    if iszero(x)
        root = compute_root_zero()
    elseif !iswide(x)
        root = compute_root(x) # In this case x never overlaps zero
    elseif xᵤ < ϵ
        root = Arb((compute_root(ϵ), compute_root_zero()))
    elseif xₗ < eps(Arb)
        root = Arb((compute_root(xᵤ), compute_root_zero()))
    else
        root = Arb((compute_root(xᵤ), compute_root(xₗ)))
    end

    return root
end

"""
    T01(u0::FractionalKdVAnsatz{Arb}, ::Ball; δ1, δ2)

Returns a function such that `T01(u0, Ball(); δ1, δ2)(x)` computes the
integral ``T_{0,1}`` from the paper.

The integral is split into three parts, one on `[0, δ1]`, one on `[δ1,
1 - δ2]` and one on `[1 - δ2, 1]`. These are These are given by
[`T011`](@ref), [`T012`](@ref) and [`T013`](@ref).
"""
function T01(
    u0::FractionalKdVAnsatz{Arb},
    evaltype::Ball;
    δ0::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    δ1::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    skip_div_u0 = false,
)
    f = T011(u0, evaltype, skip_div_u0 = true; δ0)
    g = T012(u0, evaltype, skip_div_u0 = true; δ0, δ1)
    h = T013(u0, evaltype, skip_div_u0 = true; δ1)

    return x -> begin
        asymptotic_part = f(x) + h(x)

        # Short circuit on a non-finite result
        isfinite(asymptotic_part) || return asymptotic_part

        non_asymptotic_part = g(x)

        isfinite(non_asymptotic_part) || return non_asymptotic_part

        if skip_div_u0
            return asymptotic_part + non_asymptotic_part
        else
            return (asymptotic_part + non_asymptotic_part) / u0(x)
        end
    end
end

"""
    T01(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M = 3, ϵ = one(Arb))

Returns a function for computing an **upper bound** of the integral of
`T01(u0)`, using an evaluation strategy that works asymptotically as
`x` goes to zero.

# Arguments
- `M::Integer` determines the number of terms in the expansions.
- `ϵ::Arb` determines the interval ``[-ϵ, ϵ]`` on which the expansion
  is valid.

# Implementation
It first splits the function as
```
T01(u0)(x) = inv(π) * inv(u0(x) / x^-α) * (U01(x) / x^(-α + p))
```
where `α = u0.α`, `p = u0.p` and
```
U01(x) = x^(1 + p) * ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t^u0.p
```
integrated from `0` to `1`. The factor `inv(u0(x) / x^-α)` is computed
using [`inv_u0_normalised`](@ref) so the remaining work is in bounding
the `U01 / x^(-α + p)` factor.

From the lemma in the paper we have
```
U01 / x^(-α + p) <= c + d * x^(2 + α - p)
```
with
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
```
and
```
d = 4sum((-1)^m * zeta(-α - 2m) * (2ϵ)^(2m - 2) / factorial(2m) for m = 1:Inf)
```
Here `r` is the unique root of
```
(1 - t)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)
```
on ``[0, 1]``and is computed by [`_integrand_compute_root`](@ref).

To compute an upper bound of `d` we sum the first `M - 1` terms and
for the tail we rewrite it as
```
4sum((-1)^m * zeta(-α - 2m) * (2ϵ)^(2m - 2) / factorial(2m) for m = M:Inf)
= 4(2ϵ)^(2M - 2) * (sum((-1)^m * zeta(-α - 2m) * (2x)^(2m) / factorial(2m) for m = M:Inf) / (2ϵ)^(2M))
```
Using [`clausenc_expansion_remainder`](@ref) we can compute an upper
bound of the sum divided by `(2ϵ)^2M`.

- **IMPROVE:** Include more details here, most of them are in the
  paper.
- **IMPROVE:** We could get a better bound for `d` by summing the true
  form of the first few terms instead of the simplified version used
  now. See the paper for how the true form look like.
"""
function T01(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M::Integer = 3, ϵ::Arb = Arb(1))
    α = u0.α
    p = u0.p

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    r = _integrand_compute_root(u0, Arb(0))
    c =
        gamma(1 + α) *
        sinpi(-α / 2) *
        (
            2 / (α - p) +
            gamma(-α) * gamma(1 + p) / gamma(1 - α + p) +
            hypgeom_2f1(1 + α, 1 + p, 2 + p, -one(α)) / (1 + p) -
            2r^p * (
                2r^-α / (α - p) +
                r * hypgeom_2f1(1 + α, 1 + p, 2 + p, -r) / (1 + p) +
                r * hypgeom_2f1(1 + α, 1 + p, 2 + p, r) / (1 + p)
            )
        )

    # Sum first M - 1 terms
    d = 4sum((-1)^m * zeta(-α - 2m) * (2ϵ)^(2m - 2) / factorial(2m) for m = 1:M-1)
    # Enclose remainder
    d += 4(2ϵ)^(2M - 2) * clausenc_expansion_remainder(2ϵ, -α, M)

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && Arblib.ref(x, 0) <= ϵ)

        return inv_u0(x) * (c + d * abspow(x, 3 + α)) / π
    end
end

"""
    T011(u0::FractionalKdVAnstaz{Arb}, evaltype = Ball(); δ0)

Returns a function such that `T011(u0; δ0)(x)` computes the integral
``T_{0,1,1}`` from the paper.

The strategy for evaluation is to compute an expansion for the
integrand which is then integrated termwise on the interval. The first
two terms inside the absolute value are analytic and their Taylor
expansions of degree `N - 1` are computed at `t = 0` and the error
term is enclosed. The last term inside the absolute value is not
analytic and instead we use the expansion for Clausians from the paper
and also here enclose the error term.

The value inside the absolute value has constant sign so we can remove
it. Switching integration and summation gives us terms of the form
```
∫_0^δ0 t^s * t^p dt
```
where `s` depends on the term and `p = u0.p`. This is easily
calculated to be
```
δ0^(s + p + 1 / (s + p + 1)
```
The errors are handled as constant values which are just multiplied by
the length of the interval.

**Prove:** that the expression inside the absolute value of the
integrand is of constant sign and determine the sign.
"""
function T011(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball = Ball();
    δ0::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    N::Integer = 3,
    skip_div_u0 = false,
)
    Γ = gamma
    α = u0.α
    δ0 = Arb(δ0)

    M = N ÷ 2 + 1

    return x -> begin
        # Analytic terms
        (P1, P1_E) = taylor_with_error(zero(α), union(zero(α), δ0), N) do t
            clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α)
        end
        P1_restterm = Arblib.add_error!(zero(α), P1_E * δ0^N)

        # Singular term
        (C, e, P2, P2_E) = clausenc_expansion(x * δ0, -α, M)
        C *= x^e
        for m = 1:M-1
            P2[2m] *= x^(2m)
        end
        P2_restterm = P2_E * (x * δ0)^(2M)

        # Compute the integral

        # Integrate the singular term
        singular_term = -2C * δ0^(e + u0.p + 1) / (e + u0.p + 1)

        # Integrate the analytic terms
        analytic_term = zero(α)
        full_series = P1 - 2P2
        for i = 0:N-1
            analytic_term += full_series[i] * δ0^(i + u0.p + 1) / (i + u0.p + 1)
        end

        res = singular_term + analytic_term

        # Add the error term
        res += δ0 * (P1_restterm - P2_restterm)

        if skip_div_u0
            return -res * x / π
        else
            return -res * x / (π * u0(x))
        end
    end
end

"""
    T012(u0::FractionalKdVAnsatz{Arb}, evaltype = Ball(); δ0, δ1)

Returns a function such that `T012(u0; δ0, δ1)(x)` computes the
integral ``T_{0,1,2}`` from the paper. This is the integral
```
x / (π * u0(x)) * ∫abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t^p dt
```
from `δ0` to `1 - δ1`.

To speed up the integration we want to avoid the absolute value inside
the integrand. The expression inside the absolute value is given by
```
f(t) = clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
and it has a unique root on the interval `[0, 1]` which we can isolate
using [`_integrand_compute_root`](@ref). By isolating this root we can
skip the absolute value for the parts of the interval where there is
no roots. For the part enclosing the root we keep the absolute value.
"""
function T012(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball = Ball();
    δ0::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    δ1::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    skip_div_u0 = false,
)
    mα = -u0.α
    cp = Acb(u0.p)

    a = Acb(δ0)
    b = 1 - Acb(δ1)

    return x -> begin
        integrand(t; analytic) = begin
            # Check that the real part of t is strictly between 0 and
            # 1 or return an indeterminate result
            Arblib.ispositive(Arblib.realref(t)) && Arblib.realref(t) < 1 ||
            return Arblib.indeterminate!(zero(t))

            if isreal(t)
                rt = Arblib.realref(t)

                res = Acb(
                    clausenc(x * (1 - rt), mα) + clausenc(x * (1 + rt), mα) -
                    2clausenc(x * rt, mα),
                )
            else
                res =
                    clausenc(x * (1 - t), mα) + clausenc(x * (1 + t), mα) -
                    2clausenc(x * t, mα)
            end

            return Arblib.real_abs!(res, res, analytic) * t^cp
        end

        # The same integrand as above but without the absolute value.
        # It is used for the parts of the interval where the integrand
        # is known to be non-zero.
        integrand_no_abs(t) = begin
            # Check that the real part of t is strictly between 0 and
            # 1 or return an indeterminate result
            Arblib.ispositive(Arblib.realref(t)) && Arblib.realref(t) < 1 ||
            return Arblib.indeterminate!(zero(t))

            if isreal(t)
                rt = Arblib.realref(t)

                return Acb(
                    (
                        clausenc(x * (1 - rt), mα) + clausenc(x * (1 + rt), mα) -
                        2clausenc(x * rt, mα)
                    ) * rt^u0.p,
                )

            else
                return (
                    clausenc(x * (1 - t), mα) + clausenc(x * (1 + t), mα) -
                    2clausenc(x * t, mα)
                ) * t^cp
            end
        end

        if false
            # Attempt to isolate the root of clausenc(x * (1 - t), mα) +
            # clausenc(x * (1 + t), mα) - 2clausenc(x * t, mα) using that
            # it is lower bounded by 0.5.
            f(t) =
                clausenc(x * (1 - t), mα) + clausenc(x * (1 + t), mα) - 2clausenc(x * t, mα)

            roots, flags = ArbExtras.isolate_roots(f, Arf(0.5), ubound(real(b)))

            if length(flags) == 1 && flags[1]
                # Refine the unique root
                root = ArbExtras.refine_root(f, Arb(only(roots)))
                # Get lower and upper bounds for the root
                root_lower, root_upper = getinterval(root)

                isolated_root = true
            else
                # Get lower and upper bounds for possible roots
                root_lower = roots[1][1]
                root_upper = roots[end][2]

                isolated_root = false
            end
        else
            root = _integrand_compute_root(u0, x)

            root_lower, root_upper = getinterval(root)
        end

        res1 = abs(
            real(
                Arblib.integrate(
                    integrand_no_abs,
                    a,
                    root_lower,
                    check_analytic = false,
                    rtol = 2e-5,
                    atol = 2e-5,
                    warn_on_no_convergence = false,
                    opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
                ),
            ),
        )

        res2 = real(
            Arblib.integrate(
                integrand,
                root_lower,
                root_upper,
                check_analytic = true,
                rtol = 2e-5,
                atol = 2e-5,
                warn_on_no_convergence = false,
                opts = Arblib.calc_integrate_opt_struct(0, 100, 0, 0, 0),
            ),
        )

        res3 = abs(
            real(
                Arblib.integrate(
                    integrand_no_abs,
                    root_upper,
                    b,
                    check_analytic = false,
                    rtol = 2e-5,
                    atol = 2e-5,
                    warn_on_no_convergence = false,
                    opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
                ),
            ),
        )

        res = res1 + res2 + res3

        if skip_div_u0
            return res * x / π
        else
            return res * x / (π * u0(x))
        end
    end
end

"""
    T013(u0::FractionalKdVAnstaz{Arb}, evaltype = Ball(); δ1)

Returns a function such that `T013(u0; δ1)(x)` computes the integral
``T_{0,1,3}`` from the paper.

The strategy for evaluation is the same as for [`T011`](@ref) except
that the first term is singular and the last two are analytic and
their Taylor expansion is computed at `t = 1`.

The integral that needs to be computed in this case is
```
∫_(1 - δ1)^1 (1 - t)^s * t^p dt
```
which is given by
```
Γ(1 + s) * Γ(1 + p) / Γ(2 + s + p) - B(1 + p, 1 + s; 1 - δ1)
```
where `B(a, b; z)` is the incomplete Beta-function.

If `x` is equal or very close to π (determined by `ϵ`) then the Taylor
expansion gives a very poor approximation for `clausenc(x * (t + 1),
-α)`. In this case we make use of the fact that it's 2π periodic and
even, so that
```
clausenc(x * (t + 1), -α) = clausenc(x * (t + 1) - 2π, -α) = clausenc(2π - x * (t + 1), -α)
```
, to be able to use the asymptotic expansion instead. That gives us the
integral
```
∫_(1 - δ1)^1 (2π - x * (t + 1))^s * t^p dt
```
which is given by
```
(2π - x)^(1 + p + s) * x^(-1 - p) * (B(1 + p, 1 + s, x / (2π - x)) - B(1 + p, 1 + s, (x - δ1 * x) / (2π - x)))
```
The value of `x / (2π - x)` will always be less than or equal to 1 for
`x` less than or equal to π, however due to overestimation the
enclosing ball might contain values greater than one, we therefore
have to use [`beta_inc_zeroone`](@ref) to be able to get finite
results in that case.

**Prove:** that the expression inside the absolute value of the
integrand is of constant sign and determine the sign.
"""
function T013(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball = Ball();
    δ1::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    ϵ::Arb = Arb(1e-2),
    N::Integer = 3,
    skip_div_u0 = false,
)
    Γ = gamma
    α = u0.α
    δ1 = Arb(δ1)
    π = Arb(pi)

    M = N ÷ 2 + 1

    return x -> begin
        # Determine if the asymptotic expansion or the Taylor
        # expansion should be used for the second term
        use_asymptotic = π - x < ϵ

        # Analytic terms
        (P1, P1_E) = taylor_with_error(one(α), union(1 - δ1, one(α)), N) do t
            if !use_asymptotic
                return clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
            else
                return -2clausenc(x * t, -α)
            end
        end
        P1_restterm = Arblib.add_error!(zero(α), P1_E * δ1^N)

        # Singular term
        (C, e, P2, P2_E) = clausenc_expansion(x * δ1, -α, M)
        C *= x^e
        for m = 1:M-1
            P2[2m] *= x^(2m)
        end
        P2_restterm = P2_E * (x * δ1)^(2M)

        # Compute the integral

        # Integrate the singular term
        # Using ∫_(1 - δ1)^1 |t - 1|^s*t^(p) dt = ∫_(1 - δ1)^1 (1 - t)^s*t^(p) dt
        singular_term =
            C * (
                Γ(1 + e) * Γ(1 + u0.p) / Γ(2 + e + u0.p) -
                beta_inc(1 + u0.p, 1 + e, 1 - δ1)
            )

        # Integrate the analytic terms
        # Using ∫_(1-δ1)^1 (t-1)^i*t^(p) dt = (-1)^i ∫_(1-δ1)^1 (t-1)^i*t^(p) dt
        analytic_term = zero(α)
        full_series = P1 + P2
        for i = 0:N-1
            analytic_term +=
                full_series[i] *
                (-1)^i *
                (
                    Γ(Arb(1 + i)) * Γ(1 + u0.p) / Γ(2 + i + u0.p) -
                    beta_inc(1 + u0.p, Arb(1 + i), 1 - δ1)
                )
        end

        res = singular_term + analytic_term

        # Add the error term
        res += δ1 * (P1_restterm + P2_restterm)

        if use_asymptotic
            # Handle asymptotic expansion of clausenc(x*(t + 1), -α)
            (C, e, P3, P3_E) = clausenc_expansion(2π - x * (2 - δ1), -α, M)
            P3_restterm = P3_E * (2π - x * (2 - δ1))^(2M)

            # Add the singular part to the integral
            singular_term_2 =
                C *
                (2π - x)^(1 + u0.p + e) *
                x^(-1 - u0.p) *
                (
                    beta_inc_zeroone(1 + u0.p, 1 + e, x / (2π - x)) -
                    beta_inc_zeroone(1 + u0.p, 1 + e, (x - δ1 * x) / (2π - x))
                )

            for i = 0:2:N-1
                # Only even terms
                singular_term_2 +=
                    P3[i] *
                    (2π - x)^(1 + u0.p + i) *
                    x^(-1 - u0.p) *
                    (
                        beta_inc_zeroone(1 + u0.p, Arb(1 + i), x / (2π - x)) -
                        beta_inc_zeroone(1 + u0.p, Arb(1 + i), (x - δ1 * x) / (2π - x))
                    )
            end

            # Add rest term
            singular_term_2 += δ1 * P3_restterm

            # Add to res
            res += singular_term_2
        end

        if skip_div_u0
            return res * x / π
        else
            return res * x / (π * u0(x))
        end
    end
end
