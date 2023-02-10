"""
    _integrand_I_hat(x, t, α)

Compute
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
which is part of the integrand of `T01`.

**IMPROVE:** Optimize this for performance and enclosure in `x` and
  `α`. It seems to be that the integrand is monotone in `x` for `t`
  less than the root with `x = 0`.
"""
_integrand_I_hat(x, t, α) =
    clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)

"""
    _integrand_I_hat_series(x::Arb, t::Arb, α::Arb; M = 5)
    _integrand_I_hat_series(x::Arb, t::Arb, α::Arb, C::Arb, P::ArbSeries, E::Arb)

Compute
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
which is part of the integrand of `T01`, by expanding the Clausen
functions at zero.

It only supports `t > 0`.

In general this gives better enclosures than
[`_integrand_I_hat`](@ref) when `x` is small.

The version with `M` automatically computes the expansion, using `M`
terms. The version with `C, P, E` can be used with a precomputed
expansion. In that case they should be computed as
```
C, _, P, E = clausenc_expansion(x * (1 + t), -α, M, skip_constant = true)
```
Alternatively you can precompute `C` and `P`, which don't depend on
`x` and `t` as above and then compute
```
E = clausenc_expansion_remainder(x * (1 + t), -α, M)
```

**IMPROVE:** Look at using this as the default version for small
values of `x`.
"""
function _integrand_I_hat_series(x::Arb, t::Arb, α::Arb; M = 5)
    C, _, P, E = clausenc_expansion(x * (1 + t), -α, M, skip_constant = true)

    return _integrand_I_hat_series(x, t, α, C, P, E)
end

function _integrand_I_hat_series(x::Arb, t::Arb, α::Arb, C::Arb, P::ArbSeries, E::Arb)
    M = Arblib.degree(P) ÷ 2 + 1

    # t1 = abs(1 - t)
    t1 = 1 - t
    Arblib.abs!(t1, t1)
    # t2 = 1 + t
    t2 = 1 + t

    res = C * ArbExtras.enclosure_series(-α - 1, degree = 2) do e
        x^e * (t1^e + t2^e - 2t^e)
    end

    term = zero(x)
    tmp = zero(x)

    x2 = x^2

    # t1 = abs(1 - t)^2
    Arblib.pow!(t1, t1, UInt(2))
    # t2 = (1 + t)^2
    Arblib.pow!(t2, t2, UInt(2))
    # t3 = t^2
    t3 = t^2

    for m = 1:M-1
        # term = (1 - t)^2m + (1 + t)^2m - 2t^2m
        Arblib.pow!(term, t1, UInt(m))
        Arblib.pow!(tmp, t2, UInt(m))
        Arblib.add!(term, term, tmp)
        Arblib.pow!(tmp, t3, UInt(m))
        Arblib.submul!(term, tmp, 2)

        # term *= x^2m
        Arblib.pow!(tmp, x2, UInt(m))
        Arblib.mul!(term, term, tmp)

        Arblib.addmul!(res, Arblib.ref(P, 2m), term)
    end

    # term = (1 - t)^2M + (1 + t)^2M - 2t^2M
    Arblib.pow!(term, t1, UInt(M))
    Arblib.pow!(tmp, t2, UInt(M))
    Arblib.add!(term, term, tmp)
    Arblib.pow!(tmp, t3, UInt(M))
    Arblib.submul!(term, tmp, 2)

    # term *= x^2M
    Arblib.pow!(tmp, x2, UInt(M))
    Arblib.mul!(term, term, tmp)

    Arblib.addmul!(res, E, term)

    return res
end

"""
    _integrand_I_hat_dt(x, t, α)

Compute the function [`_integrand_I_hat`](@ref) differentiated once
with respect to `t`.
"""
_integrand_I_hat_dt(x, t, α) =
    -x * (
        -clausens(x * (1 - t), -α - 1) + clausens(x * (1 + t), -α - 1) -
        2clausens(x * t, -α - 1)
    )

"""
    _integrand_compute_root(::Type{<:Fractionalkdvansatz}, x::Arb, α::Arb)

Compute the unique root of
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
in `t` on the interval `[0, 1]`. It assumes that `0 <= x <= π`.

The existence and uniqueness of the root is based on lemma
[`lemma_integrand_1`](@ref).

It uses that the root is decreasing in `x` to better handle wide input
intervals.

In practice the root is also decreasing in `α`, though this is harder
to prove in the general case. Instead we do the computations assuming
it is decreasing in `α` and then afterward verify that this indeed was
the case. If this verification fails it falls back to not assuming
that the root is decreasing in `α`.

The verification is done in the following way. Consider `αₗ, αᵤ =
getinterval(Arb, α)` and let `rₗ` be the root computed with `αᵤ` and
`rᵤ` the root computed with `αₗ`. If the derivative of
```
_integrand_I_hat(x, rₗ, α)
```
is positive then, since `_integrand_I_hat(x, rₗ, αᵤ) = 0`, we have
that `_integrand_I_hat(x, rₗ, α) <= 0`. Since the function is
increasing in `t` this ensures that the root is at `rₗ` or to the left
of it. Similarly we have that if the derivative of
```
_integrand_I_hat(x, rᵤ, α)
```
is positive then `_integrand_I_hat(x, rₗ, α) >= 0` and the root must
lie at or to the right of `rᵤ`.

If the lower bound of `x` is zero or close to zero (smaller than
`eps(Arb)`) it computes an upper bound of the root by considering the
limiting case as `x` goes to zero. The limiting root can then be bound
by computing the root of
```
(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
```

If the upper bound of `x` is close to zero, smaller than `eps(Arb)`,
we compute the root at `eps(Arb)` and use that as a lower bound. This
avoids computing with very small values of `x`.
"""
function _integrand_compute_root(::Type{<:FractionalKdVAnsatz}, x::Arb, α::Arb)
    # Compute the root at x = 0
    root_zero = let
        # Function we are computing root of
        f = t -> (1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)

        roots, flags = ArbExtras.isolate_roots(f, Arf(0.5), Arf(0.9))
        length(flags) == 1 && flags[1] || error("could not isolate root for x = 0")

        ArbExtras.refine_root(f, Arb(only(roots)))
    end

    # If x or α is wide we don't need to compute the root to very high
    # precision. We therefore take the absolute tolerance to depend on
    # their radius.
    atol = Arblib.mul_2exp!(Mag(), max(radius(x), radius(α)), -20)

    compute_root(x::Arb, α::Arb) =
        let
            # Function we are computing the root of
            f = t -> _integrand_I_hat(x, t, α)
            df = t -> _integrand_I_hat_dt(x, t, α)

            # The root is lower bounded by 1 / 2
            r_lower = Arf(0.5)

            # Find a crude upper bound for the root
            δ = root_zero - r_lower # root_lower + δ gives upper bound of root
            if Arblib.midref(x) < 1e-5 && !Arblib.ispositive(f(r_lower + δ))
                # If x is very small then the sign might not be
                # determined at root_zero. In that case it is
                # beneficial to slightly increase the bound for
                # isolate_roots to be able to isolate the root.
                δ += sqrt(eps(δ))
            end
            while Arblib.ispositive(f(r_lower + δ / 2))
                Arblib.mul_2exp!(δ, δ, -1)
            end
            r_upper = ubound(r_lower + δ)

            # Short circuit in case the sign can't be determined on
            # the lower endpoint, this happens when x is very close to
            # π
            Arblib.contains_zero(f(Arb(r_lower))) && return Arb((r_lower, r_upper))

            # Improve the enclosure of the root
            roots, flags = ArbExtras.isolate_roots(f, r_lower, r_upper)
            if length(flags) == 1 && flags[1]
                # Refine the unique root
                r = ArbExtras.refine_root(f, Arb(only(roots)); df, atol)
            else
                r = Arb((roots[1][1], roots[end][2]))
            end

            return r
        end

    root_is_decreasing(x::Arb, t::Arb, α::Arb) =
        let
            dfdα = ArbExtras.enclosure_series(
                ArbExtras.derivative_function(α -> _integrand_I_hat(x, t, α)),
                α,
                degree = 2,
            )

            Arblib.ispositive(dfdα)
        end

    xₗ, xᵤ = getinterval(Arb, x)
    xᵤ = min(Arb(π), xᵤ) # We assume that xᵤ <= π
    αₗ, αᵤ = getinterval(Arb, α)

    ϵ = eps(Arb)

    if iszero(x)
        root = root_zero
    elseif Arblib.overlaps(xᵤ, Arb(π))
        rootₗ = Arb(1 // 2) # Lower bound is 1 / 2

        rootᵤ = compute_root(xₗ, αₗ)
        if !root_is_decreasing(xₗ, rootᵤ, α)
            rootᵤ = compute_root(xₗ, α)
        end

        root = Arb((rootₗ, rootᵤ))
    elseif !iswide(x)
        # In this case x never overlaps zero
        if iswide(α)
            rootₗ = compute_root(x, αᵤ)
            if !root_is_decreasing(x, rootₗ, α)
                rootₗ = compute_root(x, α)
            end

            rootᵤ = compute_root(x, αₗ)
            if !root_is_decreasing(x, rootᵤ, α)
                rootᵤ = compute_root(x, α)
            end

            root = Arb((rootₗ, rootᵤ))
        else
            root = compute_root(x, α)
        end
    elseif xᵤ < ϵ
        rootₗ = compute_root(ϵ, αᵤ)
        if !root_is_decreasing(ϵ, rootₗ, α)
            rootₗ = compute_root(ϵ, α)
        end

        rootᵤ = root_zero

        root = Arb((rootₗ, rootᵤ))
    elseif xₗ < ϵ
        rootₗ = compute_root(xᵤ, αᵤ)
        if !root_is_decreasing(xᵤ, rootₗ, α)
            rootₗ = compute_root(xᵤ, α)
        end

        rootᵤ = root_zero

        root = Arb((rootₗ, rootᵤ))
    else
        rootₗ = compute_root(xᵤ, αᵤ)
        if !root_is_decreasing(xᵤ, rootₗ, α)
            rootₗ = compute_root(xᵤ, α)
        end

        rootᵤ = compute_root(xₗ, αₗ)
        if !root_is_decreasing(xₗ, rootᵤ, α)
            rootᵤ = compute_root(xₗ, α)
        end

        root = Arb((rootₗ, rootᵤ))
    end

    return root
end

"""
    T01(u0::FractionalKdVAnsatz{Arb}, ::Ball; δ1::Arb = Arb(1e-3), skip_div_u0 = false)

Return a function such that `T01(u0, Ball(); δ1)(x)` computes the
integral
```
inv(π * u0(x) * u0.w(x)) * x * ∫ abs(_integrand_I_hat(x, t, α)) * u0.w(x * t) dt
```
where the integration is taken from `0` to `1`.

The interval of integration is split into two parts, one from `0` to
`1 - δ1` and one from `1 - δ1` to `1`. The first part is handled using
[`T012`](@ref) and the second one using [`T013`](@ref).

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.
"""
function T01(
    u0::FractionalKdVAnsatz{Arb},
    evaltype::Ball;
    δ1::Arb = Arb(1e-3),
    skip_div_u0 = false,
)
    f = T012(u0, evaltype, skip_div_u0 = true; δ1)
    g = T013(u0, evaltype, skip_div_u0 = true; δ1)

    return x::Arb -> begin
        fx = f(x)

        isfinite(fx) || return indeterminate(x)

        res = fx + g(x)

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T01(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M = 5, ϵ = one(Arb))

Returns a function for computing an **upper bound** of the integral of
`T01(u0)`, using an evaluation strategy that works asymptotically as
`x` goes to zero.

# Arguments
- `M::Integer` determines the number of terms in the expansions.
- `ϵ::Arb` determines the interval ``[-ϵ, ϵ]`` on which the expansion
  is valid.
- `return_enclosure::Bool` if true it returns an enclosure instead of
  an upper bound by returning the interval between zero and the
  computer upper bound.

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
d = 2sum(1:N-1) do m
        (-1)^m * zeta(-α - 2m) * ϵ^(2m - 2) / factorial(2m) *
            sum(binomial(2m, 2k) / (2k + 1 + p) for k = 0:m-1)
    end +
    1 / ϵ^2 * sum((-1)^m * zeta(-α - 2m) * (2ϵ)^2m / factorial(2m) for m = N:Inf)
```
for any `N >= 1`. Here `r` is the unique root of
```
(1 - t)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)
```
on ``[0, 1]``and is computed by [`_integrand_compute_root`](@ref).

To compute an enclosure of the tail of `d` we note that it is the same
as the sum in [`clausenc_expansion_remainder`](@ref) with `x = 2ϵ`.
"""
function T01(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
    return_enclosure::Bool = false,
)
    @assert !u0.use_bhkdv # Assumes the weight is abs(x)^u0.p
    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    α = u0.α
    p = u0.p

    r = _integrand_compute_root(typeof(u0), zero(α), α)

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

    N = 10
    # Sum first N - 1 terms
    d =
        2sum(1:N-1) do m
            (-1)^m * zeta(-α - 2m) * ϵ^(2m - 2) / factorial(2m) *
            sum(binomial(2m, 2k) / (2k + 1 + p) for k = 0:m-1)
        end
    # Enclose remainder
    # Note that 4(2ϵ)^(2M - 2) = (2ϵ)^(2M) / ϵ^2
    d += 4(2ϵ)^(2N - 2) * clausenc_expansion_remainder(2ϵ, -α, N)

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && Arblib.ref(x, 0) <= ϵ)

        res = inv_u0(x) * (c + d * abspow(x, 3 + α)) / π

        if return_enclosure
            return union(zero(res), res)
        else
            return res
        end
    end
end

"""
    T012(u0::FractionalKdVAnsatz{Arb}, ::Ball; δ1::Arb = skip_div_u0 = false)

Return a functions such that `T012(u0; δ1)(x)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * x * ∫ abs(_integrand_I_hat(x, t, α)) * u0.w(x * t) dt
```
where the integration is taken from `0` to `1 - δ1`.

If `weightfactors(u0)` is true then `u0.w(x * t) = u0.w(x) * u0.w(t)`
and this can be simplified to
```
inv(π * u0(x)) * x * ∫ abs(_integrand_I_hat(x, t, α)) * u0.w(t) dt
```
**IMPROVE:** The split depending on if `weightfactors(u0)` is true or
not is currently quite ugly and could possibly be improved.

As a first step the unique root of [`_integrand_I_hat`](@ref) is
computed using [`_integrand_compute_root`](@ref). Outside the
enclosure of the root the function has a constant sign and the
absolute value can hence be moved outside of the integral there. At
the root we use a naive enclosure given by the diameter of the
interval times the enclosure of the integrand, in general this will be
very small.

Outside of the enclosure of the root we use [`Arblib.integrate`](@ref)
for the integration.

The integrand is not analytic at the endpoint `t = 0`. For computing
an enclosure the only problematic part of the integrand is the term
`clausenc(x * t, -α) * u0.w(t). To enclose it we compute an expansion
of `clausenc` and explicitly handle the multiplication with `u0.w(t)`.
This enclosure is only finite if `u0.p - u0.α - 1 > 0` , which is true
in all cases we consider.

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.
"""
function T012(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball = Ball();
    δ1::Arb = Arb(1e-3),
    skip_div_u0 = false,
)
    return x::Arb -> begin
        root = _integrand_compute_root(typeof(u0), x, u0.α)
        root_lower, root_upper = getinterval(Arb, root)

        # Compute expansion of clausenc, used for enclosure when t
        # overlaps zero. We only use it for t < ϵ, in which case the
        # argument is bounded by x * ϵ.
        ϵ = Arb(1e-2)
        M = 3
        C, e, P, E = clausenc_expansion(x * ϵ, -u0.α, M)

        # The integrand without the absolute value. Assuming that
        # weightfactors(u0) is true.
        integrand_no_abs_factors(t; analytic = false) = begin
            rt = real(t)

            if analytic && !(Arblib.ispositive(rt) && rt < 1)
                # Only analytic for 0 < t < 1
                return indeterminate(t)
            elseif Arblib.contains_zero(t)
                # In this case analytic is false and t is real.
                @assert !analytic && isreal(t)

                # Only compute an enclosure if t < ϵ, otherwise it is
                # not good enough anyway.
                if rt < ϵ
                    # Part of the integrand which is well behaved around t = 0
                    part1 =
                        (clausenc(x * (1 - rt), -u0.α) + clausenc(x * (1 + rt), -u0.α)) * u0.w(rt)

                    # Enclosure of clausenc(x * t, -α) * u0.w(t)
                    part2 =
                        C * abspow(x, e) * u0.wmulpow(rt, e) +
                        P(x * rt) * u0.w(rt) +
                        E * abspow(x * rt, 2M) * u0.w(rt)

                    return Acb(part1 - 2part2)
                else
                    return indeterminate(t)
                end
            else
                if isreal(t)
                    return _integrand_I_hat(x, rt, u0.α) * u0.w(t)
                else
                    return _integrand_I_hat(x, t, u0.α) * u0.w(t)
                end
            end
        end

        # The integrand without the absolute value. Assuming that
        # weightfactors(u0) is false.
        integrand_no_abs_doesnt_factor(t; analytic = false) = begin
            rt = real(t)

            if analytic && !(Arblib.ispositive(rt) && rt < 1)
                # Only analytic for 0 < t < 1
                return indeterminate(t)
            elseif Arblib.contains_zero(t)
                # In this case analytic is false and t is real.
                @assert !analytic && isreal(t)

                # Only compute an enclosure if t < ϵ, otherwise it is
                # not good enough anyway.
                if rt < ϵ
                    # Part of the integrand which is well behaved around t = 0
                    part1 =
                        (clausenc(x * (1 - rt), -u0.α) + clausenc(x * (1 + rt), -u0.α)) * u0.w(x * rt)

                    # Enclosure of clausenc(x * t, -α) * u0.w(t)
                    part2 =
                        C * u0.wmulpow(x * rt, e) +
                        P(x * rt) * u0.w(x * rt) +
                        E * abspow(x * rt, 2M) * u0.w(x * rt)

                    return Acb(part1 - 2part2)
                else
                    return indeterminate(t)
                end
            else
                if isreal(t)
                    return _integrand_I_hat(x, rt, u0.α) * u0.w(x * t)
                else
                    return _integrand_I_hat(x, t, u0.α) * u0.w(x * t)
                end
            end
        end

        integrand_no_abs = ifelse(
            weightfactors(u0),
            integrand_no_abs_factors,
            integrand_no_abs_doesnt_factor,
        )

        res1 = abs(
            real(
                Arblib.integrate(
                    integrand_no_abs,
                    0,
                    root_lower,
                    check_analytic = true,
                    rtol = 2e-5,
                    atol = 2e-5,
                    warn_on_no_convergence = false,
                    opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
                ),
            ),
        )

        # Naive enclosure
        res2 =
            (root_upper - root_lower) * abs(integrand_no_abs(Arb((root_lower, root_upper))))

        res3 = abs(
            real(
                Arblib.integrate(
                    integrand_no_abs,
                    root_upper,
                    1 - δ1,
                    check_analytic = true,
                    rtol = 2e-5,
                    atol = 2e-5,
                    warn_on_no_convergence = false,
                    opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
                ),
            ),
        )

        res = (res1 + res2 + res3) * x / π

        if !weightfactors(u0)
            res /= u0.w(x)
        end

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T013(u0::FractionalKdVAnstaz{Arb}, ::Ball; δ1::Arb = Arb(1e-3), skip_div_u0 = false,)

Return a function such that `T013(u0; δ1)(x)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * x * ∫ _integrand_I_hat(x, t, α) * u0.w(x * t) dt
```
where the integration is taken from `1 - δ1` to `1`.

Since `u0.w` is bounded on the interval of integration we can factor
it out as
```
inv(π * u0(x) * u0.w(x)) * x * u0.w(x * Arb((1 - δ1, 1))) * ∫ _integrand_I_hat(x, t, α) dt
```
and compute the integral explicitly. It is given by
```
∫ _integrand_I_hat(x, t, α) dt = inv(x) * (
    clausens(2x, 1 - α) -
    2clausens(x, 1 - α) +
    clausens(x * δ1, 1 - α) -
    clausens(x * (2 - δ1), 1 - α) +
    2clausens(x * (1 - δ1), 1 - α)
)
```

If `weightfactors(u0)` is true then `u0.w(x * t) = u0.w(x) * u0.w(t)`
and we can simplify the result to
```
inv(π * u0(x)) * u0.w(Arb((1 - δ1, 1))) * (
    clausens(2x, 1 - α) -
    2clausens(x, 1 - α) +
    clausens(x * δ1, 1 - α) -
    clausens(x * (2 - δ1), 1 - α) +
    2clausens(x * (1 - δ1), 1 - α)
)
```

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.

**IMPROVE:** An earlier version of this method worked by expanding the
integrand and integrating the expansion. This was more complicated and
in general it gave worse results. However, for `x` close to zero it
could give better results. If need be we could reintroduce this
method, though it is most likely a better idea to do that in the
`::Asymptotic` version instead.

**IMPROVE:** For `x` overlapping or very close to `π` the Clausen
functions are not differentiable and the computed enclosure is very
wide. This could possibly be improved by expanding them at zero.
"""
function T013(
    u0::FractionalKdVAnsatz{Arb},
    ::Ball = Ball();
    δ1::Arb = Arb(1e-3),
    skip_div_u0 = false,
)
    return x::Arb -> begin
        if weightfactors(u0)
            weight_factor = u0.w(Arb((1 - δ1, 1)))
        else
            weight_factor = u0.w(x * Arb((1 - δ1, 1))) / u0.w(x)
        end

        # Compute a tighter enclosure by expanding in x. If x is
        # closer to zero or π it is beneficial to use a higher degree
        # when computing the enclosure.
        degree = ifelse(0.5 < x < 3, 1, 4)
        s = 1 - u0.α
        clausen_factor = ArbExtras.enclosure_series(x; degree) do x
            clausens(2x, s) - 2clausens(x, s) + clausens(x * δ1, s) -
            clausens(x * (2 - δ1), s) + 2clausens(x * (1 - δ1), s)
        end

        res = weight_factor * clausen_factor / π

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
