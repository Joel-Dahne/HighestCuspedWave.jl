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

It uses that the root is increasing in `α` and decreasing in `x` to
better handle wide input intervals.

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
            root_lower = Arf(0.5)

            # Find a crude upper bound for the root
            # root_lower + δ gives upper bound of root
            δ = root_zero - root_lower
            while Arblib.ispositive(f(root_lower + δ / 2))
                Arblib.mul_2exp!(δ, δ, -1)
            end
            root_upper = ubound(root_lower + δ)

            # Short circuit in case the sign can't be determined on
            # the lower endpoint, this happens when x is very close to
            # π
            Arblib.contains_zero(f(Arb(root_lower))) && return Arb((root_lower, root_upper))

            # Improve the enclosure of the root
            roots, flags = ArbExtras.isolate_roots(f, root_lower, root_upper)
            if length(flags) == 1 && flags[1]
                # Refine the unique root
                root = ArbExtras.refine_root(f, Arb(only(roots)); df, atol)
            else
                root = Arb((roots[1][1], roots[end][2]))
            end

            return root
        end

    xₗ, xᵤ = getinterval(Arb, x)
    xᵤ = min(Arb(π), xᵤ) # We assume that xᵤ <= π
    αₗ, αᵤ = getinterval(Arb, α)

    ϵ = eps(Arb)

    if iszero(x)
        root = root_zero
    elseif Arblib.overlaps(xᵤ, Arb(π))
        root = Arb((1 // 2, compute_root(xₗ, αₗ))) # Lower bound is 1 / 2
    elseif !iswide(x)
        # In this case x never overlaps zero
        if iswide(α)
            root = Arb((compute_root(x, αᵤ), compute_root(x, αₗ)))
        else
            root = compute_root(x, α)
        end
    elseif xᵤ < ϵ
        root = Arb((compute_root(ϵ, αᵤ), root_zero))
    elseif xₗ < ϵ
        root = Arb((compute_root(xᵤ, αᵤ), root_zero))
    else
        root = Arb((compute_root(xᵤ, αᵤ), compute_root(xₗ, αₗ)))
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

Using that `u0.w(x * t) = u0.w(x) * u0.w(t)` this can be simplified to
```
inv(π * u0(x)) * x * ∫ abs(_integrand_I_hat(x, t, α)) * u0.w(t) dt
```

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
of `clausenc` and explicitly handle the multiplication with `u0.w(t) =
t^u0.p`. This enclosure is only finite if `u0.p - u0.α - 1 > 0` ,
which is true in all cases we consider.

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

        # The integrand without the absolute value.
        integrand_no_abs(t; analytic = false) = begin
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

                    # Enclosure of clausenc(x * t, -α) * t^u0.p
                    part2 =
                        C * abspow(x, e) * abspow(rt, e + u0.p) +
                        P(x * rt) * abspow(rt, u0.p) +
                        E * abspow(x * rt, 2M) * abspow(rt, u0.p)

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
Using that `u0.w(x * t) = u0.w(x) * u0.w(t)` we can simplify the
result to
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
        weight_factor = u0.w(Arb((1 - δ1, 1)))

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
