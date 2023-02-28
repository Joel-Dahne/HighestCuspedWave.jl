"""
    T012(u0::BHKdVAnsatz; δ, skip_div_u0)

Return a functions such that `T012(u0; δ)(x)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * x * ∫ abs(_integrand_I_hat(x, t, α)) * u0.w(x * t) dt
```
where the integration is taken from `0` to `1 - δ`.

As a first step we rewrite it as
```
inv(π * u0(x) * u0.wdivx(x)) * x * ∫ abs(_integrand_I_hat(x, t, α)) * t * u0.wdivx(x * t) dt
```

The integrand is not differentiable at the endpoint `t = 0`. For
computing an enclosure the only problematic part of the integrand is
the term `clausenc(x * t, -α) * u0.w(x * t). This is given by
```
clausenc(x * t, -α) * u0.w(x * t) = clausenc(x * t, -α) * abs(x * t)^(1 - u0.γ *
                                        (α + 1)) * log(u0.c + inv(abs(x * t)))
```
For `u0.γ == 1 // 2` and `u0.c = 2ℯ` there is a lemma in the paper
proving this to be increasing in `t` for `0 < t < t0` with
```
t0 = (-2gamma(1 + α) * sinpi(-α / 2) * x^(-α - 1) * (-α - 1 / 2) / zeta(-α))^inv(α + 1)
```
Furthermore it is zero at `t = 0`.

To evaluate `t0` we use that
```
gamma(1 + α) / zeta(-α) = gamma(2 + α) / ((1 + α) * zeta_deflated(-α, 1) - 1)
```
Giving us
```
t0 = (
        -2gamma(2 + α) * sinpi(-α / 2) * x^(-α - 1) * (-α - 1 / 2) /
        ((1 + α) * zeta_deflated(-α, 1) - 1)
    )^inv(α + 1)
```
Next we rewrite the power into `exp` and `log` as
```
t0 = exp(
    log(
        -2gamma(2 + α) * sinpi(-α / 2) * x^(-α - 1) * (-α - 1 / 2) /
        ((1 + α) * zeta_deflated(-α, 1) - 1)
    ) / (α + 1)
)
```
The division by `α + 1` has to be done taking into account the
removable singularity.
"""
function T012(u0::BHKdVAnsatz, ::Ball = Ball(); δ::Arb = Arb(1e-5), skip_div_u0 = false)
    # Enclosure of Arb((-1, -1 + u0.ϵ)) computed in a way so that the
    # lower endpoint is exactly -1
    α = -1 + Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    @assert u0.γ == 1 // 2
    @assert Arblib.overlaps(u0.c, 2Arb(ℯ))

    return (x::Arb; tol = Arb(1e-5)) -> begin
        # Compute t0
        t0 = begin
            a = TaylorModel(α, Arb(-1), degree = 2) do α
                num = -2gamma(2 + α) * sinpi(-α / 2) * x^(-α - 1) * (-α - 1 // 2)
                den = (1 + α) * zeta_deflated(-α, Arb(1)) - 1
                num / den
            end
            b = TaylorModel(α -> 1 + α, α, Arb(-1), degree = 2)

            compose(exp, div_removable(compose(log, a), b))(α)
        end

        integrand(t) = begin
            if (t isa Arb && Arblib.contains_zero(t)) ||
               (t isa ArbSeries && Arblib.contains_zero(t[0]))
                # If t is a ArbSeries only the constant term in the
                # result will be finite. We therefore only compute
                # with the constant part.
                t₀ = t isa ArbSeries ? t[0] : t

                # Check that clausenc(x * t, -α) * t * u0.wdivx(x * t)
                # is increasing
                if t₀ < t0
                    tᵤ = ubound(Arb, t₀)

                    # Use that t * u0.wdivx(x * t) is increasing in t.
                    # This follows from that the derivative has the
                    # same sign as x * t * u0.wdivx(x * t) = u0.w(x *
                    # t) and u0.w is increasing.
                    part1 =
                        clausenc(x * (1 - t₀), -α) +
                        clausenc(x * (1 + t₀), -α) * Arb((0, tᵤ * u0.wdivx(x * tᵤ)))

                    # Use that this term is increasing in t.
                    part2 = Arb((0, clausenc(x * tᵤ, -α) * tᵤ * u0.wdivx(x * tᵤ)))

                    if t isa Arb
                        return abs(part1 - 2part2)
                    elseif t isa ArbSeries
                        res = indeterminate(t)
                        res[0] = abs(part1 - 2part2)
                        return res
                    end
                else
                    return indeterminate(t)
                end
            else
                I =
                    clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) -
                    2clausenc(x * t, -α)

                return abs(I) * t * u0.wdivx(x * t)
            end
        end

        res = ArbExtras.integrate(integrand, Arb(0), 1 - δ, atol = tol, rtol = tol)

        res *= x / (π * u0.wdivx(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
