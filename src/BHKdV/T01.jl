"""
    T01(u0::BHKdVAnsatz, ::Ball; δ, skip_div_u0)

Returns a function such that `T01(u0, Ball(); δ)(x)` computes the
integral ``T_{0,1}`` from the paper.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T01(u0::BHKdVAnsatz, evaltype::Ball; δ::Arb = Arb(1e-5), skip_div_u0 = false)
    f = T012(u0, evaltype, skip_div_u0 = true; δ)
    g = T013(u0, evaltype, skip_div_u0 = true; δ)

    if skip_div_u0
        return x -> f(x) + g(x)
    else
        return x -> (f(x) + g(x)) / u0(x)
    end
end

"""
    T012(u0::BHKdVAnsatz; δ)

Returns a function such that `T012(u0; δ)(x; tol)` computes the
integral ``T_{0,1,2}`` from the paper using the prescribed tolerance
in the integration.
"""
function T012(u0::BHKdVAnsatz, ::Ball = Ball(); δ::Arb = Arb(1e-5), skip_div_u0 = false)
    # Enclosure of Arb((-1, -1 + u0.ϵ)) computed in a way so that the
    # lower endpoint is exactly -1
    α = -1 + Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    return (x::Arb; tol = Arb(1e-5)) -> begin
        # Compute critical point of clausenc(x * t, -α) * t^0.5
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
                if t isa ArbSeries
                    tt = t[0]
                else
                    tt = t
                end

                # Check that clausenc(x * t, -α) * t * u0.wdivx(x * t)
                # is increasing
                if tt < t0
                    tᵤ = ubound(Arb, tt)

                    # TODO: Prove monotonicity of weight factor
                    part1 =
                        clausenc(x * (1 - tt), -α) +
                        clausenc(x * (1 + tt), -α) * Arb((0, tᵤ * u0.wdivx(x * tᵤ)))

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

"""
    T013(u0::BHKdVAnsatz; δ)

Computes the integral ``T_{0,1,3}`` from the paper.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval.

As long as `1 - δ` lies to the right of the unique root of the
integrand the value inside the absolute value is positive so we can
remove the absolute value. Since the root is increasing in `x` and
decreasing in `α` it is enough to check that `1 - δ` is to the right
of the root for `x = 0` and `α = 1`, for which the root is
`inv(sqrt(2))`.

We are left with integrating the three Clausen terms
1. `clausenc(x * (1 - t), -α)`
2. `clausenc(x * (1 + t), -α)`
3. `2clausenc(x * t, -α)`
We have that the primitive functions for the three terms are given by
1. `-clausens(x * (1 - t), 1 - α) / x`
2. `clausens(x * (1 + t), 1 - α) / x`
3. `2clausens(x * t, 1 - α) / x`
Hence the integral from `1 - δ` to `1` is
```
inv(x) * (
    (-clausens(0, 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α)) -
    (-clausens(x * δ, 1 - α) + clausens(x * (2 - δ), 1 - α) - 2clausens(x * (1 - δ), 1 - α))
)
```
The multiplication by `inv(x)` can be cancelled by the multiplication
by `x` that is outside of the integral. Since `1 - α > 1` we have
`clausens(0, 1 - α) = 0`. If we also reorder the terms to more clearly
see which ones gives cancellations we get
```
clausens(x * δ, 1 - α) +
(clausens(2x, 1 - α) - clausens(x * (2 - δ), 1 - α)) -
2(clausens(x, 1 - α) - clausens(x * (1 - δ), 1 - α))
```

- **IMPROVE:** Could improve enclosures by better handling
  cancellations for `clausens(2x, 1 - α) - clausens(x * (2 - δ), 1 -
  α)` and `clausens(x, 1 - α) - clausens(x * (1 - δ), 1 - α)`. Though
  this might not be needed.
"""
function T013(u0::BHKdVAnsatz, ::Ball = Ball(); δ::Arb = Arb(1e-5), skip_div_u0 = false)
    # Check that 1 - δ is to the right of the root of the integrand
    inv(sqrt(Arb(2))) < 1 - δ || error("interval of integration contains root")

    return x::Arb -> begin
        weight_factor = let t = Arb((1 - δ, 1))
            t * u0.wdivx(x * t)
        end

        # s = 1 - α
        s = Arb((2 - u0.ϵ, 2))

        integral =
            clausens(x * δ, s) + (clausens(2x, s) - clausens(x * (2 - δ), s)) -
            2(clausens(x, s) - clausens(x * (1 - δ), s))

        integral *= weight_factor

        res = integral / (π * u0.wdivx(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
