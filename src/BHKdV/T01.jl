"""
    T01(u0::BHKdVAnsatz, ::Ball; δ1, δ2, skip_div_u0)

Returns a function such that `T01(u0, Ball(); δ1, δ2)(x)` computes the
integral ``T_{0,1}`` from the paper.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T01(
    u0::BHKdVAnsatz,
    evaltype::Ball;
    δ0::Arb = Arb(1e-5),
    δ1::Arb = Arb(1e-5),
    skip_div_u0 = false,
)
    f = T011(u0, evaltype, skip_div_u0 = true; δ0)
    g = T012(u0, evaltype, skip_div_u0 = true; δ0, δ1)
    h = T013(u0, evaltype, skip_div_u0 = true; δ1)

    if skip_div_u0
        return x -> f(x) + g(x) + h(x)
    else
        return x -> (f(x) + g(x) + h(x)) / u0(x)
    end
end

"""
    T011(u0::BHKdVAnsatz; δ0)

Computes the integral ``T_{0,1,1}`` from the paper.

It uses the fact that the integrand is strictly increasing on the
interval `[0, 0.05]` for every value of `x` and 0 at `x = 0`. This
allows us to enclose the integrand on the interval which then easily
gives an enclosure of the integral by multiplying with the size of the
interval.

- **PROVE**: That the integrand indeed is increasing on the said
  interval.
"""
function T011(u0::BHKdVAnsatz, ::Ball = Ball(); δ0::Arb = Arb(1e-5), skip_div_u0 = false)
    δ0 < 0.05 || Throw(ArgumentError("δ0 must be less than 0.05, got $δ0"))

    #Enclosure of Arb((-1, -1 + u0.ϵ)) computed in a way so that
    #the lower endpoint is exactly -1
    α = -1 + Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    return x::Arb -> begin
        integrand(t) =
            abs(
                clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) -
                2clausenc(x * t, -α),
            ) *
            t *
            u0.wdivx(x * t)

        integral = δ0 * Arb((0, integrand(δ0)))

        res = integral * x / (π * u0.wdivx(x))
        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T012(u0::BHKdVAnsatz; δ0, δ1)

Returns a function such that `T012(u0; δ0, δ1)(x; tol)` computes the
integral ``T_{0,1,2}`` from the paper using the prescribed tolerance
in the integration.
"""
function T012(
    u0::BHKdVAnsatz,
    ::Ball = Ball();
    δ0::Arb = Arb(1e-5),
    δ1::Arb = Arb(1e-5),
    skip_div_u0 = false,
)
    # Enclosure of -α so that the upper bound is exactly 1
    mα = 1 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    # Integration limits
    a = δ0
    b = 1 - δ1

    return (x::Arb; tol = Arb(1e-5)) -> begin
        integrand(t) = begin
            xt = x * t

            term = abs(
                clausenc(x * (1 - t), mα) + clausenc(x * (1 + t), mα) - 2clausenc(xt, mα),
            )

            return term * t * u0.wdivx(xt)
        end

        res = ArbExtras.integrate(integrand, a, b, atol = tol, rtol = tol)

        res *= x / (π * u0.wdivx(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T013(u0::BHKdVAnsatz; δ1)

Computes the integral ``T_{0,1,3}`` from the paper.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval. We can also notice that the value inside the absolute value
is negative so we can remove the absolute value by putting a minus
sign, which we can bake in to the weight factor.

We are left with integrating the three Clausen terms
1. `clausenc(x * (1 - t), -α)`
2. `clausenc(x * (1 + t), -α)`
3. `2clausenc(x * t, -α)`
We have that the primitive functions for the three terms are given by
1. `-clausens(x * (1 - t), 1 - α) / x`
2. `clausens(x * (1 + t), 1 - α) / x`
3. `2clausens(x * t, 1 - α) / x`
Hence the integral from `1 - δ1` to `1` is
```
inv(x) * (
    (-clausens(0, 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α)) -
    (-clausens(x * δ1, 1 - α) + clausens(x * (2 - δ1), 1 - α) - 2clausens(x * (1 - δ1), 1 - α))
)
```
The multiplication by `inv(x)` can be cancelled by the multiplication
by `x` that is outside of the integral. Since `1 - α > 1` we have
`clausens(0, 1 - α) = 0`. If we also reorder the terms to more clearly
see which ones gives cancellations we get
```
clausens(x * δ1, 1 - α) +
(clausens(2x, 1 - α) - clausens(x * (2 - δ1), 1 - α)) -
2(clausens(x, 1 - α) - clausens(x * (1 - δ1), 1 - α))
```

In the case that `x` overlaps with `π` we get issues when evaluating
`(clausens(2x, 1 - α)` since it doesn't have a good implementation in
that case. We could subtract `2π` from the argument since it is `2π`
periodic, but that doesn't solve the issue since `clausens` currently
doesn't support evaluation on intervals containing zero.
- **FIXME:** Currently we do this by assuming that `clausens` is
  monotonic in `s`. In practice this is true for small enough
  arguments but not in general. If `x` is sufficiently close to `π`
  this will thus give a correct result, but it is not rigorously
  proved.

- **TODO:** Could improve enclosures by better handling cancellations
  for `clausens(2x, 1 - α) - clausens(x * (2 - δ1), 1 - α)` and
  `clausens(x, 1 - α) - clausens(x * (1 - δ1), 1 - α)`. Though this
  might not be needed.
"""
function T013(u0::BHKdVAnsatz, ::Ball = Ball(); δ1::Arb = Arb(1e-5), skip_div_u0 = false)
    return x -> begin
        x = convert(Arb, x)

        weight_factor = let t = Arb((1 - δ1, 1))
            -t * u0.wdivx(x * t)
        end

        # s = 1 - α
        s = Arb((2 - u0.ϵ, 2))

        integral = clausens(x * δ1, s) - 2(clausens(x, s) - clausens(x * (1 - δ1), s))
        if Arblib.overlaps(x, Arb(π))
            # FIXME: This assumes that clausens is monotonic on the
            # interval. In practice this is true for small enough
            # argument. But it is not true in general.

            # Compute an enclosure of clausens(2(x - π), s) on the
            # symmetric interval [-abs(2(x - π)), abs(2(x - π))] using
            # the oddness and assuming that the maximum is attained at
            # the endpoint.
            term = Arblib.add_error!(zero(x), clausens(abs_ubound(Arb, 2(x - Arb(π))), s))
            integral += term

            # Compute an enclosure of clausens(x * (2 - δ1) - 2π, s)
            # on the symmetric interval [-abs(x * (2 - δ1) - 2π),
            # abs(x * (2 - δ1) - 2π)] using the oddness and assuming
            # the maximum is attained at the endpoint.
            term = Arblib.add_error!(
                zero(x),
                clausens(abs_ubound(Arb, x * (2 - δ1) - 2Arb(π)), s),
            )
            integral -= term
        else
            integral += clausens(2x, s) - clausens(x * (2 - δ1), s)
        end
        integral *= weight_factor

        res = integral / (π * u0.wdivx(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
