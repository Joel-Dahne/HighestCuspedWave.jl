"""
    T02(u0::BHKdVAnsatz; δ2, skip_div_u0)

Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral \$T_{0,2}\$ from the paper.

The interval of integration is `[x, π]`. Since the integrand is
singular at `y = x` we split the interval into two parts, `[x, a]` and
`[a, π]`. In principle we want to take `a = x + δ2`, however this
gives issues if `x` is a wide ball. Instead we take `a` to always be a
thin ball between `x` and `π`.

In general we take `a = ubound(x + δ2)`. However if `x` is wide then
it's beneficial to take a larger value than `δ2`, depending on the
radius of `x`. In practice we take `8radius(x)` in that case.

If `x` is very close to `π`, so that `a` would be larger than `π`,
then we use the asymptotic version for the whole interval `[x, π]`.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T02(u0::BHKdVAnsatz, evaltype::Ball; δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    f = T021(u0, evaltype, skip_div_u0 = true; δ2)
    g = T022(u0, evaltype, skip_div_u0 = true; δ2)

    return x -> begin
        x = convert(Arb, x)
        δ2 = max(δ2, 8Arblib.radius(Arb, x))
        a = Arblib.ubound(Arb, x + δ2)

        if !(a < π)
            return f(x, π)
        end

        res = f(x, a) + g(x, a)

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T021(u0::BHKdVAnsatz)

Computes the (not yet existing) integral \$T_{0,2,1}\$ from the paper.

The interval of integration is `[x, a]`. Both `x` and `a` are assumed
to be less than or equal to `π`, if they are balls which overlap `π`
anything above `π` will be ignored.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval. We can also notice that the value inside the absolute value
is positive so we can remove the absolute value.

We are left with integrating the three Clausen terms
1. `clausenc(x - y, -α)`
2. `clausenc(x + y, -α)`
3. `2clausenc(y, -α)`
We have that the primitive functions for the three terms are given by
1. `-clausens(x - y, 1 - α)`
2. `clausens(x + y, 1 - α)`
3. `2clausens(y, 1 - α)`
Hence the integral from `x` to `a` is
```
(-clausens(x - a, 1 - α) + clausens(x + a, 1 - α) - 2clausens(a, 1 - α)) -
(-clausens(0, 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α))
```
Since `1 - α > 1` we have
`clausens(0, 1 - α) = 0`. If we also reorder the terms to more clearly
see which ones gives cancellations we get
```
-clausens(x - a, 1 - α) +
(clausens(x + a, 1 - α) - clausens(2x, 1 - α)) -
2(clausens(a, 1 - α) - clausens(x, 1 - α))
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
  for `clausens(x + a, 1 - α) - clausens(2x, 1 - α)` and `clausens(a,
  1 - α) - clausens(x, 1 - α)`. Though this might not be needed.
"""
function T021(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)
        δ = a - x

        interval = union(x, a)

        weight_factor = u0.w(interval)

        # s = 1 - α
        s = Arb((2 - u0.ϵ, 2))

        integral = -clausens(x - a, s) - 2(clausens(a, s) - clausens(x, s))
        if Arblib.overlaps(x, Arb(π))
            # FIXME: This assumes that clausens is monotonic on the
            # interval. In practice this is true for small enough
            # argument. But it is not true in general.

            # Compute an enclosure of clausens(x + a, s) on the
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
                clausens(abs_ubound(Arb, x * (2 - δ2) - 2Arb(π)), s),
            )
            integral -= term
        else
            integral += clausens(x + a, s) - clausens(2x, s)
        end
        integral *= weight_factor

        res = integral / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T022(u0::BHKdVAnsatz)

Computes the (not yet existing) integral \$T_{0,2,2}\$ from the paper.

The interval of integration is given by `[a, π]`. In practice `a`
should be a thin ball to not give problems with the integration.

This is done by directly computing the integral with the integrator in
Arb.

Notice that the expression inside the absolute value is always
positive, so we can remove the absolute value.

**FIXME:** This currently assumes that
```
clausenc(y - x, s) + clausenc(y + x, s) - 2clausenc(y, s)
```
and its derivatives up to the fourth one are monotonic in `s`. This is
true for most of the interval but there are some points where it
doesn't hold. One solution would be to prove that this only happens at
some places, isolate them and handle them separately. This might be
tedious though since we would have to do it for all required
derivatives. The point where it happens does depend on `x`.
"""
function T022(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * log(10 + inv(x)))
    end

    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)

        # FIXME: Currently we assume monotonicity in s, including for
        # all derivatives.
        s_l = 1 - u0.ϵ
        s_u = one(Arb)
        integrand(y) = begin
            term_l = clausenc(y - x, s_l) + clausenc(y + x, s_l) - 2clausenc(y, s_l)
            term_u = clausenc(y - x, s_u) + clausenc(y + x, s_u) - 2clausenc(y, s_u)

            if y isa ArbSeries
                coefficients = union.(Arblib.coeffs(term_l), Arblib.coeffs(term_u))
                term_union = ArbSeries(coefficients)
            else
                term_union = union(term_l, term_u)
            end

            return term_union * y * log(10 + inv(y))
        end

        res = ArbExtras.integrate(integrand, a, Arb(π), atol = 1e-5, rtol = 1e-5)

        res = res / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
