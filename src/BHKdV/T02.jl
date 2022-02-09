"""
    T02(u0::BHKdVAnsatz; δ2, skip_div_u0)

Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral ``T_{0,2}`` from the paper.

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

Computes the (not yet existing) integral ``T_{0,2,1}`` from the paper.

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
- **TODO:** Could improve enclosures by better handling cancellations
  for `clausens(x + a, 1 - α) - clausens(2x, 1 - α)` and `clausens(a,
  1 - α) - clausens(x, 1 - α)`. Though this might not be needed.
"""
function T021(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    # s = 1 - α computed such that the upper bound is exactly 2
    s = 2 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    return (x::Arb, a::Arb = x + δ2) -> begin
        integral =
            -clausens(x - a, s) + (clausens(x + a, s) - clausens(2x, s)) -
            2(clausens(a, s) - clausens(x, s))

        # Multiply by weight inside the integral that was factored out
        integral *= u0.w(union(x, a))

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

Computes the (not yet existing) integral ``T_{0,2,2}`` from the paper.

It returns a function `f` such that `f(x, a; tol)` computes the
integral on `[a, π]` for the given value of `x` and it uses the
prescribed tolerance for the integration. In practice `a` should be a
thin ball to not give problems with the integration.

This is done by directly computing the integral with the integrator in
Arb.

Notice that due to lemma [`lemma_integrand_2`](@ref) the expression
inside the absolute value is always positive, so we can remove the
absolute value.
"""
function T022(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    # Lower and upper bounds of s = -α
    s_l = 1 - u0.ϵ
    s_u = one(Arb)

    # Enclosure of -α so that the upper bound is exactly 1
    mα = 1 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    # Upper integration limit
    b = Arb(π)

    return (x::Arb, a::Arb = x + δ2; tol = Arb(1e-5)) -> begin
        integrand(y) =
            (clausenc(y - x, mα) + clausenc(y + x, mα) - 2clausenc(y, mα)) * u0.w(y)

        res = ArbExtras.integrate(integrand, a, b, atol = tol, rtol = tol)

        res /= (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
