"""
    T02(u0::BHAnsatz; δ2, skip_div_u0)
Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral T_{0,2} from the paper.

The interval of integration is `[x, π]`. Since the integrand is
singular at `y = x` we split the interval into two parts, `[x, a]` and
`[a, π]`. In principle we want to take `a = x = δ2`, however this
gives issues if `x` is a wide ball. Instead we take `a` to always be a
thin ball between `x` and `π`.

In general we take `a = ubound(x + δ2)`. However if `x` is wide then
it's beneficial to take a larger value than `δ2`, depending on the
radius of `x`. In practice we take `8radius(x)` in that case.

If `x` is very close to `π`, so that `a` would be larger than `π`,
then we use the asymptotic version for the whole interval `[x, π]`.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T02(u0::BHAnsatz, evaltype::Ball; δ2::Arb = Arb(1e-10), skip_div_u0 = false)
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
    T021(u0::BHAnsatz)
Computes the (not yet existing) integral T_{0,2,1} from the paper.

The interval of integration is `[x, a]`. Both `x` and `a` are assumed
to be less than or equal to `π`, if they are balls which overlap `π`
anything above `π` will be ignored.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval. We can also notice that the value inside the absolute value
is negative so we can remove the absolute value by putting a minus
sign, which we can bake in to the weight factor.

We are left with integrating the log-term. This allows us to split the
integrand into three terms
1. `log(-sin((x - y) / 2)) = log(sin((y - x) / 2))`
2. `log(sin((x + y) / 2))`
3. `-2log(sin(y / 2))`

The third term is well behaved no matter the value of `x` and we can
enclose that integral by enclosing the integrand on the interval and
multiplying with the intervals size.

For the first two terms we have to differentiate between the case when
`x` overlaps with `π` and when it doesn't.

For the case when `x` doesn't overlap with `π` we have that the second
term is well behaved and we can bound its integral in the same way as
we handle the third term. For the first term we use the inequality
```
c * (y - x) / 2 <= sin((y - x) / 2) <= (y - x) / 2
```
which holds for `c = sin(δ / 2) / (δ / 2)` on `0 <= x <= π` and `x <=
y <= a`. This gives us
```
log(c * (y - x) / 2) <= log(sin((y - x) / 2)) <= log((y - x) / 2)
```
The same inequality holds after integration from `x` to `a` and
gives us
```
δ * (log(c * δ / 2) - 1) <= ∫log(sin((y - x) / 2)) <= δ * (log(δ / 2) - 1)
```
where `δ = a - x`.

If `x` overlaps with `π` then we assume that `a = π` since this is the
only valid case. We have that both the first and the second term are
singular. We can get a direct upper bound for both of them by using
that `log(sin(t)) <= 0` for all `t`, hence they are both upper bounded
by `0`.

For lower bounding the first term we use that `(y - x) / 2 < π / 2`
and hence `sin((y - x) / 2) >= 2 / π * (y - x) / 2. We can thus lower
bound the integral by integrating `log(2 / π * (y - x) / 2)` from `x`
to `π`. This integral is given by `(π - x) * (-1 + log(2 / π * (π - x)
/ 2))`, which is strictly increasing in `x` and a lower bound is hence
given by evaluating it at `Arblib.lbound(x)`.

For the second term we use that `sin((x + y) / 2) = sin(π - (x + y) /
2)` and as long as `π - (x + y) / 2 < π / 2` this is lower bounded by
`2 / π * (π - (x + y) / 2)`. The inequality `π - (x + y) / 2 < π / 2`
holds on the interval as long as `x <= π / 2`. The integral of `2 / π
* (π - (x + y) / 2)` from `x` to `π` is given by `(π - x) * (-1 +
log(4 / π * (π - x)))`, which, similar to above, is strictly
increasing in `x` and we can hence evaluate it at `Arblib.lbound(x)`
to get a lower bound.
"""
function T021(u0::BHAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-10), skip_div_u0 = false)
    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)
        δ = a - x

        interval = union(x, a)

        weight_factor = -u0.w(interval)

        if !(x < π)
            @assert Arblib.overlaps(a, Arb(π)) # Take π to be the upper integration limit

            x >= Arb(π) / 2 ||
                throw(ArgumentError("we require that x >= π / 2, got x = $x"))

            x_lower = Arblib.lbound(Arb, x)

            x_lower < π || throw(
                ArgumentError("we require that the lower bound of x is less than π"),
            )

            part1 = begin
                part1_lower = (π - x_lower) * (-1 + log(2 / Arb(π) * (π - x_lower) / 2))
                part1_upper = zero(x)

                Arb((part1_lower, part1_upper))
            end

            part2 = begin
                part2_lower = (π - x_lower) * (-1 + log(4 / Arb(π) * (π - x_lower)))
                part2_upper = zero(x)

                Arb((part2_lower, part2_upper))
            end
        else
            part1 = let c = sin(δ / 2) / (δ / 2)
                part1_lower = δ * (log(c * δ / 2) - 1)
                part1_upper = δ * (log(δ / 2) - 1)

                Arb((part1_lower, part1_upper))
            end

            part2 = log(sin((x + interval) / 2)) * δ
        end

        part3 = -2log(sin(interval / 2)) * δ

        integral = weight_factor * (part1 + part2 + part3)

        res = integral / (π * u0.w(x))
        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T022(u0::BHAnsatz)
Computes the (not yet existing) integral T_{0,2,2} from the paper.

The interval of integration is given by `[a, π]`. In practice `a`
should be a thin ball to not give problems with the integration.

This is done by directly computing the integral with the integrator in
Arb. Accounting for the fact that the integrand is non-analytic at `y
= x`.

TODO: This doesn't work when `x` is close to `π`.
"""
function T022(u0::BHAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-10), skip_div_u0 = false)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)

        # PROVE: That there are no branch cuts that interact with the
        # integral
        integrand(y; analytic::Bool) = begin
            res = log(-sin((x - y) / 2) * sin((x + y) / 2) / sin(y / 2)^2)
            Arblib.real_abs!(res, res, analytic)
            weight = y * Arblib.sqrt_analytic!(zero(y), log((y + 1) / y), analytic)
            return res * weight
        end

        res = Arblib.integrate(
            integrand,
            a,
            π,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
            warn_on_no_convergence = false,
        )
        @assert !isfinite(res) || isreal(res)
        res = real(res)

        res = res / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
