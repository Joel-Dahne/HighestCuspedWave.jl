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
function T02(u0::BHAnsatz, evaltype::Ball; δ2::Arb = Arb(1e-5), skip_div_u0 = false)
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
    T02(u0::BHAnsatz, ::Asymptotic; ϵ = Arb(2e-1))
Returns a function such that `T02(u0, Asymptotic())(x)` computes an
**upper bound** of the integral T_{0,2} from the paper using an
evaluation strategy that works asymptotically as `x` goes to 0.

It precomputes the expansions of `u0` and for that reason a number `ϵ`
has to be given, the resulting expansion will be valid for all `x <
ϵ`. The value of `ϵ` has to be less than `1`.

TODO: This is currently only partially implemented for full asymptotic
evaluation. It computes the integral after the change of variables `t
= y / x` using normal ball arithmetic. The factor outside the integral
is computed using an asymptotic approach.
"""
function T02(
    u0::BHAnsatz,
    ::Asymptotic;
    non_asymptotic_u0 = false,
    ϵ = Arb(2e-1),
)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    ϵ = convert(Arb, ϵ)
    u0_expansion = u0(ϵ, AsymptoticExpansion())

    u0_expansion_div_xlogx = empty(u0_expansion)
    for ((i, m, k, l), value) in u0_expansion
        u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
    end

    factor(x) = begin
        # PROVE: That this is monotonically decreasing on [0, 1]
        w(x) = 1 / (log(x) * sqrt(log((x + 1) / x)))

        @assert x <= ϵ

        if iszero(x)
            weight = zero(x)
        elseif Arblib.contains_zero(x) && x < 1
            weight = union(zero(x), w(Arblib.ubound(Arb, x)))
        else
            weight = w(x)
        end

        return weight / (π * eval_expansion(u0, u0_expansion_div_xlogx, x))
    end

    return x -> begin
        x = convert(Arb, x)

        # Variables for storing temporary values during integration
        x_complex = convert(Acb, x)
        xdiv2 = x_complex / 2
        tmp = zero(x_complex)
        tx = zero(x_complex)

        # TEMPORARY: This currently computes a non-asymptotic version
        # for testing purposes
        integrand!(res, t; analytic::Bool) = begin
            # The code below is an inplace version of the following code
            #res = -log(sin(x * (t - 1) / 2) * sin(x * (1 + t) / 2) / sin(t * x / 2)^2)
            #weight = t * Arblib.sqrt_analytic!(zero(t), log((t * x + 1) / (t * x)), analytic)
            #return res * weight

            Arblib.mul!(tx, t, x_complex)

            # res = sin((t - 1) * x / 2)
            Arblib.sub!(tmp, t, 1)
            Arblib.mul!(tmp, tmp, xdiv2)
            Arblib.sin!(res, tmp)

            # res *= sin((1 + t) * x / 2)
            Arblib.add!(tmp, t, 1)
            Arblib.mul!(tmp, tmp, xdiv2)
            Arblib.sin!(tmp, tmp)
            Arblib.mul!(res, res, tmp)

            # res /= sin(t * x / 2)^2
            Arblib.mul_2exp!(tmp, tx, -1)
            Arblib.sin!(tmp, tmp)
            Arblib.sqr!(tmp, tmp)
            Arblib.div!(res, res, tmp)

            Arblib.log!(res, res)
            Arblib.neg!(res, res)

            # tmp = t * sqrt(log((tx + 1) / tx))
            Arblib.add!(tmp, tx, 1)
            Arblib.div!(tmp, tmp, tx)
            Arblib.log!(tmp, tmp)
            Arblib.sqrt_analytic!(tmp, tmp, analytic)
            Arblib.mul!(tmp, tmp, t)

            Arblib.mul!(res, res, tmp)

            return
        end

        δ = Arb(1e-10)

        res = Arblib.integrate!(
            integrand!,
            zero(x_complex),
            one(x_complex) + δ,
            π / x_complex,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
            warn_on_no_convergence = false,
        )

        @assert !isfinite(res) || isreal(res)
        res = real(res)

        return res * factor(x)
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
y <= a`. In particular it also holds for `c = sin(ubound(δ) / 2) /
(ubound(δ) / 2)` due to monotonicity. This gives us
```
log(c * (y - x) / 2) <= log(sin((y - x) / 2)) <= log((y - x) / 2)
```
The same inequality holds after integration from `x` to `a` and
gives us
```
δ * (log(c * δ / 2) - 1) <= ∫log(sin((y - x) / 2)) <= δ * (log(δ / 2) - 1)
```
where `δ = a - x`.

If `x` overlaps with `π` or is very close to `π` (less than
`sqrt(eps(x)))` then we take `a = π` as the upper integration limit
since even if a smaller `a` is given this will only give a small
overestimation. We have that both the first and the second term are
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
function T021(u0::BHAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)
        δ = a - x

        interval = union(x, a)

        weight_factor = -u0.w(interval)

        if !(x < π - sqrt(eps(x)))
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
            part1 = let c = sin(Arblib.ubound(Arb, δ) / 2) / (Arblib.ubound(δ) / 2)
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

Notice that the expression inside the absolute value is always
negative, so we can replace the absolute value with a negation.
"""
function T022(u0::BHAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)

        # Variables for storing temporary values during integration
        x_complex = convert(Acb, x)
        tmp = zero(x_complex)

        # PROVE: That there are no branch cuts that interact with the
        # integral
        integrand!(res, y; analytic::Bool) = begin
            # The code below is an inplace version of the following code
            #res = -log(sin((y - x) / 2) * sin((x + y) / 2) / sin(y / 2)^2)
            #weight = y * Arblib.sqrt_analytic!(zero(y), log((y + 1) / y), analytic)
            #return res * weight

            # res = sin((y - x) / 2)
            Arblib.sub!(tmp, y, x_complex)
            Arblib.mul_2exp!(tmp, tmp, -1)
            Arblib.sin!(res, tmp)

            # res *= sin((x + y) / 2)
            Arblib.add!(tmp, x_complex, y)
            Arblib.mul_2exp!(tmp, tmp, -1)
            Arblib.sin!(tmp, tmp)
            Arblib.mul!(res, res, tmp)

            # res /= sin(y / 2)^2
            Arblib.mul_2exp!(tmp, y, -1)
            Arblib.sin!(tmp, tmp)
            Arblib.sqr!(tmp, tmp)
            Arblib.div!(res, res, tmp)

            Arblib.log!(res, res)

            Arblib.neg!(res, res)

            # tmp = y * sqrt(log((y + 1) / y))
            Arblib.add!(tmp, y, 1)
            Arblib.div!(tmp, tmp, y)
            Arblib.log!(tmp, tmp)
            Arblib.sqrt_analytic!(tmp, tmp, analytic)
            Arblib.mul!(tmp, tmp, y)

            Arblib.mul!(res, res, tmp)

            return
        end

        res = Arblib.integrate!(
            integrand!,
            zero(x_complex),
            a,
            π,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
            warn_on_no_convergence = false,
            #opts = Arblib.calc_integrate_opt_struct(0, 0, 0, 0, 1),
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
