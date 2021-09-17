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
ϵ`. The value of `ϵ` has to be less than `1 // 2`.

To begin with the factor `x * log(x) / (π * u0(x))` is factored out
from the whole expression and multiplied back in the end. Notice that
this factor is bounded in `x`.

What we are left with computing is
```
W(x) * I
```
where `W(x) = 1 / (x^2 * log(x) * sqrt(log((x + 1) / x)))` and `I`
is given by the integral
```
∫-(log(sin((y - x) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2))) * y * sqrt(log((y + 1) / y)) dy
```
for `x` to `π`.

We first get rid of the `sin` inside the `log` terms by adding and subtracting
```
(log((y - x) / 2) + log((y + x) / 2) - 2log(y / 2)) * y * sqrt(log((y + 1) / y))
```
This gives us two integrals
```
J = -∫(log((y - x) / 2) + log((y + x) / 2) - 2log(y / 2)) * y * sqrt(log((y + 1) / y)) dy
```
and
```
E = -∫((log(sin((y - x) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2))) - (log((y - x) / 2) + log((y + x) / 2) - 2log(y / 2))) * y * sqrt(log((y + 1) / y)) dy
```
with `I = J + E`. We can simplify `E` to
```
E = -∫(log(2sin((y - x) / 2) / (y - x)) + log(2sin((y + x) / 2) / (y + x)) - 2log(2sin(y / 2) / y)) * y * sqrt(log((y + 1) / y)) dy
```

For `J` we can simplify it to
```
J = -∫(log(y - x) + log(y + x) - 2log(y)) * y * sqrt(log((y + 1) / y)) dy
```
by taking out the divisions by `2`. We can simplify it further by
rewriting `log(y - x) = log(y) + log(1 - x / y)` and `log(y + x) =
log(y) + log(1 + x / y)`, giving us
```
J = -∫(log(1 - x / y) + log(1 + x / y)) * y * sqrt(log((y + 1) / y)) dy
```
Finally we can combine the two logs
```
J = -∫log(1 - (x / y)^2) * y * sqrt(log((y + 1) / y)) dy
```
Recall that this integral is taken from `x` to `π`. Next we do a
change of variables, `t = y / x`, giving us the integral
```
J = -x^2 * ∫log(1 - 1 / t^2) * t * sqrt(log((t * x + 1) / (t * x))) dt
```
from `1` to `π / x`. Assuming that `x < 1 // 2` we split this into
three parts, one from `1` to `2`, one from `2` to `1 / x` and one from
`1 / x` to `π / x`, which we call `J1`, `J2` and `J3` respectively.

For `J1` we first notice that
```
sqrt(log((t * x + 1) / (t * x))) = sqrt(log((t * x + 1) / (t * x + t)) + log((x + 1) / x))
```
Here we have that `log((t * x + 1) / (t * x + t))` is bounded for `0 <
x < 1 / 2` and `1 < t < 2` and an enclosure can easily be computed
with interval arithmetic, call this enclosure `C`. This allows us to
simplify `J1` as
```
J1 = -x^2 * sqrt(C + log((x + 1) / x)) * ∫log(1 - 1 / t^2) * t dt
```
where the integral can now be explicitly computed to be `log(3 *
sqrt(3) / 16)`. If we also consider the `W(x)` factor we find that
```
W(x) * J1 = -x^2 * sqrt(C + log((x + 1) / x)) * log(3 * sqrt(3) / 16) / (x^2 * log(x) * sqrt(log((x + 1) / x)))
```
which simplifies to
```
W(x) * J1 = -log(3 * sqrt(3) / 16) * sqrt(C + log((x + 1) / x)) / (log(x) * sqrt(log((x + 1) / x)))
```
PROVE: that `sqrt(C + log((x + 1) / x)) / (log(x) * sqrt(log((x +
1) / x)))` is zero at `x = 0` and monotonically decreasing.

For `J2` we use the Taylor expansion
```
log(1 - 1 / t^2) = - 1 / t^2 + D / t^4,
```
where `D` is the rest term, together with
```
sqrt(log((t * x + 1) / (t * x))) = sqrt(log(1 / (t * x))) + (sqrt(log((t * x + 1) / (t * x))) - sqrt(log(1 / (t * x))))
```
to split the integral into four terms
```
J21 = x^2 * ∫ sqrt(log(1 / (t * x))) / t dt
```
```
J22 = x^2 * ∫ (sqrt(log((t * x + 1) / (t * x))) - sqrt(log(1 / (t * x)))) / t dt
```
```
J23 = -D * x^2 * ∫ sqrt(log(1 / (t * x))) / t^3 dt
```
```
J24 = -D * x^2 ∫ (sqrt(log((t * x + 1) / (t * x))) - sqrt(log(1 / (t * x)))) / t^3 dt
```
We can bound `D` directly using `ArbSeries`. The integrals `J21` and
`J23` can be computed explicitly, giving us
```
J21 = 2 // 3 * x^2 * log(1 / 2x)^(3 // 2)
```
and
```
J23 = -1 // 8 * D * x^2 * (sqrt(-log(2x)) - sqrt(2π) * x^2 * erfi(sqrt(-log(4x^2))))
```
For `J22` and `J24` we will use the fact that `sqrt(log((t * x + 1) /
(t * x))) - sqrt(log(1 / (t * x)))` is bounded on `[2, 1 / x]` to
factor it out as a constant. If we let `G` be a ball containing the
range of it we get
```
J22 = G * x^2 * ∫1 / t dt = -G * x^2 * log(2x)
```
```
J24 = -D * G * x^2 * ∫1 / t^3 dt = 1 // 8 * D * G * x^2 * (1 - 4x^2)
```
TODO: Figure out how to handle the `erfi` function in a good way.
TODO: Bound G on the interval `[2, 1 / x]`.


For `J3` we notice that
```
sqrt(log((π + 1) / π)) < sqrt(log((t * x + 1) / (t * x))) < sqrt(log(2))
```
for `1 / x < t < π / x`. If we let `D` be a ball containing this
interval we have
```
J3 = -x^2 * D * ∫log(1 - 1 / t^2) * t dt
```
where the integral can be computed to be
```
∫log(1 - 1 / t^2) * t dt = 1 / (2x^2) * ((x^2 - 1) * log(1 / x^2 - 1) + (π^2 - x^2) * log(π^2 / x^2 - 1) - 2(π^2 * log(π) + log(x) - π^2 * log(x)))
```
which we will call `P(x)`. Including the `W(x)` factor we get
```
W(x) * J3 = -x^2 * D * P(x) / (x^2 * log(x) * sqrt(log((x + 1) / x)))
```
which simplifies to
```
W(x) * J3 = -D * P(x) / (log(x) * sqrt(log((x + 1) / x)))
```
PROVE: monotonicity for of P(x).

For `E` we start by taking out a factor `x^2` from the integrand,
giving us
```
E = -x^2 * ∫(log(2sin((y - x) / 2) / (y - x)) + log(2sin((y + x) / 2) / (y + x)) - 2log(2sin(y / 2) / y)) * y * sqrt(log((y + 1) / y)) / x^2 dy
```
which is to be integrated from `x` to `π`. We notice that this
integrand is uniformly bounded on the interval and we can thus enclose
the integral for all values of `x`.
TODO: Compute this bound.

"""
function T02(u0::BHAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    @warn "T02(u0, Asymptotic()) doesn't bound the error term completely yet"

    ϵ = convert(Arb, ϵ)
    @assert ϵ < 0.5

    # Compute expansion for u0 and set up W(x) method
    u0_expansion = u0(ϵ, AsymptoticExpansion())
    u0_expansion_div_xlogx = empty(u0_expansion)
    for ((i, m, k, l), value) in u0_expansion
        u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
    end

    # This gives the factor x * log(x) / (π * u0(x))
    factor(x) = inv(π * eval_expansion(u0, u0_expansion_div_xlogx, x))

    W(x) = 1 / (x^2 * log(x) * sqrt(log((x + 1) / x)))

    # Setting this to true makes the code do a few sanity checks on
    # the results, they are proved to hold but it's good to double
    # check.
    check_results = true
    return x -> begin
        x = convert(Arb, x)

        # Includes the W(x) factor
        J1 = let t = Arb((1, 2))
            C = log((t * x + 1) / (t * x + t))

            if iszero(x)
                # PROVE: That this is zero at x = 0
                fact = zero(x)
            elseif Arblib.contains_zero(x)
                # PROVE: That this is monotone
                xᵤ = ubound(Arb, x)
                fact = union(
                    zero(x),
                    sqrt(C + log((xᵤ + 1) / xᵤ)) / (log(xᵤ) * sqrt(log((xᵤ + 1) / xᵤ))),
                )
            else
                fact = sqrt(C + log((x + 1) / x)) / (log(x) * sqrt(log((x + 1) / x)))
            end

            -log(3 * sqrt(oftype(x, 3)) / 16) * fact
        end

        if check_results
            J1_check =
                -x^2 *
                W(x) *
                real(
                    Arblib.integrate(
                        (t; analytic) ->
                            log(1 - 1 / t^2) *
                            t *
                            Arblib.sqrt_analytic!(
                                zero(t),
                                log((t * x + 1) / (t * x)),
                                analytic,
                            ),
                        1 + 1e-10,
                        2,
                        check_analytic = true,
                    ),
                )

            @assert Arblib.overlaps(J1, J1_check)
        end

        # Includes the W(x) factor
        J2 = let
            # Bound the constant D using ArbSeries. We compute the
            # error term for log(1 - a) for a in [0, 1 // 4]
            D = log(1 - ArbSeries([Arb((0, 1 // 4)), 1, 0]))[2]

            if check_results
                for t in range(Arb(2), 1 / x, length = 1000)
                    @assert Arblib.overlaps(log(1 - 1 / t^2), -1 / t^2 + D / t^4)
                end
            end

            # FIXME: Using this as a temporary enclosure for the error
            # term, this should be updated with a proved enclosure.
            G = Arb((0, 1))

            J21 = 2 // 3 * log(1 / 2x)^(3 // 2)

            J22 = -G * log(2x)

            J23 =
                1 // 8 *
                D *
                (
                    sqrt(-log(2x)) -
                    sqrt(2oftype(x, π)) * x^2 * SpecialFunctions.erfi(sqrt(-log(4x^2)))
                )

            J24 = -1 // 8 * D * G * (1 - 4x^2)

            (J21 + J22 + J23 + J24) / (log(x) * sqrt(log((x + 1) / x)))
        end

        if check_results
            J2_check =
                -x^2 *
                W(x) *
                real(
                    Arblib.integrate(
                        (t; analytic) ->
                            log(1 - 1 / t^2) *
                            t *
                            Arblib.sqrt_analytic!(
                                zero(t),
                                log((t * x + 1) / (t * x)),
                                analytic,
                            ),
                        2,
                        1 / x,
                        check_analytic = true,
                    ),
                )

            @assert Arblib.overlaps(J2, J2_check)
        end

        # Includes the W(x) factor
        J3 = let
            D = Arb((sqrt(log((Arb(π) + 1) / Arb(π))), sqrt(log(Arb(2)))))

            # FIXME: Compute a proper enclosure of P(x)
            P = let π = Arb(π)
                1 / (2x^2) * (
                    (x^2 - 1) * log(1 / x^2 - 1) + (π^2 - x^2) * log(π^2 / x^2 - 1) -
                    2(π^2 * log(π) + log(x) - π^2 * log(x))
                )
            end

            -D * P / (log(x) * sqrt(log((x + 1) / x)))
        end

        if check_results
            J3_check =
                -x^2 *
                W(x) *
                real(
                    Arblib.integrate(
                        (t; analytic) ->
                            log(1 - 1 / t^2) *
                            t *
                            Arblib.sqrt_analytic!(
                                zero(t),
                                log((t * x + 1) / (t * x)),
                                analytic,
                            ),
                        1 / x,
                        π / x,
                        check_analytic = true,
                    ),
                )

            @assert Arblib.overlaps(J3, J3_check)
        end

        # Includes the W(x) factor
        E = let
            # Enclosure of the integrand
            # FIXME: Using this as a temporary enclosure for the error
            # term, this should be updated with a proved enclosure.
            D = Arb((-1, 0))

            # This gives us an enclosure of the integral
            -π * D / (log(x) * sqrt(log((x + 1) / x)))
        end

        return factor(x) * (J1 + J2 + J3 + E)

        return J1, J2, J3, E
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
`sqrt(eps(x))`) then we take `a = π` as the upper integration limit
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
