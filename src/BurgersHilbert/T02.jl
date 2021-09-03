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

To begin with the factor `x * log(x) / (π * u0(x))` is factored out
from the whole expression and added back in the end. Notice that this
factor is bounded.

What we are left with computing is
```
W(x) * I
```
where `W(x) = 1 / (x^2 * log(x) * sqrt(log((x + 1) / x)))` and `I`
given by the integral
```
∫(log(sin((y - x) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2))) * y * sqrt(log((y + 1) / y)) dy
```
for `x` to `π`.

Expanding the `log(sin)` terms at `x = 0` we find that
```
log(sin((y - x) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2)) = x^2 / 4 * (1 + cot(y / 2)^2) + O(x^4)
```
and we can therefore split the problem into computing the two
integrals
```
J = ∫(1 + cot(y / 2)^2) * y * sqrt(log((y + 1) / y)) dy
```
and
```
E = ∫(-(log(sin((y - x) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2))) - x^2 / 4 * (1 + cot(y / 2)^2)) * y * sqrt(log((y + 1) / y)) dy
```
both from `x` to `π`.

We will split `J` into several parts, one parts which grows in `x` and
other parts which are all bounded. To begin with we split the integral
domain into two intervals, one from `x` to `ϵ` which we will call `J1`
and one from `ϵ` to `π` which we will call `J2`. The integral `J2` is
bounded and an enclosure can be computed directly using
`Arblib.integrate`.

For `J1` we expand the integrand at `y = 0`, giving us
```
(1 + cot(y / 2)^2) * y * sqrt(log((y + 1) / y)) = 4sqrt(-log(y)) / y - 2sqrt(-log(y)) / log(y) + O(y)
```
We therefore split `J1` into three parts
```
J11 = ∫4sqrt(-log(y)) / y dy
```
```
J12 = ∫-2sqrt(-log(y)) / log(y) dy
```
```
J13 = ∫(1 + cot(y / 2)^2) * y * sqrt(log((y + 1) / y)) - 4sqrt(-log(y)) / y - (-2sqrt(-log(y)) / log(y)) dy
```
all with the limits `x` to `ϵ`. Notice that both `J12` and `J13` are
both bounded in `x`. The first two we can compute explicitly, giving
us
```
J11 = 8 // 3 * (sqrt(-log(ϵ)) * log(ϵ) - (-log(x))^(3 // 2))
```
and
```
J12 = 2sqrt(π) * (erf(sqrt(-log(x))) - erf(sqrt(-log(ϵ))))
```
Finally we split `J11` into the constant and the non-constant part,
`J111 = 8 // 3 * sqrt(-log(ϵ)) * log(ϵ)` and `J112 = -8 // 3 *
(-log(x))^(3 // 2)`.

The full integral is now given by
```
I = x^2 / 4 * (J111 + J112 + J12 + J13 + J2) + E
```
And we are interested in `W(x) * I`. To be able to properly cancel all
occurrences of `x` we split it into three parts
```
part1 = -2sqrt(-log(x)) / 3sqrt(log((x + 1) / x))
```
```
part2 = (J111 + J12 + J13 + J2) / 4(log(x) * sqrt(log((x + 1) / x)))
```
and
```
part3 = W(x) * E
```

For computing `J13` we need to be able to enclose the integrand for
`y` overlapping `0`. This is done by expanding the integral and
explicitly bounding all the constants, it requires quite a bit of
work. Recall that the integrand in this case is given by
```
(1 + cot(y / 2)^2) * y * sqrt(log((y + 1) / y)) - 4sqrt(-log(y)) / y - (-2sqrt(-log(y)) / log(y))
```
As a first step we take out the term `y * sqrt(log((y + 1) / y))`,
this can be checked to be bounded and monotone and is hence easy to
enclose. Next we rewrite `sqrt(log((y + 1) / y))` as
```
sqrt(log(y + 1) - log(y)) = sqrt(-log(y)) * sqrt(1 - log(y + 1) / log(y))
```
We expect the integrand to behave like `y * sqrt(-log(y))` so we
factor that out, giving us
```
y * sqrt(-log(y)) * (cot(y / 2)^2 * sqrt(1 - log(y + 1) / log(y)) - 4 / y^2 + 2 / (y * log(y)))
```
To handle the large parenthesis we put everything in one fraction,
giving us (skipping the `y * sqrt(-log(y))` in front and rewriting
`cot(y / 2) = 1 / tan(y / 2)`)
```
(y^3 * log(y) * sqrt(1 - log(y + 1) / log(y)) - 4y * log(y) * tan(y / 2)^2 + 2y^2 * tan(y / 2)^2) /
(y^3 * log(y) * tan(y / 2)^2)
```
Now consider the three expansions
```
tan(y / 2)^2 = y^2 / 4 + C1 * y^4
```
```
sqrt(1 - log(y + 1) / log(y)) = 1 - 1 // 2 * log(y + 1) / log(y) + C2 * y^2
```
and
```
log(1 + y) = y + C3 * y^2
```
where `C1`, `C2` and `C3` can be **explicitly** enclosed by balls.
Inserting this in the expression above gives us for the denominator
```
(y^3 * log(y) * tan(y / 2)^2) = log(y) * (y^5 / 4 + C1 * y^7)
```
and for the numerator we get
```
(y^3 * log(y) - 1 // 2 * y^3 * log(y + 1) + C2 * y^5 * log(y)) -
(y^3 * log(y) + 4C1 * y^5 * log(y)) +
(1 // 2 * y^4 + 2C1 * y^6)
```
We see that the two occurrences of `y^3 * log(y)` cancel out, also
expanding the `log(1 + y)` term gives
```
(-1 // 2 * y^4 - 1 // 2 * C3 * y^5 + C2 * y^5 * log(y)) -
4C1 * y^5 * log(y) + 1 // 2 * y^4 + 2C1 * y^6
```
and the `-1 // 2 * y^4` terms cancel, leaving us with
```
-1 // 2 * C3 * y^5 + C2 * y^5 * log(y) - 4C1 * y^5 * log(y) + 2C1 * y^6 =
(C2 - 4C1) * y^5 * log(y) - 1 // 2 * C3 * y^5 + 2C1 * y^6
```
Now putting the numerator and denominator together and canceling the
`y^5 * log(y)` factor we get
```
((C2 - 4C1) - 1 // 2 * C3 / log(y) + 2C1 * y^2 / log(y)) / (1 // 4 + C1 * y^5)
```
This can be explicitly bounded.

For enclosing `C1` and `C3` we can use `ArbSeries directly.

TODO: How do we enclose `C2`? It seems like the second derivative is
not even bounded?
TODO: Handle asymptotic evaluation of `part3`, notice that it should
behave like `O(x^4)`.
"""
function T02(u0::BHAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    ϵ = convert(Arb, ϵ)
    @assert ϵ < 1

    # Compute expansion for u0 and set up W(x) method
    u0_expansion = u0(ϵ, AsymptoticExpansion())
    u0_expansion_div_xlogx = empty(u0_expansion)
    for ((i, m, k, l), value) in u0_expansion
        u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
    end

    # Enclose J2
    integrand_J2!(res, y; analytic::Bool) = begin
        # The code below is an inplace version of the following code
        #res = (1 + cot(y / 2)^2)
        #weight = y * Arblib.sqrt_analytic!(zero(t), log((y + 1) / y), analytic)
        #return res * weight

        tmp = zero(res)

        # res = 1 + cot(y / 2)^2
        Arblib.mul_2exp!(res, y, -1)
        Arblib.cot!(res, res)
        Arblib.sqr!(res, res)
        Arblib.add!(res, res, 1)

        # res *= y
        Arblib.mul!(res, res, y)

        # tmp = sqrt(log((y + 1) / y))
        Arblib.add!(tmp, y, 1)
        Arblib.div!(tmp, tmp, y)
        Arblib.log!(tmp, tmp)
        Arblib.sqrt_analytic!(tmp, tmp, analytic)

        # res *= tmp
        Arblib.mul!(res, res, tmp)

        return
    end

    J2 = Arblib.integrate!(
        integrand_J2!,
        Acb(prec = precision(ϵ)),
        ϵ,
        π,
        check_analytic = true,
        rtol = 1e-10,
        atol = 1e-10,
        warn_on_no_convergence = false,
    )

    @assert !isfinite(J2) || isreal(J2)
    J2 = real(J2)

    integrand_J13!(res, y; analytic::Bool) = begin
        # The code below is an inplace version of the following code
        #res = (1 + cot(y / 2)^2)
        #weight = y * Arblib.sqrt_analytic!(zero(t), log((y + 1) / y), analytic)
        #return res * weight - 4sqrt(-log(y)) / y - (-2sqrt(-log(y)) / log(y))

        if Arblib.contains_zero(y)
            # The integrand is not analytic at y = 0 and our bounds
            # only work for real y
            if analytic || !isreal(y)
                Arblib.indeterminate!(res)
                return
            end

            # Get the upper bound of the real part of y
            yᵤ = ubound(Arb, Arblib.realref(y))

            # Bound the term y * sqrt(log((y + 1) / y)) which is
            # monotonically increasing and 0 at y = 0
            term1 = union(zero(yᵤ), yᵤ * sqrt(log((yᵤ + 1) / yᵤ)))

            # Bound the more complicated term cot(y / 2)^2 * y *
            # sqrt(log((y + 1) / y)) - 4sqrt(-log(y)) / y -
            # (-2sqrt(-log(y)) / log(y))

            # Compute C1, C2 and C3 as described in the documentation above

            # C1 and C3 can be computed directly with ArbSeries
            C1, C2, C3 = let y = ArbSeries([real(y), 1], degree = 4)
                C1 = (tan(y / 2)^2)[4]
                # FIXME: This value works in practice on [0, 1/2] but
                # is not proved
                C2 = union(-one(yᵤ), zero(yᵤ))
                C3 = log(1 + y)[2]
                C1, C2, C3
            end

            # Evaluate ((C2 - 4C1) - 1 // 2 * C3 / log(y) + 2C1 * y^2 / log(y)) / (1 // 4 + C1 * y^5)
            # We can use that 1 / log(y), y^2 / log(y) and y^5 are all
            # zero at y = 0 and monotone on [0, 1]
            fraction =
                (
                    (C2 - 4C1) - 1 // 2 * C3 * union(zero(yᵤ), 1 / log(yᵤ)) +
                    2C1 * union(zero(yᵤ), yᵤ^2 / log(yᵤ))
                ) / (1 // 4 + C1 * union(zero(yᵤ), yᵤ^5))

            # Now term2 = y * sqrt(-log(y)) * constant and is monotone on [0, ϵ]
            term2 = union(zero(yᵤ), yᵤ * sqrt(-log(yᵤ)) * fraction)

            Arblib.add!(res, res, term1 + term2)

            return
        end

        tmp = zero(res)

        # res = 1 + cot(y / 2)^2
        Arblib.mul_2exp!(res, y, -1)
        Arblib.cot!(res, res)
        Arblib.sqr!(res, res)
        Arblib.add!(res, res, 1)

        # res *= y
        Arblib.mul!(res, res, y)

        # tmp = sqrt(log((y + 1) / y))
        Arblib.add!(tmp, y, 1)
        Arblib.div!(tmp, tmp, y)
        Arblib.log!(tmp, tmp)
        Arblib.sqrt_analytic!(tmp, tmp, analytic)

        # res *= tmp
        Arblib.mul!(res, res, tmp)

        # res -= 4sqrt(-log(y)) / y
        Arblib.sub!(res, res, 4Arblib.sqrt_analytic!(zero(tmp), -log(y), analytic) / y)

        # res -= -2sqrt(-log(y)) / log(y)
        Arblib.sub!(
            res,
            res,
            -2Arblib.sqrt_analytic!(zero(tmp), (-log(y)), analytic) / log(y),
        )

        return
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

        return weight
    end

    W(x) = 1 / (x^2 * log(x) * sqrt(log((x + 1) / x)))

    return x -> begin
        x = convert(Arb, x)

        J111 = 8 // 3 * sqrt(-log(ϵ)) * log(ϵ)
        #J112 = 8 // 3 * (-log(x))^(3 // 2) # We don't actually use this value in the end
        J12 = begin
            # Use that erf(sqrt(-log(x))) is monotonically increasing and bounded by 1
            if iszero(x)
                erfsqrtlogx = one(x)
            elseif Arblib.contains_zero(x)
                erfsqrtlogx = union(
                    SpecialFunctions.erf(sqrt(-log(Arblib.ubound(Arb, x)))),
                    one(x),
                )
            else
                erfsqrtlogx = SpecialFunctions.erf(sqrt(-log(x)))
            end

            2sqrt(oftype(x, π)) * (erfsqrtlogx - SpecialFunctions.erf(sqrt(-log(ϵ))))
        end
        J13 = real(
            Arblib.integrate!(
                integrand_J13!,
                Acb(prec = precision(x)),
                x,
                ϵ,
                check_analytic = true,
                rtol = 1e-10,
                atol = 1e-10,
                warn_on_no_convergence = false,
            ),
        )

        #J = J111 + J112 + J12 + J13 + J2

        # Variables for storing temporary values during integration
        x_complex = convert(Acb, x)
        tmp = zero(x_complex)

        integrand_E!(res, y; analytic::Bool) = begin
            # The code below is an inplace version of the following code
            #res = -log(sin((y - x) / 2) * sin((x + y) / 2) / sin(y / 2)^2)
            #res -= x^2 / 4 * (1 + cot(y / 2)^2)
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

            # res -= x^2 / 4 * (1 + cot(y / 2)^2)
            Arblib.mul_2exp!(tmp, y, -1)
            Arblib.cot!(tmp, tmp)
            Arblib.sqr!(tmp, tmp)
            Arblib.add!(tmp, tmp, 1)
            Arblib.mul!(tmp, tmp, x_complex)
            Arblib.mul!(tmp, tmp, x_complex)
            Arblib.div!(tmp, tmp, 4)
            Arblib.sub!(res, res, tmp)

            # tmp = y * sqrt(log((y + 1) / y))
            Arblib.add!(tmp, y, 1)
            Arblib.div!(tmp, tmp, y)
            Arblib.log!(tmp, tmp)
            Arblib.sqrt_analytic!(tmp, tmp, analytic)
            Arblib.mul!(tmp, tmp, y)

            Arblib.mul!(res, res, tmp)

            return
        end

        # TODO: Allow asymptotic expansion
        E = real(
            Arblib.integrate!(
                integrand_E!,
                zero(x_complex),
                x + Arb(1e-20),
                π,
                check_analytic = true,
                rtol = 1e-20,
                atol = 1e-20,
                warn_on_no_convergence = false,
            ),
        )

        part1 = begin
            # PROVE: that this is monotonically decreasing in x and 1 at x = 0
            if iszero(x)
                value = one(x)
            elseif Arblib.contains_zero(x)
                xᵤ = ubound(Arb, x)
                value = union(one(x), sqrt(-log(xᵤ)) / sqrt(log((xᵤ + 1) / xᵤ)))
            else
                value = sqrt(-log(x)) / sqrt(log((x + 1) / x))
            end
            -2 // 3 * value
        end
        part2 = (J111 + J12 + J13 + J2) / 4(log(x) * sqrt(log((x + 1) / x)))
        part3 = W(x) * E

        res = part1 + part2 + part3

        return res / (π * eval_expansion(u0, u0_expansion_div_xlogx, x))
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
