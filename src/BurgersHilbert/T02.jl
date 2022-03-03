"""
    T02(u0::BHAnsatz; δ2, skip_div_u0)
Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral T_{0,2} from the paper.

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
enclosure of the integral `T02` using an evaluation strategy that
works asymptotically as `x` goes to 0.

It precomputes the expansions of `u0` and for that reason a number `ϵ`
has to be given, the resulting expansion will be valid for all `x <
ϵ`. The value of `ϵ` has to be less than `1 // 2`.

To begin with the factor `x * log(inv(x)) / (π * u0(x))` is factored
out from the whole expression and multiplied back in the end. Notice
that this factor is positive and bounded in `x`.

What we are left with computing is
```
W(x) * I
```
where `W(x) = 1 / (x^2 * log(inv(x)) * sqrt(log(1 + inv(x))))` and `I`
is given by the integral
```
∫-(log(sin((y - x) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2))) * y * sqrt(log(1 + inv(y))) dy
```
for `x` to `π`.

The change of coordinates `t = y / x` transforms the integral into
```
I = x^2 * ∫-(log(sin(x * (t - 1) / 2)) + log(sin(x * (t + 1) / 2)) - 2log(sin(x * t / 2))) * t * sqrt(log(1 + inv(x * t))) dt
```
for `1` to `π / x`.

Using that
```
log(sin(x * (t - 1) / 2)) = log(x * (t - 1) / 2) + log(sinc(x * (t - 1) / 2π))
log(sin(x * (1 + t) / 2)) = log(x * (1 + t) / 2) + log(sinc(x * (1 + t) / 2π))
log(sin(x * t / 2)) = log(x * t / 2) + log(sinc(x * t / 2π))
```
where we have used the Julia convention that `sinc(x) = sinpi(x) / (π
* x)`. we can split `I` as
```
I <= x^2 * (I1 + I2)
```
where
```
I1 = ∫ -(log(x * (t - 1) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2)) * t * sqrt(log(1 + inv(x * t))) dt
I2 = ∫ -(log(sinc(x * (t - 1) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π))) * t * sqrt(log(1 + inv(x * t))) dt
```
We now handle `I1` and `I2` separately.

# Handling `I1`
We start by noticing that
```
log(x * (t - 1) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2) =
    log(1 - 1 / t^2)
```
giving us
```
I1 = ∫ -log(1 - 1 / t^2) * t * sqrt(log(1 + inv(x * t))) dt
```
We now split the interval of integration, `[1, π / x]`, into three
parts, `[1, 2], [2, inv(x)], [inv(x), π / x]`, and treat each interval
separately. Let `I11, I12, I13` correspond to the integrals on the
separate intervals.

## Handling `I11`
We begin by noting that
```
sqrt(log(1 + inv(x * t))) = sqrt(log(x / (x + 1) + 1 / (t * (x + 1))) + log(1 + inv(x)))
```
and that `log(x / (x + 1) + 1 / (t * (x + 1)))` is decreasing in `t`,
an upper bound is thus given by its value at `t = 1` where it is zero.
Hence
```
sqrt(log(1 + inv(x * t))) = sqrt(log(x / (x + 1) + 1 / (t * (x + 1))) + log(1 + inv(x)))
    <= sqrt(log(1 + inv(x)))
```
and we have
```
I11 <= sqrt(log(1 + inv(x))) * ∫_1^2 -log(1 - 1 / t^2) * t dt
```
where the integral can now be explicitly computed to be `log(16 / (3 *
sqrt(3)))`, giving us
```
I11 <= sqrt(log(1 + inv(x))) * log(16 / (3 * sqrt(3)))
```

## Handling `I12`
Taylor expanding `-log(1 - 1 / t^2)` at `t = ∞` gives us
```
-log(1 - 1 / t^2) = 1 / t^2 + D2 / t^4,
```
where `D2` is a bound for the remainder term, which we can compute
directly using `ArbSeries`. This together with
```
sqrt(log(1 + inv(x * t))) =
    sqrt(log(1 + x * t) + log(inv(x * t))) <=
    sqrt(log(1 + x * t)) + sqrt(log(inv(x * t)))
```
and `sqrt(log(1 + x * t)) < sqrt(log(2))` gives us
```
I12 <=
    ∫ sqrt(log(inv(x * t))) / t dt
    + D2 * ∫ sqrt(log(inv(x * t))) / t^3 dt
    + sqrt(log(2)) ∫ 1 / t dt
    + D2 * sqrt(log(2)) * ∫ 1 / t^3 dt
```
For the integrals we have
```
∫ sqrt(log(inv(x * t))) / t dt = 2 // 3 * log(inv(2x))^(3 // 2)

∫ sqrt(log(inv(x * t))) / t^3 dt = (sqrt(2π) * x^2 * erfi(sqrt(log(inv(4x^2)))) - sqrt(log(inv(2x)))) / 8

∫ 1 / t dt = log(inv(x)) - log(2)

∫ 1 / t^3 dt = (1 - 4x^2) / 8
```
- **PROVE:** That these integrals are correct, in particular the
  second one.
To avoid having to use the `erfi` function we use that
```
x^2 * erfi(sqrt(log(inv(4x^2))))
= 2 / sqrt(π) * x^2 * ∫_0^sqrt(log(inv(4x^2))) exp(t^2) dt
<= 2 / sqrt(π) * x^2 (exp(sqrt(log(inv(4x^2)))^2) - 1) / sqrt(log(inv(4x^2)))
= 2 / sqrt(π) * x^2 (inv(4x^2) - 1) / sqrt(2log(inv(2x)))
= 1 / 2sqrt(2π) * (1 - 4x^2) / sqrt(log(inv(2x)))
```
This gives us the following upper bound for `I12`
```
I12 <=
    2 // 3 * log(inv(2x))^(3 // 2)
    + D2 * (1 - 4x^2) / 16sqrt(log(inv(2x)))
    - D2 * sqrt(log(inv(2x))) / 8
    + sqrt(log(2)) * log(inv(x))
    - log(2)^(3 // 2)
    + D2 * sqrt(log(2)) * (1 - 4x^2) / 8
```

## Handling `I13`
We use that `sqrt(log(1 + inv(x * t))) < sqrt(log(2))` for `1 / x < t
< π / x`, this gives us
```
I13 <= sqrt(log(2)) * ∫ -log(1 - 1 / t^2) * t dt
```
where the integral can be computed to be
```
∫-log(1 - 1 / t^2) * t dt =
    1 / (2x^2) * ((1 - x^2) * log(1 / x^2 - 1) + (x^2 - π^2) * log(π^2 / x^2 - 1) + 2(π^2 * log(π) + log(x) - π^2 * log(x)))
```
With some work this can be simplified to
```
log((π^2 - x^2) / (1 - x^2)) + log(1 - x^2) / 2x^2 - π^2 * log(1 - x^2 / π^2) / 2x^2
```
and hence
```
I13 <= sqrt(log(2)) * (log((π^2 - x^2) / (1 - x^2)) + log(1 - x^2) / 2x^2 - π^2 * log(1 - x^2 / π^2))
```

## Putting `I11`, `I12` and `I13` together
With all of this we get the following upper bound for `I1`
```
I1 <= sqrt(log(1 + inv(x))) * log(16 / (3 * sqrt(3)))
    + (
        2 // 3 * log(inv(2x))^(3 // 2)
        + D2 * (1 - 4x^2) / 16sqrt(log(inv(2x)))
        - D2 * sqrt(log(inv(2x))) / 8
        + sqrt(log(2)) * log(inv(x))
        - log(2)^(3 // 2)
        + D2 * sqrt(log(2)) * (1 - 4x^2) / 8
    ) + sqrt(log(2)) * (
        log((π^2 - x^2) / (1 - x^2))
        + log(1 - x^2) / 2x^2
        - π^2 * log(1 - x^2 / π^2) / 2x^2
    )
```
The terms `log(1 - x^2) / 2x^2` and `log(1 - x^2 / π^2) / 2x^2` can be
bounded using [`fx_div_x`](@ref).
To be able to bound this after division by `log(inv(x)) * sqrt(log(1 +
inv(x)))` there are some terms we have to take care of, all of them
coming from `I12`.
1. For the term `sqrt(log(2)) * log(inv(x))` we can directly cancel
   `log(inv(x))`.
2. For the term `2 // 3 * log(inv(2x))^(3 // 2)` we use that it is
  monotone in `x` after division by `log(inv(x)) * sqrt(log(1 +
  inv(x)))`.
  - **PROVE:** That this function indeed is monotone in `x`
3. For the term `- D2 * sqrt(log(inv(2x))) / 8` we use that
  `sqrt(log(inv(2x))) / sqrt(log(1 + inv(x)))` is monotone in `x`.
  - **PROVE:** That this function indeed is monotone in `x`

# Handling `I2`
In this case we switch back to the coordinates given by `y = xt`,
giving us
- **IMPROVE:** Consider never switching to it...
```
I2 = inv(x^2) * ∫ -(log(sinc((y - x) / 2π)) + log(sinc((x + y) / 2π)) - 2log(sinc(y / 2π))) * y * sqrt(log(1 + inv(y))) dt
```
from `x` to `π`. Using that
```
-(log(sinc((y - x) / 2π)) + log(sinc((x + y) / 2π)) - 2log(sinc(y / 2π))) / x^2
```
is bounded in `y` on `[0, π]` uniformly in `x` we can compute an
enclosure of the range.
If we let `D3` denote this enclosure we have
```
I2 = D3 * ∫ y * sqrt(log(1 + inv(y))) dy
```
in an interval sense. The integral can be enclosed directly using
[`Arblib.integrate`](@ref) and that the integrand is zero at `y = 0`
and increasing in `y` to handle `y` close to zero.

To enclose `D3` we use that the function is increasing in `y` and it
is hence enough to evaluate it at `y = 0` and `y = π`. At `y = 0` it
simplifies to
```
-2log(sinc(x / 2π)) / x^2
```
and for `y = π` we can write it as
```
-(log(sinc((x / π - 1) / 2)) + log(sinc((x / π + 1) / 2)) - 2log(sinc(1 / 2))) / x^2
```
In both cases we have used that `sinc` is even. Both of these have a
removable singularity at `x = 0` that we have to handle.
- **PROVE:** That the expression for `D3` is monotone in `y`.
"""
function T02(u0::BHAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    ϵ = convert(Arb, ϵ)
    @assert ϵ < 0.5

    # Compute expansion for u0 / (x * log(x))
    u0_expansion = u0(ϵ, AsymptoticExpansion())
    u0_expansion_div_xlogx = empty(u0_expansion)
    for ((i, m, k, l), value) in u0_expansion
        u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
    end

    # This gives the factor x * log(inv(x)) / (π * u0(x))
    factor(x) = -inv(π * eval_expansion(u0, u0_expansion_div_xlogx, x))

    # This gives the factor inv(log(inv(x)) * sqrt(log(1 + inv(x)))) in a
    # way so that it can handle x including zero
    weight_factor(x) = begin
        if iszero(x)
            return zero(x)
        elseif Arblib.contains_zero(x)
            # Use that inv(log(x) * sqrt(log(1 + inv(x)))) is
            # monotonically increasing for 0 < x < 1.
            return let xᵤ = ubound(Arb, x)
                Arb((0, inv(log(inv(xᵤ))) * sqrt(log(1 + inv(xᵤ)))))
            end
        else
            return inv(log(inv(x)) * sqrt(log(1 + inv(x))))
        end
    end

    # Bound the constant D2 using ArbSeries. We compute the error term
    # for log(1 - y) for y in [0, 1 // 4]
    D2 = log(1 - ArbSeries((Arb((0, 1 // 4)), 1, 0)))[2]

    D3 = begin
        # Lower bound at y = 0
        lower = -2fx_div_x(Arb((0, ϵ)), 2, extra_degree = 2) do x
            if x isa ArbSeries && Arblib.contains_zero(x[0])
                # sinc(::ArbSeries) doesn't handle overlap with zero
                log(fx_div_x(sin, x / 2, extra_degree = 2))
            else
                log(sinc(x / 2Arb(π)))
            end
        end

        # Upper bound at y = π
        upper =
            -fx_div_x(Arb((0, ϵ)), 2, force = true, extra_degree = 2) do x
                log(sinc((x / π - 1) / 2)) + log(sinc((x / π + 1) / 2)) -
                2log(sinc(Arb(1 // 2)))
            end

        Arb((lower, upper))
    end

    return x -> begin
        x = convert(Arb, x)
        @assert x <= ϵ

        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invlog = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            Arb((0, -inv(log(ubound(Arb, x)))))
        else
            -inv(log(x))
        end

        # Enclosure of inv(sqrt(log(1 + inv(x))))
        invsqrtlog = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            Arb((0, inv(sqrt(log(1 + inv(ubound(Arb, x)))))))
        else
            inv(sqrt(log(1 + inv(x))))
        end

        # Enclosure of weight_factor(x) * I1
        WxI1 = let
            WxI11 = log(16 / (3sqrt(Arb(3)))) * invlog

            WxI12 = let
                # Terms which require explicit cancellations

                # Enclosure of 2log(inv(2x))^(3 // 2) / 3 * weight_factor(x)
                term1 = if Arblib.contains_zero(x)
                    # Use that
                    # log(inv(2x))^(3 // 2) / (log(inv(x)) * sqrt(log(1 + inv(x))))
                    # converges to 1 at x = 0 and is decreasing
                    lower = if iszero(x)
                        one(x)
                    else
                        let xᵤ = ubound(Arb, x)
                            log(inv(2xᵤ))^(3 // 2) / (log(inv(xᵤ)) * sqrt(log(1 + inv(xᵤ))))
                        end
                    end
                    upper = one(x)
                    2Arb((lower, upper)) / 3
                else
                    2log(inv(2x))^(3 // 2) / 3 * weight_factor(x)
                end

                # Enclosure of -D2 * sqrt(log(inv(2x))) * weight_factor(x)
                term2 =
                    -D2 *
                    invlog *
                    if Arblib.contains_zero(x)
                        # Use that log(inv(2x)) / log(1 + inv(x))
                        # converges to 1 at x = 0 and is decreasing
                        lower = if iszero(x)
                            one(x)
                        else
                            let xᵤ = ubound(Arb, x)
                                sqrt(log(inv(2xᵤ))) / sqrt(log(1 + inv(xᵤ)))
                            end
                        end
                        upper = one(x)
                        Arb((lower, upper))
                    else
                        sqrt(log(inv(2x))) * invsqrtlog
                    end

                # Enclosure of sqrt(log(Arb(2))) * log(inv(x)) * weight_factor(x)
                # Explicitly cancelling the log(inv(x))
                term3 = sqrt(log(Arb(2))) * invsqrtlog

                # Terms for which we don't have to cancel things

                # Enclosure of inv(sqrt(log(inv(2x))))
                invsqrtloginv2x = if iszero(x)
                    zero(x)
                elseif Arblib.contains_zero(x)
                    lower = zero(x)
                    upper = inv(sqrt(log(inv(2ubound(Arb, x)))))
                else
                    inv(sqrt(log(inv(2x))))
                end

                remaining =
                    (
                        D2 * (1 - 4x^2) * invsqrtloginv2x / 16 - log(Arb(2))^(3 // 2) +
                        D2 * sqrt(log(Arb(2))) * (1 - 4x^2) / 8
                    ) * weight_factor(x)

                (term1 + term2 + term3) + remaining
            end

            # Enclosure of I13
            I13 =
                sqrt(log(Arb(2))) * if Arblib.contains_zero(x)
                    log((Arb(π)^2 - x^2) / (1 - x^2)) / 2 +
                    fx_div_x(y -> log(1 - y^2), x, 2) / 2 -
                    Arb(π)^2 * fx_div_x(y -> log(1 - y^2 / Arb(π)^2), x, 2) / 2
                else
                    log((Arb(π)^2 - x^2) / (1 - x^2)) / 2 + log1p(-x^2) / 2x^2 -
                    Arb(π)^2 * log1p(-x^2 / Arb(π)^2) / 2x^2
                end

            WxI11 + WxI12 + I13 * weight_factor(x)
        end

        # Enclosure of I2 / x^2
        I2divx2 = let
            # Enclosure of ∫y * sqrt(log((y + 1) / y)) dy from x to π
            integrand_I2(y; analytic) = begin
                if Arblib.contains_zero(y)
                    analytic && return Arblib.indeterminate!(zero(y))
                    @assert isreal(y)

                    # Use monotonicity to compute enclosure
                    yᵤ = ubound(Arb, real(y))
                    return Acb((0, yᵤ * sqrt(log(1 + inv(yᵤ)))))
                else
                    return y * Arblib.real_sqrtpos!(zero(y), log(1 + inv(y)), analytic)
                end
            end
            integral_I2 =
                real(Arblib.integrate(integrand_I2, x, π, check_analytic = true))

            D3 * integral_I2
        end

        return factor(x) * (WxI1 + I2divx2 * weight_factor(x))
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
