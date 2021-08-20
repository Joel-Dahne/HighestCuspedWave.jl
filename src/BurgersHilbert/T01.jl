"""
    T01(u0::BHAnsatz, ::Ball; δ1, δ2)
Returns a function such that `T01(u0, Ball(); δ1, δ2)(x)` computes the
integral T_{0,1} from the paper.
"""
function T01(u0::BHAnsatz, evaltype::Ball; δ0::Arb = Arb(1e-10), δ1::Arb = Arb(1e-10))
    f = T011(u0, evaltype; δ0)
    g = T012(u0, evaltype; δ0, δ1)
    h = T013(u0, evaltype; δ1)
    return x -> begin
        return f(x) + g(x) + h(x)
    end
end

"""
    T011(u0::BHAnsatz; δ0)
Computes the integral T_{0,1,1} from the paper.

It uses the fact that the integrand is strictly increasing on the
interval `[0, 0.5]` for every value of `x` and 0 at 0. This allows us
to enclose the integrand on the interval which then easily gives an
enclosure of the integral by multiplying with the size of the
interval.

PROVE: That the integrand indeed is increasing on the said interval.
"""
function T011(u0::BHAnsatz, ::Ball = Ball(); δ0::Arb = Arb(1e-10))
    δ0 < 0.5 || Throw(ArgumentError("δ0 must be less than 0.5, got $δ0"))
    return x -> begin
        x = convert(Arb, x)

        integrand(t) =
            abs(log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(t * x / 2)^2)) *
            t *
            sqrt(log((t * x + 1) / (t * x)))

        integral = δ0 * Arb((0, integrand(δ0)))

        return integral * x / (π * sqrt(log((x + 1) / x)) * u0(x))
    end
end

"""
    T012(u0::BHAnsatz; δ0, δ1)
Returns a function such that `T012(u0; δ0, δ1)(x)` computes the integral
T_{0,1,2} from the paper.

This is done by directly computing the integral with the integrator in
Arb. Accounting for the fact that the integrand is non-analytic at `t
= x`.
"""
function T012(u0::BHAnsatz, ::Ball = Ball(); δ0::Arb = Arb(1e-10), δ1::Arb = Arb(1e-10))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    a = δ0
    b = 1 - δ1

    return x -> begin
        x = convert(Arb, x)
        # PROVE: That there are no branch cuts that interact with the
        # integral
        integrand(t; analytic::Bool) = begin
            tx = t * x
            res = log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(tx / 2)^2)
            Arblib.real_abs!(res, res, analytic)
            weight = t * Arblib.sqrt_analytic!(zero(t), log((tx + 1) / tx), analytic)
            return res * weight
        end

        res = Arblib.integrate(
            integrand,
            a,
            b,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
        )
        @assert !isfinite(res) || isreal(res)
        res = real(res)

        return res * x / (π * sqrt(log((x + 1) / x)) * u0(x))
    end
end

"""
    T013(u0::BHAnsatz; δ1)
Computes the integral T_{0,1,3} from the paper.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval.

We are left with integrating the log-term. By noticing that the value
inside the absolute value is negative so we can remove the absolute
value by putting a minus sign. This allows us to split the integrand
into three terms
1. `log(sin(x * (1 - t) / 2))`
2. `log(sin(x * (1 + t) / 2))`
3. `-2log(sin(x * t / 2))`

For the first term we have the inequality
```
c * x * (1 - t) / 2 <= sin(x * (1 - t) / 2) <= x * (1 - t) / 2
```
which holds for `c = sin(x * δ1 / 2) / (x * δ1 / 2)` on `0 <= x <= π`
and `1 - δ1 <= t <= 1`. This gives us
```
log(c * x * (1 - t) / 2) <= log(sin(x * (1 - t) / 2)) <= log(x * (1 - t) / 2)
```
and the same inequality holds after integration for `1 - δ1` to `1`
and gives us
```
δ1 * (log(c * x * δ1 / 2) - 1) <= ∫log(sin(x * (1 - t) / 2)) <= δ1 * (log(x * δ1 / 2) - 1)
```

The third term is well behaved for any `x` bounded away from zero and
`t` on the interval. We can thus enclose the integral by directly
enclosing the integrand on the interval and multiplying by the size of
the interval.

If `x` is bounded away from `π` then the second term is also well
behaved and we can treat it similarly to the third term. However, if
`x` overlaps with `π` then it becomes singular and we have to handle
it differently. We first do the change of variables `t → 1 - s`,
together with the identity `sin(x) = sin(π - x)`, giving us the
integral `∫ log(sin(π - x + x * s / 2)) ds` from `0` to ` δ1`. We then
take the very crude upper bound `sin(π - x + x * s / 2) <= 1` and for
the lower bound we notice that `sin(π - x + x * s / 2)` is minimized
by taking `x = π`, giving us the lower bound `sin(x * s / 2) <= sin(x
- x * s / 2), which then gives the lower bound `x * s / 2`. This means
  we have the bounds
```
x * s / 2 <= sin(π - x + x * s / 2) <= 1
```
and hence
```
log(x * s / 2) <= log(sin(π - x + x * s / 2)) <= log(1) = 0
```
After integration we thus get
```
δ1 * (log(c * x * δ1 / 2) - 1) <= ∫log(sin(π - x + x * s / 2)) <= 0
```

TODO: The upper bound for when `x` contains `π` could be improved, but
this is probably not needed.
PROVE: There are several minor details here that might need to be
proved.
"""
function T013(u0::BHAnsatz, ::Ball = Ball(); δ1::Arb = Arb(1e-10))
    return x -> begin
        x = convert(Arb, x)

        weight_factor = let t = Arb((1 - δ1, 1))
            -t * sqrt(log((t * x + 1) / (t * x)))
        end

        part1 = let c = sin(x * δ1 / 2) / (x * δ1 / 2)
            part1_lower = δ1 * (log(c * x * δ1 / 2) - 1)
            part1_upper = δ1 * (log(x * δ1 / 2) - 1)

            Arb((part1_lower, part1_upper))
        end

        part2 = let t = Arb((1 - δ1, 1))
            if x < π
                log(sin(x * (1 + t) / 2)) * δ1
            else
                c = sin(x * δ1 / 2) / (x * δ1 / 2)

                part2_lower = δ1 * (log(c * x * δ1 / 2) - 1)
                part2_upper = zero(x)

                Arb((part2_lower, part2_upper))
            end
        end

        part3 = let t = Arb((1 - δ1, 1))
            -2log(sin(x * t / 2)) * δ1
        end

        integral = weight_factor * (part1 + part2 + part3)

        return integral * x / (π * sqrt(log((x + 1) / x)) * u0(x))
    end
end
