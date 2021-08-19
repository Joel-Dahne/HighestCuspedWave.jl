"""
    T02(u0::BHAnsatz; δ2)
Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral T_{0,2} from the paper.
"""
function T02(
    u0::BHAnsatz,
    evaltype::Ball;
    δ2::Arb = Arb(1e-10),
    ϵ = 1e-2, # TODO: Use this
)
    f = T021(u0, evaltype; δ2)
    g = T022(u0, evaltype; δ2)
    return x -> begin
        return f(x) + g(x)
    end
end

"""
    T021(u0::BHAnsatz; δ2)
Computes the (not yet existing) integral T_{0,2,1} from the paper.

This method assumes that `x + δ2 < π`. Otherwise it's handled by
[`T02`](@ref) directly.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval. We can also notice that the value inside the absolute value
is negative so we can remove the absolute value by putting a minus
sign, which we can bake in to the weight factor.

We are left with integrating the log-term.  This allows us to split the integrand
into three terms
1. `log(-sin((x - y) / 2)) = log(sin((y - x) / 2))`
2. `log(sin((x + y) / 2))`
3. `-2log(sin(y / 2))`

For the first term we use the inequality
```
c * (y - x) / 2 <= sin((y - x) / 2) <= (y - x) / 2
```
which holds for `c = sin(δ2 / 2) / (δ2 / 2)` on `0 <= x <= π`
and `x <= t <= x + δ2`. This gives us
```
log(c * (y - x) / 2) <= log(sin((y - x) / 2)) <= log((y - x) / 2)
```
The same inequality holds after integration from `x` to `x + δ2` and
gives us
```
δ2 * (log(c * δ2 / 2) - 1) <= ∫log(sin((y - x) / 2)) <= δ2 * (log(δ2 / 2) - 1)
```

The two remaining terms are both well behaved when `x` is bounded away
from `0` and `π`, which we assume is the case as mentioned above. We
can thus enclose the integral by directly enclosing the integrands on
the interval and multiplying with the size of the interval.
"""
function T021(u0::BHAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-10))
    return x -> begin
        x = convert(Arb, x)

        if !(x + δ2 < π)
            @warn "we don't have x + δ2 < π as required, x + δ2 = $(x + δ2)"
            return Arb(NaN)
        end

        interval = Arb((x, x + δ2))

        weight_factor = -u0.w(interval)

        part1 = let c = sin(δ2 / 2) / (δ2 / 2)
            part1_lower = δ2 * (log(c * δ2 / 2) - 1)
            part1_upper = δ2 * (log(δ2 / 2) - 1)

            Arb((part1_lower, part1_upper))
        end

        part2 = log(sin((x + interval) / 2)) * δ2

        part3 = -2log(sin(interval / 2)) * δ2

        integral = weight_factor * (part1 + part2 + part3)

        return integral / (π * u0.w(x) * u0(x))
    end
end

"""
    T022(u0::BHAnsatz; δ2)
Computes the (not yet existing) integral T_{0,2,2} from the paper.

This is done by directly computing the integral with the integrator in
Arb. Accounting for the fact that the integrand is non-analytic at `t
= x`.

TODO: This doesn't work when `x` is close to `π`.
"""
function T022(u0::BHAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-10))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    return x -> begin
        x = convert(Arb, x)
        a = x + δ2
        b = π

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
            b,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
        )
        @assert !isfinite(res) || isreal(res)
        res = real(res)

        return res / (π * u0.w(x) * u0(x))
    end
end
