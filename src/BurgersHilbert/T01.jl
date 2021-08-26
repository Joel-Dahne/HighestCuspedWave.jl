"""
    T01(u0::BHAnsatz, ::Ball; δ1, δ2, skip_div_u0)
Returns a function such that `T01(u0, Ball(); δ1, δ2)(x)` computes the
integral T_{0,1} from the paper.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T01(
    u0::BHAnsatz,
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
    T01(u0::BHAnsatz, ::Asymptotic)
Returns a function such that `T01(u0, Asymptotic())(x)` computes an
**upper bound** of the integral T_{0,1} from the paper using an
evaluation strategy that works asymptotically as `x` goes to 0.

The integral is split into one main term and several error terms.

The integral for the main term is `∫|log(1 - t) + log(1 + t) -
2log(t)|t dt` on the interval `[0, 1]`. Which turns out to be equal to
`log(2)`

PROVE: That it equals `log(2)` (Mathematica handles it)

TODO: Handle the error terms
"""
function T01(
    u0::BHAnsatz,
    ::Asymptotic;
    non_asymptotic_u0 = false,
    precompute_u0_expansion = false,
    ϵ = nothing,
)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    @warn "T01(u0, Asymptotic()) doesn't bound the error term yet"

    if precompute_u0_expansion
        ϵ = convert(Arb, ϵ)
        u0_expansion = u0(ϵ, AsymptoticExpansion())

        u0_expansion_div_xlogx = empty(u0_expansion)
        for ((i, m, k, l), value) in u0_expansion
            u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
        end
    else
        u0_expansion_div_xlogx = nothing
    end

    factor(x) = begin
        # PROVE: That this is monotonically decreasing on [0, 1]
        w(x) = sqrt(log(1 / x)) / (log(x) * sqrt(log((x + 1) / x)))

        if iszero(x)
            weight = zero(x)
        elseif Arblib.contains_zero(x) && x < 1
            weight = union(zero(x), w(Arblib.ubound(Arb, x)))
        else
            weight = w(x)
        end

        if !precompute_u0_expansion
            u0_expansion = u0(x, AsymptoticExpansion())

            u0_expansion_div_xlogx = empty(u0_expansion)
            for ((i, m, k, l), value) in u0_expansion
                u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
            end
        end

        return weight / (π * eval_expansion(u0, u0_expansion_div_xlogx, x))
    end

    return x -> begin
        x = convert(Arb, x)

        @assert !precompute_u0_expansion || x <= ϵ

        main_term = log(Arb(2))

        error_term = zero(x) # TODO

        if non_asymptotic_u0
            return x^2 * sqrt(log(1 / x)) / (π * u0.w(x) * u0(x)) * (main_term + error_term)
        else
            return factor(x) * (main_term + error_term)
        end
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
function T011(u0::BHAnsatz, ::Ball = Ball(); δ0::Arb = Arb(1e-5), skip_div_u0 = false)
    δ0 < 0.5 || Throw(ArgumentError("δ0 must be less than 0.5, got $δ0"))
    return x -> begin
        x = convert(Arb, x)

        integrand(t) =
            abs(log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(t * x / 2)^2)) *
            t *
            sqrt(log((t * x + 1) / (t * x)))

        integral = δ0 * Arb((0, integrand(δ0)))

        res = integral * x / (π * sqrt(log((x + 1) / x)))
        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
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
function T012(
    u0::BHAnsatz,
    ::Ball = Ball();
    δ0::Arb = Arb(1e-5),
    δ1::Arb = Arb(1e-5),
    skip_div_u0 = false,
)
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

        # Variables for storing temporary values during integration
        x_complex = convert(Acb, x)
        xdiv2 = x_complex / 2
        tmp = zero(x_complex)
        tx = zero(x_complex)

        # PROVE: That there are no branch cuts that interact with the
        # integral
        integrand!(res, t; analytic::Bool) = begin
            # The code below is an inplace version of the following code
            #res = log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(t * x / 2)^2)
            #Arblib.real_abs!(res, res, analytic)
            #weight = t * Arblib.sqrt_analytic!(zero(t), log((t * x + 1) / (t * x)), analytic)
            #return res * weight

            Arblib.mul!(tx, t, x_complex)

            # res = sin((1 - t) * x / 2)
            Arblib.neg!(tmp, t)
            Arblib.add!(tmp, tmp, 1)
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

            Arblib.real_abs!(res, res, analytic)

            # tmp = t * sqrt(log((tx + 1) / tx))
            Arblib.add!(tmp, tx, 1)
            Arblib.div!(tmp, tmp, tx)
            Arblib.log!(tmp, tmp)
            Arblib.sqrt_analytic!(tmp, tmp, analytic)
            Arblib.mul!(tmp, tmp, t)

            Arblib.mul!(res, res, tmp)

            return
        end

        res = Arblib.integrate!(
            integrand!,
            zero(x_complex),
            a,
            b,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
            warn_on_no_convergence = false,
        )

        @assert !isfinite(res) || isreal(res)
        res = real(res)

        res = res * x / (π * sqrt(log((x + 1) / x)))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
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
We use this bound if `x` is less than `sqrt(eps(x))` away from `π`.

TODO: The upper bound for when `x` contains `π` could be improved, but
this is probably not needed.
PROVE: There are several minor details here that might need to be
proved.
"""
function T013(u0::BHAnsatz, ::Ball = Ball(); δ1::Arb = Arb(1e-5), skip_div_u0 = false)
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
            if x < π - sqrt(eps(x))
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

        res = integral * x / (π * sqrt(log((x + 1) / x)))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
