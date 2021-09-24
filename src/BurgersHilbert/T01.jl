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

It precomputes the expansions of `u0` and for that reason a number `ϵ`
has to be given, the resulting expansion will be valid for all `x <
ϵ`. The value of `ϵ` has to be less than `1`.

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
I = ∫-(log(sin((y - x) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2))) * y * sqrt(log((y + 1) / y)) dy
```
for `0` to `x`.

The change of coordinates `t = y / x` transforms the integral into
```
I = x^2 * ∫abs(log(sin(x * (1 - t) / 2) + sin(x * (1 + t) / 2) - 2log(sin(x * t / 2)))) * t * sqrt(log((t * x + 1) / (t * x))) dt
```
from `0` to `1`.

Now consider the expansions
```
log(sin(x * (1 - t) / 2)) = log(x * (1 - t) / 2) + R₁(x)
log(sin(x * (1 + t) / 2)) = log(x * (1 + t) / 2) + R₂(x)
log(sin(x * t / 2)) = log(x * t / 2) + R₃(x)
sqrt(log((x * t + 1) / (x * t))) = sqrt(-log(t * x)) + R₄(x)
```
Where the remained terms `Rᵢ(x)` are all `O(x^2)` This allows us to
write the integral as
```
I = x^2 * ∫abs(log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2) + R₁(x) + R₂(x) - 2R₃(x)) * t * (sqrt(-log(t * x)) + R₄(x)) dt
```
Since we are only looking for an upper bound of the norm we can split
the absolute value in two and also split the `sqrt(-log(t * x)) +
R₄(x)` term in two. This gives us 4 integrals where we also skip the
`x^2` factor
```
I₁ = x^2 * ∫abs(log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2)) * t * sqrt(-log(t * x)) dt
I₂ = x^2 * ∫abs(log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2)) * t * R₄(x) dt
I₃ = x^2 * ∫abs(R₁(x) + R₂(x) - 2R₃(x)) * t * sqrt(-log(t * x)) dt
I₄ = x^2 * ∫abs(log(R₁(x) + R₂(x) - 2R₃(x)) * t *  R₄(x) dt
```
satisfying `I <= x^2 * (I₁ + I₂ + I₃ + I₄)`.

For `I₁` we can simplify the integrand further using
```
log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2) = log(1 / t^2 - 1)
```
giving us
```
I₁ = ∫abs(log(1 / t^2 - 1)) * t * sqrt(-log(t * x)) dt
```
Using that
```
sqrt(-log(t * x)) = sqrt(log(1 / x) + log(1 / t)) <= sqrt(log(1 / x)) + sqrt(log(1 / t))
```
We split this into
```
I₁₁ = sqrt(log(1 / x)) * ∫abs(log(1 / t^2 - 1)) * t dt
I₁₂ = ∫abs(log(1 / t^2 - 1)) * t * sqrt(log(1 / t)) dt
```
with `I₁ <= I₁₁ + I₁₂`. We can get rid of the absolute value by
splitting the integral at the only zero of `log(1 / t^2 - 1)` which is
given by `t = 1 / sqrt(2)`. The integral for `I₁₁` can be computed
explicitly and is given by `log(2)`. The integrand for `I₁₂` is
bounded and can be computed using `Arblib.integrate`, we just have to
be careful with handling `x` close to zero.
- **PROVE**: The value for the integral for `I₁₁`, Mathematica gives this.

For `I₂` we see the same simplification for the absolute value factor,
giving us
```
I₂ = ∫abs(log(1 / t^2 - 1)) * t * R₄(x) dt
```
Note that `
```
0 <= R₄(x) = sqrt(log((x * t + 1) / (x * t))) - sqrt(-log(t * x)) <= sqrt(log(t * x + 1))
```
where we have used that `sqrt(a + b) - sqrt(b) <= sqrt(a)`. Notice
that `sqrt(log(t * x + 1))` is monotone in both `x` and `t` and zero
for `x = 0`. This means that we can compute an enclosure `R4` of
`R₄(x)` by taking the lower endpoint to be `0` and the upper endpoint
to be `sqrt(log(x + 1))`. This gives us
```
I₂ = R4 * ∫abs(log(1 / t^2 - 1)) * t dt = x^2 * R4 * log(2)
```
where the value `log(2)` comes from noticing that it's the same
integral as for `I₁`.

For `I₃` we begin by noticing that
```
R₁(x) = log(sin(x * (1 - t) / 2)) - log(x * (1 - t) / 2) = log(sin(x * (1 - t) / 2) / (x * (1 - t) / 2)) = log(sinc(x * (1 - t) / 2π))
R₂(x) = log(sin(x * (1 + t) / 2)) - log(x * (1 + t) / 2) = log(sin(x * (1 + t) / 2) / (x * (1 + t) / 2)) = log(sinc(x * (1 + t) / 2π))
R₃(x) = log(sin(x * t / 2)) - log(x * t / 2) = log(sin(x * t / 2) / (x * t / 2)) = log(sinc(x * t / 2π))
```
where we have used the Julia convention that `sinc(x) = sinpi(x) / (π
* x)`. This gives us that
```
I₃ = ∫abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π))) * t * sqrt(-log(t * x)) dt
```
Similar to for `I₁` we use that `sqrt(-log(t * x)) <= sqrt(log(1 / x))
+ sqrt(log(1 / t))` to split it into
```
I₃₁ = sqrt(log(1 / x)) * ∫abs(log(sinc(x * (1 + t) / 2π)) + log(sinc(x * (1 - t) / 2π)) - 2log(sinc(x * t / 2π))) * t dt
I₃₂ = ∫abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π))) * t * sqrt(log(1 / t)) dt
```
The expression inside the absolute value is close to constant in `t`
on the interval and we therefore compute an enclosure of it an factor
it out of the integral. Notice that `x * (1 - t) / 2π`, `x * (1 + t) /
2π` and `x * t / 2π` all are less than `x / π`. By explicitly bounding
the second order error term, `D`, in the Taylor expansion of
`log(sinc(y))` at `y = 0` on the interval `[0, x / π]` we get that
```
log(sinc(x)) = D * x^2
```
In particular this gives us
```
log(sinc(x * (1 - t) / 2π)) = D * x^2 * (1 - t)^2 / (2π)^2
log(sinc(x * (1 + t) / 2π)) = D * x^2 * (1 + t)^2 / (2π)^2
log(sinc(x * t / 2π)) = D * x^2 * t^2 / (2π)^2
```
and hence
```
log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π)) = D / 2π^2 * x^2
```
This gives us
```
I₃₁ = x^2 * sqrt(log(1 / x)) * abs(D) / 2π^2 * ∫ t dt = x^2 * sqrt(log(1 / x)) * abs(D) / (2π)^2
```
and
```
I₃₂ = x^2 * abs(D) / 2π^2 * ∫ t * sqrt(log(1 / t)) dt
```
The remaining integral can be computed exactly and is given by `sqrt(π
/ 2) / 4`. To bound `D` we need to bound the second derivative of
`log(sinc(y))` on the interval `[0, x / π]`, the second derivative is
monotonically decreasing on the interval and can hence be enclosed by
evaluating it on the endpoints.
- **PROVE**: That the second derivative of `log(sinc(y))` is
    monotonically decreasing on `[0, x / π]`, or find another way to
    enclose `D`

For `I₄` we get, using the above computed expressions for `R₁`, `R₂`,
`R₃` and the enclosure `R4` for `R₄` to rewrite
```
I₄ = R4 * ∫abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 - t) / 2π)) - 2log(sinc(x * t / 2π))) * t dt
```
The integral is now the same as that for `I₃₁` and we hence get
```
I₄ = R4 * x^2 * abs(D) * ∫t dt = x^2 * R4 * abs(D) / 2
```

To put everything together we notice that all the different
subintegrals except `I₁₁` are bounded in `x`. We can therefore compute
an enclosure by just multiplying with the factor
```
inv(log(x) * sqrt(log((x + 1) / x)))
```
which is zero at `x = 0` and monotonically decreasing on `[0, 1]`. The
term `I₁₁ = sqrt(log(1 / x)) * log(2)` is not bounded in `x` and we
need to handle it slightly different. We instead consider
```
log(2) * sqrt(log(1 / x)) * inv(log(x) * sqrt(log((x + 1) / x)))
```
which again is zero at `x = 0` and monotonically decreasing.
- **PROVE**: That these two functions indeed are monotonically
    decreasing.
"""
function T01(u0::BHAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    ϵ = convert(Arb, ϵ)
    @assert ϵ < 0.5

    u0_expansion = u0(ϵ, AsymptoticExpansion())
    u0_expansion_div_xlogx = empty(u0_expansion)
    for ((i, m, k, l), value) in u0_expansion
        u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
    end

    # This gives the factor x * log(x) / (π * u0(x))
    factor(x) = inv(π * eval_expansion(u0, u0_expansion_div_xlogx, x))

    # This gives the factor inv(log(x) * sqrt(log((x + 1) / x))) in a
    # way so that it can handle x including zero
    weight_factor(x) = begin
        iszero(x) && return zero(x)

        if Arblib.contains_zero(x)
            # Use that inv(log(x) * sqrt(log((x + 1) / x))) is
            # monotonically decreasing for 0 < x < 1.
            xᵤ = ubound(Arb, x)

            return Arb((inv(log(xᵤ) * sqrt(log((xᵤ + 1) / xᵤ))), 0))
        end

        return inv(log(x) * sqrt(log((x + 1) / x)))
    end

    # Enclose the rest term D in the Taylor expansion of log(sinc(y))
    # on [0, ϵ / π]. Arb is not able to evaluate this directly so
    # instead we use that it is monotone so we only have to evaluate
    # it at the endpoints.
    Dᵤ = log(sinc(ArbSeries([0, 1, 0])))[2]
    Dₗ = log(sinc(ArbSeries([ϵ / π, 1, 0])))[2]
    D = Arb((Dₗ, Dᵤ))

    return x -> begin
        x = convert(Arb, x)
        @assert x <= ϵ

        # Compute the value of R4, using that it's zero at x = 0 and
        # monotone in both x and t, hence we only need to evaluate it
        # at the upper bound for x and t = 1.
        if iszero(x)
            R4 = zero(x)
        else
            xᵤ = ubound(Arb, x)
            R4 = Arb((0, sqrt(log1p(xᵤ))))
        end

        # This term includes the weight factor
        I₁₁ = begin
            if iszero(x)
                zero(x)
            elseif Arblib.contains_zero(x)
                # Use that sqrt(log(1 / x)) / (log(x) * sqrt(log((x + 1) / x)))
                # is zero at x = 0 and monotonically decreasing
                xᵤ = ubound(Arb, x)

                log(Arb(2)) *
                Arb((sqrt(log(1 / xᵤ)) / (log(xᵤ) * sqrt(log((xᵤ + 1) / xᵤ))), 0))
            else
                log(Arb(2)) * sqrt(log(1 / x)) / (log(x) * sqrt(log((x + 1) / x)))
            end
        end

        I₁₂ = begin
            integrand_I12(t; analytic) = begin
                if Arblib.contains_zero(t)
                    analytic && return Arblib.indeterminate!(zero(t))
                    @assert isreal(t)

                    # Use monotonicity of the integrand close to zero
                    # PROVE: That it is monotonic and also on which interval.
                    tᵤ = ubound(Arb, real(t))
                    tᵤ < 0.1 || return Arblib.indeterminate!(zero(t)) # FIXME: Put proper value
                    return Acb((zero(tᵤ), log(1 / tᵤ^2 - 1) * tᵤ * sqrt(log(1 / tᵤ))))
                elseif Arblib.overlaps(t, one(t))
                    analytic && return Arblib.indeterminate!(zero(t))
                    @assert isreal(t)

                    # Use monotonicity of the integrand close to one
                    # PROVE: That it is monotonic and also on which interval.
                    tₗ = lbound(Arb, real(t))
                    tₗ > 0.99 || return Arblib.indeterminate!(zero(t)) # FIXME: Put proper value
                    return Acb((log(1 / tₗ^2 - 1) * tₗ * sqrt(log(1 / tₗ)), zero(tₗ)))
                else
                    return log(1 / t^2 - 1) *
                           t *
                           Arblib.real_sqrtpos!(zero(t), log(1 / t), analytic)
                end
            end

            # Integrate from 0 to 1 / sqrt(2) and then from 1 / sqrt(2) to 1
            integral1 = real(
                Arblib.integrate(
                    integrand_I12,
                    0,
                    1 / sqrt(Arb(2)),
                    check_analytic = true,
                    atol = 1e-10,
                ),
            )
            integral2 = real(
                Arblib.integrate(
                    (t; analytic) -> -integrand_I12(t; analytic),
                    1 / sqrt(Arb(2)),
                    1,
                    check_analytic = true,
                    atol = 1e-10,
                ),
            )

            integral1 + integral2
        end

        I₂ = R4 * log(Arb(2))

        I₃₁ = let π = Arb(π)
            if iszero(x)
                xpart = zero(x)
            elseif Arblib.contains_zero(x) && x < exp(Arb(-1 // 4))
                # Use that x^2 * sqrt(log(1 / x)) is monotone for 0 < x < exp(-1 // 4)
                xᵤ = ubound(Arb, x)
                xpart = Arb((0, xᵤ^2 * sqrt(log(1 / xᵤ))))
            else
                xpart = x^2 * sqrt(log(1 / x))
            end
            xpart * abs(D) / (2π)^2
        end

        I₃₂ = let π = Arb(π)
            x^2 * abs(D) / 2π^2 * sqrt(π / 2) / 4
        end

        I₄ = x^2 * R4 * abs(D) / 2

        # Note that I₁₁ already includes the weight factor
        return factor(x) * (I₁₁ + weight_factor(x) * (I₁₂ + I₂ + I₃₁ + I₃₂ + I₄))
    end
end

"""
    T011(u0::BHAnsatz; δ0)
Computes the integral T_{0,1,1} from the paper.

It uses the fact that the integrand is strictly increasing on the
interval `[0, 0.5]` for every value of `x` and 0 at `x = 0`. This
allows us to enclose the integrand on the interval which then easily
gives an enclosure of the integral by multiplying with the size of the
interval.

**PROVE**: That the integrand indeed is increasing on the said interval.
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

- **TODO**: The upper bound for when `x` contains `π` could be
    improved, but this is probably not needed.
- **PROVE**: There are several minor details here that might need to
    be proved.
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
