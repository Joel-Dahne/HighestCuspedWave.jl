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
    T01(u0::BHAnsatz, ::Asymptotic; ϵ = Arb(2e-1))

Return a function such that `T01(u0, Asymptotic())(x)` computes an
**upper bound** of `T01(u0)` using an evaluation strategy that works
asymptotically as `x` goes to 0.

It precomputes the expansions of `u0` and for that reason a number `ϵ`
has to be given, the resulting expansion will be valid for all `x <
ϵ`. The value of `ϵ` has to be less than `1`.

To begin with the factor `x * log(inv(x) / (π * u0(x))` is factored
out from the whole expression and multiplied back in the end. Notice
that this factor is positive and bounded in `x`.

What we are left with computing is
```
W(x) * I
```
where `W(x) = 1 / (x^2 * log(inv(x)) * sqrt(log(1 + inv(x))))` and `I`
is given by the integral
```
I = ∫abs(log(sin((x - y) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2))) * y * sqrt(log(1 + inv(y))) dy
```
for `0` to `x`.

The change of coordinates `t = y / x` transforms the integral into
```
I = x^2 * ∫abs(log(sin(x * (1 - t) / 2) + sin(x * (1 + t) / 2) - 2log(sin(x * t / 2)))) * t * sqrt(log(1 + inv(x * t))) dt
```
from `0` to `1`.

Using that
```
log(sin(x * (1 - t) / 2)) = log(x * (1 - t) / 2) + log(sinc(x * (1 - t) / 2π))
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
I1 = ∫abs(log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2)) * t * sqrt(log(1 + inv(x * t))) dt
I2 = ∫abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π))) * t * sqrt(log(1 + inv(x * t))) dt
```
We now handle `I1` and `I2` separately.

# Handling `I1`
We start by noticing that
```
log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2) =
    log(1 / t^2 - 1)
```
Furthermore, since `0 <= x <= 1` and `0 <= t <= 1` we have
```
sqrt(log(1 + inv(x * t))) = sqrt(log(inv(x)) + log(inv(t)) + log(1 + x * t))
    <= sqrt(log(inv(x))) + sqrt(log(inv(t))) + sqrt(log(1 + x * t))
```
and also `sqrt(log(1 + x * t)) <= sqrt(log(1 + x))`, giving us
```
I1 <= sqrt(log(inv(x))) * ∫ abs(log(1 / t^2 - 1)) * t dt
    + ∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt
    + sqrt(log(1 + x)) * ∫ abs(log(1 / t^2 - 1)) * t dt
```
If we let
```
C1 = ∫ abs(log(1 / t^2 - 1)) * t dt

C2 = ∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt
```
this can be written as
```
I1 <= sqrt(log(inv(x))) * C1 + C2 + sqrt(log(1 + x)) * C1
```
We thus have to compute `C1` and `C2`. We can get rid of the absolute
value by splitting the integral at the only zero of `log(1 / t^2 - 1)`
which is given by `t = 1 / sqrt(2)`. The integral for `C1` can be
computed explicitly and is given by `log(2)`. The integrand for `C2`
is bounded and can be computed using `Arblib.integrate`, using that
the integrand is zero at both `t = 0` and `t = 1`, increasing at the
former and decreasing at the latter.
- **PROVE**: That we have `∫ abs(log(1 / t^2 - 1)) * t dt = log(2)`
- **TODO**: Prove that the integrand for `C2` is zero at both `t = 0`
  and `t = 1`, increasing at the former and decreasing at the latter.
  Give explicit bounds in `t` for when it is increasing respectively
  decreasing.

# Handling `I2`
For `I2` we give a uniform bound of
```
abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π)))
```
in `t` on the interval and use this to simplify the integral. Notice
that `x * (1 - t) / 2π`, `x * (1 + t) / 2π` and `x * t / 2π` all are
less than `x / π`. Using a Taylor expansion of `log(sinc(y))` at `y =
0` and explicitly enclosing the second order remainder term `D` in the
expansion
```
log(sinc(y)) = D * y^2
```
In particular this gives
```
abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π))) <=
    D * x^2 / π^2 * abs((1 - t)^2 / 4 + (1 + t)^2 / 4 - t^2 / 2)) =
    D / 2π^2 * x^2
```
Inserting this into `I2` we have
```
I2 <= D / 2π^2 * x^2 * ∫ t * sqrt(log(1 + inv(x * t))) dt
```
Similarly to for `I1` we use the inequality
```
sqrt(log(1 + inv(x * t))) = sqrt(log(inv(x)) + log(inv(t)) + log(1 + x * t))
    <= sqrt(log(inv(x))) + sqrt(log(inv(t))) + sqrt(log(1 + x * t))
```
giving us
```
I2 <= D / 2π^2 * x^2 * (
    sqrt(log(inv(x))) * ∫ t dt
    + ∫ t * sqrt(log(inv(t))) dt
    + sqrt(log(1 + x)) * ∫ t dt
)
```
We have `∫ t dt = 1 / 2` and
```
∫ t * sqrt(log(inv(t))) dt = sqrt(π / 2) / 4
```
- **PROVE:** That we have `∫ t * sqrt(log(inv(t))) dt = sqrt(π / 2) / 4`
Giving us
```
I2 <= D / 2π^2 * x^2 * (
    sqrt(log(inv(x))) / 2
    + sqrt(π / 2) / 4
    + sqrt(log(1 + x)) / 2
) = D / 4π^2 * x^2 * (
    sqrt(log(inv(x)))
    + sqrt(π / 2) / 2
    + sqrt(log(1 + x))
)
```

# Putting it together
From the above we have
```
I1 <= <= sqrt(log(inv(x))) * C1 + C2 + sqrt(log(1 + x)) * C1

I2 <= D / 4π^2 * x^2 * (sqrt(log(inv(x))) + sqrt(π / 2) / 2 + sqrt(log(1 + x)))
```
This gives us
```
W(x) * I <= inv(sqrt(log(1 + inv(x)))) * (
    C1 / sqrt(log(inv(x)))
    + (C2 + sqrt(log(1 + x)) * C1) / log(inv(x))
    + D / 4π^2 * x^2 * (
        1 / sqrt(log(inv(x)))
        + (sqrt(π / 2) / 2 + sqrt(log(1 + x))) / log(inv(x))
    )
)
```
"""
function T01(u0::BHAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ::Arb = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    ϵ = convert(Arb, ϵ)
    ϵ < 1 || throw(DomainError(ϵ, "must have ϵ < 1"))

    # Expansion for evaluating u0(x) / (x * log(x))
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
        iszero(x) && return zero(x)

        if Arblib.contains_zero(x)
            # Use that inv(log(inv(x)) * sqrt(log(1 + inv(x)))) is
            # monotonically decreasing for 0 < x < 1.
            xᵤ = ubound(Arb, x)

            return Arb((inv(log(inv(xᵤ)) * sqrt(log(1 + inv(xᵤ)))), 0))
        end

        return inv(log(inv(x)) * sqrt(log(1 + inv(x))))
    end

    # Enclosure of second derivative of log(sinc(y)) on [0, ϵ / π].
    # Giving the remainder term in the Taylor expansion.
    D = abs(log(fx_div_x(s -> sinpi(s) / π, ArbSeries((((0, ϵ / π), 1, 0)))))[2])

    # Enclosure of ∫ abs(log(1 / t^2 - 1)) * t dt
    C1 = log(Arb(2))

    # Enclosure of ∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt
    C2 = begin
        integrand_C2(t; analytic) = begin
            if Arblib.contains_zero(t)
                analytic && return Arblib.indeterminate!(zero(t))
                @assert isreal(t)

                # Use that the integrand is increasing close to zero
                # FIXME: Use the correct bounds for when it is increasing
                tᵤ = ubound(Arb, real(t))
                tᵤ < 0.1 || return Arblib.indeterminate!(zero(t))
                return Acb((zero(tᵤ), log(1 / tᵤ^2 - 1) * tᵤ * sqrt(log(inv(tᵤ)))))
            elseif Arblib.overlaps(t, one(t))
                analytic && return Arblib.indeterminate!(zero(t))
                @assert isreal(t)

                # Use that the integrand is increasing (since we
                # are skipping the absolute value in this case)
                # close to one
                # FIXME: Use the correct bounds for when it is increasing
                tₗ = lbound(Arb, real(t))
                tₗ > 0.99 || return Arblib.indeterminate!(zero(t))
                return Acb((log(1 / tₗ^2 - 1) * tₗ * sqrt(log(inv(tₗ))), zero(tₗ)))
            else
                return log(1 / t^2 - 1) *
                       t *
                       Arblib.real_sqrtpos!(zero(t), log(inv(t)), analytic)
            end
        end

        # Integrate from 0 to 1 / sqrt(2) and then from 1 / sqrt(2) to 1
        integral1 = real(
            Arblib.integrate(
                integrand_C2,
                0,
                1 / sqrt(Arb(2)),
                check_analytic = true,
                atol = 1e-10,
            ),
        )
        integral2 = real(
            Arblib.integrate(
                (t; analytic) -> -integrand_C2(t; analytic),
                1 / sqrt(Arb(2)),
                1,
                check_analytic = true,
                atol = 1e-10,
            ),
        )

        integral1 + integral2
    end

    return x::Arb -> begin
        @assert x <= ϵ

        # Enclosure of inv(sqrt(log(1 + inv(x))))
        invsqrtlog1pinvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            lower = zero(x)
            upper = let xᵤ = ubound(Arb, x)
                inv(sqrt(log(1 + inv(xᵤ))))
            end
        else
            inv(sqrt(log(1 + inv(x))))
        end

        # Enclosure of inv(log(inv(x)))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            lower = zero(x)
            upper = let xᵤ = ubound(Arb, x)
                inv(log(inv(xᵤ)))
            end
        else
            inv(log(inv(x)))
        end

        # Enclosure of sqrt(log(1 + x))
        sqrtlog1px = sqrt(Arblib.nonnegative_part!(zero(x), log1p(x)))

        res =
            invsqrtlog1pinvx * (
                C1 * sqrt(invloginvx) +
                (C2 + sqrtlog1px * C1) * invloginvx +
                D / 4Arb(π)^2 *
                x^2 *
                (sqrt(invloginvx) + (sqrt(Arb(π) / 2) / 2 + sqrtlog1px) * invloginvx)
            )

        return factor(x) * res
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

            # Check that the real part of t is strictly between 0 and
            # 1 or return an indeterminate result
            if !(Arblib.ispositive(Arblib.realref(t)) && Arblib.realref(t) < 1)
                Arblib.indeterminate!(res)
                return
            end

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
            opts = Arblib.calc_integrate_opt_struct(0, 10_000, 0, 0, 0),
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
