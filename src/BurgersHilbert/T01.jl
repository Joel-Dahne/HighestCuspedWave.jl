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

It takes as argument `ϵ` and the resulting function is valid for `x <=
ϵ`. The value of `ϵ` has to be less than `1`.

To begin with the factor `x * log(inv(x) / (π * u0(x))` is factored
out from the whole expression and multiplied back in the end. Notice
that this factor is positive and bounded in `x`.

What we are left with computing is
```
W(x) * U1(x)
```
where `W(x) = 1 / (x^2 * log(inv(x)) * sqrt(log(1 + inv(x))))`
```
U1(x) = x^2 * ∫abs(log(sin(x * (1 - t) / 2) + sin(x * (1 + t) / 2) - 2log(sin(x * t / 2)))) * t * sqrt(log(1 + inv(x * t))) dt
```
from `0` to `1`.

Using that
```
log(sin(x * (1 - t) / 2)) = log(x * (1 - t) / 2) + log(sinc(x * (1 - t) / 2π))
log(sin(x * (1 + t) / 2)) = log(x * (1 + t) / 2) + log(sinc(x * (1 + t) / 2π))
log(sin(x * t / 2)) = log(x * t / 2) + log(sinc(x * t / 2π))
```
where we have used the Julia convention that `sinc(x) = sinpi(x) / (π
* x)`. we can split `U1` as
```
U1(x) <= x^2 * (U1_m(x) + U1_r(x))
```
where
```
U1_m(x) = ∫abs(log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2)) * t * sqrt(log(1 + inv(x * t))) dt
U1_r(x) = ∫abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π))) * t * sqrt(log(1 + inv(x * t))) dt
```
We now handle `I1` and `I2` separately.

# Handling `U1_m`
We start by noticing that
```
log(x * (1 - t) / 2) + log(x * (1 + t) / 2) - 2log(x * t / 2) =
    log(1 / t^2 - 1)
```
Furthermore, since `0 <= x <= ϵ < 1` and `0 <= t <= 1` we have
```
sqrt(log(1 + inv(x * t))) = sqrt(log(inv(x)) + log(inv(t)) + log(1 + x * t))
    <= sqrt(log(inv(x))) + sqrt(log(inv(t))) + sqrt(log(1 + x * t))
```
and also `sqrt(log(1 + x * t)) <= sqrt(log(1 + x))`, giving us
```
U1_m(x) <=
    sqrt(log(inv(x))) * ∫ abs(log(1 / t^2 - 1)) * t dt
    + ∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt
    + sqrt(log(1 + x)) * ∫ abs(log(1 / t^2 - 1)) * t dt
```
We have `∫ abs(log(1 / t^2 - 1)) * t dt = log(2)` and if we let `c1 =
∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt` this can be
written as
```
I1 <= sqrt(log(inv(x))) * log(2) + c1 + sqrt(log(1 + x)) * log(2)
```

What remains is to compute an enclosure of
```
c1 = ∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt
```
The integrand is bounded so with some work an enclosure can be
computed using [`Arblib.integrate`](@ref). We can get rid of the
absolute value by splitting the integral at the only zero of `log(1 /
t^2 - 1)` given by `t = inv(sqrt(2))`. At the endpoints the integrand
is zero but not analytic, to compute enclosures we use the following
approach.

After having removed the absolute value we can split the integrand as
```
log(1 / t^2 - 1) * t * sqrt(log(inv(t)))
    = log(1 - t) * t * sqrt(log(inv(t))) +
      log(1 + t) * t * sqrt(log(inv(t))) -
      2log(t) * t * sqrt(log(inv(t)))
```
Let
```
f1(t) = log(1 - t) * sqrt(log(inv(t)))
f2(t) = log(1 + t) * sqrt(log(inv(t)))
f3(t) = log(t) * t * sqrt(log(inv(t)))
```
The integrand is then given by `t * f1(t) + t * f2(t) - 2f3(t)`. We
are thus interested in bounding `f1`, `f2` and `f3` near `t = 0` and
`t = 1`. They are all zero and `t = 0` and `t = 1`, our goal will be
to isolate the parts of the interval ``(0, 1)`` where they are not
monotone. We can then use monotonicity near the endpoints.

For `f3` we can differentiate to get
```
f3'(t) = log(inv(t)) * (2sqrt(log(inv(t)))) * (3 + 2log(t))
```
For ``t ∈ (0, 1)`` this has the unique root `t = exp(-3 / 2)`.

For `f1` and `f2` differentiation gives
```
f1'(t) = inv(2sqrt(log(inv(t)))) * (2log(t) / (1 - t) - log(1 - t) / t)

f2'(t) = -inv(2sqrt(log(inv(t)))) * (2log(t) / (1 + t) + log(1 + t) / t)
```
which is zero iff
```
2log(t) / (1 - t) - log(1 - t) / t = 0
```
or
```
2log(t) / (1 + t) + log(1 + t) / t = 0
```
- **PROVE:** That these two functions have a unique root on ``(0,
  1)``. They are both increasing in `t`, it would be enough to prove
  that.

# Handling `U1_r`
For `U1_r` we give a uniform bound of
```
abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π)))
```
in `t` on the interval and use this to simplify the integral. Notice
that `x * (1 - t) / 2π`, `x * (1 + t) / 2π` and `x * t / 2π` all lie
in the interval ``[0, x / π]``. Using a Taylor expansion of
`log(sinc(y))` at `y = 0` and explicitly enclosing the second order
remainder term `R1` in the expansion we get
```
log(sinc(y)) = R1 * y^2
```
In particular this gives
```
abs(log(sinc(x * (1 - t) / 2π)) + log(sinc(x * (1 + t) / 2π)) - 2log(sinc(x * t / 2π))) <=
    R1 * x^2 / π^2 * abs((1 - t)^2 / 4 + (1 + t)^2 / 4 - t^2 / 2)) =
    R1 / 2π^2 * x^2
```
Inserting this into `U1_r` we have
```
U1_r(x) <= R1 / 2π^2 * x^2 * ∫ t * sqrt(log(1 + inv(x * t))) dt
```
Similarly to for `U1_m` we use the inequality
```
sqrt(log(1 + inv(x * t))) = sqrt(log(inv(x)) + log(inv(t)) + log(1 + x * t))
    <= sqrt(log(inv(x))) + sqrt(log(inv(t))) + sqrt(log(1 + x * t))
```
giving us
```
U1_r(x) <= R1 / 2π^2 * x^2 * (
    sqrt(log(inv(x))) * ∫ t dt
    + ∫ t * sqrt(log(inv(t))) dt
    + sqrt(log(1 + x)) * ∫ t dt
)
```
We have `∫ t dt = 1 / 2` and
```
∫ t * sqrt(log(inv(t))) dt = sqrt(π / 2) / 4
```
Giving us
```
U1_r(x) <= R1 / 2π^2 * x^2 * (
    sqrt(log(inv(x))) / 2
    + sqrt(π / 2) / 4
    + sqrt(log(1 + x)) / 2
) = R1 / 8π^2 * x^2 * (
    2sqrt(log(inv(x)))
    + sqrt(π / 2)
    + 2sqrt(log(1 + x))
)
```

# Putting it together
From the above we have
```
U1_m(x) <= <= sqrt(log(inv(x))) * log(2) + c1 * sqrt(log(1 + x)) * log(2)

U1_r(x) <= R1 / 4π^2 * x^2 * (sqrt(log(inv(x))) + sqrt(π / 2) / 2 + sqrt(log(1 + x)))
```
This gives us
```
W(x) * U1(x) <= inv(sqrt(log(1 + inv(x)))) * (
    log(2) / sqrt(log(inv(x)))
    + (c1 + log(2) * sqrt(log(1 + x))) / log(inv(x))
    + R1 / 8π^2 * x^2 * (
        2 / sqrt(log(inv(x)))
        + (sqrt(π / 2) + 2sqrt(log(1 + x))) / log(inv(x))
    )
)
```
"""
function T01(u0::BHAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ::Arb = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert Arblib.overlaps(u0.w(x), x * sqrt(log(1 + inv(x))))
    end

    # This also checks that ϵ < 1
    inv_u0 = inv_u0_normalised(u0; ϵ)

    # This gives the factor x * log(inv(x)) / (π * u0(x))
    factor(x) = inv_u0(x) / π

    # Enclosure of ∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt
    c1 = begin
        # Unique root of derivative of log(1 - t) * sqrt(log(inv(t)))
        # on (0, 1)
        root1 = begin
            roots, flags = ArbExtras.isolate_roots(Arf(0.1), Arf(0.9)) do t
                2log(t) / (1 - t) - log(1 - t) / t
            end
            @assert only(flags)
            ArbExtras.refine_root(Arb(only(roots))) do t
                2log(t) / (1 - t) - log(1 - t) / t
            end
        end
        # Unique root of derivative of log(1 + t) * sqrt(log(inv(t)))
        # on (0, 1)
        root2 = begin
            roots, flags = ArbExtras.isolate_roots(Arf(0.1), Arf(0.9)) do t
                2log(t) / (1 + t) + log(1 + t) / t
            end
            @assert only(flags)
            ArbExtras.refine_root(Arb(only(roots))) do t
                2log(t) / (1 + t) + log(1 + t) / t
            end
        end
        # Unique root of derivative of log(t) * t *
        # sqrt(log(inv(t))) on (0, 1)
        root3 = exp(Arb(-3 // 2))

        integrand_c1(t; analytic) = begin
            if Arblib.contains_zero(t)
                analytic && return Arblib.indeterminate!(zero(t))
                @assert isreal(t)
                t = real(t)

                # Check that log(t - 1) * sqrt(log(inv(t))) is
                # monotone
                contains(t, root1) && return Arblib.indeterminate!(zero(t))

                # Check that log(t + 1) * sqrt(log(inv(t))) is
                # monotone
                contains(t, root2) && return Arblib.indeterminate!(zero(t))

                # Check that log(t) * t * sqrt(log(inv(t))) is
                # monotone
                contains(t, root3) && return Arblib.indeterminate!(zero(t))

                tᵤ = ubound(Arb, t)

                term1 = Arb((log(1 - tᵤ) * sqrt(log(inv(tᵤ))), 0))
                term2 = Arb((0, log(1 + tᵤ) * sqrt(log(inv(tᵤ)))))
                term3 = Arb((log(tᵤ) * tᵤ * sqrt(log(inv(tᵤ))), 0))

                return t * term1 + t * term2 - 2term3
            elseif Arblib.overlaps(t, one(t))
                analytic && return Arblib.indeterminate!(zero(t))
                @assert isreal(t)
                t = real(t)

                # Check that log(t - 1) * sqrt(log(inv(t))) is
                # monotone
                contains(t, root1) && return Arblib.indeterminate!(zero(t))

                # Check that log(t + 1) * sqrt(log(inv(t))) is
                # monotone
                contains(t, root2) && return Arblib.indeterminate!(zero(t))

                # Check that log(t) * t * sqrt(log(inv(t))) is
                # monotone
                contains(t, root3) && return Arblib.indeterminate!(zero(t))

                tₗ = lbound(Arb, t)

                term1 = Arb((log(1 - tₗ) * sqrt(log(inv(tₗ))), 0))
                term2 = Arb((0, log(1 + tₗ) * sqrt(log(inv(tₗ)))))
                term3 = Arb((log(tₗ) * tₗ * sqrt(log(inv(tₗ))), 0))

                return t * term1 + t * term2 - 2term3
            else
                return log(1 / t^2 - 1) *
                       t *
                       Arblib.real_sqrtpos!(zero(t), log(inv(t)), analytic)
            end
        end

        # Integrate from 0 to 1 / sqrt(2) and then from 1 / sqrt(2) to 1
        integral1 = real(
            Arblib.integrate(
                integrand_c1,
                0,
                1 / sqrt(Arb(2)),
                check_analytic = true,
                atol = 1e-10,
            ),
        )
        integral2 = real(
            Arblib.integrate(
                (t; analytic) -> -integrand_c1(t; analytic),
                1 / sqrt(Arb(2)),
                1,
                check_analytic = true,
                atol = 1e-10,
            ),
        )

        integral1 + integral2
    end

    # Enclosure of the remainder term in the Taylor expansion of
    # log(sinc(y)) on [0, ϵ / π] of degree 1.
    R1 = abs(log(fx_div_x(s -> sinpi(s) / π, ArbSeries((((0, ϵ / π), 1, 0)))))[2])

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

        # Enclosure of W(x) * U1(x)
        WU1 =
            invsqrtlog1pinvx * (
                log(Arb(2)) * sqrt(invloginvx) +
                (c1 + log(Arb(2)) * Arblib.sqrtpos(log1p(x))) * invloginvx +
                R1 / 8Arb(π)^2 *
                x^2 *
                (
                    2sqrt(invloginvx) +
                    (sqrt(Arb(π) / 2) + 2Arblib.sqrtpos(log1p(x))) * invloginvx
                )
            )

        return factor(x) * WU1
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

- **IMPROVE**: The upper bound for when `x` contains `π` could be
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
