"""
    T01(u0::BHKdVAnsatz, ::Ball; δ1, δ2, skip_div_u0)

Returns a function such that `T01(u0, Ball(); δ1, δ2)(x)` computes the
integral \$T_{0,1}\$ from the paper.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T01(
    u0::BHKdVAnsatz,
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
    T01(u0::BHKdVAnsatz, ::Asymptotic)

Returns a function such that `T01(u0, Asymptotic())(x)` computes an
**upper bound** of the integral \$T_{0,1}\$ from the paper using an
evaluation strategy that works asymptotically as `x` goes to `0`.

First of all the change of coordinates `t = y / x` leaves us with
```
x / (π * u0(x) * log(10 + inv(x))) *
    ∫ abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t * log(10 + inv(x * t)) dt
```
with the integration going from `0` to `1`.

To begin with the factor
```
x^-α * log(x) / (π * u0(x) * log(10 + inv(x)))
```
is factored out from the whole expression and multiplied back in the
end. Notice that this factor is bounded in `x`. We can bound `log(x) /
log(10 + inv(x))` by using that it is `-1` at `x = 0` and increasing.
To bound `x^-α / u0(x)` we compute the asymptotic expansion of `u0(x)`
and multiply that by `x^α` and evaluate.
- **TODO:** Deal with the fact that `u0(x) * x^α` blows up as `x ->
    0`. It's not a problem since we divide by it, but it needs to be
    handled.

What we are left with computing is
```
W(x) * I
```
where `W(x) = x^(1 + α) / log(x)` and `I` the same integral as above.

Now consider the expansions
```
clausenc(x * (1 - t), -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * (1 - t)^(-α - 1) + R(x * (1 - t))
clausenc(x * (1 + t), -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * (1 + t)^(-α - 1) + R(x * (1 + t))
clausenc(x * t, -α) = -gamma(1 + α) * sinpi(α / 2) * x^(-α - 1) * t^(-α - 1) + R(x * t)
```
where the error term `R` contains one constant term and `R(x * (1 -
t)) + R(x * (1 + t)) - 2R(x * t)` behaves like `O(x^2)`. Since we are
only looking for an upper bound we can split the absolute value and we
get the two integrals
```
I₁ = abs(gamma(1 + α) * sinpi(α / 2)) * x^(-α - 1) *
    ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(10 + inv(x * t)) dt
I₂ = ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(10 + inv(x * t)) dt
```
satisfying `I <= I₁ + I₂`. Furthermore we can split `log(10 + inv(x *
t))` as
```
log(10 + inv(x * t)) = log((10x * t + 1) / (x * t)) = log(1 + 10x * t) - log(x) - log(t)
```
Which allows us to split the two above integrals into three integrals
each.
```
I₁₁ = -log(x) * abs(gamma(1 + α) * sinpi(α / 2)) * x^(-α - 1) *
    ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
I₁₂ = -abs(gamma(1 + α) * sinpi(α / 2)) * x^(-α - 1) *
    ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
I₁₃ = abs(gamma(1 + α) * sinpi(α / 2)) * x^(-α - 1) *
    ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(1 + 10x * t) dt
I₂₁ = -log(x) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t dt
I₂₂ = -∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(t) dt
I₂₃ = ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(1 + 10x * t) dt
```

Focusing on the first three integrals the first step is to understand
the behaviour of
```
gamma(1 + α) * sinpi(α / 2) * ((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))
```
The factor `sinpi(α / 2)` converges to `-1` and can be enclosed
directly, we can hence focus our attention on
```
gamma(1 + α) * ((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1))
```
As `α -> -1` this converges to
```
-log(1 - t) - log(1 + t) + 2log(t)
```
- **TODO:** Finish this.

For the first integral we get
```
W(x) * I₁₁ = -abs(sinpi(α / 2)) *
    ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
```
As mentioned `abs(sinpi(α / 2))` can be enclosed directly. The
integral doesn't depend on `x` and can be enclosed.
- **TODO:** Enclose the integral
  `∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt`.
Next we have
```
W(x) * I₁₂ = -inv(log(x)) * abs(sinpi(α / 2)) *
    ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
```
where again the integral doesn't depend on `x` and can be enclosed.
The factor `abs(sinpi(α / 2))` can be enclosed as above and
`-inv(log(x))` can be enclosed using that it is zero at `x = 0` and
increasing.
- **TODO:** Enclose the integral
  ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt.
For the third one we have
```
W(x) * I₁₃ = inv(log(x)) * abs(sinpi(α / 2)) *
    ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(1 + 10x * t) dt
```
The factor outside the integral can be dealt with in the same way as
for `W(x) * I₁₂` The factor `log(1 + 10x * t)` inside the integral can
be enclosed using that it is zero for `x = 0` and increasing in both
`x` and `t`. We can factor out the enclosure from the integral and the
integral we are left with is the same as for `W(x) * I₁₁`.

Now for the remaining three terms the first step is to handle the term
`R₁(x) + R₂(x) - 2R₃(x)`. From the expansion of the Clausen functions
we get
```
R₁(x) + R₂(x) - 2R₃(x) = sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
Since we are interested in an upper bound of absolute value we get
```
abs(R₁(x) + R₂(x) - 2R₃(x)) <= x^2 * sum(
    abs(zeta(-α - 2m)) / factorial(2m) * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
) <= x^2 * 2sum(
    abs(zeta(-α - 2m)) * 2^2m / factorial(2m) for m = 1:Inf
)
```
The final sum doesn't depend on `x` and can be bounded, call the bound
`C`.
- **TODO:** Bound the sum
  `2sum(abs(zeta(-α - 2m)) * 2^2m / factorial(2m) for m = 1:Inf)`.

Now for the integrals we get
```
W(x) * I₂₁ = -x^(1 + α) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t dt <=
    -x^(1 + α) * x^2 * C * ∫ t dt =
    -x^(3 + α) * C / 4  <=
    -x^3 * C / 4
```
```
W(x) * I₂₂ = -x^(1 + α) / log(x) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(t) dt <=
    -x^(1 + α) / log(x) * x^2 * C * ∫ t * log(t) dt =
    -x^(3 + α) / log(x) * C * ∫ t * log(t) dt <=
    -x^3 / log(x) * C * ∫ t * log(t) dt =
    x^3 / log(x) * C / 4
```
```
W(x) * I₂₃ = x^(1 + α) / log(x) * ∫ abs(R(x * (1 - t)) + R(x * (1 + t)) - 2R(x * t)) * t * log(1 + 10x * t) dt <=
    x^(1 + α) / log(x) * x^2 * C * ∫ t * log(1 + 10x * t) dt <=
    x^(3 + α) / log(x) * C * log(1 + 10x) / 2 <=
    x^3 / log(x) * C * log(1 + 10x) / 2
```
where in the last step we have used that `log(1 + 10x * t) <= log(1 + 10x)`.
- **TODO:** Take a look at the sign for these three integrals. We
  should probably just consider the absolute value instead.

"""
function T01(u0::BHKdVAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ = Arb(2e-1))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), x * log(10 + inv(x)))
    end

    ϵ = convert(Arb, ϵ)
    @assert ϵ < 0.5

    u0_expansion = u0(ϵ, AsymptoticExpansion())
    u0_expansion_mul_xα = empty(u0_expansion)
    for ((p, q, i, j, k, l, m), value) in u0_expansion
        u0_expansion_mul_xα[(p, q, i + 1, j, k, l, m)] = value
    end

    # Enclosure of the factor
    # x^-α * log(x) / (π * u0(x) * log(10 + inv(x)))
    factor(x) = begin
        # Enclosure of log(x) / log(10 + inv(x))
        if iszero(x)
            log_factor = -one(x)
        elseif Arblib.contains_zero(x)
            # Use monotonicity, it is -1 at x = 0 and increasing
            xᵤ = ubound(Arb, x)
            log_factor = Arb((-1, log(xᵤ) / log(10 + inv(xᵤ))))
        else
            log_factor = log(x) / log(10 + inv(x))
        end

        # Enclosure of u0(x) * x^α
        if non_asymptotic_u0
            u0_factor = u0(x) * x^Arb((-1, -1 + u0.ϵ))
        else
            u0_factor = eval_expansion(u0, u0_expansion_mul_xα, x)
        end

        log_factor / (π * u0_factor)
    end

    return x -> begin
        x = convert(Arb, x)
        @assert x <= ϵ

        # Enclosure of abs(sinpi(α / 2))
        factor_I₁ = let α = Arb((-1, -1 + u0.ϵ))
            abs(sinpi(α / 2))
        end

        # Enclosure of inv(log(x))
        invlogx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            xᵤ = ubound(Arb, x)
            Arb((inv(log(xᵤ)), 0))
        else
            inv(log(x))
        end

        WxI₁₁ = begin
            # Enclosure of
            # ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t dt
            # FIXME: Rigorously enclose this
            integral_I₁₁ = Arb("0.6931477765876257")

            -factor_I₁ * integral_I₁₁
        end

        WxI₁₂ = begin
            if iszero(x)
                zero(x)
            else
                # Enclosure of
                # ∫ gamma(1 + α) * abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t * log(t) dt
                # FIXME: Rigorously enclose this
                integral_I₁₂ = Arb("-0.46668696508451124")

                -invlogx * factor_I₁ * integral_I₁₂
            end
        end

        WxI₁₃ = begin
            if iszero(x)
                zero(x)
            else
                # The integral is the same as for WxI₁₁
                integral_I₁₃ = integral_I₁₁
                # Enclosure of log(1 + 10x * t) for t ∈ [0, 1]
                factor_I₁₂ = Arb((0, log1p(10x)))

                invlogx * factor_I₁ * factor_I₁₂ * integral_I₁₃
            end
        end

        WxI₁ = WxI₁₁ + WxI₁₂ + WxI₁₃

        # Bound of 2sum((zeta(-α - 2m)) * 2^2m / factorial(2m) for m = 1:Inf)
        # FIXME: Rigorously bound this
        C = Arb("0.1725034263261214970439543405065145523")

        # Enclosure of inv(abs(gamma(1 + α) * sinpi(α / 2)))

        WxI₂₁ = x^3 * C / 2

        WxI₂₂ = x^3 * invlogx * C / 4

        WxI₂₃ = x^3 * invlogx * C * log1p(10x) / 2

        WxI₂ = WxI₂₁ + WxI₂₂ + WxI₂₃

        WxI = WxI₁ + WxI₂

        res = factor(x) * WxI

        return res
    end
end

"""
    T011(u0::BHKdVAnsatz; δ0)

Computes the integral \$T_{0,1,1}\$ from the paper.

It uses the fact that the integrand is strictly increasing on the
interval `[0, 0.05]` for every value of `x` and 0 at `x = 0`. This
allows us to enclose the integrand on the interval which then easily
gives an enclosure of the integral by multiplying with the size of the
interval.

- **FIXME:** Currently this uses `α = -1`, this should be changed to
  using the interval `[1 - u0.ϵ, 1]` for `α` once [`clausenc`](@ref)
  supports it.
- **PROVE**: That the integrand indeed is increasing on the said
  interval.
"""
function T011(u0::BHKdVAnsatz, ::Ball = Ball(); δ0::Arb = Arb(1e-5), skip_div_u0 = false)
    δ0 < 0.05 || Throw(ArgumentError("δ0 must be less than 0.05, got $δ0"))
    return x -> begin
        x = convert(Arb, x)

        # FIXME: Use the interval for s once this is supported by clausenc
        #α = Arb((-1, -1 + u0.ϵ))
        α = Arb(-1)
        integrand(t) =
            abs(
                clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) -
                2clausenc(x * t, -α),
            ) *
            t *
            log(10 + inv(x * t))

        integral = δ0 * Arb((0, integrand(δ0)))

        res = integral * x / (π * log(10 + inv(x)))
        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T012(u0::BHKdVAnsatz; δ0, δ1)

Returns a function such that `T012(u0; δ0, δ1)(x)` computes the integral
\$T_{0,1,2}\$ from the paper.

**FIXME:** This currently uses `α = -1`. We need to either compute
  with `α = [-1, -1 + u0.ϵ]` or bound the error in some other way. One
  issue is that `clausenc` doesn't allow balls overlapping integers.
  Another one is that we get very bad bounds if `α` is close to by not
  equal to `-1`. This could possibly be handled by switching to an
  integration method that only uses real values, where we have better
  bounds.
"""
function T012(
    u0::BHKdVAnsatz,
    ::Ball = Ball();
    δ0::Arb = Arb(1e-5),
    δ1::Arb = Arb(1e-5),
    skip_div_u0 = false,
)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * log(10 + inv(x)))
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
            #weight = t * log(10 + inv(tx))
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

            #s = Arb(1 - u0.ϵ)
            #s = one(Arb)
            #if isreal(t)
            #    t_real = Arblib.realref(t)
            #    Arblib.set!(res, clausenc(x * (1 - t_real), s) + clausenc(x * (1 + t_real), s) - 2clausenc(x * t_real, s))
            #else
            #    Arblib.set!(res, clausenc(x * (1 - t), s) + clausenc(x * (1 + t), s) - 2clausenc(x * t, s))
            #end

            Arblib.real_abs!(res, res, analytic)

            # tmp = t * log(10 + inv(t * x))
            Arblib.inv!(tmp, tx)
            Arblib.add!(tmp, tmp, 10)
            Arblib.log!(tmp, tmp)
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
            rtol = 1e-20,
            atol = 1e-20,
            warn_on_no_convergence = false,
            opts = Arblib.calc_integrate_opt_struct(0, 10000, 0, 1, 1),
        )

        @assert !isfinite(res) || isreal(res)
        res = real(res)

        res = res * x / (π * log(10 + inv(x)))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T013(u0::BHKdVAnsatz; δ1)

Computes the integral \$T_{0,1,3}\$ from the paper.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval. We can also notice that the value inside the absolute value
is negative so we can remove the absolute value by putting a minus
sign, which we can bake in to the weight factor.

We are left with integrating the three Clausen terms
1. `clausenc(x * (1 - t), -α)`
2. `clausenc(x * (1 + t), -α)`
3. `2clausenc(x * t, -α)`
We have that the primitive functions for the three terms are given by
1. `-clausens(x * (1 - t), 1 - α) / x`
2. `clausens(x * (1 + t), 1 - α) / x`
3. `2clausens(x * t, 1 - α) / x`
Hence the integral from `1 - δ1` to `1` is
```
inv(x) * (
    (-clausens(0, 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α)) -
    (-clausens(x * δ1, 1 - α) + clausens(x * (2 - δ1), 1 - α) - 2clausens(x * (1 - δ1), 1 - α))
)
```
Finally the multiplication by `inv(x)` can be cancelled by the
multiplication by `x` that is outside of the integral.

- **FIXME:** Currently this uses `α = -1`, this should be changed to
  using the interval `[1 - u0.ϵ, 1]` for `α` once [`clausenc`](@ref)
  supports it.
- **TODO:** This doesn't work for `x` overlapping `π` because
  `clausens` doesn't support evaluation around `2π`.
"""
function T013(u0::BHKdVAnsatz, ::Ball = Ball(); δ1::Arb = Arb(1e-5), skip_div_u0 = false)
    return x -> begin
        x = convert(Arb, x)

        weight_factor = let t = Arb((1 - δ1, 1))
            -t * log(10 + inv(x * t))
        end

        # FIXME: Use the interval for α once this is supported by clausenc
        #α = Arb((-1, -1 + u0.ϵ))
        α = Arb(-1)
        integral =
            weight_factor * (
                (-clausens(zero(Arb), 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α)) - (
                    -clausens(x * δ1, 1 - α) + clausens(x * (2 - δ1), 1 - α) -
                    2clausens(x * (1 - δ1), 1 - α)
                )
            )

        res = integral / (π * log(10 + inv(x)))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
