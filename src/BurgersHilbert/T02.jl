"""
    T02(u0::BHAnsatz; δ2, skip_div_u0)

Return a function for computing an enclosure of
```
inv(u0(x)) * U2(x) / (π * u0.w(x))
```

It computes the `inv(u0(x))` separately and then splits the rest into
the two terms
```
U21(x) / (π * u0.w(x))
```
and
```
U22(x) / (π * u0.w(x))
```
computed by [`T021`](@ref) and [`T022`](@ref) respectively.

# Arguments
- `δ2::Arb = 1e-5`: Determines the split intervals for `U21(x)` and
  `U22(x)`.
- `skip_div_u0::Bool = false`: If true it skips the factor
  `inv(u0(x))`.
"""
function T02(u0::BHAnsatz, evaltype::Ball; δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    f = T021(u0, evaltype, skip_div_u0 = true)
    g = T022(u0, evaltype, skip_div_u0 = true)

    return x::Arb -> begin
        # If x is wide it is beneficial to take a larger δ2
        δ = max(δ2, 8radius(Arb, x))
        a = ubound(Arb, x + δ)

        if a < π
            res = f(x, a) + g(x, a)
        else
            res = f(x, Arb(π))
        end

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

It takes as argument `ϵ` and the resulting function is valid for `x <=
ϵ`. The value of `ϵ` has to be less than `1 / 2`.

To begin with the factor `x * log(inv(x)) / (π * u0(x))` is factored
out from the whole expression and multiplied back in the end. Notice
that this factor is positive and bounded in `x`.

What we are left with computing is
```
W(x) * U2(x)
```
where `W(x) = 1 / (x^2 * log(inv(x)) * sqrt(log(1 + inv(x))))`
```
U2(x) = ∫-(log(sin((y - x) / 2)) + log(sin((y + x) / 2)) - 2log(sin(y / 2))) * y * sqrt(log(1 + inv(y))) dy
```
for `x` to `π`. Using that
```
log(sin((y - x) / 2)) = log((y - x) / 2) + log(sinc((y - x) / 2π))
log(sin((y + x) / 2)) = log((y + x) / 2) + log(sinc((y + x) / 2π))
log(sin(y / 2)) = log(y / 2) + log(sinc(y / 2π))
```
where we have used the Julia convention that `sinc(x) = sinpi(x) / (π
* x)`. we can split `U2` as
```
U2(x) <= U2_m(x) + U2_r(x)
```
where
```
U2_m(x) = ∫ -(log((y - x) / 2) + log((y + x) / 2) - 2log(y / 2)) * y * sqrt(log(1 + inv(y))) dy
U2_r(x) = ∫ -(log(sinc((y - x) / 2π)) + log(sinc((y + x) / 2π)) - 2log(sinc(y / 2π))) * y * sqrt(log(1 + inv(y))) dy
```
We now handle `U2_m` and `U2_r` separately.

# Handling `U2_m`
We start by noticing that
```
log((y - x) / 2) + log((y + x) / 2) - 2log(y / 2) = log(1 - (x / y)^2)
```
```
U2_m(x) = ∫ -log(1 - (x / y)^2) * y * sqrt(log(1 + inv(y))) dy
```
Due to the occurrence of `x / y` it is naturaly to switch to `t = y /
x`, giving us
```
U2_m = x^2 * ∫ -log(1 - 1 / t^2) * t * sqrt(log(1 + inv(x * t))) dt
```
We now split the interval of integration, `[1, π / x]`, into three
parts, `[1, 2], [2, inv(x)], [inv(x), π / x]`, and treat each interval
separately. Let `U2_m1, U2_m2, U2_m3` correspond to the integrals on
the separate intervals.

## Handling `U2_m1`
We begin by noting that `sqrt(log(1 + inv(x * t)))` is decreasing in
`t`, an upper bound is thus given by its value at `t = 1` where it is
`sqrt(log(1 + inv(x)))`. Hence
```
U2_m1 <= x^2 * sqrt(log(1 + inv(x))) * ∫ -log(1 - 1 / t^2) * t dt
```
where the integral can now be explicitly computed to be `log(16 / (3 *
sqrt(3)))`, giving us
```
U2_m1 <= x^2 * sqrt(log(1 + inv(x))) * log(16 / (3 * sqrt(3)))
```

## Handling `U2_m2`
Taylor expanding `-log(1 - y)` gives us
```
-log(1 - y) = y + R2 * y^2
```
where `R2` is the remainder term. If we take the remainder term on the
interval ``[0, 1 / 4]`` we get
```
-log(1 - 1 / t^2) = 1 / t^2 + R2 / t^4,
```
for ``t ∈ [2, ∞]``. This together with
```
sqrt(log(1 + inv(x * t))) =
    sqrt(log(1 + x * t) + log(inv(x * t))) <=
    sqrt(log(1 + x * t)) + sqrt(log(inv(x * t)))
```
and `sqrt(log(1 + x * t)) < sqrt(log(2))` gives us
```
U2_m2(x) <= x^2 * (
    ∫ sqrt(log(inv(x * t))) / t dt
    + R2 * ∫ sqrt(log(inv(x * t))) / t^3 dt
    + sqrt(log(2)) ∫ 1 / t dt
    + R2 * sqrt(log(2)) * ∫ 1 / t^3 dt
)
```
For the integrals we have
```
∫ sqrt(log(inv(x * t))) / t dt = 2 // 3 * log(inv(2x))^(3 // 2)

∫ sqrt(log(inv(x * t))) / t^3 dt = (sqrt(log(inv(2x))) - sqrt(2π) * x^2 * erfi(sqrt(log(inv(4x^2))))) / 8

∫ 1 / t dt = log(inv(x)) - log(2)

∫ 1 / t^3 dt = (1 - 4x^2) / 8
```
To avoid having to use the `erfi` function we use [DLMF
7.8.7]([https://dlmf.nist.gov/7.8.E7) to get
```
x^2 * erfi(sqrt(log(inv(4x^2))))
= 2 / sqrt(π) * x^2 * ∫_0^sqrt(log(inv(4x^2))) exp(t^2) dt
<= 2 / sqrt(π) * x^2 (exp(sqrt(log(inv(4x^2)))^2) - 1) / sqrt(log(inv(4x^2)))
= 2 / sqrt(π) * x^2 (inv(4x^2) - 1) / sqrt(2log(inv(2x)))
= 1 / 2sqrt(2π) * (1 - 4x^2) / sqrt(log(inv(2x)))
```
This gives us the following upper bound for `U2_m2`
```
U2_m2 <=
    2 / 3 * log(inv(2x))^(3 // 2)
    + R2 * sqrt(log(inv(2x))) / 8
    - R2 * (1 - 4x^2) / 16sqrt(log(inv(2x)))
    + sqrt(log(2)) * log(inv(x))
    - log(2)^(3 // 2)
    + R2 * sqrt(log(2)) * (1 - 4x^2) / 8
```

## Handling `U2_m3`
We use that `sqrt(log(1 + inv(x * t))) < sqrt(log(2))` for `1 / x < t
< π / x`, this gives us
```
U2_m3(x) <= sqrt(log(2)) * x^2 * ∫ -log(1 - 1 / t^2) * t dt
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
U3_m3(x) <= sqrt(log(2)) * x^2 * (log((π^2 - x^2) / (1 - x^2)) + log(1 - x^2) / 2x^2 - π^2 * log(1 - x^2 / π^2))
```

## Putting `U2_m1`, `U2_m2` and `U2_m3` together
With all of this we get the following upper bound for `U2_m`
```
U2_m(x) <= x^2 * (
    sqrt(log(1 + inv(x))) * log(16 / (3 * sqrt(3)))
    + (
        2 / 3 * log(inv(2x))^(3 // 2)
        + R2 * sqrt(log(inv(2x))) / 8
        - R2 * (1 - 4x^2) / 16sqrt(log(inv(2x)))
        + sqrt(log(2)) * log(inv(x))
        - log(2)^(3 // 2)
        + R2 * sqrt(log(2)) * (1 - 4x^2) / 8
    ) + sqrt(log(2)) * (
        log((π^2 - x^2) / (1 - x^2))
        + log(1 - x^2) / 2x^2
        - π^2 * log(1 - x^2 / π^2) / 2x^2
    )
)
```
The terms `log(1 - x^2) / 2x^2` and `log(1 - x^2 / π^2) / 2x^2` can be
bounded using [`fx_div_x`](@ref).
To be able to bound this after division by `log(inv(x)) * sqrt(log(1 +
inv(x)))` there are some terms we have to take care of, all of them
coming from `U2_m2`.
1. For the term `sqrt(log(2)) * log(inv(x))` we can directly cancel
   `log(inv(x))`.
2. For the term `2 / 3 * log(inv(2x))^(3 // 2)` we use that it is
  decreasing in `x` after division by `log(inv(x)) * sqrt(log(1 +
  inv(x)))` and tends to `2 / 3` at `x = 0`.
3. For the term `- R2 * sqrt(log(inv(2x))) / 8` we use that
  `sqrt(log(inv(2x))) / sqrt(log(1 + inv(x)))` is monotone in `x`.

# Handling `U2_r`
Recall that
```
U2_r = ∫ -(log(sinc((y - x) / 2π)) + log(sinc((x + y) / 2π)) - 2log(sinc(y / 2π))) * y * sqrt(log(1 + inv(y))) dt
```
from `x` to `π`. Using that
```
-(log(sinc((y - x) / 2π)) + log(sinc((x + y) / 2π)) - 2log(sinc(y / 2π))) / x^2
```
is bounded in `y` on `[x, π]` uniformly in `x` we can upper bound of
it. If we let `D1` denote this upper bound we have
```
U2_r(x) <= D1 * ∫ y * sqrt(log(1 + inv(y))) dy
```
The integral is taken from `x` to `π`, we can get an upper bound by
integrating from `0`. If we let
```
c2 = ∫_0^π y * sqrt(log(1 + inv(y))) dy
```
we get
```
U2_r(x) <= D1 * c2
```
The integral can be enclosed directly using [`Arblib.integrate`](@ref)
and that the integrand is zero at `y = 0` and increasing in `y` to
handle `y` close to zero.

To compute `D1` we use that the function is increasing in `y` on ``[x,
π]`` and it is hence enough to evaluate it `y = π`, where it is
```
-(log(sinc((x / π - 1) / 2)) + log(sinc((x / π + 1) / 2)) - 2log(sinc(1 / 2))) / x^2
```
This expression has a removable singularities that we handle using
[`fx_div_x`](@ref).
"""
function T02(u0::BHAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ::Arb = Arb(2e-1))
    0 < ϵ < 1 // 2 || throw(DomainError(ϵ, "must have 0 < ϵ < 1 / 2"))

    inv_u0 = inv_u0_normalised(u0; ϵ)

    # This gives the factor x * log(inv(x)) / (π * u0(x))
    factor(x) = inv_u0(x) / π

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

    # Bound the constant R2 using ArbSeries. We compute the error term
    # for -log(1 - y) for y in [0, 1 // 4]
    R2 = -log(1 - ArbSeries((Arb((0, 1 // 4)), 1, 0)))[2]

    D1 =
        -fx_div_x(Arb((0, ϵ)), 2, force = true, extra_degree = 2) do x
            log(sinc((x / π - 1) / 2)) + log(sinc((x / π + 1) / 2)) -
            2log(sinc(Arb(1 // 2)))
        end

    # Enclosure of ∫y * sqrt(log(1 + inv(y))) dy from 0 to π
    c2 = begin
        integrand_c2(y; analytic) = begin
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

        real(Arblib.integrate(integrand_c2, 0, π, check_analytic = true))
    end

    return x::Arb -> begin
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

        # Enclosure of W(x) * U2_m(x)
        WU2_m = let
            # Enclosure of W(x) * U2_m1(x)
            WU2_m1 = log(16 / (3sqrt(Arb(3)))) * invlog

            # Enclosure of W(x) * U2_m2(x)
            WU2_m2 = let
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

                # Enclosure of R2 * sqrt(log(inv(2x))) * weight_factor(x) / 8
                term2 =
                    R2 * invlog / 8 * if Arblib.contains_zero(x)
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
                        -R2 * (1 - 4x^2) * invsqrtloginv2x / 16 - log(Arb(2))^(3 // 2) +
                        R2 * sqrt(log(Arb(2))) * (1 - 4x^2) / 8
                    ) * weight_factor(x)

                (term1 + term2 + term3) + remaining
            end

            # Enclosure of U2_m3(x) / x^2
            U2_m3_div_x2 =
                sqrt(log(Arb(2))) * if Arblib.contains_zero(x)
                    log((Arb(π)^2 - x^2) / (1 - x^2)) / 2 +
                    fx_div_x(y -> log(1 - y^2), x, 2) / 2 -
                    Arb(π)^2 * fx_div_x(y -> log(1 - y^2 / Arb(π)^2), x, 2) / 2
                else
                    log((Arb(π)^2 - x^2) / (1 - x^2)) / 2 + log1p(-x^2) / 2x^2 -
                    Arb(π)^2 * log1p(-x^2 / Arb(π)^2) / 2x^2
                end

            WU2_m1 + WU2_m2 + U2_m3_div_x2 * weight_factor(x)
        end

        # Enclosure of W(x) * U2_r(x)
        WU2_r = D1 * c2 * weight_factor(x)

        return factor(x) * (WU2_m + WU2_r)
    end
end

"""
    T021(u0::BHAnsatz{Arb}; skip_div_u0)

Return a function taking the arguments `(x, a)` for computing an
enclosure of
```
inv(u0(x)) * U21(x) / (π * u0.w(x))
```
where
```
U21(x) = ∫ -log(sin((y - x) / 2) * sin((y + x) / 2) / sin(y / 2)^2) * y * sqrt(log(1 + inv(y))) dy
```
on the interval ``[x, a]``.

# Arguments
- `skip_div_u0::Bool = false`: If true it skips the factor
  `inv(u0(x))`.

# Implementation
To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval.

Using that `clausenc(x, 1) = -log(2sin(abs(x) / 2))` we are left with
```
∫ clausenc(y - x, 1) + clausenc(y + x, 1) - 2clausenc(y, 1) dy
```
on the interval ``[x, a]``. A primitive function of the integral can
be determined to be
```
clausens(y - x, 2) + clausens(y + x, 2) - 2clausens(y, 2)
```
This gives us the integral
```
(clausens(a - x, 2) + clausens(a + x, 2) - 2clausens(a, 2))
- (clausens(0, 2) + clausens(2x, 2) - 2clausens(x, 2))
```
"""
function T021(u0::BHAnsatz, ::Ball = Ball(); skip_div_u0 = false)
    return (x::Arb, a::Arb) -> begin
        K2 = let y = Arb((x, a))
            y * sqrt(log(1 + inv(y)))
        end

        integral = ArbExtras.enclosure_series(x, degree = 4) do y
            (
                (clausens(a - y, 2) + clausens(a + y, 2) - 2clausens(a, 2)) -
                (clausens(Arb(0), 2) + clausens(2y, 2) - 2clausens(y, 2))
            )
        end

        res = K2 * integral / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T022(u0::BHAnsatz{Arb}; skip_div_u0)

Return a function taking the arguments `(x, a)` for computing an
enclosure of
```
inv(u0(x)) * U22(x) / (π * u0.w(x))
```
where
```
U22(x) = ∫ -log(sin((y - x) / 2) * sin((y + x) / 2) / sin(y / 2)^2) * y * sqrt(log(1 + inv(y))) dy
```
on the interval ``[a, π]``.

This is done by numerically integrating with
[`Arblib.integrate`](@ref). To not give problems with integration `a`
should be a thin ball.

# Arguments
- `skip_div_u0::Bool = false`: If true it skips the factor
  `inv(u0(x))`.
"""
function T022(u0::BHAnsatz, ::Ball = Ball(); skip_div_u0 = false)
    return (x::Arb, a::Arb) -> begin
        # Allocate space to store temporary variables in integration
        tmp = Acb()

        integrand!(res, y; analytic::Bool) = begin
            # The code below is an inplace version of the following code
            #res = -log(sin((y - x) / 2) * sin((x + y) / 2) / sin(y / 2)^2)
            #weight = y * Arblib.sqrt_analytic!(zero(y), log(1 + inv(y)), analytic)
            #return res * weight

            # res = sin((y - x) / 2)
            Arblib.sub!(tmp, y, x)
            Arblib.mul_2exp!(tmp, tmp, -1)
            Arblib.sin!(res, tmp)

            # res *= sin((x + y) / 2)
            Arblib.add!(tmp, y, x)
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

            # tmp = y * sqrt(log(1 + inv(y)))
            Arblib.inv!(tmp, y)
            Arblib.add!(tmp, tmp, 1)
            Arblib.log!(tmp, tmp)
            Arblib.sqrt_analytic!(tmp, tmp, analytic)
            Arblib.mul!(tmp, tmp, y)

            Arblib.mul!(res, res, tmp)

            return
        end

        U22 = real(
            Arblib.integrate!(
                integrand!,
                zero(Acb),
                a,
                π,
                check_analytic = true,
                rtol = 1e-10,
                atol = 1e-10,
                warn_on_no_convergence = false,
                opts = Arblib.calc_integrate_opt_struct(0, 10_000, 0, 0, 0),
            ),
        )

        res = U22 / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
