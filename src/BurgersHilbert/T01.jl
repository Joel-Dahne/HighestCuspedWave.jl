"""
    _integrand_compute_root(u0::BHAnsatz{Arb}, x::Arb)

Compute the unique root of
```
log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(x * t / 2)^2)
```
in `t` on the interval `[0, 1]`. It assumes that `0 <= x <= π`.

For wide values of `x` it uses that the root is decreasing in `x` to
only have to evaluate at the endpoints.

If the lower bound of `x` is zero or close to zero (smaller than
`eps(x)`) it computes an upper bound of the root by considering the
limiting case as `x` goes to zero. The limiting root is then given by
`inv(sqrt(2))`, which is the root of
```
log(inv(t^2) - 1)
```

If the upper bound of `x` is close to zero, smaller than `eps(x)`, we
compute the root at `eps(x)` and use that as a lower bound. This
avoids computing with very small values of `x`.
"""
function _integrand_compute_root(u0::BHAnsatz{Arb}, x::Arb)
    # Root for x = 0
    root_zero = inv(sqrt(Arb(2)))

    compute_root(x) =
        let
            # Function we are computing the root of
            f = t -> log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(x * t / 2)^2)

            # The root is lower bounded by 1 / 2, take a value
            # slightly smaller so that we can still isolate it even if
            # it touches 1 / 2.
            root_lower = Arf(0.5) - sqrt(eps(Arf))

            # The root is upper bounded by the root at zero
            root_upper = ubound(root_zero)

            # Improve the enclosure of the root
            roots, flags = ArbExtras.isolate_roots(f, root_lower, root_upper)
            if length(flags) == 1 && flags[1]
                # Refine the unique root
                root = ArbExtras.refine_root(f, Arb(only(roots)))
            else
                root = Arb((roots[1][1], roots[end][2]))
            end

            return root
        end

    xₗ, xᵤ = getinterval(Arb, x)
    xᵤ = min(Arb(π), xᵤ) # We assume that xᵤ <= π
    ϵ = eps(x)

    if iszero(x)
        root = root_zero
    elseif !iswide(x)
        root = compute_root(x) # In this case x never overlaps zero
    elseif xᵤ < ϵ
        root = Arb((compute_root(ϵ), root_zero))
    elseif xₗ < eps(Arb)
        root = Arb((compute_root(xᵤ), root_zero))
    else
        root = Arb((compute_root(xᵤ), compute_root(xₗ)))
    end

    return root
end

"""
    T01(u0::BHAnsatz, ::Ball; δ1, skip_div_u0)

Return a function for computing an enclosure of
```
inv(u0(x)) * U1(x) / (π * x* sqrt(log(1 + inv(x))))
```

It computes the `inv(u0(x))` separately and then splits the rest into
the two terms
```
x * (U11(x) + U121(x)) / (π * sqrt(log(1 + inv(x))))
```
and
```
x * U122(x) / (π * sqrt(log(1 + inv(x))))
```
computed by [`T011`](@ref) and [`T012`](@ref) respectively.

# Arguments
- `δ1::Arb = 1e-5`: Determines the interval of integration for
  `U121(x)` and `U122(x)`.
- `skip_div_u0::Bool = false`: If true it skips the factor
  `inv(u0(x))`.
"""
function T01(u0::BHAnsatz{Arb}, evaltype::Ball; δ1::Arb = Arb(1e-5), skip_div_u0 = false)
    f = T011(u0, evaltype, skip_div_u0 = true; δ1)
    g = T012(u0, evaltype, skip_div_u0 = true; δ1)

    if skip_div_u0
        return x::Arb -> f(x) + g(x)
    else
        return x::Arb -> (f(x) + g(x)) / u0(x)
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
f3'(t) = sqrt(log(inv(t))) / 2 * (3 + 2log(t))
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
These functions can be checked to be monotone on ``[0, 1]`` (for
details see paper) with limits ``±∞`` at `t = 0` and `t = 1`, they
hence have a unique root which we can isolate. We can then use
monotonicity of `f1` and `f3` as long as we avoid these roots.

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
    R1 * x^2 / π^2 * abs((1 - t)^2 / 4 + (1 + t)^2 / 4 + t^2 / 2)) =
    R1 * x^2 / π^2 * abs(t^2 + 1 / 2) <=
    3R1 / 2π^2 * x^2
```
Inserting this into `U1_r` we have
```
U1_r(x) <= 3R1 / 2π^2 * x^2 * ∫ t * sqrt(log(1 + inv(x * t))) dt
```
Similarly to for `U1_m` we use the inequality
```
sqrt(log(1 + inv(x * t))) = sqrt(log(inv(x)) + log(inv(t)) + log(1 + x * t))
    <= sqrt(log(inv(x))) + sqrt(log(inv(t))) + sqrt(log(1 + x * t))
```
giving us
```
U1_r(x) <= 3R1 / 2π^2 * x^2 * (
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
U1_r(x) <= 3R1 / 2π^2 * x^2 * (
    sqrt(log(inv(x))) / 2
    + sqrt(π / 2) / 4
    + sqrt(log(1 + x)) / 2
) = 3R1 / 8π^2 * x^2 * (
    2sqrt(log(inv(x)))
    + sqrt(π / 2)
    + 2sqrt(log(1 + x))
)
```

# Putting it together
From the above we have
```
U1_m(x) <= <= sqrt(log(inv(x))) * log(2) + c1 * sqrt(log(1 + x)) * log(2)

U1_r(x) <= 3R1 / 4π^2 * x^2 * (sqrt(log(inv(x))) + sqrt(π / 2) / 2 + sqrt(log(1 + x)))
```
This gives us
```
W(x) * U1(x) <= inv(sqrt(log(1 + inv(x)))) * (
    log(2) / sqrt(log(inv(x)))
    + (c1 + log(2) * sqrt(log(1 + x))) / log(inv(x))
    + 3R1 / 8π^2 * x^2 * (
        2 / sqrt(log(inv(x)))
        + (sqrt(π / 2) + 2sqrt(log(1 + x))) / log(inv(x))
    )
)
```
"""
function T01(u0::BHAnsatz, ::Asymptotic; non_asymptotic_u0 = false, ϵ::Arb = Arb(2e-1))
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

                tᵤ = ubound(Arb, t)

                # Check that log(t - 1) * sqrt(log(inv(t))) is
                # monotone
                tᵤ < root1 || return Arblib.indeterminate!(zero(t))

                # Check that log(t + 1) * sqrt(log(inv(t))) is
                # monotone
                tᵤ < root2 || return Arblib.indeterminate!(zero(t))

                # Check that log(t) * t * sqrt(log(inv(t))) is
                # monotone
                tᵤ < root3 || return Arblib.indeterminate!(zero(t))

                term1 = Arb((log(1 - tᵤ) * sqrt(log(inv(tᵤ))), 0))
                term2 = Arb((0, log(1 + tᵤ) * sqrt(log(inv(tᵤ)))))
                term3 = Arb((log(tᵤ) * tᵤ * sqrt(log(inv(tᵤ))), 0))

                return t * term1 + t * term2 - 2term3
            elseif Arblib.overlaps(t, one(t))
                analytic && return Arblib.indeterminate!(zero(t))
                @assert isreal(t)
                t = real(t)

                tₗ = lbound(Arb, t)

                # Check that log(t - 1) * sqrt(log(inv(t))) is
                # monotone
                tₗ > root1 || return Arblib.indeterminate!(zero(t))

                # Check that log(t + 1) * sqrt(log(inv(t))) is
                # monotone
                tₗ > root2 || return Arblib.indeterminate!(zero(t))

                # Check that log(t) * t * sqrt(log(inv(t))) is
                # monotone
                tₗ > root3 || return Arblib.indeterminate!(zero(t))

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
                3R1 / 8Arb(π)^2 *
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
    T011(u0::BHAnsatz; δ1, skip_div_u0)

Return a function for computing an enclosure of
```
inv(u0(x)) * x * (U11(x) + U121(x)) / (π * sqrt(log(1 + inv(x))))
```
where
```
U11(x) = ∫ log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(x * t / 2)^2) * t * sqrt(log(1 + inv(x * t))) dt
U121(x) = ∫ -log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(x * t / 2)^2) * t * sqrt(log(1 + inv(x * t))) dt
```
and `U11(x)` is integrated on ``[0, r_x]`` and `U121(x)` on ``[r_x, 1
- δ1]``. Here `r_x` is the unique root of the integrand on the
interval ``[0, 1]``, as computed by [`_integrand_compute_root`](@ref).

# Arguments
- `δ1::Arb = 1e-5`: Determines the interval of integration for `U121(x)`
- `skip_div_u0::Bool = false`: If true it skips the factor
  `inv(u0(x))`.

# Implementation
This is done by directly computing the integrals using
[`Arblib.integrate`](@ref).

The only thing to take care of is that the integrand is non-analytic
at `t = 0`. In that case we make use of monotonicity to enclose the
integrand. We split it as
```
log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2)) * t * sqrt(log(1 + inv(x * t)))
- 2log(sin(x * t / 2)) * t * sqrt(log(1 + inv(x * t)))
```
For the first term we use that `t * sqrt(log(1 + inv(x * t)))` is
increasing for `t < inv(x * (exp(1 / 2) - 1))`. For the second term we
use that the whole term is decreasing for `t < inv(10x)`.
"""
function T011(u0::BHAnsatz{Arb}, ::Ball = Ball(); δ1::Arb = Arb(1e-5), skip_div_u0 = false)
    return x::Arb -> begin
        # Precompute and preallocate for integrand
        xdiv2 = Arblib.mul_2exp!(zero(x), x, -1)
        sin_factor1 = Acb()
        sin_factor2 = Acb()
        sin_factor3 = Acb()
        factor = Acb()

        integrand!(res, t; analytic) = begin
            if analytic && !Arblib.ispositive(Arblib.realref(t))
                # Not analytic at t = 0
                Arblib.indeterminate!(res)
                return
            elseif Arblib.contains_zero(t)
                @assert isreal(t)
                t = real(t)

                tᵤ = ubound(Arb, t)

                # Check that t * sqrt(log(1 + inv(x * t))) is monotone
                if !(tᵤ < inv(x * (exp(Arb(1 // 2)) - 1)))
                    Arblib.indeterminate!(res)
                    return
                end

                # Check that log(sin(x * t / 2)) * t * sqrt(log(1 + inv(x * t)))
                # is monotone
                if !(tᵤ < inv(10x))
                    Arblib.indeterminate!(res)
                    return
                end

                # Enclosure of t * sqrt(log(1 + inv(x * t)))
                tsqrt = Arb((0, tᵤ * sqrt(log(1 + inv(x * tᵤ)))))

                # Enclosure of log(sin(x * t / 2)) * t * sqrt(log(1 + inv(x * t)))
                logsin_tsqrt =
                    Arb((log(sin(tᵤ * xdiv2)) * tᵤ * sqrt(log(1 + inv(x * tᵤ))), 0))

                Arblib.set!(
                    res,
                    log(sin((1 - t) * xdiv2) * sin((1 + t) * xdiv2)) * tsqrt + logsin_tsqrt,
                )
                return
            else
                ## Enclosure of
                ## log(sin((1 - t) * x / 2) * sin((1 + t) * x / 2) / sin(t * x / 2)^2)

                # Enclosure of sin((1 - t) * x / 2)
                Arblib.neg!(sin_factor1, t)
                Arblib.add!(sin_factor1, sin_factor1, 1)
                Arblib.mul!(sin_factor1, sin_factor1, xdiv2)
                Arblib.sin!(sin_factor1, sin_factor1)

                # Enclosure of sin((1 + t) * x / 2)
                Arblib.add!(sin_factor2, t, 1)
                Arblib.mul!(sin_factor2, sin_factor2, xdiv2)
                Arblib.sin!(sin_factor2, sin_factor2)

                # Enclosure of sin(t * x / 2)^2
                Arblib.mul!(sin_factor3, t, xdiv2)
                Arblib.sin!(sin_factor3, sin_factor3)
                Arblib.sqr!(sin_factor3, sin_factor3)

                # Set res to log(sin_factor1 * sin_factor2 / sin_factor3)
                Arblib.mul!(res, sin_factor1, sin_factor2)
                Arblib.div!(res, res, sin_factor3)
                Arblib.log!(res, res)

                # Enclosure of t * sqrt(log(1 + inv(x * t)))
                Arblib.mul!(factor, t, x)
                Arblib.inv!(factor, factor)
                Arblib.add!(factor, factor, 1)
                Arblib.log!(factor, factor)
                Arblib.sqrt!(factor, factor)
                Arblib.mul!(factor, factor, t)

                Arblib.mul!(res, res, factor)
                return
            end
        end

        # -integrand!(res, t; analytic)
        mintegrand!(res, t; analytic) = begin
            integrand!(res, t; analytic)
            Arblib.neg!(res, res)
            return
        end

        # Compute root of integrand
        r_x = _integrand_compute_root(u0, x)
        @assert r_x < 1 - δ1

        rtol = 1e-10
        atol = 1e-10

        U11 = real(
            Arblib.integrate!(
                integrand!,
                zero(Acb),
                0,
                r_x,
                check_analytic = true,
                warn_on_no_convergence = false,
                opts = Arblib.calc_integrate_opt_struct(0, 10_000, 0, 0, 0);
                rtol,
                atol,
            ),
        )

        U121 = real(
            Arblib.integrate!(
                mintegrand!,
                zero(Acb),
                r_x,
                1 - δ1,
                check_analytic = true,
                warn_on_no_convergence = false,
                opts = Arblib.calc_integrate_opt_struct(0, 10_000, 0, 0, 0);
                rtol,
                atol,
            ),
        )

        res = x * (U11 + U121) / (π * sqrt(log(1 + inv(x))))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T012(u0::BHAnsatz; δ1, skip_div_u0)

Return a function for computing an enclosure of
```
inv(u0(x)) * x * U122(x) / (π * sqrt(log(1 + inv(x))))
```
where
```
U122(x) = ∫ -log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(x * t / 2)^2) * t * sqrt(log(1 + inv(x * t))) dt
```
on the interval ``[1 - δ1, 1]``.

# Arguments
- `δ1::Arb = 1e-5`: Determines the interval of integration.
- `skip_div_u0::Bool = false`: If true it skips the factor
  `inv(u0(x))`.

# Implementation
To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval.

Using that `clausenc(x, 1) = -log(2sin(abs(x) / 2))` we are left with
```
∫ clausenc(x * (1 - t), 1) + clausenc(x * (1 + t), 1) - 2clausenc(x * t, 1) dt
```
on the interval ``[1 - δ1, 1]``. A primitive function can be
determined to be
```
1 / x * (-clausens(x * (1 - t), 2) + clausens(x * (1 + t), 2) - 2clausens(x * t, 2))
```
This gives us the integral
```
1 / x * (
    (-clausens(0, 2) + clausens(2x, 2) - 2clausens(x, 2))
    - (-clausens(x * δ1, 2) + clausens(x * (2 - δ1), 2) - 2clausens(x * (1 - δ1), 2))
)
```
Since we multiply this with `x` we can explicitly cancel the `x`
factor.
"""
function T012(u0::BHAnsatz{Arb}, ::Ball = Ball(); δ1::Arb = Arb(1e-5), skip_div_u0 = false)
    return x::Arb -> begin
        # Enclosure of t * sqrt(log(1 + inv(x * t))) on the interval
        # [1 - δ1, 1]
        K1 = let t = Arb((1 - δ1, 1))
            t * sqrt(log(1 + inv(x * t)))
        end

        # Enclosure of integral of Clausen functions multiplied with x
        integral_mul_x = ArbExtras.enclosure_series(x, degree = 4) do y
            (
                (-clausens(Arb(0), 2) + clausens(2y, 2) - 2clausens(y, 2)) - (
                    -clausens(y * δ1, 2) + clausens(y * (2 - δ1), 2) -
                    2clausens(y * (1 - δ1), 2)
                )
            )
        end

        res = K1 * integral_mul_x / (π * sqrt(log((1 + inv(x)))))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
