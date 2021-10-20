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
            abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) *
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
        integral = weight_factor * (
            (-clausens(zero(Arb), 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α)) -
            (-clausens(x * δ1, 1 - α) + clausens(x * (2 - δ1), 1 - α) - 2clausens(x * (1 - δ1), 1 - α))
        )

        res = integral / (π * log(10 + inv(x)))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
