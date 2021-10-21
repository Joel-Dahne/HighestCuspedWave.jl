"""
    T02(u0::BHKdVAnsatz; δ2, skip_div_u0)

Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral \$T_{0,2}\$ from the paper.

The interval of integration is `[x, π]`. Since the integrand is
singular at `y = x` we split the interval into two parts, `[x, a]` and
`[a, π]`. In principle we want to take `a = x + δ2`, however this
gives issues if `x` is a wide ball. Instead we take `a` to always be a
thin ball between `x` and `π`.

In general we take `a = ubound(x + δ2)`. However if `x` is wide then
it's beneficial to take a larger value than `δ2`, depending on the
radius of `x`. In practice we take `8radius(x)` in that case.

If `x` is very close to `π`, so that `a` would be larger than `π`,
then we use the asymptotic version for the whole interval `[x, π]`.

If `skip_div_u0` is `true` then don't divide the integral by `u0(x)`.
"""
function T02(u0::BHKdVAnsatz, evaltype::Ball; δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    f = T021(u0, evaltype, skip_div_u0 = true; δ2)
    g = T022(u0, evaltype, skip_div_u0 = true; δ2)

    return x -> begin
        x = convert(Arb, x)
        δ2 = max(δ2, 8Arblib.radius(Arb, x))
        a = Arblib.ubound(Arb, x + δ2)

        if !(a < π)
            return f(x, π)
        end

        res = f(x, a) + g(x, a)

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T021(u0::BHKdVAnsatz)

Computes the (not yet existing) integral \$T_{0,2,1}\$ from the paper.

The interval of integration is `[x, a]`. Both `x` and `a` are assumed
to be less than or equal to `π`, if they are balls which overlap `π`
anything above `π` will be ignored.

To begin with we notice that the weight part of the integrand is well
behaved and we can just factor it out by evaluating it on the whole
interval. We can also notice that the value inside the absolute value
is positive so we can remove the absolute value.

We are left with integrating the three Clausen terms
1. `clausenc(x - y, -α)`
2. `clausenc(x + y, -α)`
3. `2clausenc(y, -α)`
We have that the primitive functions for the three terms are given by
1. `-clausens(x - y, 1 - α)`
2. `clausens(x + y, 1 - α)`
3. `2clausens(y, 1 - α)`
Hence the integral from `x` to `a` is
```
(-clausens(x - a, 1 - α) + clausens(x + a, 1 - α) - 2clausens(a, 1 - α)) -
(-clausens(0, 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α))
```

- **FIXME:** Currently this uses `α = -1`, this should be changed to
  using the interval `[1 - u0.ϵ, 1]` for `α` once [`clausenc`](@ref)
  supports it.
- **TODO:** This doesn't work for `x` overlapping `π` because
  `clausens` doesn't support evaluation around `2π`.
"""
function T021(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)
        δ = a - x

        interval = union(x, a)

        weight_factor = u0.w(interval)

        # FIXME: Use the interval for α once this is supported by clausenc
        #α = Arb((-1, -1 + u0.ϵ))
        α = Arb(-1)
        integral =
            weight_factor * (
                (-clausens(x - a, 1 - α) + clausens(x + a, 1 - α) - 2clausens(a, 1 - α)) - (
                    -clausens(zero(Arb), 1 - α) + clausens(2x, 1 - α) - 2clausens(x, 1 - α)
                )
            )

        res = integral / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end

"""
    T022(u0::BHKdVAnsatz)

Computes the (not yet existing) integral \$T_{0,2,2}\$ from the paper.

The interval of integration is given by `[a, π]`. In practice `a`
should be a thin ball to not give problems with the integration.

This is done by directly computing the integral with the integrator in
Arb.

Notice that the expression inside the absolute value is always
negative, so we can replace the absolute value with a negation.

**FIXME:** This currently uses `α = -1`. We need to either compute
  with `α = [-1, -1 + u0.ϵ]` or bound the error in some other way. One
  issue is that `clausenc` doesn't allow balls overlapping integers.
  Another one is that we get very bad bounds if `α` is close to by not
  equal to `-1`. This could possibly be handled by switching to an
  integration method that only uses real values, where we have better
  bounds.
"""
function T022(u0::BHKdVAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-5), skip_div_u0 = false)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * log(10 + inv(x)))
    end

    return (x, a = x + δ2) -> begin
        x = convert(Arb, x)
        a = convert(Arb, a)

        # Variables for storing temporary values during integration
        x_complex = convert(Acb, x)
        tmp = zero(x_complex)

        integrand!(res, y; analytic::Bool) = begin
            # The code below is an inplace version of the following code
            #res = -log(sin((y - x) / 2) * sin((x + y) / 2) / sin(y / 2)^2)
            #weight = y * log(10 + inv(y))
            #return res * weight

            # res = sin((y - x) / 2)
            Arblib.sub!(tmp, y, x_complex)
            Arblib.mul_2exp!(tmp, tmp, -1)
            Arblib.sin!(res, tmp)

            # res *= sin((x + y) / 2)
            Arblib.add!(tmp, x_complex, y)
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

            # tmp = y * log(10 + inv(y))
            Arblib.inv!(tmp, y)
            Arblib.add!(tmp, tmp, 10)
            Arblib.log!(tmp, tmp)
            Arblib.mul!(tmp, tmp, y)

            Arblib.mul!(res, res, tmp)

            return
        end

        res = Arblib.integrate!(
            integrand!,
            zero(x_complex),
            a,
            π,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
            warn_on_no_convergence = false,
            #opts = Arblib.calc_integrate_opt_struct(0, 0, 0, 0, 1),
        )
        @assert !isfinite(res) || isreal(res)
        res = real(res)

        res = res / (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
