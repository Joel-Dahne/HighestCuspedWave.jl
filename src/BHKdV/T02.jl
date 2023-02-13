"""
    T022(u0::BHKdVAnsatz)

Return a functions `f` such that `f(x, a)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * ∫ (clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * u0.w(y) dy
```
where the integration is taken from `a` to `π`.

This is done by directly computing the integral with the integrator in
Arb. In practice `a` should be a thin ball to not give problems with
the integration.

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.
"""
function T022(u0::BHKdVAnsatz, ::Ball = Ball(); δ::Arb = Arb(1e-5), skip_div_u0 = false)
    # Enclosure of -α so that the upper bound is exactly 1
    mα = 1 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    # Upper integration limit
    b = Arb(π)

    return (x::Arb, a::Arb = x + δ; tol = Arb(1e-5)) -> begin
        integrand(y) =
            (clausenc(y - x, mα) + clausenc(y + x, mα) - 2clausenc(y, mα)) * u0.w(y)

        # The integral is decreasing in x so we take the absolute
        # tolerance to depend on x
        res = ArbExtras.integrate(integrand, a, b, atol = x * tol, rtol = tol)

        res /= (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
