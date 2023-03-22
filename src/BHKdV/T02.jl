"""
    T022(u0::BHKdVAnsatz)

Returns a functions such that `T022(u0)(x, a)` computes the integral
```
inv(π * u0(x) * u0.w(x)) * ∫ abs(clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * u0.w(y) dy
```
where the integration is taken from `b` to `π`.

By [`lemma_I_positive`](@ref) the absolute value can be removed. The
integral is then computed directly using
[`ArbExtras.integrate`](@ref). In practice `b` should be a thin ball
to not give problems with the integration.

If `skip_div_u0` is true then skip the division by `u0(x)` in the
result.
"""
function T022(u0::BHKdVAnsatz, ::Ball = Ball(); δ::Arb = Arb(1e-5), skip_div_u0 = false)
    # Enclosure of -α so that the upper bound is exactly 1
    mα = 1 - Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    return (x::Arb, b::Arb = x + δ; tol = Arb(1e-5)) -> begin
        integrand(y) =
            (clausenc(y - x, mα) + clausenc(y + x, mα) - 2clausenc(y, mα)) * u0.w(y)

        # The integral is decreasing in x so we take the absolute
        # tolerance to depend on x
        res = ArbExtras.integrate(integrand, b, Arb(π), atol = x * tol, rtol = tol)

        res /= (π * u0.w(x))

        if skip_div_u0
            return res
        else
            return res / u0(x)
        end
    end
end
