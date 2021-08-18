"""
    T02(u0::BHAnsatz; δ2)
Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral T_{0,2} from the paper.
"""
function T02(
    u0::BHAnsatz,
    evaltype::Ball;
    δ2::Arb = Arb(1e-10),
    ϵ = 1e-2, # TODO: Use this
)
    f = T021(u0, evaltype; δ2)
    g = T022(u0, evaltype; δ2)
    return x -> begin
        return f(x) + g(x)
    end
end

"""
    T021(u0::BHAnsatz; δ2)
Computes the (not yet existing) integral T_{0,2,1} from the paper.

TODO: Implement this
"""
function T021(u0::BHAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-10))
    @warn "T021 not yet implemented"
    return x -> begin
        x = convert(Arb, x)

        return zero(x)
    end
end

"""
    T022(u0::BHAnsatz; δ2)
Computes the (not yet existing) integral T_{0,2,2} from the paper.

This is done by directly computing the integral with the integrator in
Arb. Accounting for the fact that the integrand is non-analytic at `t
= x`.

TODO: This doesn't work when `x` is close to `π`.
"""
function T022(u0::BHAnsatz, ::Ball = Ball(); δ2::Arb = Arb(1e-10))
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    return x -> begin
        x = convert(Arb, x)
        a = x + δ2
        b = π

        # PROVE: That there are no branch cuts that interact with the
        # integral
        integrand(y; analytic::Bool) = begin
            res = log(-sin((x - y) / 2) * sin((x + y) / 2) / sin(y / 2)^2)
            Arblib.real_abs!(res, res, analytic)
            weight = y * Arblib.sqrt_analytic!(zero(y), log((y + 1) / y), analytic)
            return res * weight
        end

        res = Arblib.integrate(
            integrand,
            a,
            b,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
        )
        @assert !isfinite(res) || isreal(res)
        res = real(res)

        return res / (π * u0.w(x) * u0(x))
    end
end
