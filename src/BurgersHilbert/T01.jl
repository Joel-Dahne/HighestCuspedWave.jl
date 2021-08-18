"""
    T01(u0::BHAnsatz, ::Ball; δ1, δ2)
Returns a function such that `T01(u0, Ball(); δ1, δ2)(x)` computes the
integral T_{0,1} from the paper.
"""
function T01(u0::BHAnsatz, evaltype::Ball; δ0 = 1e-10, δ1 = 1e-10)
    f = T011(u0, evaltype; δ0)
    g = T012(u0, evaltype; δ0, δ1)
    h = T013(u0, evaltype; δ1)
    return x -> begin
        return f(x) + g(x) + h(x)
    end
end

"""
    T011(u0::BHAnsatz; δ0)
Computes the integral T_{0,1,1} from the paper.

It uses the fact that the integrand is strictly increasing on the
interval `[0, 0.5]` for every value of `x`. The integral is thus
bounded by the length of the interval (`δ0`) times the integrands
value at the right endpoint.
"""
function T011(u0::BHAnsatz, ::Ball = Ball(); δ0 = 1e-10)
    δ0 = convert(Arb, δ0)

    return x -> begin
        x = convert(Arb, x)

        integrand(t) =
            abs(log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(t * x / 2)^2)) *
            t *
            sqrt(log((t * x + 1) / (t * x)))

        return δ0 * integrand(δ0)
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
function T012(u0::BHAnsatz, ::Ball = Ball(); δ0 = 1e-10, δ1 = 1e-10)
    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    a = convert(Arb, δ0)
    b = 1 - convert(Arb, δ1)

    return x -> begin
        x = convert(Arb, x)
        # PROVE: That there are no branch cuts that interact with the
        # integral
        integrand(t; analytic::Bool) = begin
            tx = t * x
            res = log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(tx / 2)^2)
            Arblib.real_abs!(res, res, analytic)
            weight = t * Arblib.sqrt_analytic!(zero(t), log((tx + 1) / tx), analytic)
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

        return res * x / (π * sqrt(log((x + 1) / x)) * u0(x))
    end
end

"""
    T013(u0::BHAnsatz; δ1)
Computes the integral T_{0,1,3} from the paper.

TODO: Implement this
"""
function T013(u0::BHAnsatz, ::Ball = Ball(); δ1 = 1e-10)
    @warn "T013 not yet implemented"
    return x -> begin
        x = convert(Arb, x)

        return zero(x)
    end
end
