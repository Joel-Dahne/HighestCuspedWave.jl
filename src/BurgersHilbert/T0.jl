##
## T0
##

function T0(
    u0::BHAnsatz,
    evaltype::Ball;
    δ0 = 1e-10,
    δ1 = 1e-10,
    δ2 = 1e-10,
    rtol = -1, # Not used
    atol = -1, # Not used
    show_trace = false, # Not used
)
    f = T01(u0, evaltype; δ0, δ1)
    g = T02(u0, evaltype; δ2)

    return x -> begin
        ## Integral on [0, x]
        part1 = f(x)

        if isnan(part1)
            # Short circuit on NaN
            return part1
        end

        ## Integral on [x, π]
        part2 = g(x)
        return 2(part1 + part2)
    end
end

##
## T01
##

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

# TODO: Asymptotic version of T01

"""
    T011(u0::BHAnsatz; δ0)
Computes the integral T_{0,1,1} from the paper.

TODO: Implement this
"""
function T011(u0::BHAnsatz, ::Ball = Ball(); δ0 = 1e-10)
    @warn "T011 not yet implemented"
    return x -> begin
        x = convert(Arb, x)

        return zero(x)
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
    a = convert(Arb, δ0)
    b = 1 - convert(Arb, δ1)

    return x -> begin
        x = convert(Arb, x)
        # PROVE: That there are no branch cuts that interact with the
        # integral
        integrand(t; analytic::Bool) = begin
            res = log(sin(x * (1 - t) / 2) * sin(x * (1 + t) / 2) / sin(x * t / 2)^2)
            Arblib.real_abs!(res, res, analytic)
            return res * t * x * sqrt(log((t * x + 1) / (t * x))) # TODO: Avoid hard-coding weight
        end

        res = Arblib.integrate(
            integrand,
            a,
            b,
            check_analytic = true,
            rtol = 1e-10,
            atol = 1e-10,
        )
        @assert isreal(res)
        res = real(res)

        return res * x / (2convert(Arb, π) * u0.w(x) * u0(x))
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

##
## T02
##

"""
    T02(u0::BHAnsatz; δ2)
Returns a function such that T02(u0; δ2, ϵ)(x) computes the
integral T_{0,2} from the paper.
"""
function T02(
    u0::BHAnsatz,
    evaltype::Ball;
    δ2 = 1e-2,
    ϵ = 1e-2, # TODO: Use this
)
    f = T021(u0, evaltype; δ2)
    g = T022(u0, evaltype; δ2)
    return x -> begin
        return f(x) + g(x)
    end
end

# TODO: Asymptotic version of T01

"""
    T021(u0::BHAnsatz; δ2)
Computes the (not yet existing) integral T_{0,2,1} from the paper.

TODO: Implement this
"""
function T021(u0::BHAnsatz, ::Ball = Ball(); δ2 = 1e-10)
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
function T022(u0::BHAnsatz, ::Ball = Ball(); δ2 = 1e-10)
    return x -> begin
        x = convert(Arb, x)
        a = x + δ2
        b = π

        # PROVE: That there are no branch cuts that interact with the
        # integral
        integrand(y; analytic::Bool) = begin
            res = log(-sin((x - y) / 2) * sin((x + y) / 2) / sin(y / 2)^2)
            Arblib.real_abs!(res, res, analytic)
            return res * y * sqrt(log((y + 1) / y)) # TODO: Avoid hard-coding weight
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

        return res / (2convert(Arb, π) * u0.w(x) * u0(x))
    end
end
