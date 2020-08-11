export T0

"""
    T0(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that T0(u0)(x) is the function whose supremum
on [0, π] gives C_B. The strategy for evaluation depends on type of
evaltype.
"""
T0(u0::FractionalKdVAnsatz; kwargs...) = T0(u0, Ball(); kwargs...)

function T0(u0::FractionalKdVAnsatz{arb},
            evaltype::Ball;
            rtol = -1.0,
            show_trace = false,
            )
    δ0 = parent(u0.α)("1e-6")
    δ1 = parent(u0.α)("1e-6")
    δ2 = parent(u0.α)("1e-6")

    return x -> begin
        ## Integral on [0, x] - Change to t = y/x

        # Integral on [0, δ0]
        # TODO: Implement this
        T011 = zero(u0.α)*x/(parent(u0.α)(π)*u0(x))

        # Integral on [δ1, 1]
        # TODO: Implement this
        T013 = zero(u0.α)*x/(parent(u0.α)(π)*u0(x))

        T01 = (
            T011
            + T012(u0, δ0 = δ0, δ1 = δ1, rtol = rtol, show_trace = show_trace)(x)
            + T013
        )

        ## Integral on [x, π]

        # Integral on [x, x + δ2]
        # TODO: Implement this
        T021 = zero(u0.α)/(parent(u0.α)(π)*u0.w(x)*u0(x))

        T02 = (
            T021
            + T022(u0, δ2 = δ2, rtol = rtol, show_trace = show_trace)(x)
        )

        return T01 + T02
    end
end

function T0(u0::FractionalKdVAnsatz{T}, evaltype::Asymptotic) where {T}
    error("not yet implemented")

    return x -> begin
        return zero(u0.α)
    end
end

"""
    T012(u0::FractionalKdVAnsatz{arb}, δ0, δ1, evaltype)
Returns a function such that T012(u0, δ0, δ1)(x) computes the integral
T_{0,1,2} from the paper.
"""
T012(u0::FractionalKdVAnsatz;
     δ0::arb = parent(u0.α)("1e-6"),
     δ1::arb = parent(u0.α)("1e-6"),
     rtol = -1.0,
     atol = -1.0,
     show_trace = false,
     ) = T012(u0, Ball(), δ0 = δ0, δ1 = δ1, rtol = rtol, atol = atol, show_trace = show_trace)

function T012(u0::FractionalKdVAnsatz{arb},
              evaltype::Ball;
              δ0::arb = parent(u0.α)("1e-6"),
              δ1::arb = parent(u0.α)("1e-6"),
              rtol = -1.0,
              atol = -1.0,
              show_trace = false,
              )
    return x -> begin
        CC = ComplexField(prec(parent(u0.α)))
        mα = CC(-u0.α)
        a = CC(δ0)
        b = CC(1 - δ1)

        # PROVE: That there are no branch cuts that interact with the
        # integral
        F(t, analytic = false) = ArbTools.real_abs(
            Ci(x*(1 - t), mα) + Ci(x*(1 + t), mα) - 2Ci(x*t, mα),
            analytic = analytic,
        )*t^u0.p
        res = real(ArbTools.integrate(CC, F, a, b,
                                      rel_tol = rtol,
                                      abs_tol = atol,
                                      eval_limit = 3000,
                                      verbose = Int(show_trace),
                                      ))

        return res*x/(parent(u0.α)(π)*u0(x))
    end
end

"""
    T022(u0::FractionalKdVAnsatz;
Returns a function such that T022(u0, δ0 = δ0)(x) computes the (not yet
existing) integral T_{0,2,2} from the paper. That is, the integral
from x + δ0 to π. If x + δ0 >= π the function returns zero, note that
the integral T02 might be far from zero in this case.
)
"""
T022(u0::FractionalKdVAnsatz;
     δ2::arb = parent(u0.α)("1e-6"),
     rtol = -1.0,
     atol = -1.0,
     show_trace = false,
     ) = T022(u0, Ball(), δ2 = δ2, rtol = rtol, atol = atol, show_trace = show_trace)

function T022(u0::FractionalKdVAnsatz{arb},
              evaltype::Ball;
              δ2::arb = parent(u0.α)("1e-6"),
              rtol = -1.0,
              atol = -1.0,
              show_trace = false,
              )
    return x -> begin
        if x + δ2 >= parent(u0.α)(π)
            return zero(u0.α)
        end

        CC = ComplexField(prec(parent(u0.α)))
        mα = CC(-u0.α)
        a = CC(x + δ2)
        b = CC(parent(u0.α)(π))

        # PROVE: That there are no branch cuts that interact with the
        # integral
        F(y, analytic = false) = ArbTools.real_abs(
            Ci(x - y, mα) + Ci(x + y, mα) - 2Ci(y, mα),
            analytic = analytic,
        )*y^u0.p
        res = real(ArbTools.integrate(CC, F, a, b,
                                      rel_tol = rtol,
                                      abs_tol = atol,
                                      eval_limit = 2000,
                                      verbose = Int(show_trace),
                                      ))

        return res/(parent(u0.α)(π)*u0.w(x)*u0(x))
    end
end
