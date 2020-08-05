"""
    T0(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that T0(u0)(x) is the function whose supremum
on [0, π] gives C_B. The strategy for evaluation depends on type of
evaltype.
"""
T0(u0::FractionalKdVAnsatz) = T0(u0, Ball())

function T0(u0::FractionalKdVAnsatz{arb}, evaltype::Ball)
    return x -> begin
        δ0 = parent(u0.α)("1e-6")
        δ1 = parent(u0.α)("1e-6")

        ## Integral on [0, x] - Change to t = y/x

        # Integral on [0, δ0]
        # TODO: Implement this
        T011 = zero(u0.α)*x/(parent(u0.α)(π)*u0(x))

        # Integral on [δ1, 1]
        # TODO: Implement this
        T013 = zero(u0.α)*x/(parent(u0.α)(π)*u0(x))

        T01 = (
            T011
            + T012(u0, δ0 = δ0, δ1 = δ1)(x)
            + T013
        )

        ## Integral on [x, π]
        # TODO: Implement this
        T02 = zero(u0.α)

        T02 /= parent(u0.α)(π)*u0.w(x)*u0(x)

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
     ) = T012(u0, Ball(), δ0 = δ0, δ1 = δ1)

function T012(u0::FractionalKdVAnsatz{arb},
              evaltype::Ball;
              δ0::arb = parent(u0.α)("1e-6"),
              δ1::arb = parent(u0.α)("1e-6"),
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
                                      # TODO: Fix the bug that makes this not work
                                      #rel_tol = 1e-4,
                                      eval_limit = 3000,
                                      verbose = 1,
                                      ))

        return res*x/(parent(u0.α)(π)*u0(x))
    end
end
