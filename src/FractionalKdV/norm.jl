"""
    norm(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that norm(u0)(x) computes the norm of the
operator from the paper. The strategy for evaluation depends on type
of evaltype.
"""
norm(u0::FractionalKdVAnsatz) = norm(u0, Ball())

function norm(u0::FractionalKdVAnsatz{arb}, evaltype::Ball)
    return x -> begin
        δ0 = parent(u0.α)("1e-6")
        δ1 = parent(u0.α)("1e-6")

        ## Integral on [0, x] - Change to t = y/x

        # Integral on [0, δ0]
        # TODO: Implement this
        T011 = zero(u0.α)

        # Integral on [δ0, 1 - δ1]
        CC = ComplexField(prec(parent(u0.α)))
        mα = CC(-u0.α)
        a = CC(δ0)
        b = CC(1 - δ1)

        F(t, analytic = false) = ArbTools.real_abs(
            Ci(x*(1 - t), mα) + Ci(x*(1 + t), mα) - 2Ci(x*t, mα)
        )*t^u0.p
        T012 = real(ArbTools.integrate(CC, F, a, b, eval_limit = 3000, verbose = 1))

        # Integral on [δ1, 1]
        # TODO: Implement this
        T013 = zero(u0.α)

        T01 = (T011 + T012 + T013)*x/(parent(u0.α)(π)*u0(x))

        ## Integral on [x, π]
        # TODO: Implement this
        T02 = zero(u0.α)

        T02 /= parent(u0.α)(π)*abs(x)^p*u0(x)

        return T01 + T02
    end
end

function norm(u0::FractionalKdVAnsatz{T}, evaltype::Asymptotic) where {T}
    error("not yet implemented")

    return x -> begin
        return zero(u0.α)
    end
end
