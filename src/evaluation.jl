export H, D, F0

"""
    (u0::AbstractAnsatz)(x, evaltype = Ball())
Compute `u0(x)`.

The strategy for evaluation depends on type of `evaltype`.
"""
(u0::AbstractAnsatz)(x) = u0(x, Ball())

"""
    H(u0::AbstractAnsatz, evaltype = Ball())
Returns a function such that `H(u0)(x)` computes `H[u_0](x)` from the
paper.

The strategy for evaluation depends on type of `evaltype`.
"""
H(u0::AbstractAnsatz) = H(u0, Ball())

"""
    D(u0::AbstractAnsatz, evaltype = Ball())
Returns a function such that `D(u0)(x)` computes `u_0(x)^2/2 +
H[u_0](x)`.

The strategy for evaluation depends on type of `evaltype`.
"""
D(u0::AbstractAnsatz) = D(u0, Ball())

function D(u0::AbstractAnsatz, evaltype::Ball)
    f = H(u0, evaltype)
    return x -> begin
        u0(x, evaltype)^2 / 2 + f(x)
    end
end

"""
    F0(u0::AbstractAnsatz, evaltype = Ball())
Returns a function such that `F0(u0)(x)` computes `F_0(x)` from the
paper.

The strategy for evaluation depends on type of `evaltype`.
"""
F0(u0::AbstractAnsatz) = F0(u0, Ball())

function F0(u0::AbstractAnsatz, evaltype::Ball)
    f = H(u0, evaltype)
    return x -> begin
        u0x = u0(x, evaltype)
        return ((u0x)^2 / 2 + f(x)) / (u0.w(x) * u0x)
    end
end
