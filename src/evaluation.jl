export eval_expansion, H, defect, F0

"""
    eval_expansion(u0::AbstractAnsatz, expansion, x)

Evaluate the given `expansion` at the point `x`

This is in general used for evaluation around `x = 0`. How the
expansion is represented depends on the precise ansatz used.
"""
function eval_expansion end

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
    defect(u0::AbstractAnsatz, evaltype = Ball())
Returns a function such that `defect(u0)(x)` computes `u_0(x)^2/2 +
H[u_0](x)`.

The strategy for evaluation depends on type of `evaltype`.
"""
defect(u0::AbstractAnsatz) = defect(u0, Ball())

function defect(u0::AbstractAnsatz, evaltype::Ball)
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
        # res = (u0(x)^2 / 2 + f(x)) / (u0.w(x) * u0(x))
        res = inv(u0.w(x))

        # Abort early if non-finite
        isfinite(res) || return res

        y = u0(x, evaltype)

        res /= y

        # Abort early if non-finite
        isfinite(res) || return res

        res = (y^2 / 2 + f(x)) * res

        return res
    end
end
