"""
    norm(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that norm(u0)(x) computes the norm of the
operator from the paper. The strategy for evaluation depends on type
of evaltype.
"""
norm(u0::FractionalKdVAnsatz) = norm(u0, Ball())

function norm(u0::FractionalKdVAnsatz, evaltype::Ball)
    error("not yet implemented")

    return x -> begin
        return zero(u0.α)
    end
end

function norm(u0::FractionalKdVAnsatz{T}, evaltype::Asymptotic) where {T}
    error("not yet implemented")

    return x -> begin
        return zero(u0.α)
    end
end
