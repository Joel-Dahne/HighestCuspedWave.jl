export alpha0, alpha0_estimate

"""
    alpha0(u0::FractionalKdVAnsatz; ϵ::arb = parent(u0.α)(0.1))
Upper bound the value of α₀ from the paper. Uses an asymptotic
expansion on the interval [0, ϵ] and ball arithmetic on [ϵ, π].

If ϵ = 0 then try to automatically find a suitable value for it. This
is done by first finding a small enough value so that the asymptotic
expansion is well defined and then finding one such that the error
from the asymptotic expansion is smaller than from the direct
evaluation.
"""
function alpha0(u0::FractionalKdVAnsatz{arb};
                ϵ::arb = parent(u0.α)(0),
                M::Integer = 3,
                rtol::arb = parent(u0.α)(1e-2),
                )
    # This is required for a finite value
    @assert u0.p + u0.α > 0

    if iszero(ϵ)
        ϵ = one(u0.α)

        # Make sure the asymptotic values is well defined
        while isnan(c(u0, ϵ, M = M))
            ϵ /= 2
        end

        # Make the asymptotic error smaller than the direct one
        while radius(u0(ϵ)) < radius(u0(ϵ, Asymptotic(), M = M))
            ϵ /= 2
        end
    end

    # Bound the value one [0, ϵ]
    # We have
    # abs(x)^p/u0(x) = (1 + hat(u0)(x))/a0(u0, 0)*abs(x)^(u0.p + u0.α)
    # <= (1 + c(u0, ϵ)*abs(x)^u0.p0/a0(u0, 0)*abs(x)^(u0.p + u0.α)
    # Which attains its maximum value at x = ϵ
    m1 = (1 + c(u0, ϵ)*abs(ϵ)^u0.p0)/a0(u0, 0)*abs(ϵ)^(u0.p + u0.α)

    # Bound the value on [ϵ, π] by Ball evaluation
    # TODO: So far this naive approach only works for a very limited
    # range of parameters and needs to be improved.
    f(x) = begin
        u0.w(x)/(2u0(x))
    end
    m2 = enclosemaximum(f, ϵ, parent(u0.α)(π),
                        rtol = rtol,
                        absmax = true,
                        maxevals = 10^2,
                        )

    # The maximum is in practice attained at x = π, so add that point
    m3 = f(parent(u0.α)(π))

    if m1 > m2
        @warn "We have m1 > m2 which should not happen in practice"
    end

    return max(max(m1, m2), m3)
end

"""
    alpha0_estimate(u0::FractionalKdVAnsatz)
Estimate the value of α₀ from the paper. Uses the observation that the
maximum is obtained at x = π. This always gives a lower bound of α₀.
"""
function alpha0_estimate(u0::FractionalKdVAnsatz{T}) where {T}
    π = ifelse(T == arb, parent(u0.α)(pi), pi)
    return u0.w(π)/(2u0(π))
end
