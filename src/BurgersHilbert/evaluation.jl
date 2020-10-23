function (u0::BHAnsatz)(x, ::Ball)
    res = u0.a0*(Ci(x, 2, 1) - zeta(2, d = 1))

    for n in 1:u0.N
        res += u0.b[n]*(cos(n*x) - 1)
    end

    return res
end

function H(u0::BHAnsatz, ::Ball)
    return x -> begin
        res = -u0.a0*(Ci(x, 3, 1) - zeta(3, d = 1))

        for n in 1:u0.N
            res -= u0.b[n]/n*(cos(n*x) - 1)
        end

        return res
    end
end

"""
    D(u0::BHAnsatz, xs::AbstractVector)
Returns a function such that `D(u0, xs)(b)` computes `D(u0)(x)` on the
points `x ∈ xs` with `u0.b` set to the given values. Does this in an
efficient way by precomputing as much as possible.
"""
function D(u0::BHAnsatz, xs::AbstractVector)
    u0_xs_a0_precomputed = zeros(length(xs))
    u0_xs_b_precomputed = zeros(length(xs), u0.N)
    Hu0_xs_a0_precomputed = zeros(length(xs))
    Hu0_xs_b_precomputed = zeros(length(xs), u0.N)

    for i in eachindex(xs)
        x = xs[i]

        u0_xs_a0_precomputed[i] = u0.a0*(Ci(x, 2, 1) - zeta(2, d = 1))
        Hu0_xs_a0_precomputed[i] = -u0.a0*(Ci(x, 3, 1) - zeta(3, d = 1))

        for n in 1:u0.N
            u0_xs_b_precomputed[i, n] = cos(n*x) - 1
            Hu0_xs_b_precomputed[i, n] = -(cos(n*x) - 1)/n
        end
    end

    return b -> begin
        return (
            (u0_xs_a0_precomputed .+ u0_xs_b_precomputed*b).^2 ./ 2
            .+ (Hu0_xs_a0_precomputed .+ Hu0_xs_b_precomputed*b)
        )
    end
end
