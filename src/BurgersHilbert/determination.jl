"""
    findbs!(u0::BHAnsatz)
Find values of b[n] to minimize the defect D(u0).

This is done by solving the non-linear system given by requiring that
D(u0) evaluates to zero on N collocation points.

It uses nlsolve to find the zero, however nlsolve doesn't support
arb-types so this is always done in Float64.
"""
function findbs!(u0::BHAnsatz{T}) where {T}
    if u0.N == 0
        return u0
    end

    n = u0.N
    xs = π*(1:2:2n-1)/2n

    f = D(u0, xs)

    # Initial values which originally came from following the branch
    # but are now taken as given.
    initial = T[
        -0.5332599597604509
        -0.12693649125903272
        -0.05148250768408276
        -0.025438828102209
        -0.013604301746747078
        -0.0071927373626929905
        -0.003124688630040378
        -0.14723356162361959
    ]
    initial = [initial; zeros(T, max(n - length(initial), 0))][1:n]
    sol = nlsolve(f, initial, autodiff = :forward)

    if !sol.f_converged
        @warn "Solution did not converge for u0.N = $n"
        @warn sol
    end

    copy!(u0.b, sol.zero)

    return u0
end

function findbs!(u0::BHAnsatz{arb})
    u0_float = convert(BHAnsatz{Float64}, u0)

    findbs!(u0_float)
    # TODO: Possibly perform some Newton iterations if the values are
    # required to a higher precision.

    u0.b .= parent(u0.a0).(u0_float.b)

    return u0
end

function findbs!(u0::BHAnsatz{Arb})
    u0_float = convert(BHAnsatz{Float64}, u0)

    findbs!(u0_float)
    # TODO: Possibly perform some Newton iterations if the values are
    # required to a higher precision.

    u0.b .= u0_float.b

    return u0
end
