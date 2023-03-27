"""
    findbs(u0::BHAnsatz)
Find values of `b[n]` to minimize the defect `defect(u0)`.

This is done by solving the non-linear system given by requiring that
`defect(u0)` evaluates to zero on `N` collocation points.

For performance reasons this is always done in `Float64`
"""
function findbs(u0::BHAnsatz{Float64})
    u0.N1 == 0 && return Float64[]

    n = u0.N1
    xs = π * (1:2:2n-1) / 2n

    # Function for computing defect(u0) on the points xs. It takes the
    # coefficients b[n] as an argument
    f = defect(u0, xs)

    # The initial value we use depends on if u0.v0 is set or not
    if iszero(u0.N0)
        # Initial values which originally came from following the branch
        # but are now taken as given.
        initial = [
            -0.5332599597604509
            -0.12693649125903272
            -0.05148250768408276
            -0.025438828102209
            -0.013604301746747078
            -0.0071927373626929905
            -0.003124688630040378
            -0.14723356162361959
        ]
    else
        # In this case we start zero as the initial guess
        initial = zeros()
    end
    initial = [initial; zeros(max(n - length(initial), 0))][1:n]

    sol = nlsolve(f, initial, autodiff = :forward)

    sol.f_converged || @warn "Solution did not converge for u0.N = $n" sol

    return sol.zero
end

findbs(u0::BHAnsatz{T}) where {T} =
    convert(Vector{T}, findbs(convert(BHAnsatz{Float64}, u0)))
