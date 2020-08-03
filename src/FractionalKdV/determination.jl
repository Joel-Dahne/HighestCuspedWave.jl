"""
    findp0(α)
Compute the smallest p0 ∈ (0, ∞) such that
cos(π(2α - p0)/2)*Γ(2α -p0)/(cos(π(α - p0)/2)) - Γ(1 + 2α)cos(πα)/(α*Γ(α)cos(πα/2))
is equal to zero.

We use, but don't have to prove, that p0 < 1.5(α + 1).
"""
function findp0(α)
    Γ = SpecialFunctions.gamma
    f(p) = begin
        cospi((2α - p)/2)*Γ(2α - p)/(cospi((α - p)/2)*Γ(α - p)) - Γ(1 + 2α)cospi(α)/(α*Γ(α)cospi(α/2))
    end

    n = 1000
    ps = range(0, stop = 1.51(α + 1), length = n)
    res = f.(ps)

    # Find first sign change
    i = findfirst(i -> res[i]*res[i + 1] <= 0, 1:n-1)

    p0 = first(nlsolve(p -> [f(p[1])], [ps[i]], autodiff = :forward, ftol = 1e-15).zero)

    return p0
end

function findp0(α::arb)
    Γ = Nemo.gamma
    f(p) = begin
        cospi((2α - p)/2)*Γ(2α - p)/(cospi((α - p)/2)*Γ(α - p)) - Γ(1 + 2α)cospi(α)/(α*Γ(α)cospi(α/2))
    end

    # TODO: Prove that this is the smallest positive zero
    # TODO: Refine to required precision
    p0 = ArbTools.add_error!(parent(α)(findp0(Float64(α))), parent(α)(1e-13))
    unique, _ = ArbTools.isuniqueroot(f, ArbTools.getinterval(p0)...)

    @assert unique

    # The refine_root needs to be tuned more before this is an option.
    # Also we still have the problem of not being able to evaluate f
    # on the whole interval due to 0*∞
    #zeros, flags = isolateroots(f, zero(α), 1.5(α + 1), evaltype = :taylor, refine = true)
    #@assert flags[1]
    #p0 = zeros[1]

    return p0
end

function findas!(u0::FractionalKdVAnsatz{T}) where {T}
    Γ = ifelse(T == arb, Nemo.gamma, SpecialFunctions.gamma)

    # Make a0(u0, 0)^2/2 - A0(u0, 0) zero. This is done explicitly to
    # avoid the solution a[0] = 0.
    if u0.N0 >= 0
        u0.a[0] = 2Γ(2u0.α)*sinpi((1 - 2u0.α)/2)/(Γ(u0.α)^2*sinpi((1 - u0.α)/2)^2)
    end
    push!(u0.zeroterms, (2, 0, 0))

    # The choice of p0 makes also the term a0(u0, 0)a0(u0, 1) - A0(u0,
    # 0) = 0
    if u0.N0 >= 1
        if T == arb
            @assert contains_zero(a0(u0, 0)a0(u0, 1) - A0(u0, 1))
        else
            @assert a0(u0, 0)a0(u0, 1) - A0(u0, 1) ≈ 0.0
        end
    end
    push!(u0.zeroterms, (2, 1, 0))

    # TODO: Set the rest of the values

    return u0
end

function findbs!(u0::FractionalKdVAnsatz)
    # TODO: Implement this
end
