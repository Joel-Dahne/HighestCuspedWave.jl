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
    if iswide(α)
        # PROVE: That p0 is monotone in α
        α_low, α_upp = getinterval(α)
        return ArbTools.setinterval(findp0(α_low), findp0(α_upp))
    end
    Γ = Nemo.gamma
    f(p) = begin
        cospi((2α - p)/2)*Γ(2α - p)/(cospi((α - p)/2)*Γ(α - p)) - Γ(1 + 2α)cospi(α)/(α*Γ(α)cospi(α/2))
    end

    # PROVE: That it is the smallest positive zero
    p0 = ArbTools.add_error!(parent(α)(findp0(Float64(α))), parent(α)(1e-10))

    # Check that it is a root
    unique, _ = ArbTools.isuniqueroot(f, ArbTools.getinterval(p0)...)
    @assert unique

    # Refine the root
    p0 = setinterval(ArbTools.refine_root(f, ArbTools.getinterval(p0)...)...)

    return p0
end

"""
    finda0(α)
Compute a[0] such that a0(u0, 0)^2/2 - A0(u0, 0) is zero. That is,
compute a[0] = 2Γ(2α)*sinpi((1 - 2α)/2)/(Γ(α)^2*sinpi((1 - α)/2)^2).
It makes use of the monotinicity to get good enclosures for wide
balls.
"""
function finda0(α)
    Γ = SpecialFunctions.gamma
    return 2Γ(2α)*cospi(α)/(Γ(α)^2*cospi(α/2)^2)
end

function finda0(α::arb)
    if iswide(α)
        # PROVE: That a[0] is monotone in α
        α_low, α_upp = getinterval(α)
        return ArbTools.setinterval(finda0(α_low), finda0(α_upp))
    end
    Γ = Nemo.gamma
    return 2Γ(2α)*cospi(α)/(Γ(α)^2*cospi(α/2)^2)
end

function findas!(u0::FractionalKdVAnsatz{T}) where {T}
    if u0.N0 >= 0
        u0.a[0] = finda0(u0.α)
    end
    # This makes the term (2, 0, 0) equal to zero
    push!(u0.zeroterms, (2, 0, 0))

    # The choice of p0 makes also the term (2, 1, 0), given by a0(u0,
    # 0)a0(u0, 1) - A0(u0, 0), equal to zero.
    if u0.N0 >= 1
        if T == arb
            @assert contains_zero(a0(u0, 0)a0(u0, 1) - A0(u0, 1))
        else
            @assert a0(u0, 0)a0(u0, 1) - A0(u0, 1) ≈ 0.0
        end
    end
    push!(u0.zeroterms, (2, 1, 0))

    # This corresponds to the I_3 case
    # TODO: Possibly prove and use that a[1] is monotonically
    # increasing in this case
    if u0.N0 == 1 && u0.N1 == 0
        u0.a[1] = - u0.a[0]*zeta(-1 - 2u0.α)/zeta(-1 - 2u0.α + u0.p0)
        if T == arb
            @assert contains_zero(L0(u0, 1))
        else
            @assert L0(u0, 1) ≈ 0.0
        end
    end

    # TODO: Set the rest of the values

    return u0
end

function findbs!(u0::FractionalKdVAnsatz)
    # TODO: Implement this
end
