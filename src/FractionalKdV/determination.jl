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

function findas!(u0::FractionalKdVAnsatz{T};
                 use_midpoint = true,
                 ) where {T}
    if u0.N0 >= 0
        u0.a[0] = finda0(u0.α)
    end
    # This makes the term (2, 0, 0) equal to zero
    push!(u0.zeroterms, (2, 0, 0))

    # TODO: We might not want to do it in this way, we might have to
    # order the terms by their exponent and make them zero.
    if u0.N0 >= 1
        # The choice of p0 makes also the term (2, 1, 0), given by a0(u0,
        # 0)a0(u0, 1) - A0(u0, 0), equal to zero.
        if !use_midpoint
            if T == arb
                @assert contains_zero(a0(u0, 0)a0(u0, 1) - A0(u0, 1))
            else
                @assert a0(u0, 0)a0(u0, 1) - A0(u0, 1) ≈ 0.0
            end
            push!(u0.zeroterms, (2, 1, 0))
        end

        # Compute a[1] such that L0(u0, 1) is zero.
        # TODO: Possibly make use of the monotinicity to get good
        # enclosures for wide balls.
        # TODO: This only really makes sense in the I_3 case when we
        # do not have any more a[j] terms.
        u0.a[1] = -u0.a[0]*zeta(-1 - 2u0.α)/zeta(-1 - 2u0.α + u0.p0)

        # If there are no b[n]'s this makes the term (0, 0, 1) equal
        # to zero. This corresponds to the I_3 case
        if u0.N0 == 1 && iszero(u0.N1) && !use_midpoint
            # TODO: Check if it would maybe still be more beneficial
            # to use the midpoint instead
            push!(u0.zeroterms, (0, 0, 1))
        elseif use_midpoint
            u0.a[1] = midpoint(u0.a[1])
        end
    end

    # We choose the remaining a[j]'s to make terms from the first sum
    # (the third term) in Lemma 3.4 zero.
    Γ = ifelse(T == arb, Nemo.gamma, SpecialFunctions.gamma)
    for k in 2:u0.N0
        term = zero(u0.α)
        for j in 1:div(k - 1, 2)
            term += a0(u0, j)*a0(u0, k - j)
        end
        if iseven(k)
            term += a0(u0, div(k, 2))^2/2
        end
        u0.a[k] = -term/(
            a0(u0, 0)*Γ(u0.α - k*u0.p0)*sinpi((1 - u0.α + k*u0.p0)/2)
            - Γ(2u0.α - k*u0.p0)*sinpi((1 - 2u0.α + k*u0.p0)/2)
        )

        if use_midpoint
            u0.a[k] = midpoint(u0.a[k])
        else
            push!(u0.zeroterms, (2, k, 0))
        end

    end

    return u0
end

"""
    findbs!(u0)
Find values of b[n] to minimize the defect D(u0).

This is done by solving the non-linear system given by requiring that
D(u0) evaluates to zero on N1 collocation points.

It uses nlsolve to find the zero, however nlsolve doesn't support
arb-types so this is always done in Float64.

TODO: Possibly add support for optimizing some of the the a[j]'s at
the same time,
"""
function findbs!(u0::FractionalKdVAnsatz)
    if u0.N1 == 0
        return u0
    end
    n = u0.N1
    xs = π*(1:2:2n-1)/2n

    f = D(u0, xs)
    g(b) = begin
        f(u0.a.parent, b)
    end

    initial = fill(zero(u0.α), u0.N1)
    sol = nlsolve(g, initial, autodiff = :forward)

    if !sol.f_converged
        @warn "Solution did not converge"
        @warn sol
    end

    copy!(u0.b, sol.zero)

    return u0
end

function findbs!(u0::FractionalKdVAnsatz{arb})
    T = Float64
    u0_float = FractionalKdVAnsatz(
        T(u0.α),
        T(u0.p0),
        T.(u0.a),
        T.(u0.b),
        T(u0.p),
        copy(u0.zeroterms),
    )

    findbs!(u0_float)
    # TODO: Possibly perform some Newton iterations if the values are
    # required to a higher precision.

    u0.b .= parent(u0.α).(u0_float.b)

    return u0
end
