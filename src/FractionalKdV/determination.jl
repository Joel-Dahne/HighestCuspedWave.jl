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
        cospi((2α - p)/2)*Γ(2α - p)/(cospi((α - p)/2)*Γ(α - p)) - 2Γ(2α)cospi(α)/(Γ(α)cospi(α/2))
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
    # Do the computations at a higher precision
    if iswide(α)
        # PROVE: That p0 is monotone in α
        α_low, α_upp = getinterval(α)
        return ArbTools.setunion(findp0(α_low), findp0(α_upp))
    end
    α = RealField(2prec(parent(α)))(α)
    Γ = Nemo.gamma
    C = 2Γ(2α)cospi(α)/(Γ(α)cospi(α/2))
    f(p) = begin
        cospi((2α - p)/2)*Γ(2α - p)/(cospi((α - p)/2)*Γ(α - p)) - C
    end

    # PROVE: That it is the smallest positive zero
    p0 = ArbTools.add_error!(parent(α)(findp0(Float64(α))), parent(α)(1e-12))

    found, flags = isolateroots(
        f,
        getinterval(p0)...,
        atol = 1e-15,
        rtol = 1e-15,
        maxevals = 10000,
        evaltype = :taylor,
        refine = true,
        maxfound = 1,
    )

    @assert only(flags)
    p0 = setunion(only(found)...)

    p0 = RealField(div(prec(parent(α)), 2))(p0)

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
        return ArbTools.setunion(finda0(α_low), finda0(α_upp))
    end
    Γ = Nemo.gamma
    return 2Γ(2α)*cospi(α)/(Γ(α)^2*cospi(α/2)^2)
end

function _findas(u0::FractionalKdVAnsatz)
    f = D(u0, Symbolic())
    g(a) = begin
        f(OffsetVector([u0.a[0]; a], 0:u0.N0))
    end

    initial = u0.a[1:end]
    sol = nlsolve(g, initial, autodiff = :forward, iterations = 50)

    if !sol.f_converged
        @warn "Solution did not converge for α = $(u0.α), N0 = $(u0.N0)"
        @warn sol
    end

    return sol.zero
end

"""
    findas(u0)
Find the ones of a[j] for j > 0 that makes the coefficients of the
leading terms in the asymptotic expansion zero.

This is done by solving the corresponding non-linear system.

It uses nlsolve to find the zero, however nlsolve doesn't support
arb-types so this is always done in Float64. To speed it up we start
by computing the first `n` coefficients, then `2n` and so on until we
reach `u0.N0`.
"""
function findas(u0::FractionalKdVAnsatz{T}; minstart = 16) where {T}
    if iszero(u0.N0)
        return T[]
    end
    if u0.N0 <= minstart
        return _findas(u0)
    end

    u0 = deepcopy(u0)

    start = min(max(findlast(!iszero, u0.a), 16), minstart)
    stop = u0.N0
    N0s = [start; [start*2^i for i in 1:floor(Int, log2(stop/start) - 1)]; stop]

    for i in eachindex(N0s)[2:end]
        resize!(u0.a, N0s[i] + 1)
        u0.a[N0s[i - 1]:end] .= zero(T)

        u0.a[1:end] .= _findas(u0)
    end

    return u0.a[1:end]
end

function findas(u0::FractionalKdVAnsatz{arb})
    return parent(u0.α).(findas(convert(FractionalKdVAnsatz{Float64}, u0)))
end

function findas!(u0::FractionalKdVAnsatz{T};
                 use_midpoint = true,
                 ) where {T}
    if u0.N0 >= 0
        u0.a[0] = finda0(u0.α)
    end
    # This makes the term (2, 0, 0) equal to zero
    push!(u0.zeroterms, (2, 0, 0))

    if u0.N0 == 1 && u0.N1 == 0
        # This corresponds to the I_3 case

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
        u0.a[1] = -u0.a[0]*zeta(-1 - 2u0.α)/zeta(-1 - 2u0.α + u0.p0)

        # This makes the term (0, 0, 1) equal to zero.
        if !use_midpoint
            # TODO: Check if it would maybe still be more beneficial
            # to use the midpoint instead
            push!(u0.zeroterms, (0, 0, 1))
        elseif use_midpoint && T == arb
            u0.a[1] = midpoint(u0.a[1])
        end
    else
        a = findas(u0)
        u0.a[1:end] .= a
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
    #n = u0.N1 - 1
    n = u0.N1
    xs = π*(1:2:2n-1)/2n

    f = D(u0, xs)
    g(b) = begin
        #term = zero(u0.α)
        #for n in 1:u0.N1
        #    term += n^(2 + u0.α)*b[n]
        #end

        #[0*L0(u0, 1) - term/2; f(u0.a.parent, b)]
        f(u0.a.parent, b)
    end

    initial = u0.b
    sol = nlsolve(g, initial, autodiff = :forward)

    if !sol.f_converged
        @warn "Solution did not converge for α = $(u0.α), N1 = $(u0.N1)"
        @warn sol
    end

    copy!(u0.b, sol.zero)

    return u0
end

function findbs!(u0::FractionalKdVAnsatz{arb})
    u0_float = convert(FractionalKdVAnsatz{Float64}, u0)

    findbs!(u0_float)
    # TODO: Possibly perform some Newton iterations if the values are
    # required to a higher precision.

    u0.b .= parent(u0.α).(u0_float.b)

    return u0
end
