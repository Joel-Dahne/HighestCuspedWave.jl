# Contains coefficients used in the asymptotic expansions

"""
    a0(u0::FractionalKdVAnsatz, j::Integer)
Compute a_j^0 from Lemma 3.2
"""
function a0(u0::FractionalKdVAnsatz{T}, j::Integer) where {T}
    Γ = ifelse(T == arb, Nemo.gamma, SpecialFunctions.gamma)
    return Γ(u0.α - j * u0.p0) * sinpi((1 - u0.α + j * u0.p0) / 2) * u0.a[j]
end

"""
    K0(u0::FractionalKdVAnsatz, m::Integer)
Compute K_m^0 from Lemma 3.2
"""
function K0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for j = 0:u0.N0
        res += zeta(1 - u0.α + j * u0.p0 - 2m) * u0.a[j]
    end

    res /= factorial(fmpz(2m))
    res *= (-1)^m

    return res
end

"""
    termK0(u0::FractionalKdVAnsatz, m::Integer)
Compute K_m^0 - (-1)^m/((2m)!)sum(n^(2m)b_n for n in 1:N1)
occurring in Lemma 3.2
"""
function termK0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for j = 0:u0.N0
        res += zeta(1 - u0.α + j * u0.p0 - 2m) * u0.a[j]
    end

    for n = 1:u0.N1
        res += parent(u0.α)(n)^(2m) * u0.b[n]
    end

    res /= factorial(fmpz(2m))
    res *= (-1)^m

    return res
end

"""
    A0(u0::FractionalKdVAnsatz, j::Integer)
Compute A_j^0 from Lemma 3.2
"""
function A0(u0::FractionalKdVAnsatz{T}, j::Integer) where {T}
    Γ = ifelse(T == arb, Nemo.gamma, SpecialFunctions.gamma)
    return Γ(2 * u0.α - j * u0.p0) * sinpi((1 - 2 * u0.α + j * u0.p0) / 2) * u0.a[j]
end

"""
    L0(u0::FractionalKdVAnsatz, m::Integer)
Compute L_m^0 from Lemma 3.2
"""
function L0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for j = 0:u0.N0
        res += zeta(1 - 2 * u0.α + j * u0.p0 - 2m) * u0.a[j]
    end

    res /= factorial(fmpz(2m))
    res *= (-1)^m

    return res
end

"""
    termL0(u0::FractionalKdVAnsatz, m::Integer)
Compute L_m^0 - (-1)^m/((2m)!)sum(n^(2m + α)b_n for n in 1:N1)
occurring in Lemma 3.2
"""
function termL0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for j = 0:u0.N0
        res += zeta(1 - 2 * u0.α + j * u0.p0 - 2m) * u0.a[j]
    end

    for n = 1:u0.N1
        res += parent(u0.α)(n)^(2m + u0.α) * u0.b[n]
    end

    res /= factorial(fmpz(2m))
    res *= (-1)^m

    return res
end
