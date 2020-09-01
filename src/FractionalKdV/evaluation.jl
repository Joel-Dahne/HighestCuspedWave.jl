# This file contains code for evaluation of the approximate solution
# in different ways

export H, D, F0, hat, eval_expansion

"""
    a0(u0::FractionalKdVAnsatz, j::Integer)
Compute a_j^0 from Lemma 3.2
"""
function a0(u0::FractionalKdVAnsatz{T}, j::Integer) where {T}
    Γ = ifelse(T == arb, Nemo.gamma, SpecialFunctions.gamma)
    return Γ(u0.α - j*u0.p0)*sinpi((1 - u0.α + j*u0.p0)/2)*u0.a[j]
end

"""
    K0(u0::FractionalKdVAnsatz, m::Integer)
Compute K_m^0 from Lemma 3.2
"""
function K0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for j in 0:u0.N0
        res += zeta(1 - u0.α + j*u0.p0 - 2m)*u0.a[j]
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

    for j in 0:u0.N0
        res += zeta(1 - u0.α + j*u0.p0 - 2m)*u0.a[j]
    end

    for n in 1:u0.N1
        res += parent(u0.α)(n)^(2m)*u0.b[n]
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
    return Γ(2*u0.α - j*u0.p0)*sinpi((1 - 2*u0.α + j*u0.p0)/2)*u0.a[j]
end

"""
    L0(u0::FractionalKdVAnsatz, m::Integer)
Compute L_m^0 from Lemma 3.2
"""
function L0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for j in 0:u0.N0
        res += zeta(1 - 2*u0.α + j*u0.p0 - 2m)*u0.a[j]
    end

    res /= factorial(2m)
    res *= (-1)^m

    return res
end

function L0(u0::FractionalKdVAnsatz{arb}, m::Integer)
    res = zero(u0.α)

    for j in 0:u0.N0
        res += zeta(1 - 2*u0.α + j*u0.p0 - 2m)*u0.a[j]
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

    for j in 0:u0.N0
        res += zeta(1 - 2*u0.α + j*u0.p0 - 2m)*u0.a[j]
    end

    for n in 1:u0.N1
        res += parent(u0.α)(n)^(2m + u0.α)*u0.b[n]
    end

    res /= factorial(fmpz(2m))
    res *= (-1)^m

    return res
end

"""
    E(u0::FractionalKdVAnsatz{T}, M::Integer)
Compute an enclosure of E_{u_0}(x) from the paper, skipping the x^(2M) factor.
"""
function E(u0::FractionalKdVAnsatz{T}, M::Integer) where {T}
    return x -> begin
        π = ifelse(T == arb, parent(u0.α)(pi), pi)

        # Compute error bounds for the Clausians
        E_bound1 = zero(u0.α)
        for j in 0:u0.N0
            E_bound1 += (2π)^(j*u0.p0)*abs(zeta(2M + u0.α - j*u0.p0)*u0.a[j])
        end
        E_bound1 *= 2(2π)^(2 - u0.α - 2M)/(4π^2 - x^2)

        # Compute error bounds for the Fourier terms
        E_bound2 = zero(u0.α)
        for n in 1:u0.N1
            E_bound2 += parent(u0.α)(n)^(2M)*abs(u0.b[n])
        end
        E_bound2 /= factorial(fmpz(2M))

        if T == arb
            return ball(zero(u0.α), E_bound1 + E_bound2)
        else
            return E_bound1 + E_bound2
        end
    end
end

"""
    EH(u0::FractionalKdVAnsatz{T}, M::Integer)
Compute an enclosure of E_{H^{-\alpha}u_0}(x) from the paper, skipping
the x^(2M) factor.
"""
function EH(u0::FractionalKdVAnsatz{T}, M::Integer) where {T}
    return x -> begin
        π = ifelse(T == arb, parent(u0.α)(pi), pi)

        # Compute error bounds for the Clausians
        E_bound1 = zero(u0.α)
        for j in 0:u0.N0
            E_bound1 += (2π)^(j*u0.p0)*abs(zeta(2M + 2u0.α - j*u0.p0)*u0.a[j])
        end
        E_bound1 *= 2*(2π)^(2 - 2u0.α - 2M)/(4π^2 - x^2)

        # Compute error bounds for the Fourier terms
        E_bound2 = zero(u0.α)
        for n in 1:u0.N1
            E_bound2 += parent(u0.α)(n)^(2M + u0.α)*abs(u0.b[n])
        end
        E_bound2 /= factorial(fmpz(2M))

        if T == arb
            return ball(zero(u0.α), E_bound1 + E_bound2)
        else
            return E_bound1 + E_bound2
        end
    end
end

"""
    (u0::FractionalKdVAnsatz, evaltype)(x)
Compute u0(x). The strategy for evaluation depends on type of evaltype.
"""
(u0::FractionalKdVAnsatz)(x) = u0(x, Ball())

function (u0::FractionalKdVAnsatz)(x, evaltype::Ball)
    res = zero(u0.α)

    for j in 0:u0.N0
        res += u0.a[j]*(Ci(x, 1 - u0.α + j*u0.p0) - zeta(1 - u0.α + j*u0.p0))
    end

    for n in 1:u0.N1
        res += u0.b[n]*(cos(n*x) - 1)
    end

    return res
end

function (u0::FractionalKdVAnsatz{T})(x,
                                      evaltype::Asymptotic;
                                      M::Integer = 3,
                                      ) where {T}
    @assert M >= 1 - (u0.α + u0.N0*u0.p0)/2

    # Compute approximation
    res = zero(u0.α)
    for j in 0:u0.N0
        res += a0(u0, j)*abs(x)^(-u0.α + j*u0.p0)
    end

    for m in 1:M-1
        res += termK0(u0, m)*x^(2m)
    end

    if T == arb
        res += E(u0, M)(x)*x^(2M)
    end

    return res
end

function (u0::FractionalKdVAnsatz{T})(x,
                                      evaltype::AsymptoticExpansion;
                                      M::Integer = 3,
                                      ) where {T}
    @assert M >= 1 - (u0.α + u0.N0*u0.p0)/2

    expansion = Dict{Tuple{Int, Int, Int}, T}()

    # Compute approximation
    for j in 0:u0.N0
        expansion[(1, j, 0)] = a0(u0, j)
    end

    for m in 1:M-1
        expansion[(0, 0, m)] = termK0(u0, m)
    end

    if T == arb
        expansion[(0, 0, M)] = E(u0, M)(x)
    end

    return expansion
end

"""
    H(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that H(u0)(x) computes H^{-\alpha}u_0(x) from
the paper. The strategy for evaluation depends on type of evaltype.
"""
H(u0::FractionalKdVAnsatz) = H(u0, Ball())

function H(u0::FractionalKdVAnsatz, evaltype::Ball)
    return x -> begin
        res = zero(u0.α)

        for j in 0:u0.N0
            res -= u0.a[j]*(Ci(x, 1 - 2u0.α + j*u0.p0) - zeta(1 - 2u0.α + j*u0.p0))
        end

        for n in 1:u0.N1
            res -= u0.b[n]*parent(u0.α)(n)^u0.α*(cos(n*x) - 1)
        end

        return res
    end
end

function H(u0::FractionalKdVAnsatz{T},
           evaltype::Asymptotic;
           M::Integer = 3,
           ) where {T}
    @assert M >= 1 - u0.α + u0.N0*u0.p0/2
    return x -> begin
        res = zero(u0.α)
        for j in 0:u0.N0
            res -= A0(u0, j)*abs(x)^(-2u0.α + j*u0.p0)
        end

        for m in 1:M-1
            res -= termL0(u0, m)*x^(2m)
        end

        if T == arb
            res += EH(u0, M)(x)*x^(2M)
        end

        return res
    end
end

function H(u0::FractionalKdVAnsatz{T},
           evaltype::AsymptoticExpansion;
           M::Integer = 3,
           ) where {T}
    @assert M >= 1 - u0.α + u0.N0*u0.p0/2
    return x -> begin
        expansion = Dict{Tuple{Int, Int, Int}, T}()

        for j in 0:u0.N0
            expansion[(2, j, 0)] = -A0(u0, j)
        end

        for m in 1:M-1
            expansion[(0, 0, m)] = -termL0(u0, m)
        end

        if T == arb
            expansion[(0, 0, M)] = EH(u0, M)(x)
        end

        return expansion
    end
end

"""
    D(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that D(u0)(x) computes u_0(x)^2/2 +
H^{-\alpha}u_0(x). The strategy for evaluation depends on type of
evaltype.
"""
#D(u0::FractionalKdVAnsatz) = D(u0, Ball())

function D(u0::FractionalKdVAnsatz, evaltype::Ball = Ball())
    return x -> begin
        u0(x, evaltype)^2/2 + H(u0, evaltype)(x)
    end
end

function D(u0::FractionalKdVAnsatz{T},
           evaltype::Asymptotic;
           M = 3,
           ) where {T}
    @assert M >= 1 - u0.α + u0.N0*u0.p0/2
    @assert M == 3 # This is the only one implemented currently
    return x -> begin
        # TODO: Handle terms that might be identically equal to zero
        res = zero(u0.α)

        # TODO: In particular these two might be identically equal to zero
        if u0.N0 >= 0
            res += (a0(u0, 0)^2/2 - A0(u0, 0))*abs(x)^(-2*u0.α)
        end
        if u0.N0 >= 1
            res += (a0(u0, 0)a0(u0, 1) - A0(u0, 1))*abs(x)^(-2*u0.α + u0.p0)
        end

        for k in 2:u0.N0
            term =  -A0(u0, k)
            for j in 0:div(k - 1, 2)
                term += a0(u0, j)*a0(u0, k - j)
            end
            if iseven(k)
                term += a0(u0, div(k, 2))^2/2
            end
            res += term*abs(x)^(-2*u0.α + k*u0.p0)
        end

        for k in (u0.N0 + 1):2u0.N0
            term = zero(u0.α)
            for j in max(k - u0.N0, 0):div(k - 1, 2)
                term += a0(u0, j)*a0(u0, k - j)
            end
            if iseven(k)
                term += a0(u0, div(k, 2))^2/2
            end
            res += term*abs(x)^(-2*u0.α + k*u0.p0)
        end

        termK01 = termK0(u0, 1)
        termK02 = termK0(u0, 2)
        termL01 = termL0(u0, 1)
        termL02 = termL0(u0, 2)

        res -= termL01*x^2
        res += (termK01^2/2 - termL02)*x^4
        res += termK01*termK02*x^6
        res += termK02^2/2*x^8

        for j in 0:u0.N0
            res += a0(u0, j)*termK01*abs(x)^(2 - u0.α + j*u0.p0)
        end

        for j in 0:u0.N0
            res += a0(u0, j)*termK02*abs(x)^(4 - u0.α + j*u0.p0)
        end

        if T == arb
            res += EH(u0, M)(x)*x^(2M)

            error = E(u0, M)(x)*x^(2M)

            term = zero(u0.α)
            for j in 0:u0.N0
                term += a0(u0, j)*abs(x)^(-u0.α + j*u0.p0)
            end

            res += (term + termK01*x^2 + termK02*x^4 + error/2)*error
        end

        return res
    end
end

function D(u0::FractionalKdVAnsatz{T},
           evaltype::AsymptoticExpansion;
           M::Integer = 3,
           ) where {T}
    @assert M >= 1 - u0.α + u0.N0*u0.p0/2
    return x -> begin
        # TODO: Handle terms that might be identically equal to zero
        expansion = Dict{Tuple{Int, Int, Int}, T}()

        # u0^2/2 term
        # TODO: We can reduce the number of computations but using
        # that multiplication is commutative
        expansion1 = u0(x, evaltype, M = M)
        for ((i1, j1, m1), y1) in expansion1
            for ((i2, j2, m2), y2) in expansion1
                key = (i1 + i2, j1 + j2, m1 + m2)
                expansion[key] = get(expansion, key, zero(u0.α)) + y1*y2/2
            end
        end

        # H term
        expansion2 = H(u0, evaltype, M = M)(x)
        for (key, y) in expansion2
            expansion[key] = get(expansion, key, zero(u0.α)) + y
        end

        # Terms in u0.zeroterms are supposed to be identically equal
        # to zero
        for key in u0.zeroterms
            expansion[key] = zero(u0.α)
        end

        return expansion
    end
end

"""
    F0(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that F0(u0)(x) computes F_0(x) from the paper.
The strategy for evaluation depends on type of evaltype.
"""
F0(u0::FractionalKdVAnsatz) = F0(u0, Ball())

function F0(u0::FractionalKdVAnsatz, evaltype::Ball)
    return x -> begin
        return D(u0, evaltype)(x)/(u0.w(x)*u0(x, evaltype))
    end
end

function F0(u0::FractionalKdVAnsatz{T},
            evaltype::Asymptotic;
            M::Integer = 3,
            ) where {T}
    return x -> begin
        res = zero(u0.α)

        expansion = D(u0, AsymptoticExpansion(), M = M)(x)

        for ((i, j, m), y) in expansion
            if !iszero(y)
                res += y*abspow(x, -(i - 1)*u0.α + j*u0.p0 + 2m - u0.p)
            end
        end

        ϵ = ifelse(T == arb, ArbTools.abs_ubound(x), x)
        res *= ball(parent(u0.α)(1), c(u0, ϵ)*abspow(x, u0.p0))/a0(u0, 0)

        return res
    end
end

"""
    hat(u0::FractionalKdVAnsatz)
Returns a function such that hat(u0)(x) computes û(x) from the paper.
"""
function hat(u0::FractionalKdVAnsatz, evaltype::Ball = Ball())
    return x -> begin
        (a0(u0, 0)*abs(x)^(-u0.α) - u0(x))/u0(x)
    end
end


"""
    c(u0::FractionalKdVAnsatz{T}, ϵ)
Compute the constant c_{ϵ,û₀} from Lemma 3.3.
"""
function c(u0::FractionalKdVAnsatz{T},
           ϵ;
           M::Integer = 3,
           ) where {T}
    if iszero(ϵ)
        return zero(u0.α)
    end
    @assert ϵ > 0
    @assert M >= 1 - (u0.α + u0.N0*u0.p0)/2

    numerator = zero(u0.α)
    for j in 1:u0.N0
        numerator += abs(a0(u0, j))*abs(ϵ)^((j - 1)*u0.p0)
    end

    for m in 1:M-1
        numerator += abs(termK0(u0, m))*abs(ϵ)^(2m + u0.α - u0.p0)
    end

    E_lower, E_upper = ArbTools.getinterval(abs(E(u0, M)(ϵ)*ϵ^(2M)))
    numerator += E_upper*abs(ϵ)^(u0.α - u0.p0)

    denominator = abs(a0(u0, 0))
    for j in 1:u0.N0
        denominator -= abs(a0(u0, j))*abs(ϵ)^(j*u0.p0)
    end

    for m in 1:M-1
        denominator -= abs(termK0(u0, m))*abs(ϵ)^(2m + u0.α)
    end

    denominator -= E_upper*abs(ϵ)^(u0.α)

    # TODO: Check this check
    if denominator < 0
        return parent(u0.α)(NaN)
    end

    return numerator/denominator
end

"""
    eval_expansion(u0, expansion, x)
Evaluate the given expansion. The term ((i, j, m), y) is evaluated to
y*abs(x)^(-i*u0.α + j*u0.p0 + 2m) and then they are all summed.

In general x needs to be given both when computing the expansion and
when evaluating it.
"""
function eval_expansion(u0::FractionalKdVAnsatz{T},
                        expansion::Dict{Tuple{Int, Int, Int}, T},
                        x) where {T}
    res = zero(u0.α)

    for ((i, j, m), y) in expansion
        res += y*abspow(x, -i*u0.α + j*u0.p0 + 2m)
    end

    return res
end

"""
    D(u0::FractionalKdVAnsatz, xs::AbstractVector)
Returns a function such that D(u0, xs)(a, b) computes D(u0)(x) on the
points x ∈ xs with u0.a and u0.b set to the given values. Does this in
an efficient way by precomputing as much as possible.
"""
function D(u0::FractionalKdVAnsatz, xs::AbstractVector)
    u0_xs_a_precomputed = zeros(length(xs), u0.N0 + 1)
    u0_xs_b_precomputed = zeros(length(xs), u0.N1)
    Hu0_xs_a_precomputed = zeros(length(xs), u0.N0 + 1)
    Hu0_xs_b_precomputed = zeros(length(xs), u0.N1)

    for i in eachindex(xs)
        x = xs[i]
        for j in 0:u0.N0
            u0_xs_a_precomputed[i, j + 1] = Ci(x, 1 - u0.α + j*u0.p0) - zeta(1 - u0.α + j*u0.p0)
            Hu0_xs_a_precomputed[i, j + 1] = -(Ci(x, 1 - 2u0.α + j*u0.p0) - zeta(1 - 2u0.α + j*u0.p0))
        end
        for n in 1:u0.N1
            u0_xs_b_precomputed[i, n] = cos(n*x) - 1
            Hu0_xs_b_precomputed[i, n] = -parent(u0.α)(n)^u0.α*(cos(n*x) - 1)
        end
    end

    return (a, b) -> begin
        return (
            (u0_xs_a_precomputed*a .+ u0_xs_b_precomputed*b).^2 ./ 2
            .+ (Hu0_xs_a_precomputed*a .+ Hu0_xs_b_precomputed*b)
        )
    end
end

"""
    D(u0::FractionalKdVAnsatz, evaltype::Symbolic, n::Integer)
Returns a function such that D(u0, evaltype, N)(a) computes the
coefficients of the first `u0.N0 + 1` terms in the asymptotic
expansion using the values of `a`. Does this in an efficient way by
precomputing as much as possible.
"""
function D(u0::FractionalKdVAnsatz{T},
           evaltype::Symbolic;
           M::Integer = 10
           ) where {T}
    Γ = ifelse(T == arb, Nemo.gamma, SpecialFunctions.gamma)

    # Precompute for u0
    u0_precomputed = Dict{NTuple{3, Int}, Dict{Int, T}}()
    for j in 0:u0.N0
        s = u0.α - j*u0.p0
        u0_precomputed[(1, j, 0)] = Dict(j => Γ(s)*sinpi((1 - s)/2))
    end

    for m in 1:M-1
        u0_precomputed[(0, 0, m)] = Dict(
            j => (-1)^m*zeta(1 - u0.α + j*u0.p0 - 2m)/factorial(2m)
            for j in 0:u0.N0
        )
    end

    u0_precomputed[(0, 0, M)] = Dict()

    # Precompute H(u0)
    Hu0_precomputed = Dict{NTuple{3, Int}, Dict{Int, T}}()

    for j in 0:u0.N0
        s = 2u0.α - j*u0.p0
        Hu0_precomputed[(2, j, 0)] = Dict(j => -Γ(s)*sinpi((1 - s)/2))
    end

    for m in 1:M-1
        Hu0_precomputed[(0, 0, m)] = Dict(
            j => -(-1)^m*zeta(1 - 2u0.α + j*u0.p0 - 2m)/factorial(2m)
            for j in 0:u0.N0
        )
    end

    Hu0_precomputed[(0, 0, M)] = Dict()

    ## TODO: Check that we do not encounter the error terms. This
    ## should be fine with M = 10 though.

    return a -> begin
        S = promote_type(T, typeof(a))

        # Compute u0
        u0_res = Dict{NTuple{3, Int}, S}()
        for (key, dict) in u0_precomputed
            for (j, v) in dict
                u0_res[key] = get(u0_res, key, zero(u0.α)) + v*a[j]
            end
        end

        # Compute H(u0)
        Hu0_res = Dict{NTuple{3, Int}, S}()
        for (key, dict) in Hu0_precomputed
            for (j, v) in dict
                Hu0_res[key] = get(Hu0_res, key, zero(u0.α)) + v*a[j]
            end
        end

        # Compute u0^2/2
        res = Dict{NTuple{3, Int}, S}()
        for ((i1, j1, m1), y1) in u0_res
            for ((i2, j2, m2), y2) in u0_res
                key = (i1 + i2, j1 + j2, m1 + m2)
                res[key] = get(res, key, zero(u0.α)) + y1*y2/2
            end
        end

        # Compute u0^2/2 + H(u0)
        merge!(+, res, Hu0_res)

        return getindex.(
            sort(
                [((i, j, m), -i*u0.α + j*u0.p0 + 2m, y) for ((i, j, m), y) in res],
                by = x -> Float64(getindex(x, 2)),
            )[3:u0.N0 + 2],
            3,
        )
    end
end
