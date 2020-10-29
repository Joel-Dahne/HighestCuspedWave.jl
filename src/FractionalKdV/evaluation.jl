# This file contains code for evaluation of the approximate solution
# in different ways

export hat, eval_expansion

"""
        eval_expansion(u0::FractionalKdVAnsatz, expansion, x)
    Evaluate the given expansion. The term ((i, j, m), y) is evaluated to
    y*abs(x)^(-i*u0.α + j*u0.p0 + m) and then they are all summed.

    In general x needs to be given both when computing the expansion and
    when evaluating it.
    """
function eval_expansion(u0::FractionalKdVAnsatz{T},
                        expansion::AbstractDict{NTuple{3, Int}, T},
                        x,
                        ) where {T}
    res = zero(u0.α)

    for ((i, j, m), y) in expansion
        res += y*abspow(x, -i*u0.α + j*u0.p0 + m)
    end

    return res
end

function (u0::FractionalKdVAnsatz)(x, ::Ball)
    res = zero(u0.α)

    for j in 0:u0.N0
        res += u0.a[j]*(Ci(x, 1 - u0.α + j*u0.p0) - zeta(1 - u0.α + j*u0.p0))
    end

    for n in 1:u0.N1
        res += u0.b[n]*(cos(n*x) - 1)
    end

    return res
end

function (u0::FractionalKdVAnsatz)(x, ::Asymptotic; M::Integer = 3)
    return eval_expansion(u0, u0(x, AsymptoticExpansion(); M), x)
end

function (u0::FractionalKdVAnsatz{T})(x, ::AsymptoticExpansion; M::Integer = 3) where {T}
    @assert M >= 1 - (u0.α + u0.N0*u0.p0)/2

    expansion = OrderedDict{NTuple{3, Int}, T}()

    for j in 0:u0.N0
        expansion[(1, j, 0)] = a0(u0, j)
    end

    for m in 1:M-1
        expansion[(0, 0, 2m)] = termK0(u0, m)
    end

    if T == arb
        expansion[(0, 0, 2M)] = E(u0, M)(x)
    end

    return expansion
end

function H(u0::FractionalKdVAnsatz, ::Ball)
    return x -> begin
        res = zero(u0.α)

        for j in 0:u0.N0
            res -= u0.a[j]*(Ci(x, 1 - 2u0.α + j*u0.p0) - zeta(1 - 2u0.α + j*u0.p0))
        end

        for n in 1:u0.N1
            res -= u0.b[n]*(n*one(u0.α))^u0.α*(cos(n*x) - 1)
        end

        return res
    end
end

function H(u0::FractionalKdVAnsatz, ::Asymptotic; M::Integer = 3)
    f = H(u0, AsymptoticExpansion(); M)
    return x -> eval_expansion(u0, f(x), x)
end

function H(u0::FractionalKdVAnsatz{T}, ::AsymptoticExpansion; M::Integer = 3) where {T}
    @assert M >= 1 - u0.α + u0.N0*u0.p0/2
    return x -> begin
        expansion = OrderedDict{NTuple{3, Int}, T}()

        for j in 0:u0.N0
            expansion[(2, j, 0)] = -A0(u0, j)
        end

        for m in 1:M-1
            expansion[(0, 0, 2m)] = -termL0(u0, m)
        end

        if T == arb
            expansion[(0, 0, 2M)] = EH(u0, M)(x)
        end

        return expansion
    end
end

function D(u0::FractionalKdVAnsatz, ::Asymptotic; M::Integer = 3)
    f = D(u0, AsymptoticExpansion(); M)
    return x -> eval_expansion(u0, f(x), x)
end

function D(u0::FractionalKdVAnsatz{T},
           evaltype::AsymptoticExpansion;
           M::Integer = 3,
           ) where {T}
    f = H(u0, evaltype; M)
    return x -> begin
        expansion1 = u0(x, evaltype; M)
        expansion2 = f(x)

        expansion = empty(expansion1)

        # u0^2/2 term
        # TODO: We can reduce the number of computations but using
        # that multiplication is commutative
        for ((i1, j1, m1), y1) in expansion1
            for ((i2, j2, m2), y2) in expansion1
                key = (i1 + i2, j1 + j2, m1 + m2)
                expansion[key] = get(expansion, key, zero(u0.α)) + y1*y2/2
            end
        end

        # H term
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

function F0(u0::FractionalKdVAnsatz{T}, ::Asymptotic; M::Integer = 3) where {T}
    f = D(u0, AsymptoticExpansion(); M)
    return x -> begin
        res = zero(u0.α)

        expansion = f(x)

        for ((i, j, m), y) in expansion
            if !iszero(y)
                res += y*abspow(x, -(i - 1)*u0.α + j*u0.p0 + m - u0.p)
            end
        end

        ϵ = ifelse(T == arb, ArbTools.abs_ubound(x), x)
        res *= ball(parent(u0.α)(1), c(u0, ϵ)*abspow(x, u0.p0))/a0(u0, 0)

        return res
    end
end

"""
    F0(u0, ::AsymptoticExpansion)

# NOTE
The terms in the dictionary should be interpreted as `((i, j, m), y) →
y⋅x^(-iα + jp₀ + m - p)`, which is different from most other methods.
"""
function F0(u0::FractionalKdVAnsatz{T},
            ::AsymptoticExpansion;
            M::Integer = 3,
            ) where {T}
    return x -> begin
        expansion = D(u0, AsymptoticExpansion(); M)(x)
        res = empty(expansion)

        ϵ = ifelse(T == arb, ArbTools.abs_ubound(x), x)
        C = ball(one(u0.α), c(u0, ϵ)*abspow(x, u0.p0))/a0(u0, 0)

        for ((i, j, m), y) in expansion
            res[(i - 1, j, m)] = C*y
        end

        return res
    end
end

"""
    hat(u0::FractionalKdVAnsatz)
Returns a function such that hat(u0)(x) computes û(x) from the paper.
"""
function hat(u0::FractionalKdVAnsatz, ::Ball = Ball())
    return x -> begin
        (a0(u0, 0)*abs(x)^(-u0.α) - u0(x))/u0(x)
    end
end

"""
    c(u0::FractionalKdVAnsatz{T}, ϵ)
Compute the constant c_{ϵ,û₀} from Lemma 3.3.
"""
function c(u0::FractionalKdVAnsatz, ϵ; M::Integer = 3)
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
function D(u0::FractionalKdVAnsatz{T}, ::Symbolic; M::Integer = 10) where {T}
    Γ = ifelse(T == arb, Nemo.gamma, SpecialFunctions.gamma)

    # Precompute for u0
    u0_precomputed = OrderedDict{NTuple{3, Int}, OrderedDict{Int, T}}()
    for j in 0:u0.N0
        s = u0.α - j*u0.p0
        u0_precomputed[(1, j, 0)] = OrderedDict(j => Γ(s)*sinpi((1 - s)/2))
    end

    for m in 1:M-1
        u0_precomputed[(0, 0, 2m)] = OrderedDict(
            j => (-1)^m*zeta(1 - u0.α + j*u0.p0 - 2m)/factorial(2m)
            for j in 0:u0.N0
        )
    end

    u0_precomputed[(0, 0, 2M)] = OrderedDict()

    # Precompute H(u0)
    Hu0_precomputed = OrderedDict{NTuple{3, Int}, OrderedDict{Int, T}}()

    for j in 0:u0.N0
        s = 2u0.α - j*u0.p0
        Hu0_precomputed[(2, j, 0)] = OrderedDict(j => -Γ(s)*sinpi((1 - s)/2))
    end

    for m in 1:M-1
        Hu0_precomputed[(0, 0, 2m)] = OrderedDict(
            j => -(-1)^m*zeta(1 - 2u0.α + j*u0.p0 - 2m)/factorial(2m)
            for j in 0:u0.N0
        )
    end

    Hu0_precomputed[(0, 0, 2M)] = OrderedDict()

    ## TODO: Check that we do not encounter the error terms. This
    ## should be fine with M = 10 though.

    return a -> begin
        S = promote_type(T, typeof(a))

        # Compute u0
        u0_res = OrderedDict{NTuple{3, Int}, S}()
        for (key, dict) in u0_precomputed
            for (j, v) in dict
                u0_res[key] = get(u0_res, key, zero(u0.α)) + v*a[j]
            end
        end

        # Compute H(u0)
        Hu0_res = OrderedDict{NTuple{3, Int}, S}()
        for (key, dict) in Hu0_precomputed
            for (j, v) in dict
                Hu0_res[key] = get(Hu0_res, key, zero(u0.α)) + v*a[j]
            end
        end

        # Compute u0^2/2
        res = OrderedDict{NTuple{3, Int}, S}()
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
                [((i, j, m), -i*u0.α + j*u0.p0 + m, y) for ((i, j, m), y) in res],
                by = x -> Float64(getindex(x, 2)),
            )[3:u0.N0 + 2],
            3,
        )
    end
end

function print_asymptotic_expansion_D(u0::FractionalKdVAnsatz, expansion)
    get_exponent(i, j, m) = -i*u0.α + j*u0.p0 + m
    expansion = sort(expansion, by = x -> Float64(get_exponent(x...)))
    for ((i, j, m), c) in expansion
        println("$c ⋅ x^$(get_exponent(i, j, m))")
    end

    return expansion
end

function print_asymptotic_expansion_F0(u0::FractionalKdVAnsatz, expansion)
    get_exponent(i, j, m) = -i*u0.α + j*u0.p0 + m - u0.p
    expansion = sort(expansion, by = x -> Float64(get_exponent(x...)))
    for ((i, j, m), c) in expansion
        println("$c ⋅ x^$(get_exponent(i, j, m))   $((i, j, m))")
    end

    return expansion
end
