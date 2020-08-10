export FractionalKdVAnsatz

"""
    FractionalKdVAnsatz{T}
α - Corresponds to α in the paper
p0 - Corresponds to p₀ in the paper
a - Vector corresponding to the aⱼ's in the paper, notice that the
indexing here starts at 0 and not 1.
b - Vector corresponding to the bⱼ's in the paper
p - Exponent for the weight, which is |x|^p
zeroterms - List of terms of D(u0) which guaranteed to be identically
equal to zero. It is a set of tuples on the form (i, j, m), see
eval_expansion for more information.
"""
struct FractionalKdVAnsatz{T}
    α::T
    p0::T
    a::OffsetVector{T, Vector{T}}
    b::Vector{T}
    p::T
    zeroterms::Set{Tuple{Int, Int, Int}}
end

"""
    FractionalKdVAnsatz(α::T, N0, N1, p = one(α); kwargs...)
Construct a FractionalKdVAnsatz with the given `α` value and `N0` aⱼ's
and `N1` bₙ's. It sets the weight to be used to |x|^p.

If `use_midpoint` is true then and `T` is `arb` then use only the
midpoint values of `p0`, `a`, `b` and `p` with the exception of `a[0]`
for which the proper enclosure is used.
"""
function FractionalKdVAnsatz(α::T, N0, N1, p = one(α);
                             use_midpoint = true,
                             ) where {T}
    # Using the midpoint only makes sense for T == arb
    use_midpoint = use_midpoint && T == arb

    if use_midpoint
        p = midpoint(p)
    end

    p0 = findp0(α)
    if use_midpoint
        p0 = midpoint(p0)
    end

    a = OffsetVector(fill(zero(α), N0 + 1), 0:N0)
    b = fill(zero(α), N1)

    u0 = FractionalKdVAnsatz(α, p0, a, b, p, Set{Tuple{Int, Int, Int}}())

    findas!(u0)
    findbs!(u0)

    return u0
end

"""
    FractionalKdVAnsatz(α::T)
Construct a FractionalKdVAnsatz using a heuristic choice of
parameters.
"""
function FractionalKdVAnsatz(α::T) where {T}
    ϵ1 = 0.09
    ϵ2 = 0.09

    if α < -1 + ϵ1
        # α ∈ I_1

        u0 = FractionalKdVAnsatz(α, 12, 512, (1 - α)/2, use_midpoint = true)
    elseif α < -ϵ2
        # α ∈ I_2

        if α < -0.75
            N0 = 6
            N1 = 128
        elseif α < -0.6
            N0 = 5
            N1 = 64
        elseif α < -0.35
            N0 = 4
            N1 = 16
        else
            N0 = 3
            N1 = 8
        end

        u0 = FractionalKdVAnsatz(α, N0, N1, (1 - α)/2, use_midpoint = true)
    else
        # α ∈ I_3
        u0 = FractionalKdVAnsatz(α, 1, 0, 1, use_midpoint = false)
    end

    return u0
end


function Base.getproperty(u0::FractionalKdVAnsatz, name::Symbol)
    if name == :N0
        return length(u0.a) - 1
    elseif name == :N1
        return length(u0.b)
    elseif name == :w
        return x -> abs(x)^u0.p
    else
        return getfield(u0, name)
    end
end
