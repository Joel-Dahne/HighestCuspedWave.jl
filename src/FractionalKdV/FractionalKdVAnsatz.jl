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

function FractionalKdVAnsatz(α::T, N0, N1, p = one(α)) where {T}
    p0 = findp0(α)

    a = OffsetVector(fill(zero(α), N0 + 1), 0:N0)
    b = fill(zero(α), N1)

    u0 = FractionalKdVAnsatz(α, p0, a, b, p, Set{Tuple{Int, Int, Int}}())

    findas!(u0)
    findbs!(u0)

    # Temporary for testing
    if T == arb
        u0.a[1:end] = parent(α).(rand(N0))
        u0.b .= parent(α).(rand(N1))
    else
        u0.a[1:end] = rand(N0)
        u0.b .= rand(N1)
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
