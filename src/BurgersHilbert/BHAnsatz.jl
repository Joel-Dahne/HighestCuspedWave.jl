export BHAnsatz

"""
    BHAnsatz{T}
a0 - Corresponds to the a₀ coefficient in the paper
a1 - Coefficient in front of C₂, doesn't exist in the paper yet
b - Vector corresponding to the bₙ's in the paper
"""
struct BHAnsatz{T} <: AbstractAnsatz{T}
    a0::T
    a1::T
    b::Vector{T}
end

function BHAnsatz{T}(N::Integer = 0; a1 = zero(T)) where {T}
    a0 = 2/convert(T, π)^2
    b = zeros(T, N)

    u0 = BHAnsatz{T}(a0, a1, b)

    findbs!(u0)

    return u0
end

function BHAnsatz(RR::RealField, N::Integer = 0; a1 = RR(0))
    a0 = 2/RR(π)^2
    b = fill(zero(a0), N)

    u0 = BHAnsatz{arb}(a0, a1, b)

    findbs!(u0)

    return u0
end

function Base.getproperty(u0::BHAnsatz, name::Symbol)
    if name == :N
        return length(u0.b)
    elseif name == :w
        return x -> abs(x)*sqrt(log((abs(x) + 1)/abs(x)))
    elseif name == :parent
        return parent(u0.a0)
    else
        return getfield(u0, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u0::BHAnsatz{T}) where {T}
    print(io, "BHAnsatz{$T}: N = $(u0.N)")
end

function Base.show(io::IO, ::MIME"text/plain", u0::BHAnsatz{arb})
    print(io, "BHAnsatz{arb}: N = $(u0.N), prec = $(prec(parent(u0.a0)))")
end

function Base.convert(::Type{BHAnsatz{T}}, u0::BHAnsatz) where {T}
    return BHAnsatz{T}(
        convert(T, u0.a0),
        convert(T, u0.a0),
        convert(Vector{T}, u0.b),
    )
end
