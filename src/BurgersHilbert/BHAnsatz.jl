export BHAnsatz

"""
    BHAnsatz{T}
a0 - Corresponds to the a₀ coefficient in the paper
b - Vector corresponding to the bₙ's in the paper
v0 - Tail part coming from a FractionalKdVAnsatz
"""
struct BHAnsatz{T} <: AbstractAnsatz{T}
    a0::T
    b::Vector{T}
    v0::Union{Nothing,FractionalKdVAnsatz{T}}
end

function BHAnsatz{T}(
    N::Integer = 0;
    v0::Union{Nothing,FractionalKdVAnsatz} = nothing,
) where {T}
    a0 = 2 / convert(T, π)^2
    b = zeros(T, N)

    u0 = BHAnsatz{T}(a0, b, v0)

    findbs!(u0)

    return u0
end

function Base.getproperty(u0::BHAnsatz, name::Symbol)
    if name == :N
        return length(u0.b)
    elseif name == :w
        return x -> abs(x) * sqrt(log((abs(x) + 1) / abs(x)))
    else
        return getfield(u0, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u0::BHAnsatz{T}) where {T}
    print(io, "BHAnsatz{$T}: N = $(u0.N)")
    isnothing(u0.v0) || print(io, ", v0 is set")
end

function Base.convert(::Type{BHAnsatz{T}}, u0::BHAnsatz) where {T}
    return BHAnsatz{T}(
        convert(T, u0.a0),
        convert(Vector{T}, u0.b),
        isnothing(u0.v0) ? nothing : convert(FractionalKdVAnsatz{T}, u0.v0),
    )
end
