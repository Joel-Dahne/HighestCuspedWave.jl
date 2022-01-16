export BHAnsatz

"""
    BHAnsatz{T}(a0::T, b::Vector{T}, v0::Union{Nothing,FractionalKdVAnsatz{T}})
    BHAnsatz{T}(N::Integer = 0; v0::Union{Nothing,FractionalKdVAnsatz{T}} = nothing)

Construct an ansatz for use with the Burgers-Hilbert equation. It
stores three parameters
- `a0` which is the coefficient in front of the leading Clausen
  function. Corresponds to ``a_0`` in the paper.
- `b` which is a vector of Fourier coefficients. Corresponds to ``b_n``
  in the paper.
- `v0` which can either be `nothing` or a `FractionalKdVAnsatz`. If it
  is the latter then the non-leading Clausen terms from that are added
  to ansatz.
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

Base.convert(::Type{BHAnsatz{T}}, u0::BHAnsatz) where {T} = BHAnsatz{T}(
    convert(T, u0.a0),
    convert(Vector{T}, u0.b),
    isnothing(u0.v0) ? nothing : convert(FractionalKdVAnsatz{T}, u0.v0),
)
