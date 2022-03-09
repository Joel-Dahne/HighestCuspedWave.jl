export BHAnsatz

"""
    BHAnsatz{T}(a0::T, α::T, p0::T, a::Vector{T}, b::Vector{T})
    BHAnsatz{T}(α = -0.9997, N0 = 1929, N1::Integer = 16)

Construct an ansatz for use with the Burgers-Hilbert equation. It
stores three parameters
- `a0`: The coefficient in front of the leading Clausen function. This
  should be equal to `2 / π^2`.
- `α`, `p0` and `a`: stores the parameters for the Clausen terms, the
  have the same meaning as in [`FractionalKdVAnsatz`](@ref), except
  that it starts from `a[1]`.
- `b`: The vector of Fourier coefficients.
"""
struct BHAnsatz{T} <: AbstractAnsatz{T}
    a0::T
    α::T
    p0::T
    a::Vector{T}
    b::Vector{T}
end

function BHAnsatz{T}(
    α = convert(T, -0.9997),
    N0::Integer = 1929,
    N1::Integer = 16,
) where {T}
    a0 = 2 / convert(T, π)^2

    if N0 > 0
        # Construct a FractionalKdVansatz to determine p0 and a
        v0 = FractionalKdVAnsatz(α, N0, 0)
        v0.a[1] = v0.a[0] + v0.a[1]
        if T == Arb
            v0.a[1] = midpoint(Arb, v0.a[1])
        end

        p0 = v0.p0
        a = v0.a[1:end]
    else
        α = convert(T, NaN)
        p0 = convert(T, Nan)
        a = T[]
    end

    b = zeros(T, N1)

    u0 = BHAnsatz{T}(a0, α, p0, a, b)

    copy!(u0.b, findbs(u0))

    return u0
end

function Base.getproperty(u0::BHAnsatz, name::Symbol)
    if name == :N0
        return length(u0.a)
    elseif name == :N1
        return length(u0.b)
    elseif name == :w
        return x -> let absx = abs(x)
            absx * sqrt(log(1 + inv(absx)))
        end
    else
        return getfield(u0, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u0::BHAnsatz{T}) where {T}
    print(io, "BHAnsatz{$T}: α = $(u0.α), N0 = $(u0.N0), N1 = $(u0.N1)")
end

Base.convert(::Type{BHAnsatz{T}}, u0::BHAnsatz) where {T} = BHAnsatz{T}(
    convert(T, u0.a0),
    convert(T, u0.α),
    convert(T, u0.p0),
    convert(Vector{T}, u0.a),
    convert(Vector{T}, u0.b),
)
