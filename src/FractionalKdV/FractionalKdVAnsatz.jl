export FractionalKdVAnsatz

struct FractionalKdVAnsatz{T}
    α::T
    p0::T
    a::OffsetVector{T, Vector{T}}
    b::Vector{T}
    p::T # The exponent for the weight
end

function FractionalKdVAnsatz(α::T, N0, N1, p = one(α)) where {T}
    p0 = findp0(α)

    a = OffsetVector(fill(zero(α), N0 + 1), 0:N0)
    b = fill(zero(α), N1)

    u0 = FractionalKdVAnsatz(α, p0, a, b, p)

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
