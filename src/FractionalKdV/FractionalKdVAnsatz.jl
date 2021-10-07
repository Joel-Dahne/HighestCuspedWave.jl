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
struct FractionalKdVAnsatz{T} <: AbstractAnsatz{T}
    α::T
    p0::T
    a::OffsetVector{T,Vector{T}}
    b::Vector{T}
    p::T
    zeroterms::Set{Tuple{Int,Int,Int}}
end

"""
    FractionalKdVAnsatz(α::T, N0, N1, p = one(α); kwargs...)
Construct a FractionalKdVAnsatz with the given `α` value and `N0` aⱼ's
and `N1` bₙ's. It sets the weight to be used to |x|^p.

If `use_midpoint` is true then and `T` is `arb` then use only the
midpoint values of `p0`, `a`, `b` and `p` with the exception of `a[0]`
for which the proper enclosure is used.

If `initial_a` or `initial_b` is given then use these as initial
values for `a` and `b` when solving the system giving the
coefficients. For `a` you should have `initial_a = [a1, a2, …]` and
skip the first `a0`. For `b` you give all of them. It still respects
the values of `N0` and `N1` so it either uses only some of the
coefficients or fills up with zeros.

"""
function FractionalKdVAnsatz(
    α::T,
    N0,
    N1,
    p = one(α);
    use_midpoint = true,
    initial_a::Vector{T} = T[],
    initial_b::Vector{T} = T[],
) where {T}
    # Using the midpoint only makes sense for T == arb
    use_midpoint = use_midpoint && T == arb

    if use_midpoint
        # TODO: We have to remember that we take the midpoint of p!!!
        p = midpoint(p)
    end

    p0 = findp0(α)
    if use_midpoint
        p0 = midpoint(p0)
    end

    a = OffsetVector(
        [
            zero(α)
            collect(Iterators.take(initial_a, N0))
            fill(zero(α), max(N0 - length(initial_a), 0))
        ],
        0:N0,
    )
    b = [
        collect(Iterators.take(initial_b, N1))
        fill(zero(α), max(N1 - length(initial_b), 0))
    ]

    u0 = FractionalKdVAnsatz(α, p0, a, b, p, Set{Tuple{Int,Int,Int}}())

    findas!(u0)
    findbs!(u0)

    return u0
end

"""
    FractionalKdVAnsatz(α::T)
Construct a FractionalKdVAnsatz using a heuristic choice of
parameters.
"""
function FractionalKdVAnsatz(α::T; pp = nothing) where {T}
    ϵ1 = 0.09
    ϵ2 = 0.09

    if α < -1 + ϵ1
        # α ∈ I_1

        u0 = FractionalKdVAnsatz(
            α,
            12,
            512,
            ifelse(isnothing(pp), (1 - α) / 2, pp),
            use_midpoint = true,
        )
    elseif α < -ϵ2
        # α ∈ I_2

        if α <= -0.885
            N0 = 11
            N1 = 512
            p = (1 - α) / 2
        elseif α <= -0.87
            N0 = 10
            N1 = 512
            p = (1 - α) / 2
        elseif α <= -0.85
            N0 = 9
            N1 = 256
            p = (1 - α) / 2
        elseif α <= -0.8
            N0 = 7
            N1 = 256
            p = (1 - α) / 2
        elseif α <= -0.75
            N0 = 6
            N1 = 128
            p = (1 - α) / 2
        elseif α <= -0.61
            N0 = 5
            N1 = 64
            p = (1 - α) / 2
        elseif α <= -0.36
            N0 = 4
            N1 = 16
            p = (1 - α) / 2
        elseif α <= -0.2
            N0 = 3
            N1 = 8
            p = (1 - α) / 2
        else
            N0 = 3
            N1 = 8
            p = (2 - α) / 3
        end

        u0 = FractionalKdVAnsatz(
            α,
            N0,
            N1,
            ifelse(isnothing(pp), p, pp),
            use_midpoint = true,
        )
    else
        # α ∈ I_3
        u0 = FractionalKdVAnsatz(
            α,
            1,
            0,
            ifelse(isnothing(pp), one(α), pp),
            use_midpoint = false,
        )
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
    elseif name == :parent
        return parent(u0.α)
    else
        return getfield(u0, name)
    end
end

function update_alpha(u0::FractionalKdVAnsatz{T}, α::T) where {T}
    p0 = u0.p0
    a = copy(u0.a)
    b = copy(u0.b)
    p = u0.p
    # We can no longer guarantee that any terms are zero
    zeroterms = Set{Tuple{Int,Int,Int}}()

    u0_new = FractionalKdVAnsatz(α, p0, a, b, c, p, zeroterms)

    # We have to update u0.a[0] to make the term (2, 0, 0) zero
    if u0_new.N0 >= 0
        u0_new.a[0] = finda0(u0_new.α)
        push!(u0_new.zeroterms, (2, 0, 0))
    end

    return u0_new
end

function Base.show(io::IO, ::MIME"text/plain", u0::FractionalKdVAnsatz{T}) where {T}
    println(io, "FractionalKdVAnsatz{$T} N₀ = $(u0.N0), N₁ = $(u0.N1)")
    print(io, "α = $(u0.α), p = $(u0.p)")
end

function Base.show(io::IO, ::MIME"text/plain", u0::FractionalKdVAnsatz{arb})
    println(
        io,
        "FractionalKdVAnsatz{arb} N₀ = $(u0.N0), N₁ = $(u0.N1), prec = $(precision(parent(u0.α)))",
    )
    print(io, "α = $(u0.α), p = $(u0.p)")
end

function Base.convert(::Type{FractionalKdVAnsatz{T}}, u0::FractionalKdVAnsatz) where {T}
    return FractionalKdVAnsatz{T}(
        convert(T, u0.α),
        convert(T, u0.p0),
        convert.(T, u0.a),
        convert(Vector{T}, u0.b),
        convert(T, u0.p),
        copy(u0.zeroterms),
    )
end

function Base.convert(::Type{FractionalKdVAnsatz{arb}}, u0::FractionalKdVAnsatz)
    RR = RealField(precision(BigFloat))
    return FractionalKdVAnsatz{arb}(
        RR(u0.α),
        RR(u0.p0),
        RR.(u0.a),
        RR.(u0.b),
        RR(u0.p),
        copy(u0.zeroterms),
    )
end
