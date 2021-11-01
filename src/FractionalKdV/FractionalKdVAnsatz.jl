export FractionalKdVAnsatz

"""
    FractionalKdVAnsatz{T}

Represents an ansatz for a an approximate solution for the fractional
KdV equations with `α ∈ (-1, 0)`. It has the following parameters

- `α`: the parameter defining the equation.
- `p0`: a numerical approximation of \$p_0\$.
- `a`: vector with the coefficients in front of the Clausen functions,
  corresponding to \$a_j\$ in the paper. Notice that this vector is
  indexed starting from `0` and not `1`.
- `b`: vector with the Fourier coefficients, corresponding to \$b_n\$
  in the paper.
- `p`: exponent for the weight, which is given by `abs(x)^p`.
- `zeroterms`: list of terms in the asymptotic expansion of `D(u0)`
  which guaranteed to be identically equal to zero. It is a set of
  tuples on the form `(i, j, m)`, see [`eval_expansion`](@ref) for
  more information.

There are requirements on the parameters needed to give the correct
asymptotic values. These are
- `a[0]` which needs to be that given by `finda0(α)`
- `zeroterms` needs to contain at least the element `(2, 0, 0)` to
  ensure that the defect is finite at zero.
"""
struct FractionalKdVAnsatz{T} <: AbstractAnsatz{T}
    α::T
    p0::T
    a::OffsetVector{T,Vector{T}}
    b::Vector{T}
    p::T
    zeroterms::Set{NTuple{3,Int}}

    function FractionalKdVAnsatz{T}(
        α::T,
        p0::T,
        a::OffsetVector{T,Vector{T}},
        b::Vector{T},
        p::T,
        zeroterms::Set{NTuple{3,Int}},
    ) where {T}
        u0 = new(α, p0, a, b, p, zeroterms)

        if T == Arb
            # Check requirements that we use in the code. We only do
            # this if T == Arb since that's the only case for which we
            # need rigorous values.

            Arblib.overlaps(a[0], finda0(α)) || @warn "a[0] doesn't overlap finda0(α)"

            (2, 0, 0) ∈ zeroterms || @warn "zeroterms doesn't contain the key (2, 0, 0)"
        end

        return u0
    end
end

"""
    FractionalKdVAnsatz(α::T, N0, N1, p = one(α); kwargs...)

Construct a `FractionalKdVAnsatz` with the given `α` value and `N0`
aⱼ's and `N1` bₙ's. It sets the weight to be used to `|x|^p`.

If `use_midpoint` is true and `T` is `Arb` then use only the midpoint
values of `p0`, `a`, `b` and `p` with the exception of `a[0]` for
which the proper enclosure is used.

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
    # Using the midpoint only makes sense for T == Arb
    use_midpoint = use_midpoint && T == Arb

    if use_midpoint
        # TODO: We have to remember that we take the midpoint of p!
        p = Arblib.midpoint(Arb, p)
    end

    p0 = findp0(α)
    if use_midpoint
        # TODO: We have to remember that we take the midpoint of p0!
        p0 = Arblib.midpoint(Arb, p0)
    end

    # Initiate the a vector with a[0] = finda0(α) and the rest of the
    # coefficients zero
    a = OffsetVector([finda0(α); zeros(α, N0)], 0:N0)

    # Initiate the b vector with zeros
    b = zeros(α, N1)

    # Set the initial values of a and b according to initial_a and
    # initial_b
    a[1:min(N0, length(initial_a))] .= initial_a[1:min(N0, length(initial_a))]
    b[1:min(N1, length(initial_b))] .= initial_a[1:min(N1, length(initial_b))]

    # To begin with the only guaranteed zero term is (2, 0, 0), coming
    # from the value of a[0]
    zeroterms = Set{NTuple{3,Int}}([(2, 0, 0)])

    # Create the ansatz
    u0 = FractionalKdVAnsatz{T}(α, p0, a, b, p, zeroterms)

    # Compute values for u0.a[1:end] and u0.b
    if u0.N0 == 1 && α < -0.9
        # If we are close to α = -1 and we only have one Clausen
        # function we take a[1] so that the sum of the first and
        # second term converge towards the leading Clausian for α =
        # -1.
        # TODO: In the end we probably don't want to use this, one
        # Clausen function is not enough. But it is convenient for
        # testing at the moment.
        u0.a[1] = -u0.a[0]
    else
        u0.a[1:end] .= findas(u0)
    end
    u0.b .= findbs(u0)

    if use_midpoint
        u0.a[1:end] .= midpoint.(Arb, u0.a[1:end])
    end

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
        elseif α <= -0.848
            N0 = 9
            N1 = 512
            p = (1 - α) / 2
        elseif α <= -0.825
            N0 = 8
            N1 = 512
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
    zeroterms = Set{NTuple{3,Int}}()

    # We have to update a[0] to make the term (2, 0, 0) zero
    if !isempty(a)
        a[0] = finda0(α)
        push!(zeroterms, (2, 0, 0))
    end

    return FractionalKdVAnsatz{T}(α, p0, a, b, p, zeroterms)
end

function Base.show(io::IO, ::MIME"text/plain", u0::FractionalKdVAnsatz{T}) where {T}
    println(io, "FractionalKdVAnsatz{$T} N₀ = $(u0.N0), N₁ = $(u0.N1)")
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
