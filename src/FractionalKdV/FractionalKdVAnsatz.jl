export FractionalKdVAnsatz

"""
    FractionalKdVAnsatz{T}

Represents an ansatz for a an approximate solution for the fractional
KdV equations with `α ∈ (-1, 0)`. It has the following parameters

- `α`: the parameter defining the equation.
- `p0`: a numerical approximation of ``p_0``.
- `a`: vector with the coefficients in front of the Clausen functions,
  corresponding to ``a_j`` in the paper. Notice that this vector is
  indexed starting from `0` and not `1`.
- `b`: vector with the Fourier coefficients, corresponding to ``b_n``
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

            # We don't warn for a[0] = 0, this is used in BHKdVansatz
            # and is unlikely to be a mistake.
            iszero(a[0]) ||
                Arblib.overlaps(a[0], finda0(α)) ||
                @warn "a[0] doesn't overlap finda0(α)"

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

If `auto_N0_bound` is set to a positive value then don't use the given
`N0` value but instead find the value which gives the smallest defect
using [`find_good_as`](@ref). The maximum `N0` value considered is
then given by `auto_N0_bound`.

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
    N0s::StepRange{Int,Int} = 0:1:-1,
    initial_a::Vector{T} = T[],
    initial_b::Vector{T} = T[],
    verbose = false,
) where {T}
    # Using the midpoint only makes sense for T == Arb
    use_midpoint = use_midpoint && T == Arb

    if use_midpoint
        p = Arblib.midpoint(Arb, p)
    end

    if use_midpoint
        p0 = midpoint(Arb, findp0(midpoint(Arb, α)))
    else
        p0 = findp0(α)
    end

    # Initiate the a vector with only a[0] = finda0(α)
    a = OffsetVector([finda0(α)], 0:0)

    # Initiate the b to be empty
    b = zeros(T, 0)

    # To begin with the only guaranteed zero term is (2, 0, 0), coming
    # from the value of a[0]
    zeroterms = Set{NTuple{3,Int}}([(2, 0, 0)])

    # Create the ansatz
    u0 = FractionalKdVAnsatz{T}(α, p0, a, b, p, zeroterms)

    # Compute values for u0.a[1:end]
    if !isempty(N0s)
        as = find_good_as(u0, N0s; verbose)
        resize!(u0.a, length(as) + 1)
        u0.a[1:end] .= as
    elseif N0 == 1 && α < -0.9
        # If we are close to α = -1 and we only have one Clausen
        # function we take a[1] = -a[0], so that the sum of the first
        # and second term converge towards the leading Clausian for α
        # = -1. This is mostly for testing.
        resize!(u0.a, N0 + 1)
        u0.a[1] = -u0.a[0]
    elseif N0 >= 1
        # Make room for coefficients and set initial values according
        # to initial_a
        resize!(u0.a, N0 + 1)
        u0.a[1:end] .= zero(T)
        u0.a[1:min(N0, length(initial_a))] .= initial_a[1:min(N0, length(initial_a))]

        u0.a[1:end] .= findas(u0)
    end

    if use_midpoint
        u0.a[1:end] .= midpoint.(Arb, u0.a[1:end])
    end

    # Make room for coefficients and set initial values according
    # to initial_b
    resize!(u0.b, N1)
    u0.b[1:end] .= zero(T)
    b[1:min(N1, length(initial_b))] .= initial_a[1:min(N1, length(initial_b))]

    u0.b .= findbs(u0)

    return u0
end

"""
    pick_parameters(::Type{FractionalKdVAnsatz{T}}, α::T) where {T}

Helper function to choose parameters for construction a
`FractionalKdVAnsatz{T}` depending on `α`. It returns `N0s, N1, p`,
which are the parameters to use.

**IMPROVE:** Improve heuristic parameters for `α` closer to `-1`.
"""
function pick_parameters(::Type{FractionalKdVAnsatz{T}}, α::T;) where {T}
    # Old version of parameters:
    # Store heuristic parameters, they are stored as (upper, N0, N1,
    # p) where upper is an upper bound for the α to use these values
    # for and N0, N1 and p are the corresponding parameters for the
    # ansatz.
    #parameters = [
    #    (-0.885, 6, 32, (1 - α) / 2),
    #    (-0.87, 5, 32, (1 - α) / 2),
    #    (-0.848, 9, 32, (1 - α) / 2),
    #    (-0.825, 8, 32, (1 - α) / 2),
    #    (-0.78, 7, 32, (1 - α) / 2),
    #    (-0.75, 6, 32, (1 - α) / 2),
    #    (-0.61, 5, 32, (1 - α) / 2),
    #    (-0.50, 4, 16, (1 - α) / 2),
    #    (-0.36, 4, 16, T(3 // 4)),
    #    (-1 // 3, 3, 8, T(3 // 4)),
    #    (-0.2, 3, 8, one(α)),
    #    (-0.1, 3, 8, one(α)),
    #    (-0.01, 2, 8, one(α)),
    #    (-0.0, 2, 0, one(α)),
    #]

    # Store heuristic parameters, they are stored as (upper, N0s,
    # N1, p) where upper is an upper bound for the α to use these
    # values for. N0s is the values of N0 to consider and N1 and p
    # are the corresponding parameters for the ansatz.
    parameters = [
        (-1.0, 500:1:2000, 4, (1 - α) / 2), # This is never used
        (-0.999, 500:1:2000, 4, (1 - α) / 2), # TODO: Split this up into more steps
        (-0.997, 175:1:600, 4, (1 - α) / 2),
        (-0.995, 100:1:200, 4, (1 - α) / 2),
        (-0.99, 50:1:125, 4, (1 - α) / 2),
        (-0.95, 0:1:75, 8, (1 - α) / 2),
        (-0.5, 0:1:20, 16, (1 - α) / 2),
        (-0.33, 0:1:10, 16, T(3 // 4)),
        (-0.01, 0:1:5, 8, one(α)),
        (0, 0:1:5, 0, one(α)),
    ]

    # Find the first element in parameters which α is not greater
    # than.
    i = findfirst(value -> !(α > value[1]), parameters)

    _, N0s, N1, p = parameters[i]

    return N0s, N1, p
end

"""
    FractionalKdVAnsatz(α::T; N0s = nothing, N1 = nothing, p = nothing, verbose = false)

Construct a `FractionalKdVAnsatz` with the given `α` value using a
heuristic choice of parameters.

The default values for `N0s`, `N1` and `p` can be overridden by
setting the corresponding arguments to a non-nothing value.
"""
function FractionalKdVAnsatz(
    α::T;
    N0s = nothing,
    N1 = nothing,
    p = nothing,
    verbose = false,
) where {T}
    N0s_default, N1_default, p_default = pick_parameters(FractionalKdVAnsatz{T}, α)

    if isnothing(N0s)
        N0s = N0s_default
    end
    if isnothing(N1)
        N1 = N1_default
    end
    if isnothing(p)
        p = p_default
    end

    u0 = FractionalKdVAnsatz(α, 0, N1, p, use_midpoint = true; N0s, verbose)

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
