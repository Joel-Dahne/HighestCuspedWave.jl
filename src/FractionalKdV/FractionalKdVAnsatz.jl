export FractionalKdVAnsatz

"""
    FractionalKdVAnsatz{T}

Represents an ansatz for a an approximate solution for the fractional
KdV equations with `α ∈ (-1, 0)`. It has the following parameters

- `α`: the parameter defining the equation.
- `p0`: a numerical approximation of ``p_α``.
- `a`: vector with the coefficients in front of the Clausen functions,
  corresponding to ``a_{α,j}`` in the paper. Notice that this vector is
  indexed starting from `0` and not `1`.
- `b`: vector with the Fourier coefficients, corresponding to
  ``b_{α,n}`` in the paper.
- `p`: exponent for the weight, which is given by `abs(x)^p`.
- `use_bhkdv`: When `α` is close to `-1` we it is beneficial to use a
  slightly different version for the approximation. This is
  essentially a mix of the standard version and some parts from
  [`BHKdVAnsatz`](@ref). It corresponds to the hybrid case near `α =
  -1` discussed in the paper. Setting this value to true means that
  this different version should be used. See below for details about
  what this changes.

The code assumes that `a[0]` is given by `finda0(α)`. This is required
for the leading term in the expansion of `defect(u0)` to be exactly
zero.

# Details about `use_bhkdv`
This flag changes two things in particular, the leading term and the
weight.

With `use_bhkdv = true` the leading term is given by
```
a[0] * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
Compared to
```
a[0] * clausencmzeta(x, 1 - α)
```
when the flag is not used. In practice this means that the coefficient
`a[1]` needs to take different values depending on if the flag is set
or not.

With `use_bhkdv = true` the weight is given by
```
w(x) = abs(x)^p * log(2ℯ + inv(abs(x)))
```
compared to
```
w(x) = abs(x)^p
```
otherwise. One big difference with the modified weight is that it
doesn't satisfy `w(x * y) = w(x) * w(y)`. This means that a number of
simplifications are not supported for that version of the weight, see
also [`weightfactors`](@ref).
"""
struct FractionalKdVAnsatz{T} <: AbstractAnsatz{T}
    α::T
    p0::T
    a::OffsetVector{T,Vector{T}}
    b::Vector{T}
    p::T
    use_bhkdv::Bool

    function FractionalKdVAnsatz{T}(
        α::T,
        p0::T,
        a::OffsetVector{T,Vector{T}},
        b::Vector{T},
        p::T,
        use_bhkdv::Bool = false,
    ) where {T}
        u0 = new(α, p0, a, b, p, use_bhkdv)

        if T == Arb
            # Check requirements that we use in the code. We only do
            # this if T == Arb since that's the only case for which we
            # need rigorous values.

            # We don't throw an error for a[0] = 0, this is used in
            # BHKdVansatz and is unlikely to be a mistake.
            iszero(a[0]) ||
                Arblib.overlaps(a[0], finda0(α)) ||
                error("expected a[0] to overlap finda0(α)")

        end

        return u0
    end
end

"""
    FractionalKdVAnsatz(α::T, N0, N1, p = one(α), use_midpoint = true, use_bhkdv = false)

Construct a `FractionalKdVAnsatz` such to `u0.N0 == N0` and `u0.N1 ==
N1` and with the given `p`.

If `use_midpoint` is true and `T` is `Arb` then use only the midpoint
values of `p0`, `a`, `b` and `p` with the exception of `a[0]` for
which the proper enclosure is used.
"""
function FractionalKdVAnsatz(
    α::T,
    N0::Integer,
    N1::Integer,
    p = one(α);
    use_midpoint = true,
    use_bhkdv = false,
) where {T}
    # Using the midpoint only makes sense for T == Arb
    use_midpoint = use_midpoint && T == Arb

    if use_midpoint
        p = midpoint(Arb, p)
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

    # Create the ansatz
    # use_bhkdv is always false at this stage
    u0 = FractionalKdVAnsatz{T}(α, p0, a, b, p, false)

    # Compute values for u0.a[1:end]
    if N0 >= 1
        # Make room for coefficients and set them to zero
        resize!(u0.a, N0 + 1)
        u0.a[1:end] .= 0
        u0.a[1:end] .= findas(u0)
    end

    if use_midpoint
        u0.a[1:end] .= midpoint.(Arb, u0.a[1:end])
    end

    # Make room for coefficients and set them to zero
    resize!(u0.b, N1)
    u0.b[:] .= 0
    u0.b .= findbs(u0)

    if use_bhkdv
        # Recreate ansatz with use_bhkdv set to true and update
        # u0.a[1] accordingly
        u0 = FractionalKdVAnsatz{T}(u0.α, u0.p0, u0.a, u0.b, u0.p, true)
        u0.a[1] += finda0(midpoint(Arb, u0.α))
    end

    return u0
end


"""
    FractionalKdVAnsatz(α::T, N0s, N1, p = one(α); , use_midpoint = true, use_bhkdv = false, kwargs...)

Construct a `FractionalKdVAnsatz` such to `u0.N1 == N1` and with the
given `p`. The value of `u0.N0` is determined using
[`find_good_as`](@ref), which is given `N0s` as argument.

If `use_midpoint` is true and `T` is `Arb` then use only the midpoint
values of `p0`, `a`, `b` and `p` with the exception of `a[0]` for
which the proper enclosure is used.

It also takes the keyword arguments `try_all_combinations = false`,
`threaded = true` and `verbose = false`, which are given to
[`find_good_as`](@ref).
"""
function FractionalKdVAnsatz(
    α::T,
    N0s::StepRange{Int,Int},
    N1::Integer,
    p = one(α);
    use_midpoint = true,
    try_all_combinations = false,
    use_bhkdv = false,
    threaded = true,
    verbose = false,
) where {T}
    # Using the midpoint only makes sense for T == Arb
    use_midpoint = use_midpoint && T == Arb

    if use_midpoint
        p = midpoint(Arb, p)
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

    # Create the ansatz
    # use_bhkdv is always false at this stage
    u0 = FractionalKdVAnsatz{T}(α, p0, a, b, p, false)

    # Compute values for u0.a[1:end]
    as = find_good_as(u0, N0s; try_all_combinations, threaded, verbose)
    resize!(u0.a, length(as) + 1)
    u0.a[1:end] .= as

    if use_midpoint
        u0.a[1:end] .= midpoint.(Arb, u0.a[1:end])
    end

    # Make room for coefficients and set them to zero
    resize!(u0.b, N1)
    u0.b[:] .= 0
    u0.b .= findbs(u0)

    if use_bhkdv
        # Recreate ansatz with use_bhkdv set to true and update
        # u0.a[1] accordingly
        u0 = FractionalKdVAnsatz{T}(u0.α, u0.p0, u0.a, u0.b, u0.p, true)
        u0.a[1] += finda0(midpoint(Arb, u0.α))
    end

    return u0
end

"""
    pick_parameters(::Type{FractionalKdVAnsatz{T}}, α::T) where {T}

Helper function to choose parameters for construction a
`FractionalKdVAnsatz{T}` depending on `α`. It returns `N0s, N1, p,
try_all_combinations, use_bhkdv`, which are the parameters to use.
"""
function pick_parameters(::Type{FractionalKdVAnsatz{T}}, α::T;) where {T}
    # Store heuristic parameters, they are stored as (upper, N0s,
    # N1, p) where upper is an upper bound for the α to use these
    # values for. N0s is the values of N0 to consider and N1 and p
    # are the corresponding parameters for the ansatz.

    # Compute value to use for N0s close to α = -1
    N0_approx = round(Int, Float64(0.58 / (1 + α)))
    N0s_start = round(Int, 0.95 * N0_approx)
    N0s_stop = round(Int, 1.05 * N0_approx)
    N0s_step = max((N0s_stop - N0s_start) ÷ 10, 1)

    parameters = [
        (-1.0, N0s_start:N0s_step:N0s_stop, 4, (1 - α) / 2, true), # This is never used
        (-0.999, N0_approx:1:N0_approx, 4, (1 - α) / 2, true),
        (-0.997, N0s_start:N0s_step:N0s_stop, 4, (1 - α) / 2, true),
        (-0.995, 100:1:200, 4, (1 - α) / 2, true),
        (-0.99, 50:1:125, 4, (1 - α) / 2, true),
        (-0.95, 5:1:75, 8, (1 - α) / 2, true),
        (-0.85, 0:1:50, 8, (1 - α) / 2, false),
        (-0.5, 0:1:20, 16, (1 - α) / 2, false),
        (-0.33, 0:1:10, 16, T(3 // 4), false),
        (-0.01, 0:1:5, 2, one(α), false),
        (0, 0:1:5, 0, one(α), false),
    ]

    # For Arb we pick the parameter depending on the midpoint of α.
    # This gives more consistent behaviour near the break points.
    if α isa Arb
        α = midpoint(Arb, α)
    end

    # Find the first element in parameters which is greater than or
    # equal to α
    i = findfirst(value -> α <= value[1], parameters)

    _, N0s, N1, p, use_bhkdv = parameters[i]

    try_all_combinations = α > -0.95

    return N0s, N1, p, try_all_combinations, use_bhkdv
end

"""
    FractionalKdVAnsatz(α::T; N0s = nothing, N1 = nothing, p = nothing, try_all_combinations = nothing, use_bhkdv = nothing, verbose = false)

Construct a `FractionalKdVAnsatz` with the given `α` value using a
heuristic choice of parameters.

The default values for `N0s, N1, p, try_all_combinations, use_bhkdv`
can be overridden by setting the corresponding arguments to a
non-nothing value.
"""
function FractionalKdVAnsatz(
    α::T;
    N0s = nothing,
    N1 = nothing,
    p = nothing,
    try_all_combinations = nothing,
    use_bhkdv = nothing,
    threaded = true,
    verbose = false,
) where {T}
    N0s_default, N1_default, p_default, try_all_combinations_default, use_bhkdv_default =
        pick_parameters(FractionalKdVAnsatz{T}, α)

    if isnothing(N0s)
        N0s = N0s_default
    end
    if isnothing(N1)
        N1 = N1_default
    end
    if isnothing(p)
        p = p_default
    end
    if isnothing(try_all_combinations)
        try_all_combinations = try_all_combinations_default
    end
    if isnothing(use_bhkdv)
        use_bhkdv = use_bhkdv_default
    end

    return FractionalKdVAnsatz(
        α,
        N0s,
        N1,
        p,
        use_midpoint = true;
        try_all_combinations,
        use_bhkdv,
        threaded,
        verbose,
    )
end

"""
    weightisx(u0::FractionalKdVAnsatz)

Return true if the weight of `u0` is `abs(x)`.
"""
weightisx(u0::FractionalKdVAnsatz) = isone(u0.p) && !u0.use_bhkdv

"""
    weightfactors(u0::FractionalKdVAnsatz)

Return true if the weight of `u0.w` satisfies that `u0.w(x * y) =
u0.w(x) * u0.w(y)`.
"""
weightfactors(u0::FractionalKdVAnsatz) = !u0.use_bhkdv

"""
    _modified_logabspow(x, c, y)

Helper function for computing `log(c + inv(abs(x))) * abs(x)^y` in a
way that works for `x` overlapping zero.

For `x` overlapping zero it rewrites it as
```
log(c + inv(abs(x))) * abs(x)^y = (log(1 + c * abs(x)) - log(abs(x))) * abs(x)^y
= log(1 + c * abs(x)) * abspow(x, y) - logabspow(x, 1, y)
```
"""
_modified_logabspow(x, c, y) =
    iszero(x) ? abspow(x, y) * log(1 + c * abs(x)) - logabspow(x, 1, y) :
    abspow(x, y) * log(c + inv(abs(x)))

_modified_logabspow(x::Arb, c, y) =
    Arblib.contains_zero(x) ? abspow(x, y) * log(1 + c * abs(x)) - logabspow(x, 1, y) :
    abspow(x, y) * log(c + inv(abs(x)))

_modified_logabspow(x::ArbSeries, c, y) =
    Arblib.contains_zero(x[0]) ? abspow(x, y) * log(1 + c * abs(x)) - logabspow(x, 1, y) :
    abspow(x, y) * log(c + inv(abs(x)))

"""
    _xpwinvw_bhkdv(x)

Helper function for computing `inv(log(2ℯ + inv(abs(x))))` in a way
that works for `x` overlapping zero.
"""
_xpinvw_bhkdv(x) = iszero(x) ? zero(x) : inv(log(2oftype(x, ℯ) + inv(abs(x))))

_xpinvw_bhkdv(x::Arb) =
    if iszero(x)
        return zero(x)
    elseif Arblib.contains_zero(x)
        return Arb((0, inv(log(2Arb(ℯ) + inv(abs_ubound(Arb, x))))))
    else
        return inv(log(2Arb(ℯ) + inv(abs(x))))
    end

_xpinvw_bhkdv(x::ArbSeries) =
    if Arblib.contains_zero(x[0])
        res = indeterminate(x)
        x[0] = _xpinvw_bhkdv(x[0])
        return x
    else
        return inv(log(2Arb(ℯ) + inv(abs(x))))
    end

function Base.getproperty(u0::FractionalKdVAnsatz, name::Symbol)
    if name == :N0
        return (length(u0.a) - 1)::Int
    elseif name == :N1
        return length(u0.b)::Int
    elseif name == :w
        if u0.use_bhkdv
            return x -> _modified_logabspow(x, 2oftype(u0.α, ℯ), u0.p)
        else
            return x -> abspow(x, u0.p)
        end
    elseif name == :wmulpow
        # Compute u0.w(x) * abs(x)^q in a way that works when q < 0
        # and x overlaps zero
        if u0.use_bhkdv
            return (x, q) -> _modified_logabspow(x, 2oftype(u0.α, ℯ), u0.p + q)
        else
            return (x, q) -> abspow(x, u0.p + q)
        end
    elseif name == :xpdivw
        # Compute abs(x)^u0.p / u0.w(x)
        if u0.use_bhkdv
            # abs(x)^u0.p / u0.w(x) = = inv(log(2Arb(ℯ) + inv(abs(x))))
            return _xpinvw_bhkdv
        else
            # abs(x)^u0.p / u0.w(x) = 1
            return one
        end
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

    # We have to update a[0] to make the term (2, 0, 0) zero
    if !isempty(a)
        a[0] = finda0(α)
    end

    return FractionalKdVAnsatz{T}(α, p0, a, b, p, u0.use_bhkdv)
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
        u0.use_bhkdv,
    )
end
