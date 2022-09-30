export n0_estimate, delta0_estimate, D0_estimate

"""
    n0_estimate(u0::AbstractAnsatz)

Estimate the value of `n0(u0)`.

Uses the observation that the maximum is typically obtained at `x =
π`.

**IMPROVE:** This is not necessarily true in all cases for all choices
of weights. But since it's only an estimate it doesn't really matter.
"""
n0_estimate(u0::AbstractAnsatz{T}) where {T} = u0.w(convert(T, π)) / (2u0(convert(T, π)))

"""
    delta0_estimate(u0::AbstractAnsatz; n::Integer = 100, ϵ = 0, <keyword arguments>)

Estimate the value of `delta0(u0)`.

Does this by evaluating `F0` on `n` linearly spaced points on the
interval ``[0, π]``. Uses an asymptotic expansion on the interval `[0,
ϵ]` and ball arithmetic on `[ϵ, π]`.

# Arguments
- `M::Integer`: determines the number of terms in the asymptotic
  expansions.
- `include_asymptotic = false`: If true then take include `n` x-values
  logarithmically spaced on the interval ``[1e-20, 1e-1]`` in addition
  to the normal x-values. If `ϵ` is zero then set it to `1e-1`.
- `return_values = false`: If true then also return the points and the
  computed values.
- `threaded = false`: If true then use multiple threads when
  evaluating.
"""
function delta0_estimate(
    u0::AbstractAnsatz{T};
    n::Integer = 100,
    ϵ = 0,
    M::Integer = 3,
    include_asymptotic = false,
    return_values = false,
    threaded = true,
) where {T}
    xs = collect(range(zero(T), π, length = n + 1)[2:end])
    if include_asymptotic
        @assert T == Arb # Only supported for Arb
        xs_asym = exp.(range(Arb(log(1e-20)), log(1e-1), length = n))
        xs = sort([xs_asym; xs], by = midpoint)

        if iszero(ϵ)
            ϵ = 1e-1
        end
    end
    res = similar(xs)

    # Asymptotic version might not be defined
    f1 = !iszero(ϵ) ? F0(u0, Asymptotic(); M) : missing
    f2 = F0(u0, Ball())

    if threaded
        Threads.@threads for i in eachindex(xs)
            x = xs[i]
            if x < ϵ
                res[i] = f1(x)
            else
                res[i] = f2(x)
            end
        end
    else
        for i in eachindex(xs)
            x = xs[i]
            if x < ϵ
                res[i] = f1(x)
            else
                res[i] = f2(x)
            end
        end
    end

    m = maximum(abs.(res))

    return_values && return m, xs, res
    return m
end

"""
    D0_estimate(u0::AbstractAnsatz{T}; n::Integer = 20, ϵ = 1, <keyword arguments>)

Estimate the value of `D0(u0)`.

Does this by evaluating the norm on `n` linearly spaced points on the
interval ``[0, π]``. On the interval ``[ϵ, π]`` it uses ball
arithmetic. On the interval ``[0, ϵ]`` it computes with both ball
arithmetic and an asymptotic expansion and takes the best result.

Notice that the default value for `n` is rather low since evaluation
of the integral is costly. Therefore it gives a rather poor estimate
of the norm.

# Arguments
- `M::Integer`: determines the number of terms in the asymptotic
  expansions.
- `x_error = zero(T)`: If `T == Arb` add this error to the points
  before evaluation. This is useful if you want to assure that the
  bounds holds even if the value for `x` is not a tight ball.
- `include_zero = true`: If true then also include the point `x = 0`,
  where asymptotic evaluation is used.
- `return_values = false`: If true then also return the points and the
  computed values.
- `threaded = false`: If true then use multiple threads when
  evaluating.
"""
function D0_estimate(
    u0::AbstractAnsatz{T};
    n::Integer = 20,
    ϵ = 1,
    M::Integer = 3,
    x_error = zero(T),
    include_zero = true,
    return_values = false,
    threaded = true,
) where {T}
    xs = collect(range(zero(T), π, length = n + 1)[2:end])
    res = similar(xs)

    if !iszero(x_error)
        @assert T == Arb
        Arblib.add_error!.(xs, x_error)
        # The error for the last element should not go past π
        xs[end] = min(xs[end], Arb(π))
    end

    f1 = !iszero(ϵ) ? T0(u0, Asymptotic(), ϵ = Arb(ϵ); M) : missing
    f2 = T0(u0, Ball())

    if threaded
        Threads.@threads for i in eachindex(xs)
            x = xs[i]
            if x < ϵ
                res1 = f1(x)
                res2 = f2(x)
                # Use asymptotic result if it is better, otherwise use
                # non-asymptotic result
                res[i] = res1 < ubound(res2) ? res1 : res2
            else
                res[i] = f2(x)
            end
        end
    else
        for i in eachindex(xs)
            x = xs[i]
            if x < ϵ
                res1 = f1(x)
                res2 = f2(x)
                # Use asymptotic result if it is better, otherwise use
                # non-asymptotic result
                res[i] = res1 < ubound(res2) ? res1 : res2
            else
                res[i] = f2(x)
            end
        end
    end

    if include_zero
        pushfirst!(xs, zero(T))
        pushfirst!(res, T0(u0, Asymptotic(); M)(zero(T)))
    end

    m = maximum(abs.(res))

    return_values && return m, xs, res
    return m
end
