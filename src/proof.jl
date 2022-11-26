"""
    prove(α::Arb; M, only_estimate_D0, D0_maxevals, threaded, verbose, extra_verbose)

Prove the inequality for showing existence of fixed-point. It returns
the result together with a lot of metadata.

Depending on the value of `α` it constructs either a
[`Fractionalkdvansatz`](@ref) or a [`KdVZeroAnsatz`](@ref) and then
uses the [`prove`](@ref) methods corresponding to them.

# Arguments
- `M::Integer = 10`: the number of terms to use in asymptotic expansions
  during the computations, only used for
  [`FractionalKdVAnsatz`](@ref).
- `only_estimate_D0::Bool = false`: if true it doesn't attempt to
  prove the bound for `D₀` but only uses an estimate. This doesn't
  give a rigorous proof but is useful if you only want to determine of
  the bound seems to hold.
- `D0_maxevals::Integer = 4000`: The maximum number of evaluations
  when bounding `D₀`. For [`KdVZeroAnsatz`](@ref) this number is
  multiplied by `3` since it in some cases requires more subdivisions.
- `threaded::Bool = true`: determines if it uses multiple threads for the
  computations or only a single thread.
- `verbose::Bool = false`: if true it prints information about the
  progress.
- `extra_verbose::Bool = false`: if true it prints even more
  information about the progress.
"""
function prove(
    α::Arb;
    M = 10,
    only_estimate_D0 = false,
    D0_maxevals = 4000,
    threaded = true,
    verbose = false,
    extra_verbose = false,
)
    if α < -0.16
        u0_time = @elapsed u0 = FractionalKdVAnsatz(α; threaded)
        p = u0.p
    else
        u0_time = @elapsed u0 = KdVZeroAnsatz(α, midpoint(Arb, α))
        p = one(α)
    end

    verbose && @info "Constructed u0" u0 u0_time

    if u0 isa KdVZeroAnsatz
        D0_maxevals *= 3 # This case usually requires more subdivision
    end

    proof_data =
        prove(u0; M, only_estimate_D0, D0_maxevals, threaded, verbose, extra_verbose)

    return (α = α, p = p, proof_data..., u0_time = u0_time, prec = precision(α))
end

"""
    prove(αs::Vector{Arb}; M, only_estimate_D0, D0_maxevals, executor, threaded, verbose, extra_verbose)

Apply `prove(α)` to all elements of `αs` and return result as a
dataframe.

The mapping over `αs` is done using [`Folds.map`](@ref) and can
therefore be run on several threads or processes, see the `executor`
argument.

It reports on progress using [`Transducers.withprogress`](@ref). This
is unfortunately not supported using the [`Folds`](@ref) interface and
we therefore manually call the underlying transducers.

# Arguments
- `M`, `only_estimate_D0`, `D0_maxevals`, `threaded`, `verbose`,
  `extra_verbos`: Same as for method accepting a single `α`.
- `executor = ThreadedEx(basesize = 1)`: Executor to use for
  [`Folds.map`](@ref). The default value parallelizes over all
  available threads, in which case `threaded` should be false. Using
  `DistributedEx()` parallelizes over all available processes and
  threads, in this case `threaded` should be false. It is possible to
  mix coarse and fine grained parallelism by using
  `DistributedEx(basesize = 1)` and setting `threaded` to true, then
  each process runs one `α` and uses all of its threads for that.
  This option using a large number of processes, each with 2-8 threads
  us often the best from a performance perspective.
"""
function prove(
    αs::Vector{Arb};
    M = 10,
    only_estimate_D0 = true,
    D0_maxevals = 4000,
    executor = ThreadedEx(basesize = 1),
    threaded = false,
    verbose = false,
    extra_verbose = false,
)
    executor == ThreadedEx(basesize = 1) &&
        threaded &&
        @warn "Using threaded executor with threading enabled"

    xf = Map() do α
        HighestCuspedWave.prove(
            α;
            M,
            only_estimate_D0,
            D0_maxevals,
            threaded,
            verbose,
            extra_verbose,
        )
    end

    itr = withprogress(αs, interval = 1)

    # Manually pick the collect to use depending on the executor
    if executor isa DistributedEx
        res = dcollect(xf, itr, basesize = executor.kwargs[:basesize])
    elseif executor isa ThreadedEx
        res = tcollect(xf, itr, basesize = executor.kwargs[:basesize])
    else
        res = collect(xf, itr)
    end

    return DataFrame(res)
end

"""
    proof_interval_subdivisions()
    proof_interval_subdivisions(i::Integer)
    proof_interval_subdivisions(α::AbstractFloat)

Return the subdivision used for the interval ``[1 - δ₁, -δ₂]``.

The interval ``(1, 1 - δ₁)`` is proved using [`BHKdVAnsatz`](@ref) and
``(-δ₂, 0)`` using [`KdVZeroAnsatz`](@ref). The remaining part of the
interval is split into several small balls and handled using a
combination of [`FractionalKdVAnsatz`](@ref) and
[`KdVZeroAnsatz`](@ref). This function handles the splitting into
several smaller intervals.

With no arguments it returns a vector with elements of the form
```
((αₗ, αᵤ), n)
```
This means that the interval ``[αₗ, αᵤ]`` should be split into `n`
equally sized balls. The elements of the vector are ordered with
respect to the endpoints of the intervals, with the right endpoint of
one element being the left endpoint of the next element. That this
indeed is the case can be checked with
[`proof_interval_subdivisions_check`](@ref).

If the argument `i` is given it returns element number `i` in the
vector. If the argument `α` is given return the first element whose
upper bound is not smaller than `α`.

The endpoints of the intervals are stored as `Float64`. This works
well since the endpoints can be picked arbitrarily and in particular
we can choose them to be exact `Float64` numbers.

See also [`proof_interval_subdivisions_mince`](@ref) for getting a
subinterval and splitting it into balls directly.
"""
function proof_interval_subdivisions()
    return [
        ((-0.9999, -0.9998), 2000),
        ((-0.9998, -0.9996), 1000),
        ((-0.9996, -0.9993), 1000),
        ((-0.9993, -0.999), 1000),
        ((-0.999, -0.998), 1000),
        ((-0.998, -0.996), 1000),
        ((-0.996, -0.993), 1000),
        ((-0.993, -0.99), 1000),
        ((-0.99, -0.95), 2000),
        ((-0.95, -0.935), 4000),
        ((-0.935, -0.9), 8000),
        ((-0.9, -0.85), 8000),
        ((-0.85, -0.8), 4000),
        ((-0.8, -0.7), 2000),
        ((-0.7, -0.6), 1000),
        ((-0.6, -0.5), 1500),
        ((-0.5, -0.45), 1000),
        ((-0.45, -0.41), 1000),
        ((-0.41, -0.37), 1000),
        ((-0.37, -0.33), 1000),
        ((-0.33, -0.16), 25000),
        ((-0.16, -0.1), 1000),
        ((-0.1, -0.05), 1000),
        ((-0.05, -0.025), 1000),
        ((-0.025, -0.0125), 2000),
        ((-0.0125, -0.00625), 4000),
        ((-0.00625, -0.003125), 8000),
        ((-0.003125, -0.0015625), 16000),
        ((-0.0015625, -0.0012), 16000),
    ]
end

proof_interval_subdivisions(i::Integer) = proof_interval_subdivisions()[i]

proof_interval_subdivisions(α::AbstractFloat) = proof_interval_subdivisions(
    findfirst(x -> !(x[1][2] < α), proof_interval_subdivisions()),
)

"""
    proof_interval_subdivisions_mince(i, m = nothing; thin = false)

Use [`mince`](@ref) to split the interval given by
`proof_interval_subdivisions(i)` into a number of balls given by
`proof_interval_subdivisions(i). It then returns `m` of these balls.

If `m = nothing` or `m` is larger than the number of balls it returns
all of them.
"""
function proof_interval_subdivisions_mince(i, m = nothing; thin = false)
    ((a, b), n) = proof_interval_subdivisions(i)

    if isnothing(m)
        m = n
    else
        m = min(m, n)
    end

    if m == 1
        indices = [1]
    else
        indices = round.(Int, range(1, n, m))
    end

    αs = mince(Arb((a, b)), n)[indices]

    if thin
        return midpoint.(Arb, αs)
    else
        return αs
    end
end

"""
    proof_interval_subdivisions_check(subdivision = proof_interval_subdivisions())

Check that the subdivisions are valid. Throw an error on failure.

It checks that the interval for each subdivision is a valid interval
and that the right endpoint of each interval is equal to the left
endpoint of the next interval. It doesn't check the left endpoint
of the first interval or the right endpoint of the last interval.
"""
function proof_interval_subdivisions_check(subdivision = proof_interval_subdivisions())
    # Check that the endpoints define an interval
    for ((a, b), _) in subdivision
        ArbExtras.check_interval(Arf(a), Arf(b))
    end

    # Check that all subintervals are consecutive
    for i = 1:length(subdivision)-1
        @assert subdivision[i][1][2] == subdivision[i+1][1][1]
    end
end
