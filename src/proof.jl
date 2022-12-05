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
    if midpoint(α) < -0.16
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
    prove(αs::Vector{Arb}; parallelization, M, only_estimate_D0, D0_maxevals, threaded, verbose, extra_verbose)

Apply `prove(α)` to all elements of `αs` and return result as a
dataframe.

The mapping of `αs` can be parallelized in three different ways by
setting `parallelization` accordingly. It can be either `:sequential`,
`:threaded` or `:distributed`. If set to `:sequential` then one value
is processed at a time, possibly using threading for internally. If
set to `:threaded` then threading is used to map over the values of
`αs`, in this case internal threading is always disabled. Finally
`:distributed` can be used to parallelize it over several processes,
see [`_prove_distributed`](@ref).

It reports progress using [`ProgressLogging`](@ref). For this to work
properly you need to enable some form of progress meter. For example
using
```
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
```

# Arguments
- `M`, `only_estimate_D0`, `D0_maxevals`, `threaded`, `verbose`,
  `extra_verbos`: Same as for method accepting a single `α`.
"""
function prove(
    αs::Vector{Arb};
    parallelization = :threaded,
    M = 10,
    only_estimate_D0 = false,
    D0_maxevals = 4000,
    threaded = parallelization != :threaded,
    verbose = false,
    extra_verbose = false,
)
    if parallelization == :threaded && threaded
        @warn "Using threaded parallelization with threading enabled. Disabling threading"
        threaded = false
    end

    if parallelization == :sequential
        res = similar(αs, Any)

        ProgressLogging.@progress for i in eachindex(αs)
            res[i] = HighestCuspedWave.prove(
                αs[i];
                M,
                only_estimate_D0,
                D0_maxevals,
                threaded,
                verbose,
                extra_verbose,
            )
        end

        res = DataFrame(res)
    elseif parallelization == :threaded
        res = similar(αs, Any)

        j = Threads.Atomic{Int}(0)
        ProgressLogging.@withprogress Threads.@threads for i in eachindex(αs)
            res[i] = HighestCuspedWave.prove(
                αs[i];
                M,
                only_estimate_D0,
                D0_maxevals,
                threaded,
                verbose,
                extra_verbose,
            )
            Threads.atomic_add!(j, 1)
            ProgressLogging.@logprogress j[] / length(αs)
        end

        res = DataFrame(res)
    elseif parallelization == :distributed
        res = _prove_distributed(
            αs,
            Distributed.WorkerPool(Distributed.workers());
            M,
            only_estimate_D0,
            D0_maxevals,
            threaded,
            verbose,
            extra_verbose,
        )
    else
        throw(ArgumentError("unknown parallelization type $parallelization"))
    end

    return res
end

# Helper function for _prove_distributed
_prove_distributed_f(
    α,
    M,
    only_estimate_D0,
    D0_maxevals,
    threaded,
    verbose,
    extra_verbose,
) = HighestCuspedWave.prove(
    α;
    M,
    only_estimate_D0,
    D0_maxevals,
    threaded,
    verbose,
    extra_verbose,
)

"""
    _prove_distributed(αs::Vector{Arb}, pool::Distributed.AbstractWorkerPool; kwargs)

Map [`prove`](@ref) over `αs` using the workers in the give worker
pool.

The arguments are keyword arguments are given directly to
[`prove`](@ref).

It reports progress using [`ProgressLogging`](@ref). For this to work
properly you need to enable some form of progress meter. For example
using
```
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
```
Note that all progress logging is done one the main process.
"""
function _prove_distributed(
    αs::Vector{Arb},
    pool::Distributed.AbstractWorkerPool = Distributed.WorkerPool(Distributed.workers());
    M = 10,
    only_estimate_D0 = false,
    D0_maxevals = 4000,
    threaded = true,
    verbose = false,
    extra_verbose = false,
)
    res = similar(αs, Any)

    tasks = map(αs) do α
        @async Distributed.remotecall_fetch(
            _prove_distributed_f,
            pool,
            α,
            M,
            only_estimate_D0,
            D0_maxevals,
            threaded,
            verbose,
            extra_verbose,
        )
    end

    ProgressLogging.@progress res = [fetch(task) for task in tasks]

    DataFrame(res)
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
        ((-0.16, -0.1), 500),
        ((-0.1, -0.05), 500),
        ((-0.05, -0.025), 500),
        ((-0.025, -0.0125), 1000),
        ((-0.0125, -0.00625), 2000),
        ((-0.00625, -0.003125), 4000),
        ((-0.003125, -0.0015625), 8000),
        ((-0.0015625, -0.0012), 4000),
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
