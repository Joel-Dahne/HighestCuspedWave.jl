using ClusterManagers, CSV, DataFrames, Distributed, Folds

include("helper.jl")

# Run for one α
function run_single(α; only_estimate_D0 = false)
    proof_data = HighestCuspedWave.prove(
        α,
        M = 10,
        threaded = true,
        verbose = true,
        extra_verbose = false;
        only_estimate_D0,
    )

    (α = α, proof_data..., prec = precision(α))
end

# Run on all αs
function run_map(
    αs,
    executor = ThreadedEx(basesize = 1);
    only_estimate_D0 = true,
    threaded = false,
)
    HighestCuspedWave.prove(αs; executor, only_estimate_D0, threaded)
end

# Benchmark run_single
function benchmark_single(α; only_estimate_D0 = false)
    @time "Compile run_single" run_single(α; only_estimate_D0)
    flush(stdout)
    full_time = @elapsed @time "Run run_single" data = run_single(α; only_estimate_D0)
    flush(stdout)
    return data
end

# Benchmark run_map
function benchmark_map(
    αs,
    executor = ThreadedEx(basesize = 1);
    only_estimate_D0 = true,
    threaded = false,
)
    @time "Compile run_map" run_map(αs, executor; only_estimate_D0, threaded)
    flush(stdout)
    full_time = @elapsed @time "Run run_map" data =
        run_map(αs, executor; only_estimate_D0, threaded)
    flush(stdout)
    return data
end

pool = create_workers(verbose = true)
@everywhere begin
    using Pkg
    Pkg.activate(".")
end
@everywhere begin
    using Arblib, HighestCuspedWave
    setprecision(Arb, 100)
end

i = if length(ARGS) > 0
    parse(Int, ARGS[1])
else
    6
end
m = if length(ARGS) > 1
    parse(Int, ARGS[2])
else
    64
end
αs = HighestCuspedWave.proof_interval_subdivisions_mince(i, m)

@show i m
println("Starting computations")
flush(stdout)

#benchmark_single(αs[2])

#data = benchmark_map(αs, DistributedEx(basesize = 1), threaded = false)
data = benchmark_map(
    αs,
    DistributedEx(basesize = 1),
    threaded = true,
    only_estimate_D0 = false,
)
#data = benchmark_map(αs, DistributedEx())

#@time data = run_map(αs, DistributedEx(basesize = 1), threaded = true, only_estimate_D0 = false)

HighestCuspedWave.write_proof_data(
    "PDC/data/benchmark.csv",
    HighestCuspedWave.add_rounded_data(data),
)
