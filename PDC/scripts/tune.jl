using ClusterManagers, Dates, Distributed, OrderedCollections, Folds

include("helper.jl")

function run(
    i,
    m;
    executor = DistributedEx(basesize = 1),
    only_estimate_D0 = false,
    threaded = true,
    save = true,
    dirname = "PDC/data/tune",
)
    αs = HighestCuspedWave.proof_interval_subdivisions_mince(i, m)

    t1 = @elapsed data = HighestCuspedWave.prove(αs; executor, only_estimate_D0, threaded)

    # The elapsed time will be the time for the slower process to
    # finish. Using a higher m will lead to better load balancing and
    # it is therefore also interesting to look at the average time for
    # each computation.
    t2 =
        (
            sum(data.n₀_time) +
            sum(data.δ₀_time) +
            sum(filter(isfinite, data.D₀_time)) +
            sum(data.D₀_estimate_time) +
            sum(data.u0_time)
        ) / m

    if save
        (αₗ, αᵤ), n = HighestCuspedWave.proof_interval_subdivisions(i)

        filename = "tune_$(αₗ)_$(αᵤ)_$(m)_$(n).csv"

        HighestCuspedWave.write_proof_data(
            joinpath(dirname, filename),
            HighestCuspedWave.add_rounded_data(data),
        )
    end

    return t1, t2
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

start, stop, m = read_args()
if isnothing(m)
    m = nworkers() # We want a different default value for m
end

dirname = "PDC/data/tune/tune-$(round(Dates.now(), Second))"

@show start stop m dirname
println("Starting computations")
flush(stdout)

mkpath(dirname)

t1s = OrderedDict{Int,Float64}()
t2s = empty(t1s)
t1s_estimate = empty(t1s)
t2s_estimate = empty(t1s)
for subdivision = start:stop
    t1s[subdivision], t2s[subdivision] = run(subdivision, m; dirname)
    t1s_estimate[subdivision] =
        HighestCuspedWave.proof_interval_subdivisions()[subdivision][2] / m *
        t1s[subdivision]
    t2s_estimate[subdivision] =
        HighestCuspedWave.proof_interval_subdivisions()[subdivision][2] / m *
        t2s[subdivision]
    @show subdivision t1s[subdivision] t1s_estimate[subdivision] t2s[subdivision] t2s_estimate[subdivision]
    flush(stdout)
end

display(t1s)
display(t2s)
display(t1s_estimate)
display(t2s_estimate)
