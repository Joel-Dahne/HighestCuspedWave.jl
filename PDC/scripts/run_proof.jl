using ClusterManagers, Dates, Distributed, Folds

include("helper.jl")

function run_proof(
    i,
    m;
    executor = DistributedEx(basesize = 1),
    only_estimate_D0 = false,
    threaded = true,
    save = true,
    dirname = "PDC/data/proof",
)
    αs = HighestCuspedWave.proof_interval_subdivisions_mince(i, m)

    data = HighestCuspedWave.prove(αs; executor, only_estimate_D0, threaded)

    if save
        (αₗ, αᵤ), n = HighestCuspedWave.proof_interval_subdivisions(i)

        m = ifelse(isnothing(m), n, m)

        filename = "proof_$(αₗ)_$(αᵤ)_$(m)_$(n).csv"

        HighestCuspedWave.write_proof_data(
            joinpath(dirname, filename),
            HighestCuspedWave.add_rounded_data(data),
        )
    end
end

pool = create_workers(verbose = true)
@everywhere begin
    using Arblib, HighestCuspedWave
    setprecision(Arb, 100)
end

start, stop, m = read_args()

dirname = "PDC/data/proof/proof-$(round(Dates.now(), Second))"

@show start stop m dirname
println("Starting computations")
flush(stdout)

mkpath(dirname)

for subdivision = start:stop
    @time "subdivision = $subdivision" run_proof(subdivision, m; dirname)
    flush(stdout)
end

println("Finished computations")
flush(stdout)
