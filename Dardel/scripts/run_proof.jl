using Dates

include("helper.jl")

function run_proof(
    i,
    m;
    parallelization = :distributed,
    only_estimate_D0 = false,
    threaded = true,
    verbose = false,
    extra_verbose = false,
    save = true,
    dirname = "Dardel/data/proof",
)
    αs = HighestCuspedWave.proof_interval_subdivisions_mince(i, m)

    data = HighestCuspedWave.prove(
        αs;
        parallelization,
        only_estimate_D0,
        threaded,
        verbose,
        extra_verbose,
    )

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

    # Set logging to always flush
    using Logging: global_logger
    using TerminalLoggers: TerminalLogger
    global_logger(TerminalLogger(always_flush = true))
end

start, stop, m = read_args()

dirname = "Dardel/data/proof/proof-$(round(Dates.now(), Second))"
logfile = joinpath(dirname, "log_$(start)_$(stop)_$(m)")

@info "Determined arguments" start stop m dirname logfile
@info "Starting computations"

mkpath(dirname)
open(logfile, truncate = true) do io
    println(io, "$(round(Dates.now(), Second)): Starting computations")
end

for subdivision = start:stop
    time =
        @elapsed run_proof(subdivision, m; dirname, verbose = false, extra_verbose = false)

    @info "Subdivision $subdivision took $time seconds"

    estimate_total_time = if isnothing(m)
        # Running whole interval so the time is the estimated time
        time
    else
        time * HighestCuspedWave.proof_interval_subdivisions(subdivision)[2] / m
    end


    open(logfile, append = true) do io
        println(
            io,
            "$(round(Dates.now(), Second)): Finished $subdivision in $(time)s" *
            " - estimate total time $(estimate_total_time)s",
        )
    end
end

@info "Finished computations"
open(logfile, append = true) do io
    println(io, "$(round(Dates.now(), Second)): Finished computations")
end
