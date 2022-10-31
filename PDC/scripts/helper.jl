function create_workers(
    nw = nw = parse(Int, get(ENV, "SLURM_NTASKS", "1")),
    nt = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1"));
    use_slurm = haskey(ENV, "SLURM_JOB_ID"),
    verbose = false,
)

    if verbose
        println("nw = $nw")
        println("nt = $nt")
        flush(stdout)
    end

    # Launch worker processes
    ENV["JULIA_PROJECT"] = Base.active_project()
    ENV["JULIA_NUM_THREADS"] = nt
    ENV["JULIA_WORKER_TIMEOUT"] = 300

    if use_slurm
        # Give the current sysimage explicitly. The SlurmManager
        # doesn't handle it by itself.
        sysimage = unsafe_string(Base.JLOptions().image_file)
        addprocs(SlurmManager(nw), exeflags = "--sysimage=$sysimage")
    else
        addprocs(nw)
    end

    if verbose
        println("Number of processes: ", nprocs())
        println("Number of workers: ", nworkers())
        flush(stdout)

        # For each worker gets its id, process id, hostname and number of threads
        calls = Dict(
            i => @spawnat(i, (myid(), getpid(), gethostname(), Threads.nthreads())) for
            i in procs()
        )
        for i in procs()
            id, pid, host, threads = fetch(calls[i])
            println(id, " ", pid, " ", host, " ", threads)
            flush(stdout)
        end
    end

    return workers()
end

function read_args()
    start = if length(ARGS) > 0
        parse(Int, ARGS[1])
    else
        1
    end

    stop = if length(ARGS) > 1
        parse(Int, ARGS[2])
    else
        length(HighestCuspedWave.proof_interval_subdivisions())
    end

    m = if length(ARGS) > 2
        parse(Int, ARGS[3])
    else
        nothing # This means take all
    end

    return start, stop, m
end
