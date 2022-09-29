using Distributed, ClusterManagers

function print_slurm_details()
    display(filter(startswith("SLURM") ∘ first, ENV))
    #println("SLURM_JOB_NUM_NODES = " * get(ENV, "SLURM_JOB_NUM_NODES", ""))
    #println("SLURM_JOB_NODELIST = " * get(ENV, "SLURM_JOB_NODELIST", ""))
    #println("SLURM_JOB_CPUS_PER_NODE = " * get(ENV, "SLURM_JOB_CPUS_PER_NODE", ""))
    #println("SLURM_MEM_PER_NODE = " * get(ENV, "SLURM_MEM_PER_NODE", ""))
    #println("SLURM_MEN_PER_CPU = " * get(ENV, "SLURM_MEM_PER_CPU", ""))
    #println("SLURM_NTASKS = " * get(ENV, "SLURM_NTASKS", ""))
    #println("SLURM_NTASKS_PER_NODE = " * get(ENV, "SLURM_NTASKS_PER_NODE", ""))
    #println("SLURM_NODEFILE = " * get(ENV, "SLURM_NODEFILE", ""))
end

function test_workers()
    ENV["JULIA_NUM_THREADS"] =
        parse(Int, split(get(ENV, "SLURM_JOB_CPUS_PER_NODE", "1(x1)"), "(")[1])

    # launch worker processes
    # Look at https://stackoverflow.com/questions/48631411/julia-and-slurm-setup ?
    addprocs(parse(Int, get(ENV, "SLURM_JOB_NUM_NODES", "1")) - 1)

    println("Number of processes: ", nprocs())
    println("Number of workers: ", nworkers())

    # each worker gets its id, process id and hostname
    for i in procs()
        id, pid, host, threads, wtimeout = fetch(
            @spawnat i (
                myid(),
                getpid(),
                gethostname(),
                Threads.nthreads(),
                Distributed.worker_timeout(),
            )
        )
        println(id, " ", pid, " ", host, " ", threads, " ", wtimeout)
    end
end

function test_workers_slurm()
    #nt = parse(Int, split(get(ENV, "SLURM_JOB_CPUS_PER_NODE", "1(x1)"), "(")[1])
    nt = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1"))
    ENV["JULIA_NUM_THREADS"] = nt

    println("Number of threads: ", ENV["JULIA_NUM_THREADS"])

    ENV["JULIA_WORKER_TIMEOUT"] = 300

    # launch worker processes
    # Look at https://stackoverflow.com/questions/48631411/julia-and-slurm-setup ?
    #nw = parse(Int, get(ENV, "SLURM_JOB_NUM_NODES", "1"))
    nw = parse(Int, get(ENV, "SLURM_NTASKS", "1"))
    addprocs(SlurmManager(nw), topology = :master_worker)

    println("Number of processes: ", nprocs())
    println("Number of workers: ", nworkers())

    # each worker gets its id, process id and hostname
    for i in procs()
        id, pid, host, threads, wtimeout = fetch(
            @spawnat i (
                myid(),
                getpid(),
                gethostname(),
                Threads.nthreads(),
                Distributed.worker_timeout(),
            )
        )
        println(id, " ", pid, " ", host, " ", threads, " ", wtimeout)
    end
end

function test_output()
    for i = 1:5
        println("Hej $i")
        flush(stdout)
        sleep(1)
        println("Då $i")
        flush(stdout)
    end
end

print_slurm_details()
test_workers_slurm()
test_output()
