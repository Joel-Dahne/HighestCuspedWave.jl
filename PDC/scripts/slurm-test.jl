using Distributed, ClusterManagers

include("helper.jl")

print_slurm_details() = display(filter(startswith("SLURM") ∘ first, ENV))

function test_output()
    for i = 1:5
        println("Hej $i")
        flush(stdout)
        sleep(1)
        println("Då $i")
        flush(stdout)
    end
end

function test_progress_logging(N = 1000)
    xf = Map() do x
        sleep(10 / N)
        x
    end

    dcollect(xf, withprogress(1:N; interval = 1), basesize = 1)
end

#print_slurm_details()

create_workers(verbose = true)

#test_output()

test_progress_logging()
