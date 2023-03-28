# File used when precompiling a sysimage using compile.jl

using HighestCuspedWave

αs = [
    HighestCuspedWave.proof_interval_subdivisions_mince(15, thin = true)[1]
    HighestCuspedWave.proof_interval_subdivisions_mince(20, thin = true)[1]
]

HighestCuspedWave.prove(αs, only_estimate_D0 = false)

# This avoid compilation during construction for α close to -1
FractionalKdVAnsatz(HighestCuspedWave.proof_interval_subdivisions_mince(3)[end])
