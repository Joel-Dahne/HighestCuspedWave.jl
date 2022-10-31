# File used when precompiling a sysimage using compile.jl

using HighestCuspedWave

αs = [
    HighestCuspedWave.proof_interval_subdivisions_mince(15, thin = true)[1]
    HighestCuspedWave.proof_interval_subdivisions_mince(20, thin = true)[1]
]

HighestCuspedWave.prove(αs, only_estimate_D0 = false)
