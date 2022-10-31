# Compile a sysimage to reduce compilation costs later
using PackageCompiler

create_sysimage(
    sysimage_path = joinpath(ENV["HOME"], ".julia/sysimages/HighestCuspedWave.so"),
    precompile_execution_file = "PDC/scripts/to_compile.jl",
)
