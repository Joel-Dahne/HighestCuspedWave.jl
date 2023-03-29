# Compile a sysimage to reduce compilation costs later
using PackageCompiler

sysimage_path = joinpath(ENV["HOME"], ".julia/sysimages/HighestCuspedWave.so")
precompile_execution_file = "Dardel/scripts/to_compile.jl"

mkpath(sysimage_path)

create_sysimage(; sysimage_path, precompile_execution_file)
