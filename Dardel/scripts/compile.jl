# Compile a sysimage to reduce compilation costs later
using PackageCompiler

sysimage_dir = joinpath(ENV["HOME"], ".julia/sysimages/")
sysimage_path = joinpath(sysimage_dir, "HighestCuspedWave.so")
precompile_execution_file = "Dardel/scripts/to_compile.jl"

mkpath(sysimage_dir)

create_sysimage(; sysimage_path, precompile_execution_file)
