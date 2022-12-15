"""
    round_for_publishing(n₀::Arb, δ₀::Arb, D₀::Arb; sigdigits = nothing)

Convert `n₀, δ₀, D₀` to `Float64`, rounding up to the prescribed
number of significant digits, and check that the inequality `δ₀ <= (1
- D₀)^2 / 4n₀` holds for the rounded values as well.

This is used to get upper bounds of the values in a simpler format
than the `Arb` type.

The rounding to the prescribed number of significant digits is done
using [`round`](@ref) with `RoundUp`. However the implementation in
Julia does not take into account rounding errors and the result is
hence not guaranteed to be an upper bound. We therefore try to round
and then check if the result indeed is an upper bound. If it is not an
upper bound we return the non-rounded value.
"""
function round_for_publishing(n₀::Arb, δ₀::Arb, D₀::Arb; sigdigits = nothing)
    inequality_holds = D₀ < 1 && δ₀ < (1 - D₀)^2 / 4n₀

    n₀_float = Arblib.get_d(ubound(n₀), RoundUp)
    δ₀_float = Arblib.get_d(ubound(δ₀), RoundUp)
    D₀_float = Arblib.get_d(ubound(D₀), RoundUp)

    # Try to round up, if it fails return input without rounding
    try_round = x -> begin
        y = round(x, RoundUp; sigdigits)
        return x <= y ? y : x
    end

    n₀_rounded = try_round(n₀_float)
    δ₀_rounded = try_round(δ₀_float)
    D₀_rounded = try_round(D₀_float)

    @assert isnan(n₀) || n₀ <= n₀_rounded
    @assert isnan(δ₀) || δ₀ <= δ₀_rounded
    @assert isnan(D₀) || D₀ <= D₀_rounded

    # Check that the inequality holds before rounding. Conversion to
    # Float64 loses precision so this is not guaranteed.
    inequality_holds_float =
        D₀_float < 1 && Arb(δ₀_float) < (1 - Arb(D₀_float))^2 / 4Arb(n₀_float)

    # Check that the inequality holds after rounding.
    inequality_holds_rounded =
        D₀_rounded < 1 && Arb(δ₀_rounded) < (1 - Arb(D₀_rounded))^2 / 4Arb(n₀_rounded)

    if inequality_holds
        if !inequality_holds_float
            @warn "Inequality holds before but not after conversion to Float64" n₀_float,
            δ₀_float,
            D₀_float
        elseif !inequality_holds_rounded
            @warn "Inequality holds after conversion but not after rounding" n₀_rounded,
            δ₀_rounded,
            D₀_rounded
        end
    end

    return inequality_holds_rounded, n₀_float, δ₀_float, D₀_float
end

"""
    add_rounded_data(data)

Take the data returned from [`prove`](@ref) for `-1 < α < 0` and add
rounded bounds as given by [`round_for_publishing`](@ref).
"""
function add_rounded_data(
    data::@NamedTuple{
        α::Arb,
        p::Arb,
        proved::Bool,
        proved_estimate::Bool,
        n₀::Arb,
        δ₀::Arb,
        D₀_estimate::Arb,
        D₀::Arb,
        n₀_time::Float64,
        δ₀_time::Float64,
        D₀_estimate_time::Float64,
        D₀_time::Float64,
        u0_N0::Int,
        u0_N1::Int,
        u0_time::Float64,
        prec::Int,
    };
    sigdigits = 15,
)

    inequality_holds_rounded, n₀_rounded, δ₀_rounded, D₀_rounded =
        round_for_publishing(data.n₀, data.δ₀, data.D₀; sigdigits)

    # It is proved for the rounded values if it is proved for the
    # non-rounded values and the inequality holds after rounding.
    proved_rounded = data.proved && inequality_holds_rounded

    return (;
        data.α,
        data.p,
        data.proved,
        data.proved_estimate,
        proved_rounded,
        data.n₀,
        n₀_rounded,
        data.δ₀,
        δ₀_rounded,
        data.D₀_estimate,
        data.D₀,
        D₀_rounded,
        data.n₀_time,
        data.δ₀_time,
        data.D₀_estimate_time,
        data.D₀_time,
        data.u0_N0,
        data.u0_N1,
        data.u0_time,
        data.prec,
    )
end

function add_rounded_data(data::DataFrame; sigdigits = 15)
    DataFrame(map(row -> add_rounded_data(row; sigdigits), Tables.rowtable(data)))
end

"""
    write_proof_data(filename, data)

Write the data produced by [`add_rounded_data`](@ref) in a csv-file.

Most of the data is just directly saved in the file. The exception is
the `Arb` data for which we store a string representation, but since
this is lossy we also store a dump so that we can exactly recover it.

The data can be loaded using [`load_data`](@ref).
"""
function write_proof_data(filename, data::DataFrame)
    data = copy(data)

    for col_name in names(data)
        col = data[!, col_name]
        if eltype(col) == Arb
            insertcols!(
                data,
                col_name,
                col_name * "_dump" => Arblib.dump_string.(col),
                after = true,
            )
        end
    end

    CSV.write(filename, data)
end

"""
    read_proof_data(filename)

Read data from a csv-file stored using [`save_data`](@ref).
"""
function read_proof_data(filename; check = true)
    # We give the types explicitly for two reasons.
    # 1. When the Arb values happen to be exact they are parsed as
    # floating points instead of as strings
    # 2. There is a bug in CSV.jl which gives the wrong column type in
    # some cases https://github.com/JuliaData/CSV.jl/issues/1010
    types = [
        String,
        String,
        String,
        String,
        Bool,
        Bool,
        Bool,
        String,
        String,
        Float64,
        String,
        String,
        Float64,
        String,
        String,
        String,
        String,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Int64,
        Int64,
        Float64,
        Int64,
    ]

    data = CSV.read(filename, DataFrame; types)

    for col_name in names(data)
        if endswith(col_name, "_dump")
            target_col_name = replace(col_name, "_dump" => "")

            col = data[!, col_name]
            data[!, target_col_name] =
                Arblib.load_string!.([Arb(; prec) for prec in data.prec], col)
            select!(data, Not(col_name))
        end
    end

    if check
        check_proof_data(data)
    end

    return data
end

"""
    read_proof_data_dir(dirname)
    read_proof_data_dir(dirnames::AbstractVector)

Read data from all csv-files in the given directory using
[`read_proof_data_dir`](@ref) and concatenate them into one dataframe.
The rows are sorted by `α`.

If a vector of directory names are given then read from all files in
all of those directories.
"""
function read_proof_data_dir(dirname; check = true)
    filenames = filter(endswith(".csv"), readdir(dirname, join = true))
    datas = read_proof_data.(filenames; check)
    data = vcat(datas...)

    sort!(data, :α, by = α -> midpoint(α))

    return data
end

function read_proof_data_dir(dirnames::AbstractVector; check = true)
    datas = read_proof_data_dir.(dirnames; check)
    data = vcat(datas...)

    sort!(data, :α, by = α -> midpoint(α))

    return data
end

"""
    check_proof_data(data::DataFrame)

Perform a number of sanity checks on the proof data. Throws an error
in case something fails.
"""
function check_proof_data(data::DataFrame)
    for row in eachrow(data)
        # Sanity checks on data
        -1 < row.α < 0 || error("-1 < α < 0 not satisfied")
        0 < row.p <= 1 || error("0 < p <= 1 not satisfied")

        # Check Arb data
        n₀, δ₀, D₀, D₀_e = row.n₀, row.δ₀, row.D₀, row.D₀_estimate

        # The estimate should always give a lower bound
        D₀_e > D₀ && error("D₀ estimate larger than bound")

        if row.proved
            if !(D₀ < 1 && δ₀ < (1 - D₀)^2 / 4n₀)
                error("proved is set but inequality doesn't hold")
            end

            row.proved_estimate || error("proved is set but no proved_estimate")
        end

        # Check rounded data
        if row.proved && !row.proved_rounded
            @warn "proved is set but not proved_rounded"
        end
        if !row.proved && row.proved_rounded
            error("proved is not set but proved_rounded is")
        end

        n₀_r, δ₀_r, D₀_r = Arb(row.n₀_rounded), Arb(row.δ₀_rounded), Arb(row.D₀_rounded)
        isnan(n₀) || n₀ <= n₀_r || error("n₀ is not bounded by n₀_rounded")
        isnan(δ₀) || δ₀ <= δ₀_r || error("δ₀ is not bounded by δ₀_rounded")
        isnan(D₀) || D₀ <= D₀_r || error("D₀ is not bounded by D₀_rounded")
        if row.proved_rounded
            if !(D₀_r < 1 && δ₀_r < (1 - D₀_r)^2 / 4n₀_r)
                error("proved_rounded is set but inequality doesn't hold")
            end
        end
    end
end
