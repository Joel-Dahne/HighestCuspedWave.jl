# Very basic implementation of fmpz_struct used for a few calls to Arb
# functions

mutable struct fmpz_struct
    d::Int

    function fmpz_struct()
        z = new()
        ccall(Arblib.@libflint(fmpz_init), Nothing, (Ref{fmpz_struct},), z)
        finalizer(fmpz_clear!, z)
        return z
    end

    function fmpz_struct(x::Int)
        z = new()
        ccall(Arblib.@libflint(fmpz_init_set_si), Nothing, (Ref{fmpz_struct}, Int), z, x)
        finalizer(fmpz_clear!, z)
        return z
    end
end

fmpz_clear!(x::fmpz_struct) =
    ccall(Arblib.@libflint(fmpz_clear), Nothing, (Ref{fmpz_struct},), x)

function Int(x::fmpz_struct)
    # Check x <= typemax(Int)
    ccall(Arblib.@libflint(fmpz_cmp_si), Cint, (Ref{fmpz_struct}, Int), x, typemax(Int)) <=
    0 || throw(InexactError(nameof(Int), Int, x))
    # Check x >= typemin(Int)
    ccall(Arblib.@libflint(fmpz_cmp_si), Cint, (Ref{fmpz_struct}, Int), x, typemin(Int)) >=
    0 || throw(InexactError(nameof(Int), Int, x))
    return ccall(Arblib.@libflint(fmpz_get_si), Int, (Ref{fmpz_struct},), x)
end
