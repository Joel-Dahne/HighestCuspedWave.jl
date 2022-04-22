"""
    prove(α::Arb; only_estimate_CB, threaded, verbose, extra_verbose)

Prove the inequality for showing existence of fixed-point.

Depending on the value of `α` it constructs either a
[`Fractionalkdvansatz`](@ref) or a [`KdVZeroAnsatz`](@ref) and then
uses the [`prove`](@ref) methods corresponding to them.

# Arguments
- `M::Integer`: the number of terms to use in asymptotic expansions
  during the computations, only used for
  [`Fractionalkdvansatz`](@ref).
- `only_estimate_CB::Bool = false`: if true it doesn't attempt to
  prove the bound for `C_B` but only uses an estimate. This doesn't
  give a rigorous proof but is useful if you only want to determine of
  the bound seems to hold.
- `threaded::Bool = true`: determines if it uses multiple threads for the
  computations or only a single thread.
- `verbose::Bool = false`: if true it prints information about the
  progress.
- `extra_verbose::Bool = false`: if true it prints even more
  information about the progress.
"""
function prove(
    α::Arb;
    M = 5,
    only_estimate_CB = false,
    threaded = true,
    verbose = false,
    extra_verbose = false,
)
    if α < -1 // 6
        u0_time = @elapsed u0 = FractionalKdVAnsatz(α)
        p = u0.p
    else
        u0_time = @elapsed u0 = KdVZeroAnsatz(α, midpoint(Arb, α))
        p = one(α)
    end

    verbose && @info "Constructed u0" u0 u0_time

    proof_data = prove(u0; M, only_estimate_CB, threaded, verbose, extra_verbose)

    return (p = p, proof_data..., u0_time = u0_time)
end

"""
    round_for_publishing(n₀, δ₀, C_B; sigdigits = 10)

Convert `n₀, δ₀, C_B` to `Float64`, rounding up to the prescribed
number of significant digits, and check that the inequality `δ₀ <= (1
- C_B)^2 / 4n₀` holds for the rounded values as well.

This is used to get upper bounds of the values in a simpler format
than the `Arb` type.
"""
function round_for_publishing(n₀::Arb, δ₀::Arb, C_B::Arb; sigdigits = 10)
    n₀_float = Arblib.get_d(ubound(n₀), RoundUp)
    δ₀_float = Arblib.get_d(ubound(δ₀), RoundUp)
    C_B_float = Arblib.get_d(ubound(C_B), RoundUp)

    # Check that the inequality holds before rounding. Conversion to
    # Float64 loses precision so this is not guaranteed.
    inequality_holds = Arb(δ₀_float) < (1 - Arb(C_B_float))^2 / 4Arb(n₀_float)

    if !inequality_holds
        @warn "Inequality doesn't hold after conversion" n₀_float, δ₀_float, C_B_float
        return false, n₀_float, δ₀_float, C_B_float
    end

    n₀_float_rounded = round(n₀_float, RoundUp; sigdigits)
    δ₀_float_rounded = round(δ₀_float, RoundUp; sigdigits)
    C_B_float_rounded = round(C_B_float, RoundUp; sigdigits)

    @assert n₀ <= n₀_float_rounded
    @assert δ₀ <= δ₀_float_rounded
    @assert C_B <= C_B_float_rounded

    # Check that the inequality holds after rounding
    inequality_holds =
        Arb(δ₀_float_rounded) < (1 - Arb(C_B_float_rounded))^2 / 4Arb(n₀_float_rounded)

    return inequality_holds, n₀_float_rounded, δ₀_float_rounded, C_B_float_rounded
end
