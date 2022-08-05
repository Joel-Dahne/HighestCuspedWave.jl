"""
    prove(α::Arb; only_estimate_D0, threaded, verbose, extra_verbose)

Prove the inequality for showing existence of fixed-point.

Depending on the value of `α` it constructs either a
[`Fractionalkdvansatz`](@ref) or a [`KdVZeroAnsatz`](@ref) and then
uses the [`prove`](@ref) methods corresponding to them.

# Arguments
- `M::Integer`: the number of terms to use in asymptotic expansions
  during the computations, only used for
  [`Fractionalkdvansatz`](@ref).
- `only_estimate_D0::Bool = false`: if true it doesn't attempt to
  prove the bound for `D₀` but only uses an estimate. This doesn't
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
    only_estimate_D0 = false,
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

    proof_data = prove(u0; M, only_estimate_D0, threaded, verbose, extra_verbose)

    return (p = p, proof_data..., u0_time = u0_time)
end

"""
    round_for_publishing(n₀, δ₀, D₀; sigdigits = 10)

Convert `n₀, δ₀, D₀` to `Float64`, rounding up to the prescribed
number of significant digits, and check that the inequality `δ₀ <= (1
- D₀)^2 / 4n₀` holds for the rounded values as well.

This is used to get upper bounds of the values in a simpler format
than the `Arb` type.
"""
function round_for_publishing(n₀::Arb, δ₀::Arb, D₀::Arb; sigdigits = 10)
    inequality_holds_before = D₀ < 1 && δ₀ < (1 - D₀)^2 / 4n₀

    n₀_float = Arblib.get_d(ubound(n₀), RoundUp)
    δ₀_float = Arblib.get_d(ubound(δ₀), RoundUp)
    D₀_float = Arblib.get_d(ubound(D₀), RoundUp)

    # Check that the inequality holds before rounding. Conversion to
    # Float64 loses precision so this is not guaranteed.
    inequality_holds =
        D₀_float < 1 && Arb(δ₀_float) < (1 - Arb(D₀_float))^2 / 4Arb(n₀_float)

    if !inequality_holds
        inequality_holds_before &&
            @warn "Inequality holds before but not after conversion" n₀_float,
            δ₀_float,
            D₀_float
        return false, n₀_float, δ₀_float, D₀_float
    end

    n₀_float_rounded = round(n₀_float, RoundUp; sigdigits)
    δ₀_float_rounded = round(δ₀_float, RoundUp; sigdigits)
    D₀_float_rounded = round(D₀_float, RoundUp; sigdigits)

    @assert n₀ <= n₀_float_rounded
    @assert δ₀ <= δ₀_float_rounded
    @assert D₀ <= D₀_float_rounded

    # Check that the inequality holds after rounding
    inequality_holds =
        D₀_float_rounded < 1 &&
        Arb(δ₀_float_rounded) < (1 - Arb(D₀_float_rounded))^2 / 4Arb(n₀_float_rounded)

    return inequality_holds, n₀_float_rounded, δ₀_float_rounded, D₀_float_rounded
end
