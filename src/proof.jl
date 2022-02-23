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
    format_for_publishing(α₀, δ₀, C_B)

Convert `α₀, δ₀, C_B` to `Float64` rounding up and check that the
inequality `δ₀ <= (1 - C_B)^2 / 4α₀` holds for the `Float64` values as
well.

This is used to get upper bounds of the values in a simpler format
than the `Arb` type.
"""
function format_for_publishing(α₀::Arb, δ₀::Arb, C_B::Arb)
    α₀_float = Arblib.get_d(ubound(α₀), RoundUp)
    δ₀_float = Arblib.get_d(ubound(δ₀), RoundUp)
    C_B_float = Arblib.get_d(ubound(C_B), RoundUp)

    # Check that the inequality holds
    inequality_holds = Arb(δ₀_float) <= (1 - Arb(C_B_float))^2 / 4Arb(α₀)

    return inequality_holds, α₀_float, δ₀_float, C_B_float
end
