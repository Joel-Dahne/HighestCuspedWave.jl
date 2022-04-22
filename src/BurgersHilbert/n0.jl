"""
    n0_bound(u0::BHAnsatz; M::Integer, rtol, verbose)

Enclose the value of `n₀`. This is the supremum of
```
N(x) = u0.w(x) / 2u0(x)
```
for `0 < x < π`.

It uses an asymptotic expansion with `M` terms close to zero and ball
arithmetic on the remaining.

In practice the maximum is attained at `x = π` so the function is
evaluated there first. Then starting at `ϵ = 1 / 2` we halve `ϵ` until
the bound at `x = ϵ` using the asymptotic evaluation is smaller than
the value at `π`. We the prove that the value on ``[0, ϵ]`` is bounded
by the value at `x = π` using [`bounded_by`](@ref).

On the interval `[ϵ, π]` we bound it using ball arithmetic with
[`ArbExtras.maximum_enclosure`](@ref). Notice that we do not have to
prove that the maximum is attained at `x = π`, it's only used to make
the procedure more efficient.

# Arguments
- `M::Integer = 3`: Determines the number of terms used for the asymptotic
  expansion.
- `rtol = 1e-5`: Determines the relative tolerance used on the interval ``[ϵ,
  π]``.
- `degree::Integer = 0`: Degree to use in
  [`ArbExtras.maximum_enclosure`](@ref), in practice `degree = 0` is
  enough to pick up monotonicity of the function.
- `return_details::Bool = false`: If true it returns some more details
  about the computations. It returns the value that was used for `ϵ`.
- `verbose::Bool = false`: If true it prints more information about
  the process.
"""
function n0_bound(
    u0::BHAnsatz{Arb};
    M::Integer = 3,
    rtol = 1e-5,
    degree::Integer = 0,
    return_details::Bool = false,
    verbose::Bool = false,
)
    # For non-asymptotic evaluation
    f = x -> u0.w(x) / (2u0(x))
    # For asymptotic evaluation
    inv_u0 = inv_u0_normalised(u0, ϵ = Arb(0.6); M)
    g = x -> begin
        # Enclosure of -sqrt(log(1 + inv(x))) / 2log(x)
        factor = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x) && x < 1
            lower = zero(x)
            upper = let xᵤ = ubound(Arb, x)
                -sqrt(log(1 + inv(xᵤ))) / 2log(xᵤ)
            end
            Arb((lower, upper))
        else
            -sqrt(log(1 + inv(x))) / 2log(x)
        end

        factor * inv_u0(x)
    end

    # The maximum is in practice attained at x = π
    fπ = f(Arb(π))

    # Find ϵ such that g(ϵ) < fπ
    ϵ = Arb(1 // 2)
    while !(g(ϵ) <= fπ)
        ϵ /= 2
    end

    verbose && @info "ϵ was determined to be" ϵ

    # Prove the bound on [0, ϵ]
    bounded = ArbExtras.bounded_by(
        g,
        Arf(0),
        ubound(ϵ),
        lbound(fπ),
        abs_value = true,
        degree = -1,
        threaded = true;
        verbose,
    )

    if !bounded
        verbose && @error "Could not prove bound on [0, ϵ]"
        if return_details
            return Arblib.indeterminate!(zero(Arb)), ϵ
        else
            return Arblib.indeterminate!(zero(Arb))
        end
    end

    # Bound the value on [ϵ, π]
    m = ArbExtras.maximum_enclosure(
        f,
        lbound(ϵ),
        ubound(Arb(π)),
        abs_value = true,
        point_value_max = fπ,
        threaded = true;
        degree,
        rtol,
        verbose,
    )

    if return_details
        return m, ϵ
    else
        return m
    end
end
