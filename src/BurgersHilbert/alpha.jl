"""
    alpha0(u0::BHAnsatz; M::Integer, rtol, verbose)

Enclose the value of `α₀`, given by the maximum of `u0.w(x) / 2u0(x)`
on ``[0, π]``.

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
function alpha0(
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
    # We compute the asymptotic expansion of u0 at zero and then
    # manually cancel factors between the weight and the expansion.
    # Finally care has to be take to make sure that the weight with
    # the cancelled factors can be evaluated around zero, this is done
    # using monotonicity.
    g = let u0_expansion = u0(Arb(0.6), AsymptoticExpansion(); M)
        # Enclosure of -sqrt(log(1 + inv(x))) / 2log(x)
        # PROVE: That this is monotonically increasing on [0, 1]
        factor1(x) =
            if iszero(x)
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

        # Divide all terms in the expansion by abs(x) * log(abs(x))
        u0_expansion_div_xlogx = empty(u0_expansion)
        for ((i, m, k, l), value) in u0_expansion
            u0_expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
        end

        # Enclosure of -x * log(x) / u0(x)
        factor2(x) = -inv(eval_expansion(u0, u0_expansion_div_xlogx, x))

        x -> factor1(x) * factor2(x)
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
