"""
    delta0_bound(u0::BHAnsatz{Arb}; rtol, degree, verbose)

Enclose the value of `δ₀` from the paper.

The interval `[0, π]` is first split into two intervals, `[0, ϵ]` and
`[ϵ, π]`. On the first interval we use an asymptotic expansion and on
the second we use direct ball evaluation.

The interval `[ϵ, π]` uses straight forward ball enclosures.

The interval `[0, ϵ]` requires more care because the maximum is
attained very close to the origin, around `1e-10000` and a naive
version takes a very long time. Instead we split the interval into
four smaller parts, `[0, ϵ1]`, `[ϵ1, ϵ2]`, `[ϵ2, ϵ3]` and `[ϵ3, ϵ]`.
The value of `ϵ1` is take to be small enough so that evaluating the
function on the whole interval gives a sufficiently small enclosure.
The difference between the handling of the intervals `[ϵ1, ϵ2]`, `[ϵ2,
ϵ3]` and `[ϵ3, ϵ]` are how many terms in the asymptotic expansion we
neglect (i.e. add to the error term). The value of `ϵ2` is taken
sufficiently small to be able to neglect all terms with an exponent
larger than or equal to `1e-4`, which (with the standard choice of
`u0`) skips all but 5 terms. The value of `ϵ3` is taken sufficiently
small to be able to neglect all terms with an exponent larger than or
equal to `1 // 4`. On the interval `[ϵ3, ϵ]` we keep all the
exponents.

**IMPROVE:** We could do the asymptotic part of the calculations at a
reduced precision. However this seems to only give a marginal
improvement. It would be more important to reduce the number of
allocations in that case.

# Arguments
- `rtol = 1e-3`: Determines the relative tolerance used when enclosing
  the maximum.
- `degree = 4`: Determines the degree used when enclosing the maximum.
- `return_subresults`: If true it returns not only the maximum on the
  interval but the maximum on each of the intervals `[0, ϵ1]`, `[ϵ1,
  ϵ2]`, `[ϵ2, ϵ3]`, `[ϵ3, ϵ]` and `[ϵ, π]`.
- `verbose = false`: If true it prints more information about the
  progress.
"""
function delta0_bound(
    u0::BHAnsatz{Arb};
    rtol = Arb(1e-3),
    degree = 4,
    return_subresults = false,
    verbose = false,
)
    ϵ = Arblib.midpoint(Arb("1e-1"))

    # Bound the value on [0, ϵ] using an asymptotic expansion

    # Split the interval into four pieces, [0, ϵ1], [ϵ1, ϵ2], [ϵ2, ϵ3]
    # and [ϵ3, ϵ]
    ϵ1 = Arblib.midpoint(Arb("1e-10000000"))
    ϵ2 = Arblib.midpoint(Arb("1e-100000"))
    ϵ3 = Arblib.midpoint(Arb("1e-100"))

    # Define different functions for the different intervals, f1 is
    # used in the first 2, f2 in the third and f3 in the fourth.
    f1 = F0(u0, Asymptotic(), exponent_limit = Arb(1e-4), ϵ = 2Arb(ϵ2))
    f2 = F0(u0, Asymptotic(), exponent_limit = Arb(1 // 4), ϵ = 2Arb(ϵ3))
    f3 = F0(u0, Asymptotic(), exponent_limit = nothing, ϵ = 2Arb(ϵ), M = 5)

    verbose &&
        @info "Terms in expansion on [0, ϵ2]" length(f1.u0_expansion_div_xlogx) length(
            f1.Du0_expansion_div_x2logx,
        )
    verbose &&
        @info "Terms in expansion on [ϵ2, ϵ3]" length(f2.u0_expansion_div_xlogx) length(
            f2.Du0_expansion_div_x2logx,
        )
    verbose &&
        @info "Terms in expansion on [ϵ3, ϵ]" length(f3.u0_expansion_div_xlogx) length(
            f3.Du0_expansion_div_x2logx,
        )

    m11 = f1(Arb((0, ϵ1)))

    verbose && @info "Enclosure on [0, ϵ1]" m11

    m12 = ArbExtras.maximum_enclosure(
        f1,
        ϵ1,
        ϵ2,
        abs_value = true,
        log_bisection = true,
        depth_start = 24,
        threaded = true,
        maxevals = typemax(Int),
        depth = typemax(Int);
        degree,
        rtol,
        verbose,
    )

    verbose && @info "Enclosure on [ϵ1, ϵ2]" m12

    m13 = ArbExtras.maximum_enclosure(
        f2,
        ϵ2,
        ϵ3,
        abs_value = true,
        log_bisection = true,
        depth_start = 17,
        threaded = true,
        maxevals = typemax(Int),
        depth = typemax(Int);
        degree,
        rtol,
        verbose,
    )

    verbose && @info "Enclosure on [ϵ2, ϵ3]" m13

    m14 = ArbExtras.maximum_enclosure(
        f3,
        ϵ3,
        ϵ,
        abs_value = true,
        depth_start = 10,
        log_bisection = true,
        threaded = true,
        maxevals = typemax(Int),
        depth = typemax(Int);
        degree,
        rtol,
        verbose,
    )

    verbose && @info "Enclosure on [ϵ3, ϵ]" m14

    m1 = max(m11, m12, m13, m14)

    verbose && @info "Enclosure on [0, ϵ]" m1

    # Bound the value on [ϵ, π]
    m2 = ArbExtras.maximum_enclosure(
        F0(u0),
        ϵ,
        Arblib.ubound(Arb(π)),
        abs_value = true,
        depth_start = 7,
        threaded = true,
        maxevals = 5000;
        degree,
        rtol,
        verbose,
    )

    verbose && @info "Enclosure on [ϵ, π]" m2

    if return_subresults
        return max(m1, m2), m11, m12, m13, m14, m2
    else
        return max(m1, m2)
    end
end
