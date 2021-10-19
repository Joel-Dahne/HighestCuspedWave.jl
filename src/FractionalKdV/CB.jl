"""
    CB(u0::FractionalKdVAnsatz{arb}; ϵ::arb = parent(u0.α)(0.1))
Upper bound the value of C_B from the paper. TODO: Uses an asymptotic
expansion on the interval [0, ϵ] and ball arithmetic on [ϵ, π].
"""
function CB(
    u0::FractionalKdVAnsatz{arb};
    ϵ::arb = parent(u0.α)(0.1),
    rtol::arb = parent(u0.α)(1e-2),
    show_trace = false,
)
    # Bound the value one [0, ϵ]
    # TODO: Implement this
    n1 = enclosemaximum(
        T0(u0, Asymptotic()),
        parent(u0.α)(0),
        ϵ,
        rtol = rtol,
        absmax = true,
        maxevals = 10^3,
        show_trace = show_trace,
    )

    # Bound the value one [ϵ, π]
    # TODO: This does not fully work yet
    tol = 1e-4 * Float64(rtol)
    n2 = enclosemaximum(
        T0(u0, Ball(), rtol = tol, atol = tol),
        ϵ,
        parent(u0.α)(π),
        rtol = rtol,
        absmax = true,
        maxevals = 10^3,
        show_trace = show_trace,
    )

    return max(n1, n2)
end

"""
    CB_bounded_by(u0::AbstractAnsatz{Arb}, C::Arf; ϵ = 0, ...)

Attempt to prove that \$C_B\$ is bounded by `C`. Returns `true` on
success and `false` on failure.

If `ϵ = 0` it tries to determine an optimal choice of `ϵ` by starting
with `ϵ = 1` and then iteratively decreasing it until the bound
`T0(u0, Asymptotic())(ϵ) < C` holds.

It then uses the asymptotic version on the interval `[0, ϵ]` and the
non-asymptotic version on `[ϵ, π]`.
"""
function CB_bounded_by(
    u0::FractionalKdVAnsatz{Arb},
    C::Arf;
    ϵ::Arf = zero(Arf),
    threaded = false,
    verbose = false,
)
    # Determine ϵ such that the bound holds when evaluated at x = ϵ
    # with the asymptotic version.
    f = T0(u0, Asymptotic())

    if iszero(ϵ)
        ϵ = one(ϵ)
        while !(f(Arb(ϵ)) < C)
            # We use the factor 0.8 instead of, say, 0.5 to get
            # slightly better (higher) values for ϵ.
            ϵ *= 0.8
        end

        verbose && @show ϵ
    end


    # Check that the bound holds on [0, ϵ]
    asymptotic_bound = ArbExtras.bounded_by(f, Arf(0), ϵ, C, degree = -1; threaded, verbose)

    if !asymptotic_bound
        verbose && @info "Bound doesn't hold on [0, ϵ]"
        return false
    end

    g = T0(u0, Ball(), δ2 = Arf(1e-1), skip_div_u0 = true)

    h(x) = begin
        res = g(x)

        isfinite(res) || return res

        # Compute a tighter enclosure of u0 by using evaluation with
        # ArbSeries. In practice this will pick up that u0 is monotone
        # on the interval and therefore very efficient.
        u0x =
            Arb(ArbExtras.extrema_series(u0, Arblib.getinterval(x)..., degree = 1)[1:2])

        res /= u0x

        return res
    end

    # Check that the bound holds on [ϵ, π]
    non_asymptotic_bound = ArbExtras.bounded_by(
        h,
        ϵ,
        ubound(Arb(π)),
        C,
        degree = -1,
        depth_start = 8;
        maxevals = 1000,
        threaded,
        verbose,
    )

    if !non_asymptotic_bound
        verbose && @info "Bound doesn't hold on [ϵ, π]"
        return false
    end

    return true
end
