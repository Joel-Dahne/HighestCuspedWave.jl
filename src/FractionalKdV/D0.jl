"""
    D0_bounded_by(u0::AbstractAnsatz{Arb}, C::Arf; ϵ = 0, ...)

Attempt to prove that `D₀` is bounded by `C`. Returns `true` on
success and `false` on failure.

If `ϵ = 0` it tries to determine an optimal choice of `ϵ` by starting
with `ϵ = 1` and then iteratively decreasing it until the bound
`T0(u0, Asymptotic())(ϵ) < C` holds.

It then uses the asymptotic version on the interval `[0, ϵ]` and the
non-asymptotic version on `[ϵ, π]`.
"""
function D0_bounded_by(
    u0::FractionalKdVAnsatz{Arb},
    C::Arf;
    M::Integer = 3,
    ϵ::Arf = zero(Arf),
    maxevals = 1000,
    threaded = false,
    verbose = false,
)
    # Determine ϵ such that the bound holds when evaluated at x = ϵ
    # with the asymptotic version.
    if iszero(ϵ)
        ϵ = one(ϵ)
        while !(T0(u0, Asymptotic(), ϵ = Arb(1.1ϵ); M)(Arb(ϵ)) < C)
            # We use the factor 0.8 instead of, say, 0.5 to get
            # slightly better (higher) values for ϵ.
            ϵ *= 0.8
        end

        verbose && @info "ϵ for asymptotic evaluation to satisfy bound" ϵ
    end

    f = T0(u0, Asymptotic(), ϵ = Arb(1.1ϵ), return_enclosure = true; M)

    # Check that the bound holds on [0, ϵ]
    asymptotic_bound = ArbExtras.bounded_by(f, Arf(0), ϵ, C, degree = -1; threaded, verbose)

    if !asymptotic_bound
        verbose && @info "Bound doesn't hold on [0, ϵ]"
        return false
    end

    verbose && @info "Bound holds on [0, ϵ]"

    g = T0(u0, Ball(), skip_div_u0 = true)

    h(x) = begin
        res = g(x)

        isfinite(res) || return res

        u0x = ArbExtras.enclosure_series(u0, x)

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
        depth_start = ifelse(isone(u0.p), 4, 8);
        maxevals,
        threaded,
        verbose,
    )

    if !non_asymptotic_bound
        verbose && @info "Bound doesn't hold on [ϵ, π]"
        return false
    end

    verbose && @info "Bound holds on [ϵ, π]"

    return true
end
