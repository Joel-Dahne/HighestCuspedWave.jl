"""
    D0_bounded_by(u0::AbstractAnsatz{Arb}, C::Arf; ϵ = 0, ...)

Attempt to prove that `D₀` is bounded by `C`. Returns `true` on
success and `false` on failure.

If `ϵ = 0` it tries to determine an optimal choice of `ϵ` by starting
with `ϵ = 1` and then iteratively decreasing it until the bound
`T0(u0, Asymptotic())(ϵ) < C` holds.

It then uses the asymptotic version on the interval `[0, ϵ]` and the
non-asymptotic version on `[ϵ, π]`.

It tries to find `ϵ2` such that the bound can be proved with one
evaluation on ``[0, ϵ2]``. It starts with `ϵ2 = ϵ` and then squares it
until it succeeds. On the interval ``[ϵ2, ϵ]`` it then uses
logarithmic bisection.
"""
function D0_bounded_by(
    u0::FractionalKdVAnsatz{Arb},
    C::Arf;
    M::Integer = 3,
    ϵ::Arf = zero(Arf),
    maxevals = 4000,
    threaded = true,
    verbose = false,
)
    # Determine ϵ such that the bound holds when evaluated at x = ϵ
    # with the asymptotic version.
    if iszero(ϵ)
        ϵ = one(ϵ)
        while !(T0(u0, Asymptotic(), ϵ = Arb(1.01ϵ); M)(Arb(ϵ)) < C)
            # We use the factor 0.8 instead of, say, 0.5 to get
            # slightly better (higher) values for ϵ.
            ϵ *= 0.8
        end

        verbose && @info "ϵ for asymptotic evaluation to satisfy bound" ϵ
    end

    f = T0(u0, Asymptotic(), ϵ = Arb(1.01ϵ), return_enclosure = true; M)

    # Find ϵ2 such that the bound holds
    ϵ2 = ϵ
    while !(f(Arb((0, ϵ2))) < C)
        ϵ2 = ϵ2 < 0.5 ? ϵ2^2 : ϵ2 / 2
    end
    verbose && @info "Bound holds on [0, ϵ2]" ϵ2

    # Check that the bound holds on [0, ϵ]
    asymptotic_bound = ArbExtras.bounded_by(
        f,
        ϵ2,
        ϵ,
        C,
        degree = -1,
        depth = 100,
        log_bisection = true;
        depth_start = ifelse(weightisx(u0), 5, 0),
        maxevals,
        threaded,
        verbose,
    )

    if !asymptotic_bound
        verbose && @info "Bound doesn't hold on [ϵ2, ϵ]"
        return false
    end

    verbose && @info "Bound holds on [ϵ2, ϵ]"

    g = T0(u0, Ball(), skip_div_u0 = true)

    h(x) = begin
        res = g(x)

        isfinite(res) || return res

        u0x = ArbExtras.enclosure_series(u0, x)

        res /= u0x

        return res
    end

    # For α close to -1 the value of ϵ can get very small and in some
    # cases the non-asymptotic evaluation doesn't even satisfy the
    # bound when evaluated exactly at ϵ. In this case the bound can
    # never be proved. Check for this case and return quickly in case
    # it fails.
    if !(h(Arb(ϵ)) < C)
        verbose && @info "Bound doesn't hold at x = ϵ"
        return false
    end

    # Check that the bound holds on [ϵ, π]
    non_asymptotic_bound = ArbExtras.bounded_by(
        h,
        ϵ,
        ubound(Arb(π)),
        C,
        degree = -1,
        depth_start = ifelse(iswide(u0.α), 8, 7);
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
