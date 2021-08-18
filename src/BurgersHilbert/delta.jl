function delta0(
    u0::BHAnsatz{Arb};
    M::Integer = 3,
    rtol = 1e-3,
    degree = 4,
    maxevals = 5000,
    verbose = false,
)
    ϵ = Arb(1e-2)

    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    # TODO: Bound the value on [0, ϵ] using an asymptotic expansion
    @warn "asymptotic part not implemented"
    m1 = zero(ϵ)

    # Bound the value on [ϵ, π]
    m2 = ArbExtras.maximum_enclosure(
        F0(u0),
        Arblib.lbound(ϵ),
        Arblib.ubound(Arb(π)),
        abs_value = true,
        threaded = true;
        degree,
        maxevals,
        rtol,
        verbose,
    )

    return max(m1, m2)
end
