"""
    alpha0(u0::FractionalKdVAnsatz; M::Integer, rtol, threaded, verbose)

Enclose the value of `α₀` from the paper.

It uses an asymptotic expansion with `M` terms close to zero and ball
arithmetic on the remaining.

In practice the maximum is attained at `x = π` so the function is
evaluated there first. Then starting at `ϵ = π/2` we halve `ϵ` until
the bound on `[0, ϵ]` is smaller than the value at `π`.

On the interval `[ϵ, π]` we bound it using ball arithmetic with
[`ArbExtras.maximum_enclosure`](@ref). Notice that we do not have to
prove that the maximum is attained at `x = π`, it's only used to make
the procedure more efficient.

If `threaded = true` it enables threading when calling
[`ArbExtras.maximum_enclosure`](@ref). If `verbose = true` print more
information about the process.
"""
function alpha0(
    u0::FractionalKdVAnsatz{Arb};
    M::Integer = 3,
    rtol = 1e-5,
    threaded = true,
    verbose = false,
)
    # This is required for a finite value
    @assert u0.p + u0.α > 0

    # This is the function we want to find the maximum of
    f(x) = u0.w(x) / (2u0(x))

    # The maximum is in practice attained at x = π, so evaluate at this point
    m1 = f(Arb(π))

    ϵ = Arb(π) / 2
    # Make sure ϵ is small enough so that the asymptotic value is well
    # defined
    while isnan(c(u0, ϵ; M))
        ϵ /= 2
    end

    verbose && @info "ϵ value for finite bound" ϵ

    # Continue halving ϵ until the asymptotic value is smaller
    # than m1, but don't go to far. We have
    # abs(x)^p/u0(x) = (1 + hat(u0)(x))/a0(u0, 0)*abs(x)^(u0.p + u0.α)
    # <= (1 + c(u0, ϵ)*abs(x)^u0.p0/a0(u0, 0)*abs(x)^(u0.p + u0.α)
    # Which attains its maximum value at x = ϵ
    m2 = (1 + c(u0, ϵ; M) * abs(ϵ)^u0.p0) / a0(u0, 0) * abs(ϵ)^(u0.p + u0.α)
    while !(m2 < m1) && ϵ > 1e-10
        ϵ /= 2
        m2 = (1 + c(u0, ϵ; M) * abs(ϵ)^u0.p0) / a0(u0, 0) * abs(ϵ)^(u0.p + u0.α)
    end

    verbose && @info "ϵ value for getting smaller value than at π" ϵ

    m2 > m1 && @warn "We have m1 > m2 which should not happen in practice"

    # Bound the value on [ϵ, π] by Ball evaluation
    m3 = ArbExtras.maximum_enclosure(
        f,
        Arblib.lbound(ϵ),
        Arblib.ubound(Arb(π)),
        abs_value = true,
        point_value_max = m1, # m1 is a lower bound of the maximum
        atol = 4Arblib.radius(m1); # We cannot expect to do better than m1
        rtol,
        threaded,
        verbose,
    )

    return max(m1, m2, m3)
end
