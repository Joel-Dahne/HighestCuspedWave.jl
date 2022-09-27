"""
    findp0(α)

Compute `p0` such that
```
cospi((2α - p) / 2) * gamma(2α - p) / (cospi((α - p) / 2) * gamma(α - p)) =
    2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```

The right hand side has a removable singularity at `α = 1 / 2`. To
avoid this we use that
```
gamma(2α) = gamma(2α + 2) / (2α * (2α + 1))
```
and
```
cospi(α) / (2α + 1) = sinpi(α + 1 / 2) / (2α + 1) = π / 2 * sinc(α + 1 / 2)
```
where `sinc(x) = sinpi(x) / (π * x)`, following the Julia notation, we have
```
gamma(2α) * cospi(α)
= gamma(2α + 2) * sinpi(α + 1 / 2) / (2α * (2α + 1))
= π * sinc(α + 1 / 2) * gamma(2α + 2) / 4α
```
This gives us
```
2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
= π * sinc(α + 1 / 2) * gamma(2α + 2) / (2α * gamma(α) * cospi(α / 2))
```

In practice `p0` is monotone in `α` but we have not proved it. However
we never actually compute `p0` for wide values of `α` since we usually
take the midpoint. If this change we might return to using the
monotonicity and then we have to prove it.

**TODO:** Do we need to prove that `p0` is the smallest root in the
interval? We don't strictly use it.
"""
function findp0(α)
    f(p) = begin
        cospi((2α - p) / 2) * gamma(2α - p) / (cospi((α - p) / 2) * gamma(α - p)) -
        π * sinc(α + 0.5) * gamma(2α + 2) / (2α * gamma(α) * cospi(α / 2))
    end

    # We use, but don't have to prove, that p0 < 1.5(α + 1)
    n = 1000
    ps = range(0, stop = 1.51(α + 1), length = n)[2:end]
    res = f.(ps)

    # Find first sign change
    i = findfirst(i -> res[i] * res[i+1] <= 0, 1:n-1)

    p0 = first(nlsolve(p -> [f(p[1])], [ps[i]], autodiff = :forward, ftol = 1e-15).zero)

    return p0
end

function findp0(α::Arb)
    if iswide(α)
        @warn "findp0 doesn't handle wide α values well" α
        #α_low, α_upp = getinterval(Arb, α)
        #return Arb((findp0(α_low), findp0(α_upp)))
    end

    # We do the computations at a higher precision
    α = setprecision(α, 2precision(α))

    C = π * sinc(α + 1 // 2) * gamma(2α + 2) / (2α * gamma(α) * cospi(α / 2))
    f(p) = begin
        cospi((2α - p) / 2) * gamma(2α - p) / (cospi((α - p) / 2) * gamma(α - p)) - C
    end

    # Find an approximation in Float64
    p0_approx = Arb(findp0(Float64(α)), prec = precision(α))

    # Perform some Newton iterations at higher precision
    p0_approx = let n = 3, p0 = Arb(p0_approx, prec = precision(α))
        for i = 1:n
            y = f(ArbSeries((p0, 1)))
            p0 = Arblib.midpoint(Arb, p0 - y[0] / y[1])
        end
        p0
    end

    # We want to widen the approximation p0_approx slightly so that it
    # definitely contains a root. We do this in a very simple way by
    # widening it step by step until refine_root succeeds. This is
    # definitely not the most efficient way, but good enough.
    fp0 = f(ArbSeries((p0_approx, 1)))
    step = max(eps(p0_approx), abs(fp0[0] / fp0[1]))

    p0 = let p0_approx = copy(p0_approx)
        i = 0
        while true
            Arblib.add_error!(p0_approx, 2^i * step)
            p0 = ArbExtras.refine_root(f, p0_approx, strict = true)
            isnan(p0) || break
            i += 1
            i > 100 && throw(ErrorException("could not isolate p0"))
        end
        p0
    end

    return setprecision(p0, precision(α) ÷ 2)
end

"""
    finda0(α)

Compute `a0` such that `a0(u0, 0)^2/2 - A0(u0, 0)` is zero. That is,
compute
```
a[0] = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```

We can handle the removable singularity at `α = -1 / 2` using the same
approach as in [`findp0`](@ref) to write
```
gamma(2α) * cospi(α)
= gamma(2α + 2) * sinpi(α + 1 / 2) / (2α * (2α + 1))
= π * sinc(α + 1 / 2) * gamma(2α + 2) / 4α
```
giving us
```
a[0] = π * sinc(α + 1 / 2) * gamma(2α + 2) / (2α * gamma(α)^2 * cospi(α / 2)^2)
```

For wide values of `α` it is important to compute tight enclosures. In
practice `a[0]` is increasing in `α` and to get good enclosures we
therefore first try to prove that `a[0]` indeed is increasing. This is
done using [`ArbExtras.bounded_by`](@ref) to prove that minus the
derivative is bounded by `0` from above, in which case the derivative
is bounded by `0` from below. If it is proved to be increasing we
evaluate on the endpoints. If the proof fails we evaluate it using
[`ArbExtras.enclosure_series`](@ref) and print a warning.
"""
function finda0(α)
    f(α) = π * _sinc(α + 1 // 2) * gamma(2α + 2) / (2α * gamma(α)^2 * cospi(α / 2)^2)

    if iswide(α)
        # In this case α must be an Arb

        # Try to prove that a0 is increasing on α

        # Function for evaluating derivative of -f(α)
        mdf(α::Arb) = -f(ArbSeries((α, 1)))[1]
        mdf(α::ArbSeries) =
            if iszero(Arblib.degree(α))
                ArbSeries(df(α[0]))
            else
                Arblib.derivative(-f(ArbSeries(α, degree = Arblib.degree(α) + 1)))
            end

        # Try to prove that mdf is bounded by 0
        α_low, α_upp = getinterval(α)
        is_increasing = ArbExtras.bounded_by(mdf, α_low, α_upp, Arf(0), maxevals = 10)

        if is_increasing
            return Arb((finda0(Arb(α_low)), finda0(Arb(α_upp))))
        else
            @warn "a0 not proved to be increasing on α" α
            return ArbExtras.enclosure_series(f, α, degree = 4)
        end
    end

    return f(α)
end

function _findas(u0::FractionalKdVAnsatz; verbose = true)
    f = D(u0, Symbolic())
    g(a) = f(OffsetVector([u0.a[0]; a], 0:u0.N0))

    initial = u0.a[1:end]

    # IMPROVE: This zero finding problem is not always that stable. It
    # might be helpful to look at stabilising it more for larger
    # systems. Possibly by tuning the tolerance or other parameters.
    sol = nlsolve(g, initial, autodiff = :forward, iterations = 50)

    if verbose && !sol.f_converged
        @warn "Solution for u0.a did not converge" u0.α u0.N0
    end

    return sol.zero, sol.f_converged
end

"""
    findas(u0)

Find values for `u0.a[j]` for `j > 0` that makes the coefficients of
the leading terms in the asymptotic expansion of the defect zero.

This is done by solving the corresponding non-linear system.

It uses [`nlsolve`](@ref) to find the zero, however `nlsolve` doesn't
support `Arb` so this is always done in `Float64`. To speed it up we
start by computing the first `n` coefficients, then `2n` and so on
until we reach `u0.N0`.
"""
function findas(u0::FractionalKdVAnsatz{T}; minstart = 16) where {T}
    if iszero(u0.N0)
        return T[]
    end
    if u0.N0 <= minstart
        return _findas(u0)[1]
    end

    u0 = deepcopy(u0)

    start = min(max(findlast(!iszero, u0.a), 16), minstart)
    stop = u0.N0
    N0s = [start; [start * 2^i for i = 1:floor(Int, log2(stop / start) - 1)]; stop]

    for i in eachindex(N0s)[2:end]
        resize!(u0.a, N0s[i] + 1)
        u0.a[N0s[i-1]:end] .= zero(T)

        u0.a[1:end] .= _findas(u0)[1]
    end

    return u0.a[1:end]
end

function findas(u0::FractionalKdVAnsatz{Arb})
    u0_float = convert(FractionalKdVAnsatz{Float64}, u0)

    return convert(Vector{Arb}, findas(u0_float))
end

"""
    _find_good_as!(u0::FractionalKdVAnsatz; N0_bound, return_defects, threaded, verbose)

Find `N0` and `u0.a` such that the defect is minimized.

It checks `N0` starting from `0` to `N0_bound` and computes `u0.a`
with [`_findas`](@ref). It then takes the `N0` which gave the smallest
defect.
"""
function _find_good_as!(
    u0::FractionalKdVAnsatz;
    N0_bound = 30,
    return_defects = false,
    threaded = false,
    verbose = false,
)
    @assert iszero(u0.N1)
    resize!(u0.a, 1) # Keep only a[0]

    # Vectors for storing defects and coefficients
    defects = [delta0_estimate(u0; threaded)]
    ass = [copy(u0.a[1:end])]

    N0 = 1
    while N0 <= N0_bound
        # Compute ansatz
        resize!(u0.a, N0 + 1)
        u0.a[N0] = zero(eltype(u0.a))

        as, converged = _findas(u0, verbose = false)
        u0.a[1:end] .= as
        push!(ass, as)

        # If the solution didn't converge we stop
        if !converged
            #verbose && @info "No convergence" N0
            break
        end

        # Compute defect
        push!(defects, delta0_estimate(u0; threaded))

        # If the defect is worse than with N0 = 0 we stop
        if defects[end] > defects[1]
            #verbose && @info "Defect worse than for N0 = 0" N0
            break
        end

        N0 += 1
    end

    best_N0 = findmin(defects)[2] - 1
    verbose && @info "Determined best N0" best_N0

    return ass[best_N0+1]
end

function find_good_as(
    u0::FractionalKdVAnsatz{T};
    N0_bound = 30,
    threaded = false,
    verbose = false,
) where {T}
    if T == Float64
        u0_float = deepcopy(u0)
    else
        u0_float = convert(FractionalKdVAnsatz{Float64}, u0)
    end
    empty!(u0_float.b)

    return convert(Vector{T}, _find_good_as!(u0_float; N0_bound, threaded, verbose))
end


"""
    _finda1a2(α)

Compute `a1` and `a2` so that the leading terms in the asymptotics
are exactly zero.

This only works with an expansion of exactly three terms, i.e. it
should not include `a3` or higher. It is mainly intended for testing
when `α -> 0`.

For explanation of the procedure see [`expansion_as`](@ref).
"""
function _finda1a2(α)
    a0 = finda0(α)
    p0 = findp0(α)

    z(i, j) = zeta(-1 - i * α + j * p0)

    z1 = z(2, 1) * z(1, 2) - z(2, 2) * z(1, 1)
    z2 = z(2, 2) * z(1, 0) - z(1, 2) * z(2, 0)
    z3 = z(2, 1) * z(1, 0) - z(1, 1) * z(2, 0)

    d = z1
    v1 = a0 * z2
    v2 = -a0 * z3

    a1 = v1 / d
    a2 = v2 / d

    return a1, a2
end

"""
    findbs(u0, initial)

Find values of `u0.b[n]` to minimize the defect `D(u0)`.

This is done by solving the non-linear system given by requiring that
`D(u0)` evaluates to zero on `u0.N1` collocation points.

It uses [`nlsolve`](@ref) to find the zero, however `nlsolve` doesn't
support `Arb` so this is always done in `Float64`.
"""
function findbs(u0::FractionalKdVAnsatz{T}; verbose = true) where {T}
    if u0.N1 == 0
        return T[], true
    end

    n = u0.N1
    xs = π * (1:2:2n-1) / 2n

    f = D(u0, xs)
    g(b) = f(u0.a.parent, b)

    initial = u0.b
    sol = nlsolve(g, initial, autodiff = :forward)

    if verbose && !sol.f_converged
        @warn "Solution for u0.b did not converge" u0.α u0.N1
    end

    return sol.zero, sol.f_converged
end

function findbs(u0::FractionalKdVAnsatz{Arb}; verbose = true)
    u0_float = convert(FractionalKdVAnsatz{Float64}, u0)

    return convert(Vector{Arb}, findbs(u0_float; verbose)[1])
end
