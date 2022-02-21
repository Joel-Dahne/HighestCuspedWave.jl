"""
    findp0(α)

Compute `p0` such that
```
cospi((2α - p) / 2) * gamma(2α - p) / (cospi((α - p) / 2) * gamma(α - p)) =
    2gamma(2α)cospi(α) / (gamma(α)cospi(α / 2))
```

We use, but don't have to prove, that `p0 < 1.5(α + 1)`.

When `α::Arb` it uses monotonicity of `α` to compute a tighter
interval.
- **PROVE:** That `p0` is monotone in `α`

**TODO:** Do we need to prove that `p0` is the smallest root in the
interval?
"""
function findp0(α)
    f(p) = begin
        cospi((2α - p) / 2) * gamma(2α - p) / (cospi((α - p) / 2) * gamma(α - p)) -
        2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
    end

    n = 1000
    ps = range(0, stop = 1.51(α + 1), length = n)
    res = f.(ps)

    # Find first sign change
    i = findfirst(i -> res[i] * res[i+1] <= 0, 1:n-1)

    p0 = first(nlsolve(p -> [f(p[1])], [ps[i]], autodiff = :forward, ftol = 1e-15).zero)

    return p0
end

function findp0(α::Arb)
    if iswide(α)
        α_low, α_upp = getinterval(Arb, α)
        return Arb((findp0(α_low), findp0(α_upp)))
    end

    if Float64(α) == -0.5
        # There is a removable singularity at α = -0.5 so evaluation
        # at exactly that point fails. We handle this by creating a
        # small ball around α and use the monotonicity. This is mostly
        # meant for testing, where -0.5 tends to come up frequently...
        α = Arblib.add_error!(copy(α), Arb(1e-15))
        return findp0(α)
    end

    α = setprecision(α, 2precision(α))

    C = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
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
It makes use of the monotinicity to get good enclosures for wide
balls.
- **PROVE:** That `a[0]` is monotone in `α`. Using
  [`ArbExtras.enclosure_series`](@ref) doesn't work very well.
"""
function finda0(α)
    if iswide(α)
        # In this case α must be an Arb
        α_low, α_upp = getinterval(Arb, α)
        return Arb((finda0(α_low), finda0(α_upp)))
    end

    return 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
end

function _findas(u0::FractionalKdVAnsatz)
    f = D(u0, Symbolic())
    g(a) = f(OffsetVector([u0.a[0]; a], 0:u0.N0))

    initial = u0.a[1:end]
    sol = nlsolve(g, initial, autodiff = :forward, iterations = 50)

    if !sol.f_converged
        @warn "Solution did not converge for α = $(u0.α), N0 = $(u0.N0)"
        @warn sol
    end

    return sol.zero
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
        return _findas(u0)
    end

    u0 = deepcopy(u0)

    start = min(max(findlast(!iszero, u0.a), 16), minstart)
    stop = u0.N0
    N0s = [start; [start * 2^i for i = 1:floor(Int, log2(stop / start) - 1)]; stop]

    for i in eachindex(N0s)[2:end]
        resize!(u0.a, N0s[i] + 1)
        u0.a[N0s[i-1]:end] .= zero(T)

        u0.a[1:end] .= _findas(u0)
    end

    return u0.a[1:end]
end

function findas(u0::FractionalKdVAnsatz{Arb})
    u0_float = convert(FractionalKdVAnsatz{Float64}, u0)

    return convert(Vector{Arb}, findas(u0_float))
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
function findbs(u0::FractionalKdVAnsatz{T}) where {T}
    if u0.N1 == 0
        return T[]
    end

    n = u0.N1
    xs = π * (1:2:2n-1) / 2n

    f = D(u0, xs)
    g(b) = f(u0.a.parent, b)

    initial = u0.b
    sol = nlsolve(g, initial, autodiff = :forward)

    if !sol.f_converged
        @warn "Solution did not converge for α = $(u0.α), N1 = $(u0.N1)"
        @warn sol
    end

    return sol.zero
end

function findbs(u0::FractionalKdVAnsatz{Arb})
    u0_float = convert(FractionalKdVAnsatz{Float64}, u0)

    return convert(Vector{Arb}, findbs(u0_float))
end
