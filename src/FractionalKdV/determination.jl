"""
    findp0(α)

Compute `p0` such that
```
cospi((2α - p) / 2) * gamma(2α - p) / (cospi((α - p) / 2) * gamma(α - p)) =
    2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```
See [`equation_p0`](@ref).

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

In practice the root we compute is the smallest positive root.
However, we do not make use of the fact that it indeed is the smallest
positive root and therefore we do not attempt to prove it in any way.

The version with `α::Arb` computes an enclosure of the root, though in
practice we don't use this fact in the rest of the code. The generic
version only computes an approximation of the root.
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
    iswide(α) && @warn "findp0 doesn't handle wide α values well" α

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

Compute `a0` such that the leading term in the asymptotic expansion of
```
H(u0) + u0^2 / 2
```
is zero. That is, compute
```
a0 = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```
See [`equation_a0`](@ref).

We can handle the removable singularity at `α = -1 / 2` using the same
approach as in [`findp0`](@ref) to write
```
gamma(2α) * cospi(α)
= gamma(2α + 2) * sinpi(α + 1 / 2) / (2α * (2α + 1))
= π * sinc(α + 1 / 2) * gamma(2α + 2) / 4α
```
giving us
```
a0 = π * sinc(α + 1 / 2) * gamma(2α + 2) / (2α * gamma(α)^2 * cospi(α / 2)^2)
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
        mdf = ArbExtras.derivative_function(α -> -f(α))

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

function _findas(
    u0::FractionalKdVAnsatz{T};
    use_defect2 = true,
    threaded = true,
    verbose = true,
) where {T}
    iszero(u0.N0) && return T[], true

    if use_defect2
        f = defect2(u0, Symbolic(); threaded)
    else
        f = defect(u0, Symbolic())
    end
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

Numerically find values for `u0.a[j]` for `j = 1:u0.N0` that makes the
coefficients of the leading terms in the asymptotic expansion of the
defect approximately zero.

This is done by solving the corresponding non-linear system.

It uses [`nlsolve`](@ref) to find the zero, however `nlsolve` doesn't
support `Arb` so this is always done in `Float64`. To speed it up we
start by computing the first `n` coefficients, then `2n` and so on
until we reach `u0.N0`.
"""
function findas(
    u0::FractionalKdVAnsatz{T};
    minstart = 16,
    use_defect2 = true,
    threaded = true,
    verbose = true,
) where {T}
    if iszero(u0.N0)
        return T[]
    end
    if u0.N0 <= minstart
        return _findas(u0; use_defect2, threaded, verbose)[1]
    end

    u0 = deepcopy(u0)

    start = min(max(findlast(!iszero, u0.a), 16), minstart)
    stop = u0.N0
    N0s = [start; [start * 2^i for i = 1:floor(Int, log2(stop / start) - 1)]; stop]

    for i in eachindex(N0s)[2:end]
        resize!(u0.a, N0s[i] + 1)
        u0.a[N0s[i-1]:end] .= zero(T)

        u0.a[1:end] .= _findas(u0; use_defect2, threaded, verbose)[1]
    end

    return u0.a[1:end]
end

function findas(u0::FractionalKdVAnsatz{Arb}; use_defect2 = true, verbose = true)
    u0_float = convert(FractionalKdVAnsatz{Float64}, u0)

    # Compute an accurate value of a[0]. The conversion from a
    # wide ball gives large errors
    u0_float.a[0] = Float64(finda0(Arb(u0_float.α)))

    return convert(Vector{Arb}, findas(u0_float; use_defect2, verbose))
end

"""
    _find_good_as!(u0::FractionalKdVAnsatz, N0s::StepRange{Int,Int} = 0:30; iter_use_best_as, use_defect2, return_defects, threaded, verbose)

Find `N0` and `u0.a` such that the defect is minimized. Returns the
vector `a` and the estimated defect, `N0` is implicitly given by the
length.

It checks the values of `N0` in `N0s`. For each `N0` it computes
`u0.a` with [`_findas`](@ref). It then computes an estimate of the
global defect using [`delta0_estimate`](@ref) and takes the `N0` which
gave the smallest defect.

If `iter_use_best_as` is true then use the coefficients from the best
previous iteration as starting point for the zero finding, otherwise
it uses the iteration just before. In general this gives slightly
better results, but not always.

If `use_defect2` is true then use [`defect2`](@ref) instead of
[`defect`](@ref) for computing the coefficients in the asymptotic
expansion.

If `return_defects` is true it returns `N0s, defects`, where `defects`
is a vector of defects for the different values of `N0`. If the
iteration over `N0s` stopped early it returns a truncated `N0s`.
"""
function _find_good_as!(
    u0::FractionalKdVAnsatz{T},
    N0s::StepRange{Int,Int} = 0:1:30;
    iter_use_best_as = true,
    use_defect2 = true,
    return_defects = false,
    threaded = false,
    verbose = false,
) where {T}
    @assert iszero(u0.N1)
    resize!(u0.a, 1) # Keep only a[0]

    resize_with_zero!(a::AbstractVector, n::Int) = begin
        i = lastindex(a)
        resize!(a, n)
        if i < lastindex(a)
            a[i+1:end] .= zero(eltype(a))
        end
        return a
    end

    # If N0s.start is large compute the coefficients for it using
    # findas which is faster than directly using _findas
    if N0s.start >= 256
        resize_with_zero!(u0.a, N0s.start + 1)
        as = findas(u0, verbose = false; use_defect2, threaded)
        u0.a[1:end] .= as
    end

    # Vectors for storing defects and coefficients
    defects = T[]
    ass = typeof(u0.a)[]

    for N0 in N0s
        # Compute approximation
        resize_with_zero!(u0.a, N0 + 1)
        as, converged = _findas(u0, verbose = false; use_defect2, threaded)
        u0.a[1:end] .= as
        push!(ass, as)

        # If the solution didn't converge we stop
        if !converged
            verbose && @info "No convergence" N0
            break
        end

        # Compute defect
        push!(defects, delta0_estimate(u0; threaded))

        # If the defect is worse than the first one we stop
        if defects[end] > defects[1]
            verbose &&
                @info "Defect worse than for at start" N0 N0s.start defects[begin] defects[end]
            break
        end

        if iter_use_best_as
            # Use as from best N0
            best_index = findmin(defects)[2]
            best_N0 = N0s[best_index]
            best_as = ass[best_index]
            u0.a[1:best_N0] .= best_as
            u0.a[best_N0+1:end] .= 0
        end
    end

    return_defects && return N0s[1:length(defects)], defects

    best_defect, best_index = findmin(defects)
    best_N0 = N0s[best_index]
    best_as = ass[best_index]

    verbose && @info "Determined best N0" best_N0 best_defect

    return best_as, best_defect
end

"""
    find_good_as(u0::FractionalKdVAnsatz, N0s = StepRange{Int,Int} = 0:1:30; kwargs)

This function is similar to [`findas`](@ref) in that it tries to find
values for `u0.a[j]` with `j = 1:u0.N0` that makes the coefficients of
the leading terms in the asymptotic expansion of the defect
approximately zero. It is different in that it doesn't take a fixed
value for `N0` but a range `N0s`. It then tries to find the `N0` in
this range that gives the smallest global defect.

# Arguments
See [`_find_good_as`](@ref) for the implementation and the possible
arguments.

In addition to the arguments supported by [`_find_good_as`](@ref) it
accepts the argument `try_all_combinations`. If this argument is set
to true then it runs [`_find_good_as`](@ref) several times with
different combinations of arguments and then take the best result.
More precisely it tries the four combinations of setting
`iter_use_best_as` and `use_defect2` to either true or false. If
`try_all_combinations` is set then the values of the arguments
`iter_use_best_as` and `use_defect2` are ignored and `return_defects`
is assumed to be false. This argument is useful in particular near `α
= -0.9` where the zero finding problem has turned out to be very
unstable and it is hard to find one combination of arguments that work
in all cases. It should be avoided near `α = -1` since it would take a
long time to compute the result.
"""
function find_good_as(
    u0::FractionalKdVAnsatz{T},
    N0s::StepRange{Int,Int} = 0:1:30;
    iter_use_best_as = true,
    use_defect2 = true,
    try_all_combinations = false,
    return_defects = false,
    threaded = false,
    verbose = false,
) where {T}
    if T == Float64
        u0_float = deepcopy(u0)
    else
        u0_float = convert(FractionalKdVAnsatz{Float64}, u0)
        # Compute an accurate value of a[0]. The conversion from a
        # wide ball gives large errors
        u0_float.a[0] = Float64(finda0(Arb(u0_float.α)))
    end
    empty!(u0_float.b)

    if try_all_combinations
        @assert !return_defects

        res_as = Matrix{Vector{T}}(undef, 2, 2)
        res_defect = Matrix{T}(undef, 2, 2)
        for (i, iter_use_best_as) in enumerate((false, true))
            for (j, use_defect2) in enumerate((false, true))
                verbose &&
                    @info "iter_use_best_as = $iter_use_best_as, use_defect2 = $use_defect2"
                res_as[i, j], res_defect[i, j] = _find_good_as!(
                    u0_float,
                    N0s;
                    iter_use_best_as,
                    use_defect2,
                    threaded,
                    verbose,
                )

                empty!(u0_float.b)
            end
        end
        best_defect, best_index = findmin(res_defect)

        verbose && @info "Determined best result" best_index best_defect

        res = (res_as[best_index], best_defect)
    else
        res = _find_good_as!(
            u0_float,
            N0s;
            iter_use_best_as,
            use_defect2,
            return_defects,
            threaded,
            verbose,
        )
    end

    if return_defects
        res[1], convert(Vector{T}, res[2])
    else
        return convert(Vector{T}, res[1])
    end
end


"""
    _finda1a2(α)

Compute `a1` and `a2` so that the leading terms in the asymptotics
are exactly zero.

This only works with an expansion of exactly three terms, i.e. it
should not include `a3` or higher. It is mainly intended for testing
when `α -> 0`.

For explanation of the procedure see [`expansion_as`](@ref).

This function is mostly used for testing.
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
    _findbs(u0::FractionalKdVAnsatz{T}; verbose = true)

Numerically find values of `u0.b[n]` to minimize `defect(u0)`.

This is done by solving the non-linear system given by requiring that
`defect(u0)` evaluates to zero on `u0.N1` collocation points.

It uses [`nlsolve`](@ref) to find the zero. It returns the zero
together with information about the convergence. If `verbose` is true
then it prints a warning if the convergence failed.
"""
function _findbs(u0::FractionalKdVAnsatz{T}; verbose = true) where {T}
    if u0.N1 == 0
        return T[], true
    end

    n = u0.N1
    xs = π * (1:2:2n-1) / 2n

    f = defect(u0, xs)
    g(b) = f(u0.a.parent, b)

    initial = u0.b
    sol = nlsolve(g, initial, autodiff = :forward)

    if verbose && !sol.f_converged
        @warn "Solution for u0.b did not converge" u0.α u0.N1
    end

    return sol.zero, sol.f_converged
end


"""
    findbs(u0::FractionalKdVAnsatz; verbose = false)

Numerically find values of `u0.b[n]` to minimize `defect(u0)`.

See [`_findbs`](@ref) for the implementation. The underscore method
doesn't support `Arb` so it converts `u0` to use `Float64` in that
case.
"""
findbs(u0::FractionalKdVAnsatz; verbose = true) = _findbs(u0; verbose)[1]

function findbs(u0::FractionalKdVAnsatz{Arb}; verbose = true)
    u0_float = convert(FractionalKdVAnsatz{Float64}, u0)
    # Compute an accurate value of a[0]. The conversion from a
    # wide ball gives large errors
    u0_float.a[0] = Float64(finda0(Arb(u0_float.α)))

    return convert(Vector{Arb}, _findbs(u0_float; verbose)[1])
end
