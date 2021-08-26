"""
    eval_expansion(u0::BHAnsatz, expansion, x)

Evaluate the given expansion. The term `((i, m, k, l), y)` is
evaluated to `y * log(abs(x))^i * abs(x)^(-k*u0.v0.α + l*u0.v0.p0 +
m)` and then they are all summed.

In general `x` needs to be given both when computing the expansion and
when evaluating it.

The terms `log(abs(x))^i * abs(x)^exponent` requires some extra work
to bound when `x` contains zero (and `exponent > 0`). Due to the
absolute value it's enough to bound it on the interval `[0,
abs_ubound(x)]`. For `x = 0` it's zero and for the upper bound it's
easily computed. The critical points are given by `x = 1` (if `i > 1`)
and `x = exp(-i / exponent)`. For `x = 1` the value is zero and hence
not important, the other critical point we have to take into account.

"""
function eval_expansion(
    u0::BHAnsatz{T},
    expansion::AbstractDict{NTuple{4,S},T},
    x,
) where {T,S}
    res = zero(x)

    logabsx = log(abs(x))
    if x isa Arb
        x_upper = Arblib.abs_ubound(Arb, x)
        logx_upper = log(x_upper)
    end

    for ((i, m, k, l), y) in expansion
        if iszero(k) && iszero(l)
            exponent = m
        else
            exponent = -k * u0.v0.α + l * u0.v0.p0 + m
        end

        if iszero(i)
            factor = abspow(x, exponent)
        else
            if x isa Arb && Arblib.contains_zero(x)
                if exponent > 0
                    if !iszero(x)
                        factor = union(zero(res), logx_upper^i * x_upper^exponent)

                        critical_point = exp(oftype(res, -i) / exponent)
                        if Arblib.overlaps(x, critical_point)
                            factor = union(
                                factor,
                                log(critical_point)^i * critical_point^exponent,
                            )
                        end
                    end
                elseif iszero(exponent)
                    i < 0 || throw(
                        ErrorException("zero x-exponent with non-negative log-exponent"),
                    )
                    if !iszero(x)
                        #x_upper = Arblib.abs_ubound(Arb, x)
                        x_upper < 1 || throw(ErrorException("division by log(1)"))
                        factor = union(zero(res), logx_upper^i * x_upper^exponent)
                    end
                else
                    throw(ErrorException("negative x-exponent"))
                end
            else
                factor = logabsx^i * abspow(x, exponent)
            end
        end
        res += y * factor
    end

    return res
end

"""
    _eval_expansion!(u0::BHAnsatz, expansion, x)

Fast, inplace version of [`eval_expansion`](@ref).

It only works for `Arb` and `ArbSeries`. It takes the expansion as a
vector instead of a dictionary and assumes that the `exponent`
occurring in [`eval_expansion`](@ref) is precomputed and added to the
vector instead of `(m, k, l)`.

It assumes that `x >= 0` an any negative parts of the ball `x` are
ignored. The case when `x` contains zero is only supported for the
type `Arb` and not for `ArbSeries`, for `ArbSeries` is just returns
`NaN` in that case.
"""
function _eval_expansion!(
    res::T,
    u0::BHAnsatz{Arb},
    expansion::Vector{Tuple{Int,Arb,Arb}},
    x::T,
) where {T<:Union{Arb,ArbSeries}}
    if T == Arb && Arblib.contains_nonpositive(x)
        iszero(x) && return Arblib.zero!(res)
        # Use the upper bound of x
        x = Arblib.ubound(Arb, x)
        contained_zero = true
    elseif T == ArbSeries && Arblib.contains_nonpositive(x[0])
        for i = 0:Arblib.degree(res)
            res[i] = NaN
        end
        return res
    else
        contained_zero = false
    end

    Arblib.zero!(res)

    len = length(x)
    term = zero(x)

    logx = log(x)
    invlogx = inv(logx)
    z = zero(x)

    for (i, exponent, y) in expansion
        # term = x^exponent
        if T == Arb
            Arblib.pow!(term, x, exponent)
        elseif T == ArbSeries
            Arblib.pow_arb_series!(term, x, exponent, len)
        end

        if i == -1
            # term *= inv(log(x))
            if T == Arb
                Arblib.mul!(term, term, invlogx)
            elseif T == ArbSeries
                Arblib.mullow!(term, term, invlogx, len)
            end
        elseif i != 0
            # term *= log(x)^i
            throw(ArgumentError("i = $i: this should not occur in the cases we care about"))
        end

        Arblib.mul!(term, term, y)

        # Use that each term is monotonically increasing and zero at x
        # = 0, except the constant term
        if contained_zero && !(iszero(i) && iszero(exponent))
            Arblib.union!(term, term, z)
        end

        # res += term
        if T == Arb
            Arblib.add!(res, res, term)
        elseif T == ArbSeries
            Arblib.add_series!(res, res, term, len)
        end
    end

    return res
end

function (u0::BHAnsatz{T})(x, ::Ball) where {T}
    conv = ifelse(
        T == arb,
        parent(u0.a0),
        ifelse(T == ArbSeries, a -> convert(Arb, a), a -> convert(T, a)),
    )

    res = u0.a0 * (Ci(x, 2, 1) - zeta(conv(2), d = 1))
    res += u0.a1 * (Ci(x, 2) - zeta(conv(2)))

    for n = 1:u0.N
        res += u0.b[n] * (cos(n * x) - 1)
    end

    if !isnothing(u0.v0)
        res += u0.v0(x)
    end

    return res
end

function (u0::BHAnsatz)(x, ::Asymptotic; M::Integer = 3)
    return eval_expansion(u0, u0(x, AsymptoticExpansion(); M), x)
end

"""
    (u0::AbstractAnsatz)(x, ::AsymptoticExpansion; M = 3)

Return a dictionary containing the terms in the asymptotic expansion
of `u0` which can then be evaluated with [`eval_expansion`](@ref).

The highest term, `x^2M`, an error term is which makes sure that
evaluation of the expansion gives an enclosure of the result.
"""
function (u0::BHAnsatz{Arb})(x, ::AsymptoticExpansion; M::Integer = 3)
    # TODO: Check this
    @assert M >= 3

    res = OrderedDict{NTuple{4,Int},Arb}()

    # Error term
    res[(0, 2M, 0, 0)] = 0

    # First Clausian
    γ = Arb(Irrational{:γ}())
    res[(1, 1, 0, 0)] = -Arb(π) / 2 * u0.a0
    res[(0, 1, 0, 0)] = -(γ - 1) * Arb(π) / 2 * u0.a0
    for m = 1:M-1
        res[(0, 2m, 0, 0)] = (-1)^m * zeta(Arb(2 - 2m), d = 1) / factorial(Arb(2m)) * u0.a0
    end
    # PROVE: Give an expression for the error term. The one we use
    # here is twice the coefficient for x^2M which should be
    # larger than required but this is yet to be proved.
    Arblib.add_error!(
        res[(0, 2M, 0, 0)],
        2abs(zeta(Arb(2 - 2M), d = 1) / factorial(Arb(2M))) * u0.a0,
    )

    # Second Clausian
    if !iszero(u0.a1)
        C, _, p, E = Ci_expansion(x, Arb(2), M)
        res[(0, 1, 0, 0)] += C * u0.a1
        for m = 1:M-1
            res[(0, 2m, 0, 0)] += p[2m] * u0.a1
        end
        Arblib.add_error!(res[(0, 2M, 0, 0)], E * u0.a1)
    end

    # Fourier terms
    if !iszero(u0.N)
        for m = 1:M-1
            res[(0, 2m, 0, 0)] +=
                (-1)^m * sum(Arb(n)^(2m) * u0.b[n] for n = 1:u0.N) / factorial(Arb(2m))
        end
        Arblib.add_error!(
            res[(0, 2M, 0, 0)],
            sum(Arb(n)^(2M) * abs(u0.b[n]) for n = 1:u0.N) / factorial(Arb(2M)),
        )
    end

    # Clausians coming from u0.v0
    if !isnothing(u0.v0)
        let α = u0.v0.α, p0 = u0.v0.p0
            for j = 1:u0.v0.N0
                C, _, p, E = Ci_expansion(x, 1 - α + j * p0, M)
                res[(0, 0, 1, j)] = C * u0.v0.a[j]
                for m = 1:M-1
                    res[(0, 2m, 0, 0)] += p[2m] * u0.v0.a[j]
                end
                Arblib.add_error!(res[(0, 2M, 0, 0)], E)
            end
        end
    end

    return res
end

function H(u0::BHAnsatz{T}, ::Ball) where {T}
    conv = ifelse(
        T == arb,
        parent(u0.a0),
        ifelse(T == ArbSeries, a -> convert(Arb, a), a -> convert(T, a)),
    )

    return x -> begin
        res = -u0.a0 * (Ci(x, 3, 1) - zeta(conv(3), d = 1))
        res -= u0.a1 * (Ci(x, 3) - zeta(conv(3)))

        for n = 1:u0.N
            res -= u0.b[n] / n * (cos(n * x) - 1)
        end

        # Add Clausians coming from u0.v0
        if !isnothing(u0.v0)
            let α = u0.v0.α, p0 = u0.v0.p0
                for j = 1:u0.v0.N0
                    res -= u0.v0.a[j] * (Ci(x, 2 - α + j * p0) - zeta(2 - α + j * p0))
                end
            end
        end

        return res
    end
end

function H(u0::BHAnsatz, ::Asymptotic; M::Integer = 3)
    f = H(u0, AsymptoticExpansion(); M)

    return x -> eval_expansion(u0, f(x), x)
end

"""
    H(u0::BHAnsatz, ::AsymptoticExpansion; M = 3)

Return a dictionary containing the terms in the asymptotic expansion
of `H(u0)` which can then be evaluated with [`eval_expansion`](@ref).

The highest term, `x^2M`, an error term is which makes sure that
evaluation of the expansion gives an enclosure of the result.
"""
function H(u0::BHAnsatz{Arb}, ::AsymptoticExpansion; M::Integer = 3)
    # TODO: Check this
    @assert M >= 3

    π = Arb(Irrational{:π}())
    γ = Arb(Irrational{:γ}())
    γ₁ = stieltjes(Arb, 1)

    return x -> begin
        res = OrderedDict{NTuple{4,Int},Arb}()

        # Error term
        res[(0, 2M, 0, 0)] = 0

        # First Clausian
        res[(2, 2, 0, 0)] = -1 // 4 * u0.a0
        res[(1, 2, 0, 0)] = (3 // 4 - γ / 2) * u0.a0
        for m = 1:M-1
            if m == 1
                term = (3γ - γ^2 - 2γ₁ - 7 // 2 + (π^2) / 12) / 2
            else
                term = zeta(Arb(3 - 2m), d = 1)
            end
            res[(0, 2m, 0, 0)] = -(-1)^m * term * u0.a0 / factorial(Arb(2m))
        end
        # PROVE: Give an expression for the error term. The one we use
        # here is twice the coefficient for x^2M which should be
        # larger than required but this is yet to be proved.
        Arblib.add_error!(
            res[(0, 2M, 0, 0)],
            2abs(zeta(Arb(3 - 2M), d = 1) / factorial(Arb(2M))) * u0.a0,
        )

        # Second Clausian
        if !iszero(u0.a1)
            res[(1, 2, 0, 0)] += -u0.a1 / 2
            for m = 1:M-1
                if m == 1
                    term = -Arb(3 // 2)
                else
                    term = -zeta(Arb(3 - 2m))
                end
                res[(0, 2m, 0, 0)] += (-1)^m * term * u0.a1 / factorial(Arb(2m))
            end
            Arblib.add_error!(
                res[(0, 2M, 0, 0)],
                2(2π)^(4 - 2M) * zeta(2M - 2) / (4π^2 - x^2),
            )
        end

        # Fourier terms
        if !iszero(u0.N)
            for m = 1:M-1
                res[(0, 2m, 0, 0)] -=
                    (-1)^m * sum(Arb(n)^(2m - 1) * u0.b[n] for n = 1:u0.N) /
                    factorial(Arb(2m))
            end
            Arblib.add_error!(
                res[(0, 2M, 0, 0)],
                sum(Arb(n)^(2M - 1) * abs(u0.b[n]) for n = 1:u0.N) / factorial(Arb(2M)),
            )
        end

        # Clausians coming from u0.v0
        if !isnothing(u0.v0)
            let α = u0.v0.α, p0 = u0.v0.p0
                for j = 1:u0.v0.N0
                    C, _, p, E = Ci_expansion(x, 2 - α + j * p0, M)
                    res[(0, 1, 1, j)] = -C * u0.v0.a[j]
                    for m = 1:M-1
                        res[(0, 2m, 0, 0)] -= p[2m] * u0.v0.a[j]
                    end
                    Arblib.add_error!(res[(0, 2M, 0, 0)], E)
                end
            end
        end

        return res
    end
end

function D(u0::BHAnsatz, ::Asymptotic; M::Integer = 3)
    f = D(u0, AsymptoticExpansion(); M)

    return x -> eval_expansion(u0, f(x), x)

    if isnothing(u0.v0) && false
        5
    else
        # This is more or less a temporary solution
        f = x -> u0(x, Asymptotic(), M = M)
        g = H(u0, Asymptotic(), M = M)
        return x -> begin
            f(x)^2 / 2 + g(x)
        end
    end
end

function D(u0::BHAnsatz, evaltype::AsymptoticExpansion; M::Integer = 3)
    f = x -> u0(x, evaltype; M)
    g = H(u0, evaltype; M)

    return x -> begin
        expansion1 = f(x)
        expansion2 = g(x)
        expansion = empty(expansion1)

        # u0^2/2 term
        let expansion1 = collect(expansion1)
            z = zero(x) # Avoid allocating zero multiple times
            for (i, (key1, a1)) in enumerate(expansion1)
                expansion[2 .* key1] = get(expansion, 2 .* key1, z) + a1^2 / 2
                for j = i+1:length(expansion1)
                    (key2, a2) = expansion1[j]
                    key = key1 .+ key2
                    expansion[key] = get(expansion, key, z) + a1 * a2
                end
            end
        end

        # H term
        merge!(+, expansion, expansion2)

        # The term (2, 2, 0, 0) is identically equal to zero due to
        # the choice of a0
        # TODO: Check this
        @assert Arblib.contains_zero(expansion[(2, 2, 0, 0)])
        delete!(expansion, (2, 2, 0, 0))

        return expansion
    end
end

function F0(u0::BHAnsatz, evaltype::Ball)
    f = H(u0, evaltype)
    return x -> begin
        y = u0(x, evaltype)

        return (y^2 / 2 + f(x)) / (u0.w(x) * y)
    end
end

"""
    F0(u0::BHAnsatz{Arb}, ::Asymptotic)

Returns a function such that `F0(u0, Asymptotic())(x)` is computed
accurately for small values of `x`.

This method is one of the bottlenecks of [`delta0`](@ref) and for this
reasons it is heavily optimized. The method assumes that `x >= 0` and
any negative parts of the ball `x` are ignored.

It precomputes the expansions of `u0` and `D(u0)` and for that reason
a number `ϵ` has to be given, the resulting expansion will be valid
for all `x < ϵ`. The value of `ϵ` has to be less than `1`.

If `exponent_limit` is set to some number then all terms in the
expansions with an exponent equal to or greater than this limit will
be collapsed into one term.
TODO: Be more detail about how this works and make sure that it works
with the log-factors as well.
"""
function F0(
    u0::BHAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 3,
    ϵ::Arb = Arb(1e-2),
    exponent_limit::Union{Arb,Nothing} = nothing,
)
    @assert ϵ < 1

    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), abs(x) * sqrt(log((abs(x) + 1) / abs(x))))
    end

    u0_expansion = u0(ϵ, AsymptoticExpansion(); M)
    Du0_expansion = D(u0, AsymptoticExpansion(); M)(ϵ)

    # Divide the expansion by x * log(x) and x^2 * log(x)
    # respectively, also precompute the exponents.
    u0_expansion_div_xlogx = Vector{Tuple{Int,Arb,Arb}}(undef, length(u0_expansion))
    Du0_expansion_div_x2logx = Vector{Tuple{Int,Arb,Arb}}(undef, length(Du0_expansion))
    for (index, ((i, m, k, l), value)) in enumerate(u0_expansion)
        u0_expansion_div_xlogx[index] = (i - 1, -k * u0.v0.α + l * u0.v0.p0 + m - 1, value)
    end
    for (index, ((i, m, k, l), value)) in enumerate(Du0_expansion)
        Du0_expansion_div_x2logx[index] =
            (i - 1, -k * u0.v0.α + l * u0.v0.p0 + m - 2, value)
    end

    # If an exponent limit is given, collapse all terms with an
    # exponent greater than or equal to that.
    if !isnothing(exponent_limit)
        let coeff = sum(
                getindex.(filter(x -> x[2] >= exponent_limit, u0_expansion_div_xlogx), 3),
            )
            filter!(x -> !(x[2] >= exponent_limit), u0_expansion_div_xlogx)
            push!(u0_expansion_div_xlogx, (0, exponent_limit, coeff))
        end
        let coeff = sum(
                getindex.(filter(x -> x[2] >= exponent_limit, Du0_expansion_div_x2logx), 3),
            )
            filter!(x -> !(x[2] >= exponent_limit), Du0_expansion_div_x2logx)
            push!(Du0_expansion_div_x2logx, (0, exponent_limit, coeff))
        end
    end

    # The function inv(sqrt(log((abs(x) + 1) / abs(x))))
    # PROVE: That this is monotonically increasing on [0, 1]
    w!(res, x) = Arblib.set!(res, inv(sqrt(log((x + 1) / x))))

    # Version of w! that also handles x overlapping 0
    weight!(res, x) = begin
        iszero(x) && return Arblib.zero!(res)

        if x isa Arb && Arblib.contains_zero(x)
            w!(res, Arblib.ubound(Arb, x))
            return Arblib.union!(res, res, zero(res))
        end

        return w!(res, x)
    end

    return x -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && x[0] <= ϵ)

        res = weight!(zero(x), x)

        res /= _eval_expansion!(zero(x), u0, u0_expansion_div_xlogx, x)

        res *= _eval_expansion!(zero(x), u0, Du0_expansion_div_x2logx, x)

        return res
    end
end

"""
    D(u0::BHAnsatz, xs::AbstractVector)
Returns a function such that `D(u0, xs)(b)` computes `D(u0)(x)` on the
points `x ∈ xs` with `u0.b` set to the given values. Does this in an
efficient way by precomputing as much as possible.

NOTE: This is **not** rigorous!
"""
function D(u0::BHAnsatz, xs::AbstractVector)
    b = copy(u0.b)
    u0.b .= 0

    u0_xs_a0_precomputed = zeros(length(xs))
    u0_xs_b_precomputed = zeros(length(xs), u0.N)
    Hu0_xs_a0_precomputed = zeros(length(xs))
    Hu0_xs_b_precomputed = zeros(length(xs), u0.N)

    Hu0 = H(u0)
    for (i, x) in enumerate(xs)
        u0_xs_a0_precomputed[i] = u0(x) # u0.a0 * (Ci(x, 2, 1) - zeta(2, d = 1))
        Hu0_xs_a0_precomputed[i] = Hu0(x) # -u0.a0 * (Ci(x, 3, 1) - zeta(3, d = 1))

        #if !iszero(u0.a1)
        #    u0_xs_a0_precomputed[i] += u0.a1 * (Ci(x, 2) - zeta(2))
        #    Hu0_xs_a0_precomputed[i] -= u0.a1 * (Ci(x, 3) - zeta(3))
        #end

        for n = 1:u0.N
            u0_xs_b_precomputed[i, n] = cos(n * x) - 1
            Hu0_xs_b_precomputed[i, n] = -(cos(n * x) - 1) / n
        end
    end

    copy!(u0.b, b)

    return b -> begin
        return (
            (u0_xs_a0_precomputed + u0_xs_b_precomputed * b) .^ 2 / 2 +
            (Hu0_xs_a0_precomputed + Hu0_xs_b_precomputed * b)
        )
    end
end
