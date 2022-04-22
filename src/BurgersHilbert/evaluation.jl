"""
    eval_expansion(u0::BHAnsatz, expansion, x)

Evaluate the given expansion. The term `((i, m, k, l), y)` is
evaluated to `y * log(abs(x))^i * abs(x)^(-k*u0.α + l*u0.p0 +
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
            exponent = -k * u0.α + l * u0.p0 + m
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
                    else
                        factor = zero(x)
                    end
                elseif iszero(exponent)
                    i < 0 || throw(
                        ErrorException("zero x-exponent with non-negative log-exponent"),
                    )
                    if !iszero(x)
                        x_upper < 1 || throw(ErrorException("division by log(1)"))
                        factor = union(zero(res), logx_upper^i * x_upper^exponent)
                    else
                        factor = zero(x)
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
    _eval_expansion!(u0::BHAnsatz, expansion, x, invlogx = inv(log(x)))

Fast, inplace version of [`eval_expansion`](@ref).

It only works for `Arb` and `ArbSeries`. It takes the expansion as a
vector instead of a dictionary and assumes that the `exponent`
occurring in [`eval_expansion`](@ref) is precomputed and added to the
vector instead of `(m, k, l)`.

It also optionally takes the argument `invlogx` which should be a
precomputed value of `inv(log(x))`. Notice that if `x` overlaps with
zero then care needs to be taken so that this is a finite enclosure,
to not have to implement that logic outside of this method it checks
if `invlogx` is finite in the case when `x` overlaps with zero and
otherwise it computes it.

It assumes that `x >= 0` and any negative parts of the ball `x` are
ignored. The case when `x` contains zero is only supported for the
type `Arb` and not for `ArbSeries`, for `ArbSeries` is just returns
`NaN` in that case.
"""
function _eval_expansion!(
    res::T,
    u0::BHAnsatz{Arb},
    expansion::Vector{Tuple{Int,Arb,Arb}},
    x::T,
    invlogx::T = inv(log(x)),
    term_buffer::T = zero(x),
) where {T<:Union{Arb,ArbSeries}}
    if T == Arb && Arblib.contains_nonpositive(x)
        iszero(x) && return Arblib.zero!(res)
        if !isfinite(invlogx)
            # Compute an enclosure of inv(log(x)) using that it is
            # zero at x = 0 and decreasing. We don't expect this case
            # to be important for the performance.
            invlogx = Arb((inv(log(ubound(Arb, x))), 0))
        end
        # We do the computations with the upper bound of x. For each
        # non-constant term in the expansion we then take the union
        # with zero.
        x = ubound(Arb, x)
        contained_zero = true
        z = zero(x) # Used when taking the union with zero
    elseif T == ArbSeries && Arblib.contains_nonpositive(Arblib.ref(x, 0))
        for i = 0:Arblib.degree(res)
            res[i] = NaN
        end
        return res
    else
        contained_zero = false
    end

    Arblib.zero!(res)

    len = length(x)

    for (i, exponent, y) in expansion
        # Handle the constant term directly
        if iszero(i) && iszero(exponent)
            if T == Arb
                Arblib.add!(res, res, y)
            elseif T == ArbSeries
                if iszero(res)
                    Arblib.set_coeff!(res, 0, y)
                else
                    res0 = Arblib.ref(res, 0)
                    Arblib.add!(res0, res0, y)
                    Arblib.normalise!(res)
                end
            end
            continue
        end

        # term_buffer = x^exponent
        if T == Arb
            Arblib.pow!(term_buffer, x, exponent)
        elseif T == ArbSeries
            Arblib.pow_arb_series!(term_buffer, x, exponent, len)
        end

        if i == -1
            # term_buffer *= inv(log(x))
            if T == Arb
                Arblib.mul!(term_buffer, term_buffer, invlogx)
            elseif T == ArbSeries
                Arblib.mullow!(term_buffer, term_buffer, invlogx, len)
            end
        elseif i != 0
            # term_buffer *= log(x)^i
            throw(ArgumentError("i = $i: this should not occur in the cases we care about"))
        end

        Arblib.mul!(term_buffer, term_buffer, y)

        # Use that each term is monotonically increasing and zero at x
        # = 0. This doesn't hold for the constant term but this case
        # we have already handled above.
        if contained_zero
            Arblib.union!(term_buffer, term_buffer, z)
        end

        # res += term_buffer
        if T == Arb
            Arblib.add!(res, res, term_buffer)
        elseif T == ArbSeries
            Arblib.add_series!(res, res, term_buffer, len)
        end
    end

    return res
end

function (u0::BHAnsatz{T})(x, ::Ball) where {T}
    # Main term
    res = u0.a0 * clausencmzeta(x, 2, 1)

    # Clausen terms
    for j = 1:u0.N0
        s = 1 - u0.α + j * u0.p0
        res += u0.a[j] * clausencmzeta(x, s)
    end

    # Fourier terms
    for n = 1:u0.N1
        res += u0.b[n] * (cos(n * x) - 1)
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

The highest term, `x^2M`, is a remainder term is which makes sure that
evaluation of the expansion gives an enclosure of the result.
"""
function (u0::BHAnsatz{Arb})(x, ::AsymptoticExpansion; M::Integer = 3)
    res = OrderedDict{NTuple{4,Int},Arb}()

    # Remainder term
    res[(0, 2M, 0, 0)] = 0

    # Main term
    γ = Arb(Irrational{:γ}())
    res[(1, 1, 0, 0)] = -Arb(π) / 2 * u0.a0
    res[(0, 1, 0, 0)] = -(γ - 1) * Arb(π) / 2 * u0.a0
    for m = 1:M-1
        res[(0, 2m, 0, 0)] = (-1)^m * dzeta(Arb(2 - 2m)) / factorial(2m) * u0.a0
    end
    res[(0, 2M, 0, 0)] += clausenc_expansion_remainder(x, Arb(2), 1, M)

    # Clausens terms
    for j = 1:u0.N0
        s = 1 - u0.α + j * u0.p0
        C, _, p, E = clausenc_expansion(x, s, M)
        res[(0, 0, 1, j)] = C * u0.a[j]
        for m = 1:M-1
            res[(0, 2m, 0, 0)] += p[2m] * u0.a[j]
        end
        res[(0, 2M, 0, 0)] += E * u0.a[j]
    end

    # Fourier terms
    if !iszero(u0.N1)
        for m = 1:M-1
            res[(0, 2m, 0, 0)] +=
                (-1)^m * sum(Arb(n)^(2m) * u0.b[n] for n = 1:u0.N1) / factorial(2m)
        end
        Arblib.add_error!(
            res[(0, 2M, 0, 0)],
            sum(Arb(n)^(2M) * abs(u0.b[n]) for n = 1:u0.N1) / factorial(2M),
        )
    end

    return res
end

function H(u0::BHAnsatz{T}, ::Ball) where {T}
    return x -> begin
        res = -u0.a0 * clausencmzeta(x, 3, 1)

        # Clausen terms
        for j = 1:u0.N0
            s = 2 - u0.α + j * u0.p0
            res -= u0.a[j] * clausencmzeta(x, s)
        end

        # Fourier terms
        for n = 1:u0.N1
            res -= u0.b[n] / n * (cos(n * x) - 1)
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

The highest term, `x^2M`, is a remainder term is which makes sure that
evaluation of the expansion gives an enclosure of the result.
"""
function H(u0::BHAnsatz{Arb}, ::AsymptoticExpansion; M::Integer = 3)
    γ = Arb(Irrational{:γ}())
    γ₁ = stieltjes(Arb, 1)

    return x -> begin
        res = OrderedDict{NTuple{4,Int},Arb}()

        # Remainder term
        res[(0, 2M, 0, 0)] = 0

        # Main term
        res[(2, 2, 0, 0)] = -1 // 4 * u0.a0
        res[(1, 2, 0, 0)] = (3 - 2γ) / 4 * u0.a0
        for m = 1:M-1
            if m == 1
                # Note that we divide this by factorial(2)
                term = (36γ - 12γ^2 - 24γ₁ - 42 + Arb(π)^2) / 24
            else
                term = dzeta(Arb(3 - 2m))
            end
            res[(0, 2m, 0, 0)] = -(-1)^m * term * u0.a0 / factorial(2m)
        end

        res[(0, 2M, 0, 0)] += clausenc_expansion_remainder(x, Arb(3), 1, M)

        # Clausen terms
        for j = 1:u0.N0
            C, _, p, E = clausenc_expansion(x, 2 - u0.α + j * u0.p0, M)
            res[(0, 1, 1, j)] = -C * u0.a[j]
            for m = 1:M-1
                res[(0, 2m, 0, 0)] -= p[2m] * u0.a[j]
            end
            res[(0, 2M, 0, 0)] += E * u0.a[j]
        end

        # Fourier terms
        if !iszero(u0.N1)
            for m = 1:M-1
                res[(0, 2m, 0, 0)] -=
                    (-1)^m * sum(Arb(n)^(2m - 1) * u0.b[n] for n = 1:u0.N1) / factorial(2m)
            end
            Arblib.add_error!(
                res[(0, 2M, 0, 0)],
                sum(Arb(n)^(2M - 1) * abs(u0.b[n]) for n = 1:u0.N1) / factorial(2M),
            )
        end

        return res
    end
end

function D(u0::BHAnsatz, ::Asymptotic; M::Integer = 3)
    f = D(u0, AsymptoticExpansion(); M)

    return x -> eval_expansion(u0, f(x), x)
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
        @assert Arblib.contains_zero(expansion[(2, 2, 0, 0)])
        delete!(expansion, (2, 2, 0, 0))

        return expansion
    end
end

"""
    F0(u0::BHAnsatz{Arb}, ::Asymptotic)

Returns a function such that `F0(u0, Asymptotic())(x)` is computed
accurately for small values of `x`.

This method is one of the bottlenecks of [`delta0_bound`](@ref) and
for this reasons it is heavily optimized. The method assumes that `x
>= 0` and any negative parts of the ball `x` are ignored.

It precomputes the expansions of `u0` and `D(u0)` and for that reason
a number `ϵ` has to be given, the resulting expansion will be valid
for all `x < ϵ`. The value of `ϵ` has to be less than `1`.

If `exponent_limit` is set to some number then all terms in the
expansions with an exponent equal to or greater than this limit will
be collapsed into term which is treated as an error. More precisely
any term of the form
```
y * log(abs(x))^i * abs(x)^e
```
for which `e >= exponent_limit` and `i == 0 || i == -1` will be
rewritten as `z * x^exponent_limit` where `z` is an enclosure of `y *
log(abs(x))^i * abs(x)^(e - exponent_limit)` on the interval `[0, ϵ]`.
This interval is computed by using that `abs(x)^(e - exponent_limit) ∈
[0, 1]` for `abs(x) < 1` and that `log(x)^(-1) ∈ [inv(log(ϵ)), 0]`.
Notice that for the last interval to not be huge we need `ϵ` to be
bounded away from `1`.
"""
function F0(
    u0::BHAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 3,
    ϵ::Arb = Arb(1e-2),
    exponent_limit::Union{Arb,Nothing} = nothing,
)
    @assert ϵ < 1

    u0_expansion = u0(ϵ, AsymptoticExpansion(); M)
    Du0_expansion = D(u0, AsymptoticExpansion(); M)(ϵ)

    # Divide the expansion by x * log(x) and x^2 * log(x)
    # respectively, also precompute the exponents.
    u0_expansion_div_xlogx = Vector{Tuple{Int,Arb,Arb}}(undef, length(u0_expansion))
    Du0_expansion_div_x2logx = Vector{Tuple{Int,Arb,Arb}}(undef, length(Du0_expansion))
    for (index, ((i, m, k, l), value)) in enumerate(u0_expansion)
        u0_expansion_div_xlogx[index] = (i - 1, -k * u0.α + l * u0.p0 + m - 1, value)
    end
    for (index, ((i, m, k, l), value)) in enumerate(Du0_expansion)
        Du0_expansion_div_x2logx[index] = (i - 1, -k * u0.α + l * u0.p0 + m - 2, value)
    end

    # If an exponent limit is given, collapse all terms with an
    # exponent greater than or equal to that.
    if !isnothing(exponent_limit)
        # A ball representing [0, 1]
        unit_interval = Arblib.unit_interval!(Arb())
        # A ball representing the interval [inv(log(ϵ)), 0]
        invlog_interval = Arb((inv(log(ϵ)), 0))

        # Compute the coefficient in front of the error term for u0
        coeff_u0 = zero(Arb)
        for (i, exponent, y) in u0_expansion_div_xlogx
            if exponent >= exponent_limit && (i == 0 || i == -1)
                z = y * unit_interval
                if i == -1
                    z *= invlog_interval
                end
                coeff_u0 += z
            end
        end
        # Filter out the terms included in the error term for u0
        filter!(
            ((i, exponent, _),) -> !(exponent >= exponent_limit && (i == 0 || i == -1)),
            u0_expansion_div_xlogx,
        )
        # Add error term to expansion
        push!(u0_expansion_div_xlogx, (0, exponent_limit, coeff_u0))

        # Compute the coefficient in front of the error term for Du0
        coeff_Du0 = zero(Arb)
        for (i, exponent, y) in Du0_expansion_div_x2logx
            if exponent >= exponent_limit && (i == 0 || i == -1)
                z = y * unit_interval
                if i == -1
                    z *= invlog_interval
                end
                coeff_Du0 += z
            end
        end
        # Filter out the terms included in the error term for Du0
        filter!(
            ((i, exponent, _),) -> !(exponent >= exponent_limit && (i == 0 || i == -1)),
            Du0_expansion_div_x2logx,
        )
        # Add error term to expansion
        push!(Du0_expansion_div_x2logx, (0, exponent_limit, coeff_Du0))
    end

    # The function sqrt(log(1 + inv(x)))
    w!(res::Arb, x::Arb) = begin
        Arblib.inv!(res, x)
        Arblib.add!(res, res, 1)
        Arblib.log!(res, res)
        return Arblib.sqrt!(res, res)
    end
    w!(res::ArbSeries, x::ArbSeries) = begin
        len = length(x)
        Arblib.inv_series!(res, x, len)
        Arblib.add!(res, res, 1)
        Arblib.log_series!(res, res, len)
        return Arblib.sqrt_series!(res, res, len)
    end

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && Arblib.ref(x, 0) <= ϵ)

        tmp = zero(x)
        buffer = zero(x)

        if x isa Arb
            # invlogx = inv(log(x))
            invlogx = log(x)
            invlogx = Arblib.inv!(invlogx, invlogx)

            res = zero(x)

            # res = inv(sqrt(log(1 + inv(x))))
            if iszero(x)
                Arblib.zero!(x)
            elseif Arblib.contains_zero(x)
                w!(res, ubound(Arb, x))
                Arblib.inv!(res, res)
                return Arblib.union!(res, res, zero(res))
            else
                w!(res, x)
                Arblib.inv!(res, res)
            end

            Arblib.div!(
                res,
                res,
                _eval_expansion!(tmp, u0, u0_expansion_div_xlogx, x, invlogx, buffer),
            )
            Arblib.mul!(
                res,
                res,
                _eval_expansion!(tmp, u0, Du0_expansion_div_x2logx, x, invlogx, buffer),
            )
        elseif x isa ArbSeries
            len = length(x)

            # invlogx = inv(log(x))
            invlogx = log(x)
            invlogx = Arblib.inv_series!(invlogx, invlogx, len)

            # res = sqrt(log(1 + inv(x)))
            res = w!(zero(x), x)

            Arblib.mullow!(
                res,
                res,
                _eval_expansion!(tmp, u0, u0_expansion_div_xlogx, x, invlogx, buffer),
                len,
            )
            Arblib.div_series!(
                res,
                _eval_expansion!(tmp, u0, Du0_expansion_div_x2logx, x, invlogx, buffer),
                res,
                len,
            )
        end

        return res
    end
end

"""
    inv_u0_normalised(u0::BHAnsatz{Arb}; M = 3, ϵ = 1 / 2)

Return a function for evaluation `-abs(x) * log(abs(x)) / u0(x)` for
`x` close to zero.

# Arguments
- `M::Integer = 3` determines the number of terms in the asymptotic
  expansions.
- `ϵ::Arb = 1 / 2` determines the interval ``[-ϵ, ϵ]`` on which the
  expansion is valid. Must be less than `1`.

# Implementation
It computes an expansion of `u0` at `x = 0` and explicitly handles the
cancellation with `-abs(x) * log(abs(x))`.
"""
function inv_u0_normalised(u0::BHAnsatz{Arb}; M::Integer = 3, ϵ::Arb = Arb(1 // 2))
    0 < ϵ < 1 || throw(DomainError(ϵ, "must have 0 < ϵ < 1"))

    expansion = u0(ϵ, AsymptoticExpansion(); M)

    # Divide all terms in the expansion by abs(x) * log(abs(x))
    expansion_div_xlogx = empty(expansion)
    for ((i, m, k, l), value) in expansion
        expansion_div_xlogx[(i - 1, m - 1, k, l)] = value
    end

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && abs(x) <= ϵ) || (x isa ArbSeries && abs(Arblib.ref(x, 0)) <= ϵ)

        return -inv(eval_expansion(u0, expansion_div_xlogx, x))
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
    u0_xs_b_precomputed = zeros(length(xs), u0.N1)
    Hu0_xs_a0_precomputed = zeros(length(xs))
    Hu0_xs_b_precomputed = zeros(length(xs), u0.N1)

    Hu0 = H(u0)
    for (i, x) in enumerate(xs)
        u0_xs_a0_precomputed[i] = u0(x)
        Hu0_xs_a0_precomputed[i] = Hu0(x)

        for n = 1:u0.N1
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
