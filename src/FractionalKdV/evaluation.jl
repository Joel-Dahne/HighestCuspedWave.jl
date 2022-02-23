# This file contains code for evaluation of the approximate solution
# in different ways

export hat, eval_expansion

"""
    eval_expansion(u0::FractionalKdVAnsatz, expansion, x; offset_i = 0, offset = 0)

Evaluate the given expansion. The term `((i, j,  m), y)` is evaluated
to
```
y*abs(x)^(-i * u0.α + j * u0.p0 + m)
```
and then they are all summed.

In general `x` needs to be given both when computing the expansion and
when evaluating it.

The arguments `offset_i` and `offset` can be set to adjust the
exponent, in that case the exponent will be given by
```
-(i + offset_i) * u0.α + j * u0.p0 + m + offset
```
"""
function eval_expansion(
    u0::FractionalKdVAnsatz{T},
    expansion::AbstractDict{NTuple{3,Int},T},
    x;
    offset_i::Integer = 0,
    offset = 0,
) where {T}
    res = zero(u0.α)

    for ((i, j, m), y) in expansion
        if !iszero(y)
            exponent = -(i + offset_i) * u0.α + j * u0.p0 + m + offset

            res += y * abspow(x, exponent)
        end
    end

    return res
end

function (u0::FractionalKdVAnsatz)(x, ::Ball)
    res = zero(u0.α)

    for j = 0:u0.N0
        s = 1 - u0.α + j * u0.p0
        res += u0.a[j] * clausencmzeta(x, s)
    end

    for n = 1:u0.N1
        res += u0.b[n] * (cos(n * x) - 1)
    end

    return res
end

(u0::FractionalKdVAnsatz)(x, ::Asymptotic; M::Integer = 3) =
    eval_expansion(u0, u0(x, AsymptoticExpansion(); M), x)

function (u0::FractionalKdVAnsatz{Arb})(x, ::AsymptoticExpansion; M::Integer = 3)
    res = OrderedDict{NTuple{3,Int},Arb}()

    # Initiate even powers of x
    for m = 1:M
        res[(0, 0, 2m)] = 0
    end

    # Clausen terms
    for j = 0:u0.N0
        s = 1 - u0.α + j * u0.p0
        C, _, p, E = clausenc_expansion(x, s, M)

        # Check for the special case when s overlaps with an odd
        # integer. Notice that we don't need to do anything special if
        # s happens to overlap with two integers, we will just get NaN
        # as a result.
        contains_int, n = unique_integer(s)
        if contains_int && isodd(n)
            # The term corresponding to C and p[2((n - 1) ÷ 2)]
            # coincides and diverge so are handled separately. The
            # rest we treat normally.
            @assert !isfinite(C) && !isfinite(p[2((n-1)÷2)])
            for m = 1:M-1
                m == (n - 1) ÷ 2 && continue # Skip this term
                res[(0, 0, 2m)] += p[2m] * u0.a[j]
            end

            # TODO: The approach below gives good bounds for x close
            # to 0 bound much worse than required if x is not that
            # small. It might be that we need to improve it in the end.

            # The term has an x-factor like x^(s - 1), we factor out
            # x^(-i * u0.α + 1) from this where i is as large as
            # possible but so that we still have -i * u0.α + 1 < s -
            # 1, and enclose the rest.
            i = findfirst(i -> !(-i * u0.α + 1 < s - 1), 1:10) - 1
            D = clausenc_expansion_odd_s_singular(x, s, -i * u0.α + 1)
            res[(i, 0, 1)] = get(res, (1, 0, 1), zero(x)) + D

            @info "Encountered a term with s overlapping an integer in expansion for u0" s i (
                -i * u0.α + 1
            ) D
        else
            res[(1, j, 0)] = C * u0.a[j]
            for m = 1:M-1
                res[(0, 0, 2m)] += p[2m] * u0.a[j]
            end
        end

        Arblib.add_error!(res[(0, 0, 2M)], E)
    end

    # Fourier terms
    if !iszero(u0.N1)
        for m = 1:M-1
            res[(0, 0, 2m)] +=
                (-1)^m * sum(Arb(n)^(2m) * u0.b[n] for n = 1:u0.N1) / factorial(2m)
        end
        Arblib.add_error!(
            res[(0, 0, 2M)],
            sum(Arb(n)^(2M) * abs(u0.b[n]) for n = 1:u0.N1) / factorial(2M),
        )
    end

    return res
end

function H(u0::FractionalKdVAnsatz, ::Ball)
    return x -> begin
        res = zero(u0.α)

        for j = 0:u0.N0
            s = 1 - 2u0.α + j * u0.p0
            res -= u0.a[j] * clausencmzeta(x, s)
        end

        for n = 1:u0.N1
            res -= u0.b[n] * n^u0.α * (cos(n * x) - 1)
        end

        return res
    end
end

function H(u0::FractionalKdVAnsatz, ::Asymptotic; M::Integer = 3)
    f = H(u0, AsymptoticExpansion(); M)
    return x -> eval_expansion(u0, f(x), x)
end

function H(u0::FractionalKdVAnsatz{T}, ::AsymptoticExpansion; M::Integer = 3) where {T}
    return x -> begin
        res = OrderedDict{NTuple{3,Int},Arb}()

        # Initiate even powers of x
        for m = 1:M
            res[(0, 0, 2m)] = 0
        end

        # Clausen terms
        for j = 0:u0.N0
            s = 1 - 2u0.α + j * u0.p0
            C, _, p, E = clausenc_expansion(x, s, M)

            # Check for the special case when s overlaps with an odd
            # integer.
            contains_int, n = unique_integer(s)
            if contains_int && isodd(n)
                # TODO: It seems like we should never encounter this.
                # Maybe we do and in that case we have to deal with
                # it, by it seems like we might be able to avoid it :)
                @error "Encountered a term with s overlapping an odd integer" *
                       "in expansion for H(u0), we don't expect this to happen" s
            end

            res[(2, j, 0)] = -C * u0.a[j]
            for m = 1:M-1
                res[(0, 0, 2m)] -= p[2m] * u0.a[j]
            end
            res[(0, 0, 2M)] += E * u0.a[j]
        end

        # Fourier terms
        if !iszero(u0.N1)
            for m = 1:M-1
                res[(0, 0, 2m)] -=
                    (-1)^m * sum(Arb(n)^(2m + u0.α) * u0.b[n] for n = 1:u0.N1) /
                    factorial(2m)
            end
            Arblib.add_error!(
                res[(0, 0, 2M)],
                sum(Arb(n)^(2M + u0.α) * abs(u0.b[n]) for n = 1:u0.N1) / factorial(2M),
            )
        end

        return res
    end
end

function D(u0::FractionalKdVAnsatz, ::Asymptotic; M::Integer = 3)
    f = D(u0, AsymptoticExpansion(); M)
    return x -> eval_expansion(u0, f(x), x)
end

function D(
    u0::FractionalKdVAnsatz{T},
    evaltype::AsymptoticExpansion;
    M::Integer = 3,
) where {T}
    f = H(u0, evaltype; M)
    return x -> begin
        expansion1 = u0(x, evaltype; M)
        expansion2 = f(x)

        expansion = empty(expansion1)

        # u0^2/2 term
        let expansion1 = collect(expansion1)
            z = zero(u0.α) # Avoid allocating zero multiple times
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

        # Terms in u0.zeroterms are supposed to be identically equal
        # to zero
        for (i, j, m) in u0.zeroterms
            expansion[(i, j, m)] = zero(u0.α)
        end

        return expansion
    end
end

"""
    F0(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M = 3, ϵ = one(Arb))

Return a function for evaluating `F0(u0)(x)` accurately for small
values of `x`.

It splits `F0(u0)` as
```
inv(u0(x) / x^-u0.α) * (D(u0)(x) / x^(u0.p - u0.α))
```
It computes `inv(u0(x) / x^-u0.α)` using [`inv_u0_normalised`](@ref).
For the other factor it computes the expansion of `D(u0)(x)` and
explicitly cancels the division by `x^(u0.p - u0.α)`.

# Arguments
- `M::Integer` determines the number of terms in the asymptotic
  expansions.
- `ϵ::Arb` determines the interval ``[-ϵ, ϵ]`` on which the expansion
  is valid.
"""
function F0(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M::Integer = 3, ϵ::Arb = Arb(1))
    Du0_expansion = D(u0, AsymptoticExpansion(); M)(ϵ)

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && Arblib.ref(x, 0) <= ϵ)

        res = eval_expansion(u0, Du0_expansion, x, offset = -u0.p, offset_i = -1)

        return res * inv_u0(x)
    end
end

"""
    inv_u0_normalised(u0::FractionalKdVAnsatz{Arb}; M = 3, ϵ = one(Arb))

Return a function for evaluation `x^-u0.α / u0(x)` for `x` close to
zero.

It computes an expansion of `u0` at `x = 0` and explicitly handles the
cancellation with `x^-u0.α`.

# Arguments
- `M::Integer` determines the number of terms in the asymptotic
  expansions.
- `ϵ::Arb` determines the interval ``[-ϵ, ϵ]`` on which the expansion
  is valid.
"""
function inv_u0_normalised(u0::FractionalKdVAnsatz{Arb}; M::Integer = 3, ϵ::Arb = one(Arb))
    expansion = u0(ϵ, AsymptoticExpansion(); M)

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && abs(x) <= ϵ) || (x isa ArbSeries && abs(Arblib.ref(x, 0)) <= ϵ)

        return inv(eval_expansion(u0, expansion, x, offset_i = -1))
    end
end

"""
    hat(u0::FractionalKdVAnsatz)
Returns a function such that hat(u0)(x) computes û(x) from the paper.
"""
function hat(u0::FractionalKdVAnsatz, ::Ball = Ball())
    return x -> begin
        (a0(u0, 0) * abs(x)^(-u0.α) - u0(x)) / u0(x)
    end
end

"""
    c(u0::FractionalKdVAnsatz{T}, ϵ)

Compute the constant ``c_{\\epsilon,\\hat{u}_0}`` from Lemma 3.3.

This constant satisfies that for all `abs(x) < ϵ` we have
```
abs(hat(u0)(x)) <= c(u0, ϵ) * abs(x)^u0.p0
```
where
```
hat(u0)(x) = (a0(u0, 0) * abs(x)^-u0.α - u0(x)) / u0(x)
```

To compute `c` we first compute the asymptotic expansion of `u0` and
from this we compute an asymptotic expansion for the numerator and
denominator respectively. More precisely we rewrite it as
```
((a0(u0, 0) * abs(x)^-u0.α - u0(x)) / abs(x)^(-u0.α + u0.p0)) / (u0(x) / abs(x)^-u0.α) * abs(x)^u0.p0
```
and it's enough that we compute an upper bound for
```
(a0(u0, 0) * abs(x)^-u0.α - u0(x)) / abs(x)^(-u0.α + u0.p0)
```
and a lower bound for
```
u0(x) / abs(x)^-u0.α
```
to get `c`.

We can compute the expansion for
```
(a0(u0, 0) * abs(x)^-u0.α - u0(x)) / abs(x)^(-u0.α + u0.p0)
```
by taking the expansion of `u0`, removing the term `(1, 0, 0)`
corresponding to `a0(u0, 0) * abs(x)^-u0.α` subtracting `-u0.α +
u0.p0` from the exponents and negating all of them. To get an upper
bound we can take the absolute value of all coefficients. The
resulting expansion gives an upper bound and is increasing for `x > 0`
so can be bounded on `[-ϵ, ϵ]` by evaluating it at `x = ϵ`.

To get the expansion of
```
u0(x) / abs(x)^-u0.α
```
we only have to subtract `-u0.α` from all the exponents. To get a
lower bound we notice that the leading term is a constant. By taking
the absolute value of the leading term and minus the absolute value of
all the other terms we get something which gives a lower bound and is
decreasing for `x > 0` so can be lower bounded on `[-ϵ, ϵ]` by
evaluating it at `x = ϵ`.
"""
function c(u0::FractionalKdVAnsatz{Arb}, ϵ::Arb; M::Integer = 3)
    iszero(ϵ) && return zero(ϵ)

    @assert ϵ > 0

    expansion = u0(ϵ, AsymptoticExpansion(); M)

    # Set up expansion for numerator. Skip the term (1, 0, 0) and for
    # the other terms subtract -α + p0 from the exponent (by
    # subtracting 1 from i and j) and take the absolute value of y.
    expansion_numerator = empty(expansion)
    for ((i, j, m), y) in expansion
        (i, j, m) == (1, 0, 0) && continue
        expansion_numerator[(i - 1, j - 1, m)] = abs(y)
    end

    # Set up the expansion for the numerator. Subtract -α from the
    # exponent (by subtracting 1 from i) and set the coefficients to
    # minus the absolute value except for the constant term (0, 0, 0)
    # which is set to the absolute value.
    expansion_denominator = empty(expansion)
    for ((i, j, m), y) in expansion
        expansion_denominator[(i - 1, j, m)] = -abs(y)
    end
    expansion_denominator[(0, 0, 0)] = abs(expansion_denominator[(0, 0, 0)])

    # Evaluate the expansion at x = ϵ
    numerator = eval_expansion(u0, expansion_numerator, ϵ)
    denominator = eval_expansion(u0, expansion_denominator, ϵ)

    Arblib.isnegative(denominator) &&
        error("didn't expected denominator to be negative, got $denominator")

    return numerator / denominator
end

"""
    D(u0::FractionalKdVAnsatz, xs::AbstractVector)
Returns a function such that D(u0, xs)(a, b) computes D(u0)(x) on the
points x ∈ xs with u0.a and u0.b set to the given values. Does this in
an efficient way by precomputing as much as possible.
"""
function D(u0::FractionalKdVAnsatz, xs::AbstractVector)
    u0_xs_a_precomputed = zeros(length(xs), u0.N0 + 1)
    u0_xs_b_precomputed = zeros(length(xs), u0.N1)
    Hu0_xs_a_precomputed = zeros(length(xs), u0.N0 + 1)
    Hu0_xs_b_precomputed = zeros(length(xs), u0.N1)

    for i in eachindex(xs)
        x = xs[i]
        for j = 0:u0.N0
            u0_xs_a_precomputed[i, j+1] = clausencmzeta(x, 1 - u0.α + j * u0.p0)
            Hu0_xs_a_precomputed[i, j+1] = -clausencmzeta(x, 1 - 2u0.α + j * u0.p0)
        end
        for n = 1:u0.N1
            u0_xs_b_precomputed[i, n] = cos(n * x) - 1
            Hu0_xs_b_precomputed[i, n] = -n^u0.α * (cos(n * x) - 1)
        end
    end

    return (a, b) -> begin
        return (
            (u0_xs_a_precomputed * a .+ u0_xs_b_precomputed * b) .^ 2 ./ 2 .+
            (Hu0_xs_a_precomputed * a .+ Hu0_xs_b_precomputed * b)
        )
    end
end

"""
    D(u0::FractionalKdVAnsatz, evaltype::Symbolic; M::Integer = 5)

Return a function such that `D(u0, evaltype, N)(a)` computes the
coefficients in the asymptotic expansion with indices `3` to `u0.N0 +
1` using the values from `a`.

This is used in [`_findas`](@ref) for numerically finding values for
`a`.
"""
function D(u0::FractionalKdVAnsatz{T}, ::Symbolic; M::Integer = 5) where {T}
    # Given a key get its exponent
    key_exponent = ((i, j, m),) -> -i * u0.α + j * u0.p0 + m

    # Precompute for u0
    u0_precomputed = OrderedDict{NTuple{3,Int},OrderedDict{Int,T}}()
    for j = 0:u0.N0
        s = u0.α - j * u0.p0
        u0_precomputed[(1, j, 0)] = OrderedDict(j => gamma(s) * sinpi((1 - s) / 2))
    end

    for m = 1:M-1
        u0_precomputed[(0, 0, 2m)] = OrderedDict(
            j => (-1)^m * zeta(1 - u0.α + j * u0.p0 - 2m) / factorial(2m) for j = 0:u0.N0
        )
    end

    # Precompute H(u0)
    Hu0_precomputed = OrderedDict{NTuple{3,Int},OrderedDict{Int,T}}()

    for j = 0:u0.N0
        s = 2u0.α - j * u0.p0
        Hu0_precomputed[(2, j, 0)] = OrderedDict(j => -gamma(s) * sinpi((1 - s) / 2))
    end

    for m = 1:M-1
        Hu0_precomputed[(0, 0, 2m)] = OrderedDict(
            j => -(-1)^m * zeta(1 - 2u0.α + j * u0.p0 - 2m) / factorial(2m) for j = 0:u0.N0
        )
    end

    # We want to find the value of the maximum exponent we return. To
    # do this we first compute ALL keys we will encounter and then
    # sort and find the value for the maximum one we return.
    all_keys = OrderedDict{NTuple{3,Int},T}()
    for key in keys(Hu0_precomputed)
        if !haskey(all_keys, key)
            all_keys[key] = key_exponent(key)
        end
    end
    for key1 in keys(u0_precomputed)
        for key2 in keys(u0_precomputed)
            key = key1 .+ key2
            if !haskey(all_keys, key)
                all_keys[key] = key_exponent(key)
            end
        end
    end
    # Find the keys we care about
    returned_keys = collect(sort(all_keys, by = Float64 ∘ key_exponent))[3:u0.N0+2]
    # Find the largest key/exponent we care about
    maxkey, maxexponent = returned_keys[end]

    # Check that M is large enough
    @assert maxexponent < 2M

    # Filter out any keys larger than the max exponent
    filter!(keyvalue -> !(key_exponent(keyvalue[1]) > maxexponent), u0_precomputed)
    filter!(keyvalue -> !(key_exponent(keyvalue[1]) > maxexponent), Hu0_precomputed)

    # Sort the dictionaries by exponent
    sort!(u0_precomputed, by = Float64 ∘ key_exponent)
    sort!(Hu0_precomputed, by = Float64 ∘ key_exponent)

    # Function to compute the dictionaries u0_res and H0_res
    sum_dict(precomputed, a, S) = begin
        res = empty(precomputed, S)
        @inbounds for (key, dict) in precomputed
            for (j, v) in dict
                res[key] = get(res, key, zero(S)) + v * a[j]
            end
        end
        return res
    end

    return a -> begin
        S = promote_type(T, eltype(a))

        u0_res = sum_dict(u0_precomputed, a, S)
        Hu0_res = sum_dict(Hu0_precomputed, a, S)

        # Compute u0^2/2
        res = empty(u0_res)
        u0_res_vector = collect(u0_res)
        @inbounds for (i, (key1, y1)) in enumerate(u0_res_vector)
            key = 2 .* key1
            key_exponent(key) > maxexponent && continue
            res[key] = get(res, key, zero(S)) + y1^2 / 2
            for j = i+1:length(u0_res_vector)
                (key2, y2) = u0_res_vector[j]
                key = key1 .+ key2
                key_exponent(key) > maxexponent && break
                res[key] = get(res, key, zero(S)) + y1 * y2
            end
        end

        # Compute u0^2/2 + H(u0)
        merge!(+, res, Hu0_res)

        return collect(values(sort(res, by = Float64 ∘ key_exponent)))[3:u0.N0+2]
    end
end

"""
    print_asymptotic_expansion_D(u0::FractionalKdVAnsatz, expansion)

Debug method for printing an asymptotic expansion in a human readable
way.
"""
function print_asymptotic_expansion_D(u0::FractionalKdVAnsatz, expansion)
    get_exponent(i, j, m) = -i * u0.α + j * u0.p0 + m
    expansion = sort(expansion, by = x -> Float64(get_exponent(x...)))
    for ((i, j, m), c) in expansion
        println("$c ⋅ x^$(get_exponent(i, j, m))   $((i, j, m))")
    end

    return expansion
end

"""
    print_asymptotic_expansion_F0(u0::FractionalKdVAnsatz, expansion)

Debug method for printing an asymptotic expansion in a human readable
way.
"""
function print_asymptotic_expansion_F0(u0::FractionalKdVAnsatz, expansion)
    get_exponent(i, j, m) = -i * u0.α + j * u0.p0 + m - u0.p
    expansion = sort(expansion, by = x -> Float64(get_exponent(x...)))
    for ((i, j, m), c) in expansion
        println("$c ⋅ x^$(get_exponent(i, j, m))   $((i, j, m))")
    end

    return expansion
end
