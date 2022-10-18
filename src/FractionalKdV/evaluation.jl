# This file contains code for evaluation of the approximate solution
# in different ways

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

If `u0_skipped_u0_main` is true as well as `u0.use_bhkdv` then in
addition to the terms in `expansion` also add the term
```
u0.a[0] * (C1 - C2 * abs(x)^u0.p0) * abs(x)^(-(1 + offset_i) * u0.α + offset)
```
where `C1` and `C2` are the leading terms in the asymptotic expansions
of `clausencmzeta(x, 1 - u0.α)` and `clausencmzeta(x, 1 - u0.α +
u0.p0)` respectively. This is mean to be used when the expansion was
computed with `u0(x, AsympototicExpansion(), bhkdv_skip_main = true)`
and gives better enclosures for `α` close to `-1`.
"""
function eval_expansion(
    u0::FractionalKdVAnsatz{T},
    expansion::AbstractDict{NTuple{3,Int},T},
    x;
    offset_i::Integer = 0,
    offset = 0,
    bhkdv_skipped_u0_main = false,
) where {T}
    res = zero(u0.α)

    if u0.use_bhkdv && bhkdv_skipped_u0_main
        s = 1 - u0.α
        # The first argument to clausenc_expansion is not used in this
        # case, so we can take it to be zero.
        C1 = clausenc_expansion(Arb(0), s, 3)[1]
        C2 = clausenc_expansion(Arb(0), s + u0.p0, 3)[1]
        exponent = -(1 + offset_i) * u0.α + offset
        res += u0.a[0] * (C1 - C2 * abspow(x, u0.p0)) * abspow(x, exponent)
    end

    for ((i, j, m), y) in expansion
        if !iszero(y)
            exponent = -(i + offset_i) * u0.α + j * u0.p0 + m + offset

            res += y * abspow(x, exponent)
        end
    end

    return res
end

"""
    clausencmzeta_diff(x, s, ϵ)

Compute
```
clausencmzeta(x, s) - clausencmzeta(x, s + ϵ)
```
in a way that works well when `ϵ` is small.

At the moment the only optimization is that for wide values of `x` it
uses a zero order approximation to get better enclosures.
"""
clausencmzeta_diff(x, s, ϵ) = clausencmzeta(x, s) - clausencmzeta(x, s + ϵ)

function clausencmzeta_diff(x::Arb, s::Arb, ϵ::Arb)
    if iswide(x)
        # Use a zero order approximation
        deriv = -(clausens(x, s - 1) - clausens(x, s + ϵ - 1))
        mid = midpoint(Arb, x)
        res =
            add_error(clausencmzeta(mid, s) - clausencmzeta(mid, s + ϵ), (x - mid) * deriv)
    else
        res = clausencmzeta(x, s) - clausencmzeta(x, s + ϵ)
    end

    return res
end

function clausencmzeta_diff(x::ArbSeries, s::Arb, ϵ::Arb)
    res = zero(x)
    x₀ = x[0]

    res[0] = clausencmzeta_diff(x₀, s, ϵ)
    if iswide(x₀)
        # Use a zero order approximation for each term
        mid = midpoint(Arb, x₀)
        for i = 1:Arblib.degree(x)
            if i % 2 == 0
                deriv = -(clausens(x₀, s - i - 1) - clausens(x₀, s + ϵ - i - 1))
                term = add_error(
                    clausenc(mid, s - i) - clausenc(mid, s + ϵ - i),
                    (x₀ - mid) * deriv,
                )
                res[i] = (-1)^(i ÷ 2) * term / factorial(i)
            else
                deriv = clausenc(x₀, s - i - 1) - clausenc(x₀, s + ϵ - i - 1)
                term = add_error(
                    clausens(mid, s - i) - clausens(mid, s + ϵ - i),
                    (x₀ - mid) * deriv,
                )
                res[i] = -(-1)^(i ÷ 2) * term / factorial(i)
            end
        end
    else
        for i = 1:Arblib.degree(x)
            if i % 2 == 0
                res[i] =
                    (-1)^(i ÷ 2) * (clausenc(x₀, s - i) - clausenc(x₀, s + ϵ - i)) /
                    factorial(i)
            else
                res[i] =
                    -(-1)^(i ÷ 2) * (clausens(x₀, s - i) - clausens(x₀, s + ϵ - i)) /
                    factorial(i)
            end
        end
    end

    # Compose the Taylor series for the result with that of the input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return res
end


function (u0::FractionalKdVAnsatz)(x, ::Ball)
    res = zero(u0.α)

    if u0.use_bhkdv
        s = 1 - u0.α
        res += u0.a[0] * clausencmzeta_diff(x, s, u0.p0)
    else
        s = 1 - u0.α
        res += u0.a[0] * clausencmzeta(x, s)
    end

    for j = 1:u0.N0
        s = 1 - u0.α + j * u0.p0
        res += u0.a[j] * clausencmzeta(x, s)
    end

    for n = 1:u0.N1
        res += u0.b[n] * (cos(n * x) - 1)
    end

    return res
end

(u0::FractionalKdVAnsatz)(x, ::Asymptotic; M::Integer = 5) =
    eval_expansion(u0, u0(x, AsymptoticExpansion(); M), x)

function (u0::FractionalKdVAnsatz{Arb})(
    x,
    ::AsymptoticExpansion;
    M::Integer = 5,
    bhkdv_skip_main = false,
)
    res = OrderedDict{NTuple{3,Int},Arb}()

    # Initiate even powers of x
    for m = 1:M
        res[(0, 0, 2m)] = 0
    end
    for j = 0:u0.N0
        res[(1, j, 0)] = 0
    end

    if u0.use_bhkdv
        s = 1 - u0.α
        C1, _, p1, E1 = clausenc_expansion(x, s, M)
        C2, _, p2, E2 = clausenc_expansion(x, s + u0.p0, M)
        if !bhkdv_skip_main
            # See the argument u0_skipped_u0_main to eval_expansion
            res[(1, 0, 0)] += C1 * u0.a[0]
            res[(1, 1, 0)] += -C2 * u0.a[0]
        end
        for m = 1:M-1
            res[(0, 0, 2m)] += (p1[2m] - p2[2m]) * u0.a[0]
        end
        Arblib.add_error!(res[(0, 0, 2M)], (E1 - E2) * u0.a[0])
    else
        s = 1 - u0.α
        C, _, p, E = clausenc_expansion(x, s, M)
        res[(1, 0, 0)] += C * u0.a[0]
        for m = 1:M-1
            res[(0, 0, 2m)] += p[2m] * u0.a[0]
        end
        Arblib.add_error!(res[(0, 0, 2M)], E * u0.a[0])
    end

    # Clausen terms
    for j = 1:u0.N0
        s = 1 - u0.α + j * u0.p0
        C, _, p, E = clausenc_expansion(x, s, M)

        # Below we handle the special case when s contains an odd
        # integer. When s is close to an odd integer we have very
        # large cancellations between two of the terms in the
        # expansion and it turns out to be beneficial to use the same
        # method as when s overlaps an odd integer to compute an
        # enclosure.
        if is_approx_integer(s) && isodd(round(Float64(s)))
            s = union(s, Arb(round(Float64(s))))
        end

        # Check for the special case when s overlaps with an odd
        # integer. Notice that we don't need to do anything special if
        # s happens to overlap with two integers, we will just get NaN
        # as a result.
        contains_int, n = unique_integer(s)
        if contains_int && isodd(n)
            # The term corresponding to C and p[2((n - 1) ÷ 2)]
            # coincides and diverge so are handled separately. The
            # rest we treat normally.
            for m = 1:M-1
                m == (n - 1) ÷ 2 && continue # Skip this term
                res[(0, 0, 2m)] += p[2m] * u0.a[j]
            end

            # IMPROVE: The approach below gives good bounds for x
            # close to 0 but much worse than required if x is not that
            # small.

            # The term has an x-factor like x^(s - 1), we factor out
            # x^(-i * u0.α + 1) from this where i is as large as
            # possible but so that we still have -i * u0.α + 1 < s -
            # 1, and enclose the rest.
            i = findfirst(i -> !(-i * u0.α + 1 < s - 1), 1:100) - 1
            D = clausenc_expansion_odd_s_singular(x, s, -i * u0.α + 1)
            res[(i, 0, 1)] = get(res, (i, 0, 1), zero(x)) + D * u0.a[j]
        else
            res[(1, j, 0)] += C * u0.a[j]
            for m = 1:M-1
                res[(0, 0, 2m)] += p[2m] * u0.a[j]
            end
        end

        Arblib.add_error!(res[(0, 0, 2M)], E * u0.a[j])
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

        if u0.use_bhkdv
            s = 1 - 2u0.α
            res -= u0.a[0] * (clausencmzeta(x, s) - clausencmzeta(x, s + u0.p0))
        else
            s = 1 - 2u0.α
            res -= u0.a[0] * clausencmzeta(x, s)
        end

        for j = 1:u0.N0
            s = 1 - 2u0.α + j * u0.p0
            res -= u0.a[j] * clausencmzeta(x, s)
        end

        for n = 1:u0.N1
            res -= u0.b[n] * n^u0.α * (cos(n * x) - 1)
        end

        return res
    end
end

function H(u0::FractionalKdVAnsatz, ::Asymptotic; M::Integer = 5)
    f = H(u0, AsymptoticExpansion(); M)
    return x -> eval_expansion(u0, f(x), x)
end

function H(u0::FractionalKdVAnsatz{T}, ::AsymptoticExpansion; M::Integer = 5) where {T}
    return x -> begin
        res = OrderedDict{NTuple{3,Int},Arb}()

        # Initiate even powers of x
        for m = 1:M
            res[(0, 0, 2m)] = 0
        end
        for j = 0:u0.N0
            res[(2, j, 0)] = 0
        end

        if u0.use_bhkdv
            s = 1 - 2u0.α
            C1, _, p1, E1 = clausenc_expansion(x, s, M)
            C2, _, p2, E2 = clausenc_expansion(x, s + u0.p0, M)
            res[(2, 0, 0)] -= C1 * u0.a[0]
            res[(2, 1, 0)] -= -C2 * u0.a[0]
            for m = 1:M-1
                res[(0, 0, 2m)] -= (p1[2m] - p2[2m]) * u0.a[0]
            end
            Arblib.add_error!(res[(0, 0, 2M)], (E1 - E2) * u0.a[0])
        else
            s = 1 - 2u0.α
            C, _, p, E = clausenc_expansion(x, s, M)
            res[(2, 0, 0)] -= C * u0.a[0]
            for m = 1:M-1
                res[(0, 0, 2m)] -= p[2m] * u0.a[0]
            end
            Arblib.add_error!(res[(0, 0, 2M)], E * u0.a[0])
        end

        # Clausen terms
        for j = 1:u0.N0
            s = 1 - 2u0.α + j * u0.p0
            C, _, p, E = clausenc_expansion(x, s, M)

            # Below we handle the special case when s contains an odd
            # integer. When s is close to an odd integer we have very
            # large cancellations between two of the terms in the
            # expansion and it turns out to be beneficial to use the same
            # method as when s overlaps an odd integer to compute an
            # enclosure.
            if is_approx_integer(s) && isodd(round(Float64(s)))
                s = union(s, Arb(round(Float64(s))))
            end

            # Check for the special case when s overlaps with an odd
            # integer.
            contains_int, n = unique_integer(s)
            if contains_int && isodd(n)
                # The term corresponding to C and p[2((n - 1) ÷ 2)]
                # coincides and diverge so are handled separately. The
                # rest we treat normally.
                for m = 1:M-1
                    m == (n - 1) ÷ 2 && continue # Skip this term
                    res[(0, 0, 2m)] -= p[2m] * u0.a[j]
                end

                # IMPROVE: The approach below gives good bounds for x
                # close to 0 but much worse than required if x is not that
                # small.

                # The term has an x-factor like x^(s - 1), we factor out
                # x^(-i * u0.α + 1) from this where i is as large as
                # possible but so that we still have -i * u0.α + 1 < s -
                # 1, and enclose the rest.
                i = findfirst(i -> !(-i * u0.α + 1 < s - 1), 1:100) - 1
                D = clausenc_expansion_odd_s_singular(x, s, -i * u0.α + 1)
                res[(i, 0, 1)] = get(res, (i, 0, 1), zero(x)) - D * u0.a[j]
            else
                res[(2, j, 0)] -= C * u0.a[j]
                for m = 1:M-1
                    res[(0, 0, 2m)] -= p[2m] * u0.a[j]
                end
            end

            Arblib.add_error!(res[(0, 0, 2M)], E * u0.a[j])
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

function D(u0::FractionalKdVAnsatz, ::Asymptotic; M::Integer = 5)
    f = D(u0, AsymptoticExpansion(); M)
    return x -> eval_expansion(u0, f(x), x)
end

function D(
    u0::FractionalKdVAnsatz{T},
    evaltype::AsymptoticExpansion;
    M::Integer = 5,
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
    F0(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M = 5, ϵ = one(Arb))

Return a function for evaluating `F0(u0)(x)` accurately for small
values of `x`.

# Arguments
- `M::Integer` determines the number of terms in the asymptotic
  expansions.
- `ϵ::Arb` determines the interval ``[-ϵ, ϵ]`` on which the expansion
  is valid.

# Implementation
It splits `F0(u0)` as
```
inv(u0(x) / x^-u0.α) * (D(u0)(x) / x^(u0.p - u0.α))
```
It computes `inv(u0(x) / x^-u0.α)` using [`inv_u0_normalised`](@ref).
For the other factor it computes the expansion of `D(u0)(x)` and
explicitly cancels the division by `x^(u0.p - u0.α)`.
"""
function F0(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M::Integer = 5, ϵ::Arb = Arb(1))
    Du0_expansion = D(u0, AsymptoticExpansion(); M)(ϵ)

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && Arblib.ref(x, 0) <= ϵ)

        res = eval_expansion(u0, Du0_expansion, x, offset = -u0.p, offset_i = -1)

        return res * inv_u0(x)
    end
end

"""
    inv_u0_normalised(u0::FractionalKdVAnsatz{Arb}; M = 5, ϵ = one(Arb))

Return a function for evaluation `x^-u0.α / u0(x)` for `x` close to
zero.

# Arguments
- `M::Integer` determines the number of terms in the asymptotic
  expansions.
- `ϵ::Arb` determines the interval ``[-ϵ, ϵ]`` on which the expansion
  is valid.

# Implementation
It computes an expansion of `u0` at `x = 0` and explicitly handles the
cancellation with `x^-u0.α`.
"""
function inv_u0_normalised(u0::FractionalKdVAnsatz{Arb}; M::Integer = 5, ϵ::Arb = one(Arb))
    bhkdv_skip_main = u0.use_bhkdv
    expansion = u0(ϵ, AsymptoticExpansion(); M, bhkdv_skip_main)

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && abs(x) <= ϵ) || (x isa ArbSeries && abs(Arblib.ref(x, 0)) <= ϵ)

        bhkdv_skipped_u0_main = bhkdv_skip_main
        return inv(eval_expansion(u0, expansion, x, offset_i = -1; bhkdv_skipped_u0_main))
    end
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
        u0_precomputed[(1, j, 0)] = OrderedDict(j => gamma(s) * cospi(s / 2))
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
        if s == -1
            # -gamma(s) * cospi(s / 2) has a removable singularity at
            # -s = -1, where it takes the value π / 2
            Hu0_precomputed[(2, j, 0)] = OrderedDict(j => π / 2)
        else
            Hu0_precomputed[(2, j, 0)] = OrderedDict(j => -gamma(s) * cospi(s / 2))
        end
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

    return a::AbstractVector -> begin
        @assert length(a) == u0.N0 + 1
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
    D2(u0::FractionalKdVAnsatz, evaltype::Symbolic; M::Integer = 5)

Return a function such that `D2(u0, evaltype, N)(a)` computes the
coefficients in the asymptotic expansion with indices `3` to `N0 + 2`
using the values from `a`.

This is an alternative implementation to [`D`](@ref).

# Implementation
The terms in the expansion have exponents of the form
```
-i * α + j * p0 + 2m
```
with `i` equal to `0`, `1` or `2`, `j` ranging from `0` to `2N0` and
`m` from `1` to infinity. We say that the terms with a given `i` are
of type `i`.

The exponents for the first two terms are given by `-2α` and `-2α +
p0`. We start by finding the exponents for term number `3` to `N0 +
2`.

In the computations we split the terms into singular and analytic
terms. The singular terms are the ones coming from the Clausen
function with argument `s` that has exponent `s - 1`. The analytic
terms are those with an exponent of the form `2m`. Since we square
`u0` we also need to take into account products of singular and
analytic ones.
"""
function D2(
    u0::FractionalKdVAnsatz{T},
    ::Symbolic;
    M::Integer = 5,
    threaded = true,
) where {T}
    # First step is to compute the exponents which we need

    # Given i, j, m get the corresponding exponent
    key_exponent = ((i, j, m),) -> -i * u0.α + j * u0.p0 + 2m

    # We never have to consider m above this value
    M_upper = ceil(Int, Float64(key_exponent((2, 2u0.N0, 0)) / 2))

    # Exponents for term of type 0, 1 and 2 represented as (i, j, m)
    exponents0 = [(0, 0, m) for m = 1:M_upper]
    exponents1 = [(1, j, m) for m = 1:M_upper, j = 0:u0.N0]
    exponents2 = [(2, j, 0) for j = 0:2u0.N0]

    # Keys of exponents in sorted order
    exponents = sort!([exponents0; exponents1[:]; exponents2], by = key_exponent)[3:u0.N0+2]

    # This is the maximum value of j and m we need
    J = maximum(key -> key[2], exponents)
    M = maximum(key -> key[3], exponents)

    # Next step is to precompute the coefficients in the expansions of
    # the Clausen functions

    u0_precomputed_singular = map(0:u0.N0) do j
        s = u0.α - j * u0.p0
        gamma(s) * cospi(s / 2)
    end
    u0_precomputed_analytic =
        [(-1)^m * zeta(1 - u0.α + j * u0.p0 - 2m) / factorial(2m) for m = 1:M, j = 0:u0.N0]

    Hu0_precomputed_singular = map(0:u0.N0) do j
        s = 2u0.α - j * u0.p0
        if s == -1
            # -gamma(s) * cospi(s / 2) has a removable singularity at
            # -s = -1, where it takes the value π / 2
            π / 2
        else
            -gamma(s) * cospi(s / 2)
        end
    end
    Hu0_precomputed_analytic = [
        -(-1)^m * zeta(1 - 2u0.α + j * u0.p0 - 2m) / factorial(2m) for m = 1:M, j = 0:u0.N0
    ]

    return a::AbstractVector -> begin
        @assert length(a) == length(u0_precomputed_singular)

        u0_res_singular = u0_precomputed_singular .* a.parent
        u0_res_analytic = u0_precomputed_analytic * a.parent
        Hu0_res_singular = Hu0_precomputed_singular .* a.parent
        Hu0_res_analytic = Hu0_precomputed_analytic * a.parent

        # Compute u0_res_singular^2 / 2
        u02_res_singular = zeros(eltype(u0_res_singular), J + 1)
        if threaded && J >= 256
            Threads.@threads for i = 1:J+1
                @inbounds for j = 1:i÷2
                    u02_res_singular[i] += u0_res_singular[j] * u0_res_singular[i-j+1]
                end
                if isodd(i)
                    @inbounds u02_res_singular[i] += u0_res_singular[i÷2+1]^2 / 2
                end
            end
        else
            @inbounds for i = 1:J+1
                for j = 1:i÷2
                    u02_res_singular[i] += u0_res_singular[j] * u0_res_singular[i-j+1]
                end
                if isodd(i)
                    u02_res_singular[i] += u0_res_singular[i÷2+1]^2 / 2
                end
            end
        end

        # Compute u0_res_analytic^2 / 2
        u02_res_analytic = zeros(eltype(u0_res_analytic), M)
        @inbounds for i = 1:M
            for j = 1:(i-1)÷2
                u02_res_analytic[i] += u0_res_analytic[j] * u0_res_analytic[i-j]
            end
            if iseven(i)
                u02_res_analytic[i] += u0_res_analytic[i÷2]^2 / 2
            end
        end

        # Compute u0_res_singular * u0_res_analytic / 2
        u02_res_singular_analytic = u0_res_singular * transpose(u0_res_analytic)

        # We don't need the above computed results to make changes
        # inplace
        res_singular = u02_res_singular
        @inbounds for k = 1:min(length(Hu0_res_singular), length(res_singular))
            res_singular[k] += Hu0_res_singular[k]
        end
        res_analytic = u02_res_analytic
        @inbounds for k = 1:min(length(Hu0_res_analytic), length(res_analytic))
            res_analytic[k] += Hu0_res_analytic[k]
        end
        res_singular_analytic = u02_res_singular_analytic

        res = Vector{eltype(res_singular)}(undef, u0.N0)
        @inbounds for k in eachindex(exponents, res)
            i, j, m = exponents[k]
            if i == 0
                res[k] = res_analytic[m]
            elseif i == 1
                res[k] = res_singular_analytic[j+1, m]
            else
                res[k] = res_singular[j+1]
            end
        end

        return res
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
