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

For `α` close to `-1` there's a lot of terms in the expansion and for
performance reasons we therefore make use of inplace calculations to
reduce the number of allocations.
"""
function eval_expansion(
    u0::FractionalKdVAnsatz{Arb},
    expansion::AbstractDict{NTuple{3,Int},Arb},
    x::Union{Arb,ArbSeries};
    offset_i::Integer = 0,
    offset::Arb = Arb(0),
    bhkdv_skipped_u0_main = false,
)
    res = zero(x)
    exponent = zero(u0.α)
    term = zero(x)

    if u0.use_bhkdv && bhkdv_skipped_u0_main
        # We don't have to be very careful with allocations here

        s = 1 - u0.α
        # The first argument to clausenc_expansion is not used in this
        # case, so we can take it to be zero.
        C1 = clausenc_expansion(Arb(0), s, 3)[1]
        C2 = clausenc_expansion(Arb(0), s + u0.p0, 3)[1]

        # exponent = -(1 + offset_i) * u0.α + offset
        Arblib.mul!(exponent, u0.α, -(1 + offset_i))
        Arblib.add!(exponent, exponent, offset)

        # term = u0.a[0] * (C1 - C2 * abspow(x, u0.p0)) * abspow(x, exponent)
        Arblib.mul!(term, u0.a[0], (C1 - C2 * abspow(x, u0.p0)))
        Arblib.mul!(term, term, abspow(x, exponent))

        # res += term
        Arblib.add!(res, res, term)
    end

    for ((i, j, m), y) in expansion
        if !iszero(y)
            #exponent = -(i + offset_i) * u0.α + j * u0.p0 + m + offset
            Arblib.add!(exponent, offset, m)
            Arblib.addmul!(exponent, u0.p0, j)
            Arblib.submul!(exponent, u0.α, i + offset_i)

            # res += y * abspow(x, exponent)
            if x isa Arb
                abspow!(term, x, exponent)
                Arblib.addmul!(res, y, term)
            elseif x isa ArbSeries
                # No abspow! for ArbSeries
                Arblib.mul!(term, abspow(x, exponent), y)
                Arblib.add!(res, res, term)
            end
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

It first uses a zero order approximation in `x` and for computing that
it uses a zero order approximation in `s`.
"""
clausencmzeta_diff(x, s, ϵ) = clausencmzeta(x, s) - clausencmzeta(x, s + ϵ)

function clausencmzeta_diff(x::Arb, s::Arb, ϵ::Arb)
    if iswide(x)
        # Use a zero order approximation in x
        deriv_x = -(clausens(x, s - 1) - clausens(x, s + ϵ - 1))
        mid_x = midpoint(Arb, x)
        res = add_error(clausencmzeta_diff(mid_x, s, ϵ), (x - mid_x) * deriv_x)
    elseif iswide(s)
        # Use a zero order approximation in s
        deriv_s = clausencmzeta(x, s, 1) - clausencmzeta(x, s + ϵ, 1)
        mid_s = midpoint(Arb, s)
        res = add_error(clausencmzeta_diff(x, mid_s, ϵ), (s - mid_s) * deriv_s)
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
            # See the argument u0_skipped_u0_main to eval_expansion as
            # well as _F0_bhkdv
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
            res -= u0.a[0] * clausencmzeta_diff(x, s, u0.p0)
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

function H(
    u0::FractionalKdVAnsatz{Arb},
    ::AsymptoticExpansion;
    M::Integer = 5,
    bhkdv_skip_main = false, # See _F0_bhkdv for use
    bhkdv_skip_singular_j_until = 0, # See _F0_bhkdv for use
)
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
            if !bhkdv_skip_main
                res[(2, 0, 0)] -= C1 * u0.a[0]
                res[(2, 1, 0)] -= -C2 * u0.a[0]
            end
            for m = ifelse(bhkdv_skip_main, 2, 1):M-1
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
            if contains_int &&
               isodd(n) &&
               !(u0.use_bhkdv && j <= bhkdv_skip_singular_j_until)
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
                if !(u0.use_bhkdv && j <= bhkdv_skip_singular_j_until)
                    res[(2, j, 0)] -= C * u0.a[j]
                    res[(0, 0, 2)] -= p[2] * u0.a[j]
                end
                for m = 2:M-1
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
inv(u0(x) / x^-u0.α) * (x^u0.p / u0.w(x)) * (D(u0)(x) / x^(u0.p - u0.α))
```
It computes `inv(u0(x) / x^-u0.α)` using [`inv_u0_normalised`](@ref)
and `x^u0.p / u0.w(x)` using `w.xpdivw`. For the third factor it
computes the expansion of `D(u0)(x)` and explicitly cancels the
division by `x^(u0.p - u0.α)`.
"""
function F0(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M::Integer = 5, ϵ::Arb = Arb(1))
    u0.use_bhkdv && return _F0_bhkdv(u0, Asymptotic(); M, ϵ)

    Du0_expansion = D(u0, AsymptoticExpansion(); M)(ϵ)

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && Arblib.ref(x, 0) <= ϵ)

        res = eval_expansion(u0, Du0_expansion, x, offset = -u0.p, offset_i = -1)

        return res * inv_u0(x) * u0.xpdivw(x)
    end
end

"""
    _F0_bhkdv(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M = 5, ϵ = one(Arb), bhkdv_skip_singular_j_until::Integer = 100,)

Implementation of [`F0`](@ref) when `u0.use_bhkdv` is true.

# Arguments
- `M::Integer` determines the number of terms in the asymptotic
  expansions.
- `ϵ::Arb` determines the interval ``[-ϵ, ϵ]`` on which the expansion
  is valid.
- `bhkdv_skip_singular_j_until::Integer = 100` is used to determine
  the number of terms in the expansion of `H(u0)` to threat
  separately. See the implementation below for details on how it is
  used.

# Implementation
Similarly to the default version it splits `F0(u0)` as
```
inv(u0(x) / x^-u0.α) * (x^u0.p / u0.w(x)) * (D(u0)(x) / x^(u0.p - u0.α))
```
It computes `inv(u0(x) / x^-u0.α)` using [`inv_u0_normalised`](@ref)
and `x^u0.p / u0.w(x)` using `w.xpdivw`. The difference is in how the
third factor is handled.

We have that
```
D(u0)(x) = u0(x)^2 / 2 + H(u0)(x)
```
To compute it we split `u0` into two parts and `H(u0)` into three
parts.

For this we let
```
c(s) = gamma(s) * cospi(s / 2)
```
Note that this has a removable singularity at `s = -1` and to compute
better enclosures it is beneficial to rewrite it as
```
c(s) = π * gamma(s + 2) * sinc((1 + s) / 2) / 2s
```

## Parts of `u0`
For `u0` the first part is the two singular terms in the expansion of
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
Given by
```
u0_part1 = a0 * (c(α) - c(α - p0) * x^p0) * x^-α
```
The second part is the remaining terms in the expansion, computed
using `bhkdv_skip_main = true`.

## Parts of `H(u0)`
For `H(u0)` the first part is the two singular terms and the `x^2`
term in the expansion of
```
-a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
```
Given by
```
Hu0_part1 = -a0 * (c(2α) - c(2α - p0) * x^p0 + K * x^(2 + 2α)) * x^(-2α)
```
with `K = -(zeta(-1 - 2α) - zeta(-1 - 2α + p0)) / 2`. The second part
is the singular term and the `x^2` term in the expansion of
```
-u0.a[j] * clausencmzeta(x, 1 - 2α + j * p0)
```
for `j` from `1` to `bhkdv_skip_singular_j_until`. Given by
```
Hu0_part2 = -sum(1:bhkdv_skip_singular_j_until) do j
    s = 1 - 2α + j * p0

    (
        gamma(1 - s) * sinpi(s / 2) * abspow(x, s - 1) -
        zeta(s - 2) / 2 * abspow(x, 2)
    ) * u0.a[j]
end
```
The third part is given by all the remaining terms in the expansion,
computed using `bhkdv_skip_main = true` and
`bhkdv_skip_singular_j_until`.

## Putting the parts together
With the above split of `u0` and `H(u0)` we can write `D(u0)` as
```
D(u0)(x) = u0_part1(x)^2 / 2 + u0_part1(x) * u0_part2(x) + u0_part2(x)^2 / 2 +
    Hu0_part1(x) + Hu0_part2(x) + Hu0_part3(x)
```
We group the terms into two parts
```
part1 = u0_part1(x)^2 / 2 + Hu0_part1(x)
part2 = u0_part1(x) * u0_part2(x) + u0_part2(x)^2 / 2 + Hu0_part2(x) + Hu0_part3(x)
```

## Computing `part1 / x^(p - α)`
We have that
```
part1 = u0_part1(x)^2 / 2 + Hu0_part1(x)
    = a0^2 * (c(α) - c(α - p0) * x^p0)^2 * x^(-2α) / 2
      - a0 * (c(2α) - c(2α - p0) * x^p0 + K * x^(2 + 2α)) * x^(-2α)
    = a0 * (
        a0 * c(α)^2 / 2 - c(2α)
        - (a0 * c(α) * c(α - p0) - c(2α - p0)) * x^p0
        + a0 * c(α - p0)^2 / 2 * x^2p0
        - K * x^(2 + 2α)
      ) * x^(-2α)
```
By construction `a0 = 2c(2α) / c(α)^2`, giving us `a0 * c(α)^2 / 2 -
c(2α) = 0` and allows us to simplify it as
```
part1 = 2c(2α) / c(α)^2 * (
    c(2α - p0) - c(2α) * c(α - p0) / c(α)
    + c(2α) * (c(α - p0) / c(α))^2 * x^p0
    + (zeta(-1 - 2α) - zeta(-1 - 2α + p0)) / 2 * x^(2 + 2α - p0)
) * x^(-2α + p0)
```
where we have inserted the value for `K`. To get better enclosures for
wide `x` it is beneficial to write it as
```
part1 = 2c(2α) / c(α)^2 * (
    c(2α - p0) - c(2α) * c(α - p0) / c(α)
    + (
        (zeta(-1 - 2α) - zeta(-1 - 2α + p0)) / 2 +
        c(2α) * (c(α - p0) / c(α))^2 * x^(-2 - 2α + 2p0)
    ) * x^(2 + 2α - p0)
) * x^(-2α + p0)
```
This does however not always work when `x` overlaps zero since in some
cases `-2 - 2α + 2p0` overlaps zero. For `x` overlapping zero we
therefore use the previous formulation.

Finally we compute a tight enclosure in `α` of this using a high
degree expansion in `α`.

## Computing `part2 / x^(p - α)`
We compute it by splitting it in the following way
```
(u0_part1(x) / x^-α) * (u0_part2(x) / x^p) +
(u0_part2(x) / x^p) * (u0_part2(x) / x^-α) / 2 +
Hu0_part2(x) / x^(p - α) +
Hu0_part3(x) / x^(p - α)
```
Both `u0_part2` and `Hu0_part3` give decent enclosures directly. For
`u0_part1` it is enough to use that `a0 = 2c(2α) / c(α)^` and use
[`ArbExtras.enclosure_series`](@ref). For `Hu0_part2` we need to do
slightly more work.

### Computing `Hu0_part2 / x^(p - α)`
Recall that it is given by
```
Hu0_part2 / x^(p - α) = -sum(1:bhkdv_skip_singular_j_until) do j
    s = 1 - 2α + j * p0

    (
        gamma(1 - s) * sinpi(s / 2) * abspow(x, s + α - 1 - p) -
        zeta(s - 2) / 2 * abspow(x, 2 + α - p)
    ) * u0.a[j]
end
```
If we let
```
f(α, j, x) =  gamma(2α - j * p0) * cospi((2α - j * p0) / 2) * x^(-α + j * p0 - p) -
    zeta(-1 - 2α + j * p0) / 2 * x^(2 + α - p)
```
We can write this as
```
Hu0_part2 / x^(p - α) = -sum(1:bhkdv_skip_singular_j_until) do j
    u0.a[j] * f(α, j, x)
end
```
To get a better enclosure in `α` of `f(α, j, x)` we first compute an
enclosure of the derivative in `α`. If it is non-zero we evaluate at
the endpoints of `α`, otherwise we just a zero order approximation.
This is however not enough to get good enclosures, we also need to
better handle the removable singularity for `f`.

As a first step we use that
```
zeta(-1 - 2α + j * p0) = zeta_deflated(-1 - 2α + j * p0) - inv(2 + 2α - j * p0)
```
and
```
gamma(2α - j * p0) = gamma(3 + 2α - j * p0) / rising(2α - j * p0, 3)
    = inv(2 + 2α - j * p0) * gamma(3 + 2α - j * p0) / rising(2α - j * p0, 2)
```
This allows us to write `f` as
```
f(α, j, x) =  inv(2 + 2α - j * p0) * (
        gamma(3 + 2α - j * p0) / rising(2α - j * p0, 2) * cospi((2α - j * p0) / 2) * x^(-α + j * p0 - p)
        + 1 / 2 * x^(2 + α - p)
    ) - zeta_deflated(-1 - 2α + j * p0) / 2 * x^(2 + α - p)
```
How to best evaluate this depends on which of the exponents `-α + j *
p0 - p` and `2 + α - p` are largest. In practice we have that for `j =
1` the last one is largest and for `j >= 2` the first one is largest.
For wide values of `α` it sometimes happen that they overlap.
For non-zero `x` the precise way of evaluation is not important, it is
only for `x` overlapping zero for which we need to take care so that
all factors are finite. We handle the cases `j = 1` and `j >= 2`
separately, we then also have a third case for when `x` overlaps zero
and the exponents overlap.

#### `j >= 2`
Factoring out `x^(2 + α - p)` we get
```
f(α, j, x) =  (
    inv(2 + 2α - j * p0) * (
        gamma(3 + 2α - j * p0) / rising(2α - j * p0, 2) * cospi((2α - j * p0) / 2) * x^(-2 - 2α + j * p0)
        + 1 / 2
    ) - zeta_deflated(-1 - 2α + j * p0) / 2
) * x^(2 + α - p)
```

The term
```
zeta_deflated(-1 - 2α + j * p0) / 2
```
we can in general evaluate directly. However, if the argument is close
to, but does not contain, the removable singularity at `1` the
implementation gives very bad enclosures. In that case it is better to
slightly widen the argument to also contain the removable singularity,
in which case a different, better, algorithm is used.

The term
```
inv(2 + 2α - j * p0) * (
    gamma(3 + 2α - j * p0) / rising(2α - j * p0, 2) * cospi((2α - j * p0) / 2) * x^(-2 - 2α + j * p0)
    + 1 / 2
)
```
can also be evaluated directly as long as `2 + 2α - j * p0` is not too
close to zero, where there is a removable singularity. To better
handle this case we let `t1 = 2 + 2α - j * p0` and write it as
```
inv(t1) * (
    gamma(t1 + 1) / rising(t1 - 2, 2) * cospi((t1 - 2) / 2) * x^(-t1)
    + 1 / 2
)
```
Inserting `t1 = 0` it is straight forward to see that we indeed have a
removable singularity. When `t1` is not too close to zero we evaluate
this directly. Otherwise we widen `t1` to include zero and use
[`fx_div_x`](@ref).

#### `j = 1`
In this case we instead write `f` as
```
f(α, j, x) =  inv(2 + 2α - j * p0) * (
        gamma(3 + 2α - j * p0) / rising(2α - j * p0, 2) * cospi((2α - j * p0) / 2)
        + 1 / 2 * x^(2 + 2α - j * p0)
    ) * x^(-α + j * p0 - p)
    - zeta_deflated(-1 - 2α + j * p0) / 2 * x^(2 + α - p)
```
and then use the same approach as for `j >= 2` to compute a good
enclosure.

#### Overlapping exponents and `x` containing zero
In this case we have that `-α + j * p0 - p` overlaps with `2 + α - p`,
which means that `2 + 2α - j * p0` overlaps with zero. We keep the
original formulation
```
f(α, j, x) =  inv(2 + 2α - j * p0) * (
        gamma(3 + 2α - j * p0) / rising(2α - j * p0, 2) * cospi((2α - j * p0) / 2) * x^(-α + j * p0 - p)
        + 1 / 2 * x^(2 + α - p)
    ) - zeta_deflated(-1 - 2α + j * p0) / 2 * x^(2 + α - p)
```
and handle the removable singularity.
"""
function _F0_bhkdv(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
    bhkdv_skip_singular_j_until::Integer = 50,
)
    @assert u0.use_bhkdv

    u0_expansion = u0(ϵ, AsymptoticExpansion(), bhkdv_skip_main = true; M)
    Hu0_expansion =
        H(u0, AsymptoticExpansion(), bhkdv_skip_main = true; M, bhkdv_skip_singular_j_until)(
            ϵ,
        )

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    # c(s) = gamma(s) * cospi(s / 2)
    c(s) =
        if (s isa ArbSeries && is_approx_integer(s[0]) && round(Float64(s[0])) == -1)
            # _sinc performs poorly for wide arguments close to zero, it
            # is better to slightly widen the argument to include zero
            t = (1 + s) / 2
            t[0] = union(t[0], Arb(0))
            π * gamma(s + 2) * _sinc(t) / 2s
        else
            π * gamma(s + 2) * _sinc((1 + s) / 2) / 2s
        end

    f(α, j, x) = begin
        # Main argument for f1
        t1 = 2 + 2α - j * u0.p0
        if t1 isa ArbSeries && is_approx_integer(t1[0])
            # Widen argument to contain 0 so that it uses the algorithm
            # that explicitly handles the removable singularity.
            t1[0] = union(t1[0], Arb(0))
        end

        # Main argument for f2
        t2 = -1 - 2α + j * u0.p0
        if t2 isa ArbSeries && is_approx_integer(t2[0])
            # Widen argument to contain 1 so that it uses the algorithm
            # that explicitly handles the removable singularity.
            t2[0] = union(t2[0], Arb(1))
        end

        f2 = zeta_deflated(t2, Arb(1)) / 2

        t1_contains_zero = Arblib.contains_zero(t1 isa ArbSeries ? t1[0] : t1)

        if t1_contains_zero && Arblib.contains_zero(x)
            # The exponents overlap so don't factor out any power of x
            f1 = fx_div_x(t1, enclosure_degree = -1, force = true) do t1
                gamma(t1 + 1) * cospi((t1 - 2) / 2) / rising(t1 - 2, 2) *
                abspow(x, (j * u0.p0 - t1) / 2 + 1 - u0.p) +
                1 // 2 * abspow(x, (j * u0.p0 + t1) / 2 + 1 - u0.p)
            end

            return f1 - f2 * abspow(x, 2 + α - u0.p)
        elseif j >= 2
            if t1_contains_zero
                # Handle the removable singularity
                f1 = fx_div_x(t1, enclosure_degree = -1) do t1
                    gamma(t1 + 1) * cospi((t1 - 2) / 2) / rising(t1 - 2, 2) * abspow(x, -t1) + 1 // 2
                end
            else
                f1 =
                    (
                        gamma(t1 + 1) * cospi((t1 - 2) / 2) / rising(t1 - 2, 2) *
                        abspow(x, -t1) + 1 // 2
                    ) / t1
            end

            return (f1 - f2) * abspow(x, 2 + α - u0.p)
        else
            if t1_contains_zero
                # Handle the removable singularity
                f1 = fx_div_x(t1, enclosure_degree = -1) do t1
                    gamma(t1 + 1) * cospi((t1 - 2) / 2) / rising(t1 - 2, 2) +
                    1 // 2 * abspow(x, t1)
                end
            else
                f1 =
                    (
                        gamma(t1 + 1) * cospi((t1 - 2) / 2) / rising(t1 - 2, 2) +
                        1 // 2 * abspow(x, t1)
                    ) / t1
            end

            return f1 * abspow(x, -α + j * u0.p0 - u0.p) - f2 * abspow(x, 2 + α - u0.p)
        end
    end

    # Derivative of f w.r.t. α
    df(α::Arb, j, x) = f(ArbSeries((α, 1)), j, x)[1]
    df(α::ArbSeries, j, x) =
        Arblib.derivative(f(ArbSeries(α, degree = Arblib.degree(α) + 1), j, x))

    return x::Arb -> begin
        @assert x <= ϵ

        # abspow(x, y::ArbSeries) only supports y of degree at most 2
        # when x overlaps with zero. We therefore lower the degree
        # used in this case.

        # part1 / x^(u0.p - u0.α)
        part1_divpα = ArbExtras.enclosure_series(
            u0.α,
            degree = ifelse(Arblib.contains_zero(x), 1, 10),
        ) do α
            if Arblib.contains_zero(x)
                2c(2α) / c(α)^2 *
                abspow(x, -α + u0.p0 - u0.p) *
                (
                    c(2α - u0.p0) - 2c(2α) * c(α - u0.p0) / c(α) + (
                        (zeta(-1 - 2α) - zeta(-1 - 2α + u0.p0)) / 2 *
                        abspow(x, 2 + 2α - u0.p0) +
                        c(2α) * (c(α - u0.p0) / c(α))^2 * abspow(x, u0.p0)
                    )
                )
            else
                2c(2α) / c(α)^2 *
                abspow(x, -α + u0.p0 - u0.p) *
                (
                    c(2α - u0.p0) - 2c(2α) * c(α - u0.p0) / c(α) +
                    (
                        (zeta(-1 - 2α) - zeta(-1 - 2α + u0.p0)) / 2 +
                        c(2α) * (c(α - u0.p0) / c(α))^2 * abspow(x, -2 - 2α + 2u0.p0)
                    ) * abspow(x, 2 + 2α - u0.p0)
                )
            end
        end

        # u0_part1 / x^-α
        u0_part1_divα = ArbExtras.enclosure_series(u0.α, degree = 1) do α
            2c(2α) / c(α)^2 * (c(α) - c(α - u0.p0) * abspow(x, u0.p0))
        end

        # u0_part2 / x^p
        u0_part2_divp = eval_expansion(u0, u0_expansion, x, offset = -u0.p)
        # u0_part2 / x^-α
        u0_part2_divα = eval_expansion(u0, u0_expansion, x, offset_i = -1)

        # Hu0_part2 / x^(p - α)
        Hu0_part2_divpα =
            -sum(1:min(bhkdv_skip_singular_j_until, u0.N0), init = zero(x)) do j
                # Compute derivative in α
                if Arblib.contains_zero(x)
                    deriv_α = df(u0.α, j, x)
                else
                    deriv_α = ArbExtras.enclosure_series(α -> df(α, j, x), u0.α, degree = 4)
                end

                if Arblib.contains_zero(deriv_α)
                    # Use a zero order approximation
                    mid_α = midpoint(Arb, u0.α)
                    term = u0.a[j] * add_error(f(mid_α, j, x), (u0.α - mid_α) * deriv_α)
                else
                    # Evaluate at the endpoints
                    term =
                        u0.a[j] * union(
                            f(ArbExtras.enclosure_lbound(u0.α), j, x),
                            f(ArbExtras.enclosure_ubound(u0.α), j, x),
                        )
                end

                term
            end

        # Hu0_part3 / x^(p - α)
        Hu0_part3_divpα =
            eval_expansion(u0, Hu0_expansion, x, offset = -u0.p, offset_i = -1)

        # part2 / x^(u0.p - u0.α)
        part2_divpα =
            u0_part1_divα * u0_part2_divp +
            u0_part2_divp * u0_part2_divα / 2 +
            Hu0_part2_divpα +
            Hu0_part3_divpα

        res = part1_divpα + part2_divpα

        return res * inv_u0(x) * u0.xpdivw(x)
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

    # In some cases the parameters are very close to singularities of
    # the functions and it is therefore beneficial to do the
    # calculations in Arb and then convert back to T
    α = convert(Arb, u0.α)
    p0 = convert(Arb, u0.p0)

    u0_precomputed_singular = map(0:u0.N0) do j
        s = α - j * p0
        convert(T, gamma(s) * cospi(s / 2))
    end
    u0_precomputed_analytic =
        T[(-1)^m * zeta(1 - α + j * p0 - 2m) / factorial(2m) for m = 1:M, j = 0:u0.N0]

    Hu0_precomputed_singular = map(0:u0.N0) do j
        s = 2α - j * p0
        if s == -1
            # -gamma(s) * cospi(s / 2) has a removable singularity at
            # -s = -1, where it takes the value π / 2
            convert(T, π) / 2
        else
            convert(T, -gamma(s) * cospi(s / 2))
        end
    end
    Hu0_precomputed_analytic =
        T[-(-1)^m * zeta(1 - 2α + j * p0 - 2m) / factorial(2m) for m = 1:M, j = 0:u0.N0]

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
