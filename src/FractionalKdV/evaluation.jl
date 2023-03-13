# This file contains code for evaluation of the approximate solution
# in different ways

"""
    eval_expansion(u0::FractionalKdVAnsatz, expansion, x; offset_i = 0, offset = 0)

Evaluate the given expansion. The term `((i, j,  m), y)` is evaluated
to
```
y * abs(x)^(-i * u0.α + j * u0.p0 + m)
```
and then they are all summed.

The arguments `offset_i` and `offset` can be set to adjust the
exponent, in that case the exponent will be given by
```
-(i + offset_i) * u0.α + j * u0.p0 + m + offset
```

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
)
    res = zero(x)
    exponent = zero(u0.α)
    term = zero(x)

    for ((i, j, m), y) in expansion
        if !iszero(y)
            # exponent = -(i + offset_i) * u0.α + j * u0.p0 + m + offset
            Arblib.add!(exponent, offset, m)
            Arblib.addmul!(exponent, u0.p0, j)
            Arblib.submul!(exponent, u0.α, i + offset_i)

            # res += y * abspow(x, exponent)
            abspow!(term, x, exponent)
            if x isa Arb
                Arblib.addmul!(res, y, term)
            elseif x isa ArbSeries
                Arblib.mul!(term, term, y)
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

    return ArbExtras.compose_zero!(res, res, x)
end


function (u0::FractionalKdVAnsatz)(x, ::Ball)
    s = 1 - u0.α
    if u0.use_bhkdv
        res = u0.a[0] * clausencmzeta_diff(x, s, u0.p0)
    else
        res = u0.a[0] * clausencmzeta(x, s)
    end

    onemα = copy(s)
    for j = 1:u0.N0
        if s isa Arb
            Arblib.fma!(s, u0.p0, j, onemα)
        else
            s = onemα + j * u0.p0
        end
        term = clausencmzeta(x, s)
        if res isa Arb
            Arblib.addmul!(res, u0.a[j], term)
        elseif res isa ArbSeries
            Arblib.mul!(term, term, u0.a[j])
            Arblib.add!(res, res, term)
        else
            res += u0.a[j] * term
        end
    end

    for n = 1:u0.N1
        res += u0.b[n] * (cos(n * x) - 1)
    end

    return res
end

(u0::FractionalKdVAnsatz)(x, ::Asymptotic; M::Integer = 5) =
    eval_expansion(u0, u0(x, AsymptoticExpansion(); M), x)

"""
    (u0::FractionalKdVAnsatz{Arb})(x::Arb, ::AsymptoticExpansion; M, skip_main)

Compute an expansion of `u0` around zero with remainder terms valid
on the interval ``[-abs(x), abs(x)]``.

# Arguments
- `M::Integer = 5`: Determines the order of the expansion. The
  remainder term is of order `x^2M`.
- `skip_main = false`: If true then don't add the `x^-α` term from the
  expansion of the leading Clausen to the expansion. If `u0.use_bhkdv`
  is true then skip the `x^-α` and `x^(-α + p0)` terms from the
  expansion of the two Clausen functions that are multiplied with
  `u0.a[0]`. These would then have to be added separately to the
  result.
"""
function (u0::FractionalKdVAnsatz{Arb})(
    x,
    ::AsymptoticExpansion;
    M::Integer = 5,
    skip_main = false,
)
    res = OrderedDict{NTuple{3,Int},Arb}()

    # Initiate the coefficients
    for m = 1:M
        res[(0, 0, 2m)] = 0
    end
    for j = 0:u0.N0
        res[(1, j, 0)] = 0
    end

    if u0.use_bhkdv
        s = 1 - u0.α
        C1, _, p1, E1 = clausenc_expansion(x, s, M, skip_constant = true)
        C2, _, p2, E2 = clausenc_expansion(x, s + u0.p0, M, skip_constant = true)
        if !skip_main
            res[(1, 0, 0)] += C1 * u0.a[0]
            res[(1, 1, 0)] += -C2 * u0.a[0]
        end
        for m = 1:M-1
            res[(0, 0, 2m)] += (p1[2m] - p2[2m]) * u0.a[0]
        end
        Arblib.add_error!(res[(0, 0, 2M)], (E1 - E2) * u0.a[0])
    else
        s = 1 - u0.α
        C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)
        if !skip_main
            res[(1, 0, 0)] += C * u0.a[0]
        end
        for m = 1:M-1
            res[(0, 0, 2m)] += p[2m] * u0.a[0]
        end
        Arblib.add_error!(res[(0, 0, 2M)], E * u0.a[0])
    end

    # Clausen terms
    for j = 1:u0.N0
        s = 1 - u0.α + j * u0.p0
        C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

        # Below we handle the special case when s contains an odd
        # integer. When s is close to an odd integer we have very
        # large cancellations between two of the terms in the
        # expansion and it turns out to be beneficial to use the same
        # method as when s contains an odd integer to compute an
        # enclosure.
        tol = s < 4 ? 0.01 : 0.001 # For large s it is worse to be close to an integer
        if is_approx_integer(s; tol) && isodd(round(Float64(s)))
            s = union(s, Arb(round(Float64(s))))
        end

        # Check for the special case when s contains an odd integer.
        # This requires extra work to get a finite result.
        contains_int, n = unique_integer(s)
        if contains_int && isodd(n)
            # The term corresponding to C and p[2((n - 1) ÷ 2)]
            # coincides and diverge so are handled separately. The
            # rest we treat normally.
            for m = 1:M-1
                m == (n - 1) ÷ 2 && continue # Skip this term
                res[(0, 0, 2m)] += p[2m] * u0.a[j]
            end

            m = n ÷ 2
            # We split the term into three, with coefficients K1,
            # K2 and K3 respectively.
            K1, K2 = clausenc_expansion_odd_s_singular_K1_K2(s, m)
            K3 = clausenc_expansion_odd_s_singular_K3(m)

            # The terms corresponding to K1 and K2 are multiplied
            # by a pure power and we can add them directly to the
            # expansion.
            res[(1, j, 0)] += K1 * u0.a[j]
            res[(0, 0, 2m)] += K2 * u0.a[j]

            # The term corresponding to K3 is multiplied by
            # x^2m * (abs(x)^(s - (2m + 1)) - 1) / (s - (2m + 1)).
            # We factor out x^(-i * u0.α + 1) from this where i is
            # as large as possible but so that we still have -i *
            # u0.α + 1 < s - 1, and enclose the rest.
            i = findfirst(i -> !(-i * u0.α + 1 < s - 1), 1:100) - 1
            # Enclosure of
            # x^(2m + i * u0.α - 1) * (abs(x)^(s - (2m + 1)) - 1) / (s - (2m + 1)).
            # on the interval [0, x].
            D = x_pow_s_x_pow_t_m1_div_t(Arb((0, x)), 2m + i * u0.α - 1, s - (2m + 1))
            res[(i, 0, 1)] = get(res, (i, 0, 1), zero(x)) - K3 * D * u0.a[j]
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
        s = 1 - 2u0.α
        if u0.use_bhkdv
            res = -u0.a[0] * clausencmzeta_diff(x, s, u0.p0)
        else
            res = -u0.a[0] * clausencmzeta(x, s)
        end

        onem2α = copy(s)
        for j = 1:u0.N0
            if s isa Arb
                Arblib.fma!(s, u0.p0, j, onem2α)
            else
                s = onem2α + j * u0.p0
            end
            term = clausencmzeta(x, s)
            if res isa Arb
                Arblib.submul!(res, u0.a[j], term)
            elseif res isa ArbSeries
                Arblib.mul!(term, term, u0.a[j])
                Arblib.sub!(res, res, term)
            else
                res -= u0.a[j] * term
            end
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

"""
    H(u0::FractionalKdVAnsatz{Arb}, ::AsymptoticExpansion; M, skip_main, skip_singular_j_until)

Return a function such that `H(u0)(x)` computes an expansion of
`H(u0)` around zero with remainder terms valid on the interval
``[-abs(x), abs(x)]``.

# Arguments
- `M::Integer = 5`: Determines the order of the expansion. The
  remainder term is of order `x^2M`.
- `skip_main = false`: If true then don't add the `x^-2α` and `x^2`
  term from the expansion of the leading Clausen to the expansion. If
  `u0.use_bhkdv` is true then skip the `x^-2α`, `x^(-2α + p0)` and
  `x^2` terms from the expansion of the two Clausen functions that are
  multiplied with `u0.a[0]`. These would then have to be added
  separately to the result.
- `skip_singular_j_until::Integer = 0`: If positive then the terms
  `x^(-2α + j * p0)` and `x^2` are treated differently in the
  expansion for the Clausen terms with `j` from `1` to
  ``skip_singular_j_until`. The terms are split in correspondence with
  [`clausenc_expansion_odd_s_singular_K1_K2`](@ref) and
  [`clausenc_expansion_odd_s_singular_K3`](@ref). The `K1` and `K2`
  terms are added to the expansion, but the `K3` term is not. You
  would then have to add
  ```
  -K3 * sum(1:skip_singular_j_until) do j
    u0.a[j] * x_pow_s_x_pow_t_m1_div_t(x, 2, s - (2m + 1))
  end
  ```
  to the result to account for the `K3` term. This option is used to
  avoid large cancellations when the arguments for the Clausen
  functions are very close to the removable singularity at `s = 3`,
  which happens for `α` close to `-1`.
"""
function H(
    u0::FractionalKdVAnsatz{Arb},
    ::AsymptoticExpansion;
    M::Integer = 5,
    skip_main = false,
    skip_singular_j_until = 0,
)
    return x -> begin
        res = OrderedDict{NTuple{3,Int},Arb}()

        # Initiate the coefficients
        for m = 1:M
            res[(0, 0, 2m)] = 0
        end
        for j = 0:u0.N0
            res[(2, j, 0)] = 0
        end

        if u0.use_bhkdv
            s = 1 - 2u0.α
            C1, _, p1, E1 = clausenc_expansion(x, s, M, skip_constant = true)
            C2, _, p2, E2 = clausenc_expansion(x, s + u0.p0, M, skip_constant = true)
            if !skip_main
                res[(2, 0, 0)] -= C1 * u0.a[0]
                res[(2, 1, 0)] -= -C2 * u0.a[0]
            end
            for m = ifelse(skip_main, 2, 1):M-1
                res[(0, 0, 2m)] -= (p1[2m] - p2[2m]) * u0.a[0]
            end
            Arblib.add_error!(res[(0, 0, 2M)], (E1 - E2) * u0.a[0])
        else
            s = 1 - 2u0.α
            C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)
            if !skip_main
                res[(2, 0, 0)] -= C * u0.a[0]
            end
            for m = ifelse(skip_main, 2, 1):M-1
                res[(0, 0, 2m)] -= p[2m] * u0.a[0]
            end
            Arblib.add_error!(res[(0, 0, 2M)], E * u0.a[0])
        end

        # Clausen terms
        for j = 1:u0.N0
            s = 1 - 2u0.α + j * u0.p0
            C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

            # Below we handle the special case when s contains an odd
            # integer. When s is close to an odd integer we have very
            # large cancellations between two of the terms in the
            # expansion and it turns out to be beneficial to use the same
            # method as when s contains an odd integer to compute an
            # enclosure.
            tol = s < 4 ? 0.01 : 0.001 # For large s it is worse to be close to integer
            if is_approx_integer(s) && isodd(round(Float64(s)))
                s = union(s, Arb(round(Float64(s))))
            end

            # Check for the special case when s contains an odd integer.
            # This requires extra work to get a finite result.
            contains_int, n = unique_integer(s)
            if contains_int && isodd(n) && j > skip_singular_j_until
                # The term corresponding to C and p[2((n - 1) ÷ 2)]
                # coincides and diverge so are handled separately. The
                # rest we treat normally.
                for m = 1:M-1
                    m == (n - 1) ÷ 2 && continue # Skip this term
                    res[(0, 0, 2m)] -= p[2m] * u0.a[j]
                end

                m = n ÷ 2
                # We split the term into three, with coefficients K1,
                # K2 and K3 respectively.
                K1, K2 = clausenc_expansion_odd_s_singular_K1_K2(s, m)
                K3 = clausenc_expansion_odd_s_singular_K3(m)

                # The terms corresponding to K1 and K2 are multiplied
                # by a pure power and we can add them directly to the
                # expansion.
                res[(2, j, 0)] -= K1 * u0.a[j]
                res[(0, 0, 2m)] -= K2 * u0.a[j]

                # The term corresponding to K3 is multiplied by
                # x^2m * (abs(x)^(s - (2m + 1)) - 1) / (s - (2m + 1)).
                # We factor out x^(-i * u0.α + 1) from this where i is
                # as large as possible but so that we still have -i *
                # u0.α + 1 < s - 1, and enclose the rest.
                i = findfirst(i -> !(-i * u0.α + 1 < s - 1), 1:100) - 1
                # Enclosure of
                # x^(2m + i * u0.α - 1) * (abs(x)^(s - (2m + 1)) - 1) / (s - (2m + 1)).
                # on the interval [0, x].
                D = x_pow_s_x_pow_t_m1_div_t(Arb((0, x)), 2m + i * u0.α - 1, s - (2m + 1))
                res[(i, 0, 1)] = get(res, (i, 0, 1), zero(x)) - K3 * D * u0.a[j]
            else
                if j <= skip_singular_j_until
                    K1, K2 = clausenc_expansion_odd_s_singular_K1_K2(s, 1)
                    res[(2, j, 0)] -= K1 * u0.a[j]
                    res[(0, 0, 2)] -= K2 * u0.a[j]
                else
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

function defect(u0::FractionalKdVAnsatz, ::Asymptotic; M::Integer = 5)
    f = defect(u0, AsymptoticExpansion(); M)
    return x -> eval_expansion(u0, f(x), x)
end

function defect(
    u0::FractionalKdVAnsatz{T},
    evaltype::AsymptoticExpansion;
    M::Integer = 5,
    skip_singular_j_until = 0,
) where {T}
    f = H(u0, evaltype; M, skip_singular_j_until)
    return x -> begin
        expansion1 = u0(x, evaltype; M)
        expansion2 = f(x)

        expansion = empty(expansion1)

        # u0^2 / 2 term
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

        # The term (2, 0, 0) is identically equal to zero since a[0] =
        # finda0(α)
        @assert Arblib.contains_zero(expansion[(2, 0, 0)])
        expansion[(2, 0, 0)] = 0

        return expansion
    end
end

"""
    F0(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M = 5, ϵ = one(Arb))

Return a function for evaluating `F0(u0)(x)` accurately for small
values of `x`.

# Arguments
- `M::Integer = 5`: Determines the number of terms in the asymptotic
  expansions.
- `ϵ::Arb = 1`: Determines the interval ``[-ϵ, ϵ]`` on which the expansion
  is valid.
- `skip_singular_j_until = nothing`: Argument given to `defect(u0)`
  when computing the expansion to skip some of the terms, they are
  then added back separately. If set to `nothing` a default value is
  used. The default value is to non-zero if `1 - 2α + p0 > 2.8` and is
  then set to the largest `j` such that `1 - 2α + j * p0 < 3.2`. The
  motivation for this is that we mainly want to avoid the singularity
  at `s = 3` for the expansions of the Clausen functions.

# Implementation
It splits `F0(u0)` as
```
inv(u0(x) / x^-u0.α) * (x^u0.p / u0.w(x)) * (defect(u0)(x) / x^(u0.p - u0.α))
```
It computes `inv(u0(x) / x^-u0.α)` using [`inv_u0_normalised`](@ref)
and `x^u0.p / u0.w(x)` using `w.xpdivw`. For the third factor it
computes the expansion of `defect(u0)(x)` and explicitly cancels the
division by `x^(u0.p - u0.α)`.
"""
function F0(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
    skip_singular_j_until = nothing,
)
    # This case as a separate implementation
    u0.use_bhkdv && return _F0_bhkdv(u0, Asymptotic(); M, ϵ, skip_singular_j_until)

    if isnothing(skip_singular_j_until)
        if 1 - 2u0.α + u0.p0 > 2.8
            skip_singular_j_until = findlast(j -> 1 - 2u0.α + j * u0.p0 < 3.2, 1:u0.N0)::Int
        else
            skip_singular_j_until = 0
        end
    end

    defect_u0_expansion = defect(u0, AsymptoticExpansion(); M, skip_singular_j_until)(ϵ)

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && Arblib.ref(x, 0) <= ϵ)

        res = eval_expansion(u0, defect_u0_expansion, x, offset = -u0.p, offset_i = -1)

        if skip_singular_j_until > 0
            # Add skipped terms
            res += let v1 = 2 + u0.α - u0.p, v2 = -2 - 2u0.α
                -clausenc_expansion_odd_s_singular_K3(1) *
                sum(1:min(skip_singular_j_until, u0.N0), init = zero(x)) do j
                    u0.a[j] * x_pow_s_x_pow_t_m1_div_t(x, v1, v2 + j * u0.p0)
                end
            end
        end

        return res * inv_u0(x) * u0.xpdivw(x)
    end
end

"""
    _F0_bhkdv(u0::FractionalKdVAnsatz{Arb}, ::Asymptotic; M = 5, ϵ, skip_singular_j_until)

Implementation of [`F0`](@ref) when `u0.use_bhkdv` is true.

# Arguments
- `M::Integer`: Determines the number of terms in the asymptotic
  expansions.
- `ϵ::Arb = 1`: Determines the interval ``[-ϵ, ϵ]`` on which the
  expansion is valid.
- `skip_singular_j_until = nothing`: Argument given to `H(u0)` when
  computing the expansion to skip some of the terms, they are then
  added back separately. If set to `nothing` a default value is used.
  The default value is set to the largest `j` such that `1 - 2α + j *
  p0 < 3.2`, but with a maximum value of `100`. The motivation for
  this is that we mainly want to avoid the singularity at `s = 3` for
  the expansions of the Clausen functions.

# Implementation
Similarly to the default version it splits `F0(u0)` as
```
inv(u0(x) / x^-u0.α) * (x^u0.p / u0.w(x)) * (defect(u0)(x) / x^(u0.p - u0.α))
```
It computes `inv(u0(x) / x^-u0.α)` using [`inv_u0_normalised`](@ref)
and `x^u0.p / u0.w(x)` using `w.xpdivw`. The difference is in how the
third factor is handled.

We have that
```
defect(u0)(x) = u0(x)^2 / 2 + H(u0)(x)
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
using `skip_main = true`.

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
is the remaining terms in the expansion, save for a few specific
terms. The expansion is computed using `skip_main = true` and
`skip_singular_j_until`. The `skip_singular_j_until` argument means we
have to also add the third part
```
-K3 * sum(1:skip_singular_j_until) do j
  u0.a[j] * x_pow_s_x_pow_t_m1_div_t(x, 2, s - (2m + 1))
end
```

## Putting the parts together
With the above split of `u0` and `H(u0)` we can write `defect(u0)` as
```
defect(u0)(x) = u0_part1(x)^2 / 2 + u0_part1(x) * u0_part2(x) + u0_part2(x)^2 / 2 +
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
    c(2α - p0) - 2c(2α) * c(α - p0) / c(α)
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
They all give decent enclosures in `α` directly. For `Hu0_part3` we
compute a better enclosure in `x` using
[`ArbExtras.enclosure_series`](@ref).
"""
function _F0_bhkdv(
    u0::FractionalKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 5,
    ϵ::Arb = Arb(1),
    skip_singular_j_until = nothing,
)
    @assert u0.use_bhkdv

    # Store in variable set to ::Int. Otherwise it fails with
    # determining the type properly.
    skip_singular_j_until_int::Int = if isnothing(skip_singular_j_until)
        min(findlast(j -> 1 - 2u0.α + j * u0.p0 < 3.2, 1:u0.N0)::Int, 100)
    else
        skip_singular_j_until
    end

    u0_expansion = u0(ϵ, AsymptoticExpansion(), skip_main = true; M)
    Hu0_expansion = H(
        u0,
        AsymptoticExpansion(),
        skip_main = true,
        skip_singular_j_until = skip_singular_j_until_int;
        M,
    )(
        ϵ,
    )

    inv_u0 = inv_u0_normalised(u0; M, ϵ)

    # c(s) = gamma(s) * cospi(s / 2)
    c(s) = begin
        sinc_sp2 = _sinc((1 + s) / 2)

        if s isa ArbSeries && round(Float64(s[0])) == -1
            # _sinc performs poorly for arguments close to zero, it is
            # often better to slightly widen the argument to include
            # zero
            t = (1 + s) / 2
            t[0] = union(t[0], Arb(0))
            sinc_sp2_wide = _sinc(t)

            for i = 0:Arblib.degree(s)
                sinc_sp2[i] = intersect(sinc_sp2[i], sinc_sp2_wide[i])
            end
        end

        return π * gamma(s + 2) * sinc_sp2 / 2s
    end

    # We make use of the values c(α), c(2α), c(α - u0.p0), c(2α -
    # u0.p0) and zeta(-1 - 2α) - zeta(-1 - 2α + u0.p0). Since they
    # don't depend on x it is possible to precompute them once and for
    # all. However since they are used inside
    # ArbExtras.enclosure_series it is a bit complicated to determine
    # exactly for which values of α we need to compute it. Instead we
    # compute it on the go and store in a dictionary.
    # IMPROVE: Dict is not fully thread-safe and in principle this
    # could lead to issues. However, in our case the value is fully
    # determined from the key and it doesn't matter if we lose any
    # updates since it will just compute the value directly instead.
    # But it could be worth looking more into exactly what could
    # happen.

    # Dicts for storing values
    c_precomputed_Arb = Dict{Arb,NTuple{5,Arb}}()
    c_precomputed_ArbSeries = Dict{ArbSeries,NTuple{5,ArbSeries}}()

    # Function to compute values for a given α
    c_precompute(α) =
        () -> (
            c(α),
            c(2α),
            c(α - u0.p0),
            c(2α - u0.p0),
            zeta(-1 - 2α) - zeta(-1 - 2α + u0.p0),
        )

    return (x::Arb) -> begin
        @assert x <= ϵ

        # part1 / x^(u0.p - u0.α)
        degree = if Arblib.contains_zero(x)
            # abspow(x, y::ArbSeries) only supports y of degree at most 2
            # when x overlaps with zero. We therefore lower the degree
            # used in this case.
            1
        elseif Float64(u0.α) < -0.99985
            # Degree 3 gives slightly better enclosures but many more
            # allocations. For α close to -1 we are worried about
            # memory usage and it seems like degree 2 is good enough
            # in that case as well.
            2
        else
            3
        end
        part1_divpα = ArbExtras.enclosure_series(u0.α; degree) do α
            cα, c2α, cαmp0, c2αmp0, z = if α isa Arb
                get!(c_precompute(α), c_precomputed_Arb, α)
            else # α isa ArbSeries
                get!(c_precompute(α), c_precomputed_ArbSeries, α)
            end
            if Arblib.contains_zero(x)
                2c2α / cα^2 *
                abspow(x, -α + u0.p0 - u0.p) *
                (
                    c2αmp0 - 2c2α * cαmp0 / cα + (
                        z / 2 * abspow(x, 2 + 2α - u0.p0) +
                        c2α * (cαmp0 / cα)^2 * abspow(x, u0.p0)
                    )
                )
            else
                2c2α / cα^2 *
                abspow(x, -α + u0.p0 - u0.p) *
                (
                    c2αmp0 - 2c2α * cαmp0 / cα +
                    (z / 2 + c2α * (cαmp0 / cα)^2 * abspow(x, -2 - 2α + 2u0.p0)) *
                    abspow(x, 2 + 2α - u0.p0)
                )
            end
        end

        # u0_part1 / x^-u0.α
        u0_part1_divα = ArbExtras.enclosure_series(u0.α, degree = 1) do α
            cα, c2α, cαmp0, _, _ = if α isa Arb
                get!(c_precompute(α), c_precomputed_Arb, α)
            else # α isa ArbSeries
                get!(c_precompute(α), c_precomputed_ArbSeries, α)
            end
            2c2α / cα^2 * (cα - cαmp0 * abspow(x, u0.p0))
        end

        # u0_part2 / x^p
        u0_part2_divp = eval_expansion(u0, u0_expansion, x, offset = -u0.p)

        # u0_part2 / x^-α
        u0_part2_divα = eval_expansion(u0, u0_expansion, x, offset_i = -1)

        # Hu0_part2 / x^(p - α)
        Hu0_part2_divpα =
            eval_expansion(u0, Hu0_expansion, x, offset = -u0.p, offset_i = -1)

        # Hu0_part3 / x^(p - α)
        Hu0_part3_divpα = let v1 = 2 + u0.α - u0.p, v2 = -2 - 2u0.α, tmp = zero(v2)
            -clausenc_expansion_odd_s_singular_K3(1) *
            ArbExtras.enclosure_series(x) do x
                S = zero(x)
                for j = 1:min(skip_singular_j_until_int, u0.N0)
                    # Compute
                    # S += u0.a[j] * x_pow_s_x_pow_t_m1_div_t(x, v1, v2 + j * p0)
                    Arblib.fma!(tmp, u0.p0, j, v2)
                    term = x_pow_s_x_pow_t_m1_div_t(x, v1, tmp)
                    if x isa Arb
                        Arblib.addmul!(S, term, u0.a[j])
                    else
                        Arblib.mul!(term, term, u0.a[j])
                        Arblib.add!(S, S, term)
                    end
                end
                S
            end
        end

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

If `u0.use_bhkdv` is true then it gives the argument `skip_main` when
computing the expansion. This means that it doesn't include the term
```
u0.a[0] * (C1 - C2 * abs(x)^u0.p0) * abs(x)^-α
```
where `C1` and `C2` are the leading terms in the asymptotic expansions
of `clausencmzeta(x, 1 - u0.α)` and `clausencmzeta(x, 1 - u0.α +
u0.p0)` respectively, in the expansion. This term is instead added
separately. This gives better enclosures for `α` close to `-1`.
"""
function inv_u0_normalised(u0::FractionalKdVAnsatz{Arb}; M::Integer = 5, ϵ::Arb = one(Arb))
    expansion = u0(ϵ, AsymptoticExpansion(), skip_main = u0.use_bhkdv; M)

    C1::Arb, C2::Arb = if u0.use_bhkdv
        let s = 1 - u0.α
            C1 = clausenc_expansion_main(s)
            C2 = clausenc_expansion_main(s + u0.p0)

            C1, C2
        end
    else
        indeterminate(Arb), indeterminate(Arb)
    end

    return x::Union{Arb,ArbSeries} -> begin
        @assert (x isa Arb && abs(x) <= ϵ) || (x isa ArbSeries && abs(Arblib.ref(x, 0)) <= ϵ)

        if u0.use_bhkdv
            term = u0.a[0] * (C1 - C2 * abspow(x, u0.p0))
        else
            term = zero(x)
        end

        return inv(term + eval_expansion(u0, expansion, x, offset_i = -1))
    end
end

"""
    defect(u0::FractionalKdVAnsatz, xs::AbstractVector)

Return a function such that defect(u0, xs)(a, b) computes
defect(u0)(x) on the points x ∈ xs with u0.a and u0.b set to `a` and
`b` respectively. Does this in an efficient way by precomputing as
much as possible.

This is used in [`_findbs`](@ref) for numerically finding values for
`b`.
"""
function defect(u0::FractionalKdVAnsatz, xs::AbstractVector)
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
    defect(u0::FractionalKdVAnsatz, evaltype::Symbolic; M::Integer = 5)

Return a function such that `defect(u0, evaltype, N)(a)` computes the
coefficients in the asymptotic expansion with indices `3` to `u0.N0 +
1` using the values from `a`.

This is used in [`_findas`](@ref) for numerically finding values for
`a`. An alternative implementation is given by [`defect2`](@ref).
"""
function defect(u0::FractionalKdVAnsatz{T}, ::Symbolic; M::Integer = 5) where {T}
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
    defect2(u0::FractionalKdVAnsatz, evaltype::Symbolic; M::Integer = 5)

Return a function such that `defect2(u0, evaltype, N)(a)` computes the
coefficients in the asymptotic expansion with indices `3` to `N0 + 2`
using the values from `a`.

This is an alternative implementation to [`defect`](@ref). In
principle they give the same result but this one is more optimized. In
practice they give slightly different results due to rounding errors
and this is sometimes important.

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
function defect2(
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

    # Check if there are any exponents of type 1
    skip_singular_analytic = !any(exponent -> exponent[1] == 1, exponents)

    # Next step is to precompute the coefficients in the expansions of
    # the Clausen functions

    # In some cases the parameters are very close to singularities of
    # the functions and it is therefore beneficial to do the
    # calculations in Arb and then convert back to T
    α = convert(Arb, u0.α)
    two_α = Arblib.mul_2exp!(zero(α), α, 1)
    p0 = convert(Arb, u0.p0)

    # For α close to -1 this method has a high memory footprint. To
    # reduce memory pressure we do most of the calculations using
    # inplace arithmetic.
    s = zero(α)
    tmp1 = zero(α)
    tmp2 = zero(α)

    # Compute gamma(s) * cospi(s / 2)
    gamma_cospi!(res, s, tmp) = begin
        Arblib.mul_2exp!(res, s, -1)
        Arblib.cos_pi!(res, res)
        Arblib.gamma!(tmp, s)
        Arblib.mul!(res, res, tmp)
    end

    # Compute (-1)^m * zeta(s - 2m) / factorial(2m)
    zeta_div_factorial!(res, s, m) = begin
        Arblib.sub!(res, s, 2m)
        Arblib.zeta!(res, res)
        Arblib.mul!(res, res, (-1)^m)
        Arblib.div!(res, res, factorial(2m))
    end

    u0_precomputed_singular = map(0:u0.N0) do j
        # s = α - j * p0
        Arblib.fma!(s, p0, -j, α)
        # gamma(s) * cospi(s / 2)
        T(gamma_cospi!(tmp1, s, tmp2))
    end

    u0_precomputed_analytic = Matrix{T}(undef, M, u0.N0 + 1)
    for j = 0:u0.N0
        # s = 1 - α + j * p0
        Arblib.fma!(s, p0, -j, α)
        Arblib.neg!(s, s)
        Arblib.add!(s, s, 1)
        for m = 1:M
            # (-1)^m * zeta(s - 2m) / factorial(2m)
            u0_precomputed_analytic[m, j+1] = T(zeta_div_factorial!(tmp1, s, m))
        end
    end

    Hu0_precomputed_singular = map(0:u0.N0) do j
        # s = 2α - j * p0
        Arblib.fma!(s, p0, -j, two_α)
        if s == -1
            # -gamma(s) * cospi(s / 2) has a removable singularity at
            # -s = -1, where it takes the value π / 2
            convert(T, π) / 2
        else
            # -gamma(s) * cospi(s / 2)
            gamma_cospi!(tmp1, s, tmp2)
            T(Arblib.neg!(tmp1, tmp1))
        end
    end

    Hu0_precomputed_analytic = Matrix{T}(undef, M, u0.N0 + 1)
    for j = 0:u0.N0
        # s = 1 - 2α + j * p0
        Arblib.fma!(s, p0, -j, two_α)
        Arblib.neg!(s, s)
        Arblib.add!(s, s, 1)
        for m = 1:M
            # -(-1)^m * zeta(s - 2m) / factorial(2m)
            zeta_div_factorial!(tmp1, s, m)
            Hu0_precomputed_analytic[m, j+1] = T(Arblib.neg!(tmp1, tmp1))
        end
    end

    return a::AbstractVector -> begin
        @assert length(a) == length(u0_precomputed_singular)

        ## Compute all singular terms

        u0_res_singular = u0_precomputed_singular .* a.parent

        # Compute u0_res_singular^2 / 2
        res_singular = zeros(eltype(u0_res_singular), J + 1)
        if threaded && J >= 256
            Threads.@threads for i = 1:J+1
                @inbounds for j = 1:i÷2
                    res_singular[i] += u0_res_singular[j] * u0_res_singular[i-j+1]
                end
                if isodd(i)
                    @inbounds res_singular[i] += u0_res_singular[i÷2+1]^2 / 2
                end
            end
        else
            @inbounds for i = 1:J+1
                for j = 1:i÷2
                    res_singular[i] += u0_res_singular[j] * u0_res_singular[i-j+1]
                end
                if isodd(i)
                    res_singular[i] += u0_res_singular[i÷2+1]^2 / 2
                end
            end
        end

        # Add singular terms from H(u0)
        @inbounds for k = 1:min(length(Hu0_precomputed_singular), length(res_singular))
            res_singular[k] += Hu0_precomputed_singular[k] * a.parent[k]
        end

        ## Compute all analytic terms

        u0_res_analytic = u0_precomputed_analytic * a.parent

        # Compute u0_res_analytic^2 / 2
        res_analytic = zeros(eltype(u0_res_analytic), M)
        @inbounds for i = 1:M
            for j = 1:(i-1)÷2
                res_analytic[i] += u0_res_analytic[j] * u0_res_analytic[i-j]
            end
            if iseven(i)
                res_analytic[i] += u0_res_analytic[i÷2]^2 / 2
            end
        end

        Hu0_res_analytic = Hu0_precomputed_analytic * a.parent
        # Add analytic terms from H(u0)
        @inbounds for k = 1:min(length(Hu0_res_analytic), length(res_analytic))
            res_analytic[k] += Hu0_res_analytic[k]
        end

        # Compute all terms which are singular terms times an analytic
        # term. In some cases there are no such terms in the expansion
        # and then we skip the computation.
        res_singular_analytic = if skip_singular_analytic
            Matrix{eltype(u0_res_singular)}(undef, 0, 0)
        else
            # Compute u0_res_singular * u0_res_analytic / 2
            u0_res_singular * transpose(u0_res_analytic)
        end

        res = resize!(u0_res_singular, u0.N0) # Reuse u0_res_singular
        @inbounds for k in eachindex(exponents, res)
            i, j, m = exponents[k]
            if i == 0
                res[k] = res_analytic[m]
            elseif i == 1
                @assert !skip_singular_analytic
                res[k] = res_singular_analytic[j+1, m]
            else
                res[k] = res_singular[j+1]
            end
        end

        return res
    end
end
