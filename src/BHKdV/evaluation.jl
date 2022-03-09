# PROVE: That all occurrences of force = true indeed are zero

"""
    eval_expansion(u0::BHKdVAnsatz, expansion, x)

Evaluate the given expansion.

It requires that `0 <= x < 1`, any negative parts of `x` are ignored.

The terms are stored as `((p, q, i, j, k, l, m), y)`. The parameters `(i, j,
k l, m)` correspond to the term
```
y * x^(i * α + j * p0 - k*u0.v0.α + l*u0.v0.p0 + m)
```
where `α ∈ (-1, -1 + u0.ϵ]` and `p0 = 1 + α + (1 + α)^2 / 2`.

The parameter `p` corresponds to multiplication by the factor
```
a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) * x^-α
```
to the power `p`, which is the part of the expansion for the main term
which is not even powers of `x`. The parameter `q` corresponds to
multiplication by the factor
```
-a0 * (
    gamma(2α) * cospi(α) * x^(-2α) -
    gamma(2α - p0) * cospi((2α - p0) / 2) * x^(-2α + p0) +
    (-zeta(1 - 2α - 2) / 2 + zeta(1 - 2α + p0 - 2) / 2) * x^2
)
```
to the power `q`, which is the part of the expansion for `H` applied
to the main term which is not on the form `x^2m` for `m >= 2`.

# Handling of `p` and `q`
The method currently does not support `q != 0` and only supports `p ==
0` and `p == 1`. These cases are handled specially in [`F0`](@ref) and
for that reason we don't bother implementing them here.

When `x != 0` the function
```
a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0)
```
has a removable singularity at `α = -1`. This removable singularity
can be handled using [`fx_div_x`](@ref). For `x = 0` this doesn't
work. For now we use that it converges to
```
(1 - γ - log(x)) / π
```
and that this gives an upper bound.
- **PROVE:** That `(1 - γ - log(x)) / π` gives an upper bound.
  Alternatively we don't use this method for `x = 0` (using
  [`inv_u0_bound`](@ref) instead).

The arguments `use_approx_p_and_q` and `alpha_interval` are used for
computing approximate values and is mainly intended for testing.
"""
function eval_expansion(
    u0::BHKdVAnsatz{Arb},
    expansion::AbstractDict{NTuple{7,Int},Arb},
    x::Arb;
    use_approx_p_and_q = false,
    alpha_interval = :full,
)
    @assert x < 1

    # We only care about the non-negative part of x
    x = Arblib.nonnegative_part!(zero(x), x)

    # Enclosure of α
    if alpha_interval == :full
        α = Arb((-1, -1 + u0.ϵ))
    elseif alpha_interval == :endpoint
        α = -1 + u0.ϵ
    else
        throw(ArgumentError("unexpected value alpha_interval = $alpha_interval"))
    end

    # In-place method for computing the exponent i * α + j * p0 -
    # k*u0.v0.α + l*u0.v0.p0 + m in a way that also accounts for
    # the dependence between α and p0
    exponent = zero(α)

    # Precomputed values and buffers for computing the exponent
    lower = zero(exponent)
    upper = zero(exponent)
    # Upper bounds of α and p0
    α_upper = -1 + u0.ϵ
    p0_upper = u0.ϵ + u0.ϵ^2 / 2
    # Enclosure of p0 - α
    p0mα = 1 + Arblib.nonnegative_part!(zero(α), (1 + α)^2) / 2
    _exponent!(exponent, i, j, k, l, m) = begin
        # Compute the part i * α + j * p0 + m
        # Note that i * α + j * p0 = 3j / 2 + (i + 2j) * α + α^2 / 2
        # is increasing in α if i + 2j is non-negative.
        if i + 2j >= 0
            # It is increasing in α, evaluated at endpoints

            # Lower bound at α = -1 can be done with integers
            # For α = -1 we get
            # i * α + j * (1 + α + (1 + α)^2 // 2) + m = m - 1
            Arblib.set!(lower, m - i)

            # Upper bound i * α_upper + j * p0_upper + m
            Arblib.set!(upper, m)
            Arblib.addmul!(upper, α_upper, i)
            Arblib.addmul!(upper, p0_upper, j)

            Arblib.union!(exponent, lower, upper)

            # If the lower bound is zero we want to avoid any spurious
            # negative parts
            iszero(lower) && Arblib.nonnegative_part!(exponent, exponent)
        else
            # IMPROVE: We could rewrite the expression to only have α
            # in one place, reducing overestimations.
            Arblib.set!(exponent, m)
            Arblib.addmul!(exponent, α, i + j)
            Arblib.addmul!(exponent, p0mα, j)
        end

        # Add - k*u0.v0.α + l*u0.v0.p0
        Arblib.submul!(exponent, u0.v0.α, k)
        Arblib.addmul!(exponent, u0.v0.p0, l)

        return exponent
    end

    # Compute enclosure of the p-coefficient
    # a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0)
    # Note that this depends on x and is only used when x doesn't contain zero
    if !Arblib.contains_zero(x)
        p_coefficient = let interval = Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))
            extra_degree = 2
            # Enclosure of rgamma(α) / (α + 1)
            rgamma1_div_α = fx_div_x(s -> rgamma(s - 1), interval; extra_degree)
            # Enclosure of rgamma(2α) / (α + 1)
            rgamma2_div_α = fx_div_x(s -> rgamma(2(s - 1)), interval; extra_degree)
            # Enclosure of cospi(α / 2) / (α + 1)
            cos_div_α = fx_div_x(s -> cospi((s - 1) / 2), interval; extra_degree)

            # Enclosure of a0 * (α + 1)
            a0_mul_α =
                2cospi(Arb((-1, -1 + u0.ϵ))) / rgamma2_div_α / (cos_div_α / rgamma1_div_α)^2

            # Enclosure of
            # (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) / (α + 1)
            coefficient_div_α = fx_div_x(interval; extra_degree, force = true) do r
                if Arblib.contains_zero(r[0])
                    # Enclosure of rgamma(α) / (α + 1)
                    rgamma1_div_α = fx_div_x(s -> rgamma(s - 1), r; extra_degree)
                    # Enclosure of rgamma(α - p0) / (α + 1)^2
                    rgamma2_div_α2 =
                        fx_div_x(s -> rgamma(-1 - s^2 / 2), r, 2; extra_degree)
                    # Enclosure of cospi(α / 2) / (α + 1)
                    cos1_div_α = fx_div_x(s -> cospi((s - 1) / 2), r; extra_degree)
                    # Enclosure of cospi((α - p0) / 2) / (α + 1)^2
                    cos2_div_α2 =
                        fx_div_x(s -> cospi((-1 - s^2 / 2) / 2), r, 2; extra_degree)

                    cos1_div_α / rgamma1_div_α -
                    cos2_div_α2 / rgamma2_div_α2 * x^(r + r^2 / 2)
                else
                    gamma(r - 1) * cospi((r - 1) / 2) -
                    gamma(-1 - r^2 / 2) * cospi((-1 - r^2 / 2) / 2) * x^(r + r^2 / 2)
                end
            end

            a0_mul_α * coefficient_div_α
        end
    end

    res = zero(x)

    for ((p, q, i, j, k, l, m), y) in expansion
        if !iszero(y)
            if iszero(p)
                _exponent!(exponent, i, j, k, l, m)

                term = abspow(x, exponent)
            elseif isone(p)
                # Add -α to the exponent coming from the p factor
                _exponent!(exponent, i - 1, j, k, l, m)

                if iszero(x)
                    if Arblib.ispositive(exponent)
                        term = zero(x)
                    else
                        term = Arblib.indeterminate!(zero(x))
                    end
                elseif Arblib.contains_zero(x)
                    if Arblib.ispositive(exponent)
                        lower = zero(x)
                    else
                        lower = Arblib.indeterminate!(zero(x))
                    end

                    # FIXME: Prove that this gives an upper bound
                    upper = let x = ubound(Arb, x)
                        (1 - Arb(Irrational{:γ}()) - log(x)) / π * abspow(x, exponent)
                    end

                    term = Arb((lower, upper))
                else
                    term = abspow(x, exponent)
                    Arblib.mul!(term, term, p_coefficient)
                end
            else
                use_approx_p_and_q ||
                    throw(ArgumentError("only p == 0 or p == 1 supported, got p = $p"))

                _exponent!(exponent, i, j, k, l, m)

                term = abspow(x, exponent)

                let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, a0 = finda0(α)
                    p_factor =
                        a0 *
                        (
                            gamma(α) * cospi(α / 2) -
                            gamma(α - p0) * cospi((α - p0) / 2) * x^p0
                        ) *
                        x^-α
                    term *= p_factor^p
                end
            end

            # We don't support evaluation of non-zero q
            if !iszero(q)
                use_approx_p_and_q ||
                    throw(ArgumentError("only q == 0 is supported, got q = $q"))
                let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, a0 = finda0(α)
                    q_factor =
                        -a0 * (
                            gamma(2α) * cospi(α) * x^(-2α) -
                            gamma(2α - p0) * cospi((2α - p0) / 2) * x^(-2α + p0) +
                            (-zeta(1 - 2α - 2) / 2 + zeta(1 - 2α + p0 - 2) / 2) * x^2
                        )
                    term *= q_factor^q
                end
            end

            #res += y * term
            Arblib.addmul!(res, y, term)
        end
    end

    return res
end

"""
    expansion_ispositive(u0::BHKdVAnsatz, expansion, ϵ)

Attempt to prove that the given expansion is positive on the interval
`(0, ϵ]`. Returns true on success and false on failure. It requires
that `0 < ϵ < 1`

The method can only handle expansions on certain simple forms. More
precisely it requires that all terms have a key of the form
```
(0, 0, 0, 0, 1, l, 0)
```
or
```
(0, 0, 0, 0, 0, 0, m)
```
with `l >= 1` and `m >= 2`

It first considers all keys of the form `(0, 0, 0, 0, 1, l, 0)`, they
correspond to terms of the form
```
y * x^^(-u0 - v0.α + l*u0.v0.p0)
```
Let `y₁` be the coefficient for `l = 1` and `S` be the sum of the
coefficients with `l > 1`. It checks that `S` is negative and that `y₁
+ S > 0`, ensuring us that this sum is positive for all `0 < x < 1`
and a lower bound is given by `(y₁ + S) * x^(-u0.v0.α +
u0.v0.p0)` for all `0 < x < 1`.

The next step is to prove that the sum of the terms with `m > 0` are
smaller than (y₁ + S) * x^(-u0.v0.α + u0.v0.p0)`. This is done
by noting that it's enough to check it for `x = ϵ`

**TODO:** Check that this is correct.
"""
function expansion_ispositive(
    u0::BHKdVAnsatz{Arb},
    expansion::AbstractDict{NTuple{7,Int},Arb},
    ϵ::Arb,
)
    @assert 0 < ϵ < 1
    @assert 0 < -u0.v0.α + u0.v0.p0 < 2

    # Isolate all keys of the from (0, 0, 0, 0, 1, l, 0)
    expansion_1 = filter(expansion) do ((p, q, i, j, k, l, m), y)
        p == q == i == j == m == 0 && k == 1 && l >= 1
    end

    # Isolate all keys of the from (0, 0, 0, 0, 0, 0, m)
    expansion_2 = filter(expansion) do ((p, q, i, j, k, l, m), y)
        p == q == i == j == k == l == 0 && m >= 2
    end

    # Check that this was all keys
    @assert length(expansion) == length(expansion_1) + length(expansion_2)

    y₁ = expansion_1[(0, 0, 0, 0, 1, 1, 0)]
    # Remove they key (0, 0, 0, 0, 1, 1, 0) from the expansion and sum
    # the rest of the values.
    delete!(expansion_1, (0, 0, 0, 0, 1, 1, 0))
    S = sum(values(expansion_1))

    S < 0 || return false
    y₁ + S > 0 || return false

    a = (y₁ + S) * ϵ^(-u0.v0.α + u0.v0.p0)
    b = eval_expansion(u0, expansion_2, ϵ)

    a > b || return false

    return true
end


"""
    (u0::BHKdVAnsatz)(x, ::Ball)

Evaluate the ansatz `u0` at the point `x`.

The tail term is evaluated directly.

To evaluate the main term, given by
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
we make use of the fact that this converges to
```
2 / π^2 * clausencmzeta(x, 2, 1)
```
, see [`lemma_bhkdv_main_term_limit`](@ref), which is the main term
for `BHAnsatz`, as `α -> -1`.

Combining the above with [`lemma_bhkdv_main_term_monotonicity`](@ref)
we notice that it is enough to evaluate the limiting expression as
well as at `α = -1 + u0.ϵ` to get an enclosure.

- **FIXME:** The above approach doesn't work directly for `ArbSeries`
because we don't have monotonicity on `α` for all the derivatives.
However for now we still compute at the endpoints and take the union,
this means that the enclosure will not be correct for all value of `x`.
"""
function (u0::BHKdVAnsatz{Arb})(x::Union{Arb,ArbSeries}, ::Ball)
    # Main term

    # Compute limiting expression and at upper bound for α
    res_lower = 2 / Arb(π)^2 * clausencmzeta(x, 2, 1)
    res_upper = let α = -1 + u0.ϵ, a0 = finda0(α), p0 = 1 + α + (1 + α)^2 / 2
        a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
    end

    if x isa Arb
        res = Arb((res_lower, res_upper))
    elseif x isa ArbSeries
        @warn "Non-rigorous bounds for main term with ArbSeries" maxlog = 1
        coefficients = union.(Arblib.coeffs(res_lower), Arblib.coeffs(res_upper))
        res = ArbSeries(coefficients)
    end

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.N0
        s = 1 - u0.v0.α + j * u0.v0.p0
        res += u0.v0.a[j] * clausencmzeta(x, s)
    end

    # Fourier terms
    for n = 1:u0.v0.N1
        res += u0.v0.b[n] * (cos(n * x) - 1)
    end

    return res
end

(u0::BHKdVAnsatz)(x, ::Asymptotic; M::Integer = 3) =
    eval_expansion(u0, u0(x, AsymptoticExpansion(); M), x)

"""
    (u0::BHKdVAnsatz)(x, ::AsymptoticExpansion; M = 3)

Return a dictionary containing the terms in the asymptotic expansion
of `u0` which can then be evaluated with [`eval_expansion`](@ref).

The highest term, `x^2M`, is an error term is which makes sure that
evaluation of the expansion gives an enclosure of the result when
evaluated at `|y| < |x|`.

For the tail term the expansions are easily computed exactly like for
`BHAnsatz`. For the main term we have to be a bit more careful.

For the main term the coefficients in front of `x^2m` is given by
```
a0 * (-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / factorial(2m)
```
It has a removable singularity at `α = -1`. To compute an enclosure we
use that
```
a0 = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))^2
```
and write it as
```
(-1)^m * 2cospi(α)
    * inv(rgamma(2α) / (1 + α))
    * (rgamma(α) / (1 + α))^2
    * inv(cospi(α / 2) / (1 + α))^2
    * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / (1 + α)
    / factorial(2m)
```
, where `rgamma` is the reciprocal gamma function given by `rgamma(s)
= inv(gamma(s))`, and use [`fx_div_x`](@ref) to compute enclosures of
the individual factors.
- **TODO:** Compute enclosure of remainder term in similar way.

The only remaining part of the expansion of the main term is
```
a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) * x^-α
```
which we don't evaluate at all yet. Instead store implicitly in the
expansion.

See [`eval_expansion`](@ref) for more details about how the
coefficients are stored.
"""
function (u0::BHKdVAnsatz{Arb})(x, ::AsymptoticExpansion; M::Integer = 3)
    @assert M >= 3

    res = OrderedDict{NTuple{7,Int},Arb}()

    # Initiate even powers of x
    for m = 1:M
        res[(0, 0, 0, 0, 0, 0, 2m)] = 0
    end

    # Main term

    # Leading term
    res[(1, 0, 0, 0, 0, 0, 0)] = 1

    # x^2m terms
    # Interval corresponding to 1 + α
    interval = Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))
    # Enclosure of rgamma(α) / (α + 1)
    rgamma1_div_α = fx_div_x(s -> rgamma(s - 1), interval, extra_degree = 2)
    # Enclosure of rgamma(2α) / (α + 1)
    rgamma2_div_α = fx_div_x(s -> rgamma(2(s - 1)), interval, extra_degree = 2)
    # Enclosure of cospi(α / 2) / (α + 1)
    cos_div_α = fx_div_x(s -> cospi((s - 1) / 2), interval, extra_degree = 2)
    for m = 1:M-1
        # Enclosure of
        # (zeta(1 - (s - 1) - 2m) - zeta(2 + (1 + (s - 1))^2 / 2 - 2m)) / (α + 1)
        zeta_div_α = if m == 1
            # zeta(x::ArbSeries) doesn't handle balls containing
            # zero but centered at a negative number well. For
            # that reason we take a symmetric interval in this
            # case.
            fx_div_x(
                s -> zeta(1 - (s - 1) - 2m) - zeta(2 + (1 + (s - 1))^2 / 2 - 2m),
                union(-interval, interval),
                extra_degree = 2,
                force = true,
            )
        else
            fx_div_x(
                s -> zeta(1 - (s - 1) - 2m) - zeta(2 + (1 + (s - 1))^2 / 2 - 2m),
                interval,
                extra_degree = 2,
            )
        end

        coefficient =
            (-1)^m *
            2cospi((Arb((-1, -1 + u0.ϵ)))) *
            inv(rgamma2_div_α) *
            rgamma1_div_α^2 *
            inv(cos_div_α)^2 *
            zeta_div_α / factorial(2m)

        res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
    end

    # Remainder term for main term
    # FIXME: Compute proper remainder term
    @warn "Non-rigorous error term implemented for main term" maxlog = 1
    Arblib.add_error!(
        res[(0, 0, 0, 0, 0, 0, 2M)],
        2abs(dzeta(Arb(2 - 2M)) / factorial(2M)) * u0.v0.a0,
    )

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.N0
        s = 1 - u0.v0.α + j * u0.v0.p0
        C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)
        res[(0, 0, 0, 0, 1, j, 0)] = C * u0.v0.a[j]
        for m = 1:M-1
            res[(0, 0, 0, 0, 0, 0, 2m)] += p[2m] * u0.v0.a[j]
        end
        res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.a[j]
    end

    # Fourier terms
    if !iszero(u0.v0.N1)
        for m = 1:M-1
            res[(0, 0, 0, 0, 0, 0, 2m)] +=
                (-1)^m * sum(Arb(n)^(2m) * u0.v0.b[n] for n = 1:u0.v0.N1) / factorial(2m)
        end
        Arblib.add_error!(
            res[(0, 0, 0, 0, 0, 0, 2M)],
            sum(Arb(n)^(2M) * abs(u0.v0.b[n]) for n = 1:u0.v0.N1) / factorial(2M),
        )
    end

    return res
end

"""
    H(u0::BHKdVAnsatz, ::Ball)

Returns a function such that `H(u0, Ball())(x)` evaluates
``H^α[u0](x)``.

The transform of the main term is given by
```
-a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
```
we make use of the fact that this converges to
```
-2 / π^2 * clausencmzeta(x, 3, 1)
```
, see [`lemma_bhkdv_main_term_H_limit`](@ref), which is `H` applied to
the the main term for `BHAnsatz`, as `α -> -1`.

Combining the above with
[`lemma_bhkdv_main_term_H_monotonicity`](@ref) we notice that it is
enough to evaluate the limiting expression as well as at `α = -1 +
u0.ϵ` to get an enclosure.

- **FIXME:** The above approach doesn't work directly for `ArbSeries`
because we don't have monotonicity on `α` for all the derivatives.
However for now we still compute at the endpoints and take the union,
this means that the enclosure will not be correct for all value of `x`.

For the tail term we need to make sure that we correctly handle the
fact that the transform depends on the value of `α`. As long as `u0.ϵ`
is sufficiently small we get good enough bounds by just using it as a
ball directly.
"""
function H(u0::BHKdVAnsatz{Arb}, ::Ball)
    return x::Union{Arb,ArbSeries} -> begin
        # Main term

        # Compute limiting expression and at upper bound for α
        res_lower = let α = -1 + u0.ϵ, a0 = finda0(α), p0 = 1 + α + (1 + α)^2 / 2
            -a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
        end
        res_upper = -2 / Arb(π)^2 * clausencmzeta(x, 3, 1)

        if x isa Arb
            res = Arb((res_lower, res_upper))
        elseif x isa ArbSeries
            @warn "Non-rigorous bounds for main term with ArbSeries" maxlog = 1
            coefficients = union.(Arblib.coeffs(res_lower), Arblib.coeffs(res_upper))
            res = ArbSeries(coefficients)
        end

        # Tail term

        let α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α
            # IMPROVE: We might be able to compute better enclosures
            # in α, allowing us to work with larger values of u0.ϵ, by
            # expanding in s. This would require a fair amount of work
            # though, we would probably have to extract the formula
            # from _clausenc_zeta.

            # Clausen terms
            for j = 1:u0.v0.N0
                s = 1 - α - u0.v0.α + j * u0.v0.p0
                res -= u0.v0.a[j] * clausencmzeta(x, s)
            end

            # Fourier terms
            for n = 1:u0.v0.N1
                res -= u0.v0.b[n] * n^α * (cos(n * x) - 1)
            end
        end

        return res
    end
end

function H(
    u0::BHKdVAnsatz,
    ::Asymptotic;
    M::Integer = 3,
    skip_singular_j_until = 0,
    use_approx_p_and_q = false,
)
    f = H(u0, AsymptoticExpansion(); M, skip_singular_j_until)

    return x -> eval_expansion(u0, f(x), x; use_approx_p_and_q)
end

"""
    H(u0::BHKdVAnsatz, ::AsymptoticExpansion; M = 3, skip_singular_j_until::Integer = 0,)

Return a dictionary containing the terms in the asymptotic expansion
of `u0` which can then be evaluated with [`eval_expansion`](@ref).

The highest term, `x^2M`, is an error term is which makes sure that
evaluation of the expansion gives an enclosure of the result when
evaluated at `|y| < |x|`.

For the main term the coefficients in front of `x^2m` is given by
```
-a0 * (-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / factorial(2m)
```
For `m >= 2` it has a removable singularity and we handle it in the
same way as in [`u0`](@ref).
- **TODO:** Compute enclosure of remainder term in similar way.

The remaining part of the expansion of `H` applied to the main term is
```
-a0 * (
    gamma(2α) * cospi(α) * abs(x)^(-2α) -
    gamma(2α - p0) * cospi((1 - 2α + p0) / 2) * abs(x)^(-2α + p0)
    (zeta(1 - 2α - 2) / 2 - zeta(1 - 2α + p0 - 2) / 2) * abs(x)^2
)
```
which we don't evaluate at all yet. Instead we store it implicitly in
the expansion.

For both the Clausen terms and the Fourier terms we let `α` be a ball.
This gives good enclosures for the Fourier terms. For the Clausen
terms it give good enclosures unless `j` is small. For small values of
`j` the two terms
```
gamma(α + u0.v0.α - j * u0.v0.p0) * cospi((α + u0.v0.α - j * u0.v0.p0) / 2) *
    x^-(α + u0.v0.α - j * u0.v0.p0)
```
and
```
-zeta(-1 - α + u0.v0.α - j * u0.v0.p0) / 2 * x^2
```
have very large cancellations.

If `skip_singular_j_until` is greater than zero then skip the two
above terms for all Clausen functions in the tail from `j = 1` to
`skip_singular_j_until`. This is used in [`F0`](@ref) where these
terms are handled separately.

If `approximate_j_one_singular` is true it only computes an
approximation of the two terms for `j = 1`. This is only for testing.

The argument `alpha_interval` can be set to either `:full` or
`:endpoint`. In the former case it uses the full interval of `α` when
computing the tail terms, in the latter it uses `α = -1 + u0.ϵ`, this
can be used for testing.

See [`eval_expansion`](@ref) for more details about how the
coefficients are stored.
"""
function H(
    u0::BHKdVAnsatz{Arb},
    ::AsymptoticExpansion;
    M::Integer = 3,
    skip_singular_j_until::Integer = 0,
    approximate_j_one_singular::Bool = false,
    alpha_interval = :full,
)
    @assert M >= 3

    skip_singular_j_until > u0.v0.N0 && throw(
        ArgumentError("can't skip more j-terms than there are, j = $j, N0 = $(u0.v0.N0)"),
    )

    # Enclosure of α
    if alpha_interval == :full
        α = Arb((-1, -1 + u0.ϵ))
    elseif alpha_interval == :endpoint
        α = -1 + u0.ϵ
    else
        throw(ArgumentError("unexpected value alpha_interval = $alpha_interval"))
    end

    return x -> begin
        res = OrderedDict{NTuple{7,Int},Arb}()

        # Initiate even powers of x
        for m = 1:M
            res[(0, 0, 0, 0, 0, 0, 2m)] = 0
        end

        # Main term

        # Three leading terms
        res[(0, 1, 0, 0, 0, 0, 0)] = 1

        # x^2m terms with m >= 2
        # Interval corresponding to 1 + α
        interval = Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))
        # Enclosure of rgamma(α) / (α + 1)
        rgamma1_div_α = fx_div_x(s -> rgamma(s - 1), interval, extra_degree = 2)
        # Enclosure of rgamma(2α) / (α + 1)
        rgamma2_div_α = fx_div_x(s -> rgamma(2(s - 1)), interval, extra_degree = 2)
        # Enclosure of cospi(α / 2) / (α + 1)
        cos_div_α = fx_div_x(s -> cospi((s - 1) / 2), interval, extra_degree = 2)
        for m = 2:M-1
            # Enclosure of
            # (zeta(1 - 2(s - 1) - 2m) - zeta(2 - (s - 1) + (1 + (s - 1))^2 / 2 - 2m)) / (α + 1)
            zeta_div_α = fx_div_x(
                s ->
                    zeta(1 - 2(s - 1) - 2m) - zeta(2 - (s - 1) + (1 + (s - 1))^2 / 2 - 2m),
                interval,
                extra_degree = 2,
                force = true,
            )

            coefficient =
                -(-1)^m *
                2cospi(Arb((-1, -1 + u0.ϵ))) *
                inv(rgamma2_div_α) *
                rgamma1_div_α^2 *
                inv(cos_div_α)^2 *
                zeta_div_α / factorial(2m)

            res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
        end

        # Remainder term for main term
        # FIXME: Compute proper remainder term
        @warn "No error bounds for error term of main term" maxlog = 1
        Arblib.add_error!(
            res[(0, 0, 0, 0, 0, 0, 2M)],
            2abs(dzeta(Arb(3 - 2M)) / factorial(2M)) * u0.v0.a0,
        )

        # Tail term

        # Clausen terms
        for j = 1:u0.v0.N0
            s = 1 - α - u0.v0.α + j * u0.v0.p0
            C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

            if j > skip_singular_j_until
                if isone(j) && approximate_j_one_singular
                    C2, _, p2, _ = let s = 1 - (-1 + u0.ϵ) - u0.v0.α + j * u0.v0.p0
                        clausenc_expansion(x, s, M, skip_constant = true)
                    end
                    res[(0, 0, -1, 0, 1, j, 0)] = -C2 * u0.v0.a[j]
                    res[(0, 0, 0, 0, 0, 0, 2)] -= p2[2] * u0.v0.a[j]
                else
                    res[(0, 0, -1, 0, 1, j, 0)] = -C * u0.v0.a[j]
                    res[(0, 0, 0, 0, 0, 0, 2)] -= p[2] * u0.v0.a[j]
                end
            end

            for m = 2:M-1
                res[(0, 0, 0, 0, 0, 0, 2m)] -= p[2m] * u0.v0.a[j]
            end
            res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.a[j]
        end

        # Fourier terms
        if !iszero(u0.v0.N1)
            for m = 1:M-1
                res[(0, 0, 0, 0, 0, 0, 2m)] -=
                    (-1)^m * sum(n^α * Arb(n)^(2m) * u0.v0.b[n] for n = 1:u0.v0.N1) /
                    factorial(2m)
            end
            Arblib.add_error!(
                res[(0, 0, 0, 0, 0, 0, 2M)],
                sum(n^α * Arb(n)^(2M) * abs(u0.v0.b[n]) for n = 1:u0.v0.N1) / factorial(2M),
            )
        end

        return res
    end
end

function D(u0::BHKdVAnsatz, ::Asymptotic; M::Integer = 3, skip_singular_j_until = 0)
    f = D(u0, AsymptoticExpansion(); M, skip_singular_j_until)

    return x -> eval_expansion(u0, f(x), x)
end

function D(
    u0::BHKdVAnsatz,
    evaltype::AsymptoticExpansion;
    M::Integer = 3,
    skip_singular_j_until = 0,
    alpha_interval = :full,
)
    f = x -> u0(x, evaltype; M)
    g = H(u0, evaltype; M, skip_singular_j_until, alpha_interval)

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

        return expansion
    end
end

"""
    F0_bound(u0::BHKdVAnsatz{Arb}, evaltype::Ball = Ball())

Return a function `f` such that the absolute value of `f(x)` is an
upper bound of the absolute value of `F0(u0)(x)`.

More precisely this computes
```
D(u0)(x) / (u0.w(x) * u0.v0(x))
```
Since `u0.v0(x)` gives a lower bound of `u0(x)`, by
[`lemma_bhkdv_monotonicity_alpha`](@ref), this gives a value which has
the same sign as `F0(x)` but is larger in magnitude. This holds as
long as `u0.v0(x)` is positive at least, which is easily checked.

**IMPROVE:** The computation of `u0.v0(x)` and `u0(x)` have many
calculations in common. We could improve performance by using this.
"""
function F0_bound(u0::BHKdVAnsatz{Arb}, evaltype::Ball = Ball())
    g = D(u0, evaltype)

    return x::Union{Arb,ArbSeries} -> begin
        invweight = inv(u0.w(x))

        isfinite(invweight) || return invweight

        invu0v0 = inv(u0.v0(x))

        isfinite(invu0v0) || return invu0v0

        if (x isa Arb && !Arblib.ispositive(invu0v0)) ||
           (x isa ArbSeries && !Arblib.ispositive(Arblib.ref(invu0v0, 0)))
            error("expected u0.v0(x) to be positive, got inv(u0.v0(x)) = $invu0v0")
        end

        return g(x) * invu0v0 * invweight
    end
end


"""
    F0(u0::BHKdVAnsatz{Arb}, ::Asymptotic)

Returns a function such that an **upper bound** of `F0(u0)(x)` is
computed accurately for small values of `x`.

It requires that `0 <= x < 1`, any negative parts of `x` are ignored.

Recall that the expression we are interested in bounding is
```
(u0(x)^2 / 2 + H(u0)(x)) / (u0.w(x) * u0(x))
```
with
```
u0.w(x) = x^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x))
```

# Split into three factors
As a first step we split the expression into three factors which are
all bounded as `x -> 0` that we bound separately.
We write it as
```
(gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x))
* (log(inv(x)) / log(u0.c + inv(x)))
* ((u0(x)^2 / 2 + H(u0)(x)) / (gamma(1 + α) * log(inv(x)) * x^(1 - u0.γ * (1 + α) - α) * (1 - x^p0)))
```
The first factor we is bounded using [`inv_u0_bound`](@ref). The
second factor is bounded by noticing that it is `1` at `x = 0` and
decreasing in `x`.

For the third factor we can get an upper bound, in absolute value, by
taking `γ = 0`, giving us
```
F = (u0(x)^2 / 2 + H(u0)(x)) / (gamma(1 + α) * log(inv(x)) * (1 - x^p0) * x^(1 - α))
```
We now explain how to bound `F`

# Bounding `F`
Getting an accurate bound for `F` requires more work.

We use the asymptotic expansion of `u0(x)^2 / 2 + H(u0)(x)`, but some
of the terms in the expansion require extra care. We therefore split
the expansion into three parts
1. The two leading terms, with keys `(2, 0, 0, 0, 0, 0, 0)` and `(0, 1,
    0, 0, 0, 0, 0)`, which we call `P` and `Q` respectively.
2. The leading terms of the Clausen terms in the tail for small values
   of `j`.
3. The remaining terms.
This splits `F` into three terms, which we will call `T1, T2, T3`.

## Handling `T1`: `P` and `Q`
The terms `P` and `Q` are given by
```
P^2 = a0^2 / 2 * (c(α)^2 - 2c(α) * c(α - p0) * x^p0 + c(α - p0)^2 * x^2p0) * x^-2α
```
and
```
Q = -a0 * ((c(2α) - c(2α - p0) * x^p0) * x^-2α - (zeta(-2α - 1) / 2 - zeta(-2α + p0 - 1) / 2) * x^2)
```
where `c(a) = gamma(a) * cospi(a / 2)`, similar to
[`inv_u0_bound`](@ref).

By construction `a0` is such that the terms with exponent `x^-2α`
cancel out. This leaves us with
```
P + Q = a0 * (
    (c(2α - p0) - a0 * c(α) * c(α - p0)) * x^(-2α + p0) +
    (zeta(-2α - 1) - zeta(-2α + p0 - 1)) / 2 * x^2 +
    a0 * c(α - p0)^2 / 2 * x^(-2α + 2p0)
)
```
where we can note that `-2α + p0 < 2 < -2α + 2p0`. We are now
interested in bounding
```
(P + Q) / (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0))
```
If we cancel the exponents and reorder the factors slightly we get
```
a0 / gamma(1 + α) * (
    (c(2α - p0) - a0 * c(α) * c(α - p0)) * x^(-α + p0 - 1) / (log(x) * (1 - x^p0)) +
    (zeta(-2α - 1) - zeta(-2α + p0 - 1)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0)) +
    a0 * c(α - p0)^2 / 2 * x^(-α + 2p0 - 1) / (log(x) * (1 - x^p0))
)
```
The factor `a0 / gamma(1 + α)` has a removable singularity at `α = -1`
and can be enclosed in the same way as in [`inv_u0_bound`](@ref).

First we focus on the term
```
T11 = (c(2α - p0) - a0 * c(α) * c(α - p0)) * x^(-α + p0 - 1) / (log(x) * (1 - x^p0))
```
This term is small and fairly stable in `α` (positive and increasing).
To compute an enclosure we split it into three factors
```
T111 = x^(-α + p0 - 1) / log(x)
T112 = (1 + α) / (1 - x^p0)
T113 = (c(2α - p0) - a0 * c(α) * c(α - p0)) / (1 + α)
```
For the first one we can directly get an enclosure using that
```
-α + p0 - 1 = -α + (1 + α + (1 + α)^2 / 2) - 1 = (1 + α)^2 / 2
```
For the second one we note that it is increasing in `x` and takes the
value `1 + α` at `x = 0`. For non-zero `x` we can handle the removable
singularity in `α`.

For the third term, `T113`, we note that `a0` can be written as
```
a0 = 2c(2α) / c(α)^2
```
which allows us to simplify it to
```
T113 = (c(2α - p0) - 2c(2α) * c(α - p0) / c(α)) / (1 + α)
```
Again we handle the removable singularity at `α = -1` as before, this
one requires a bit more work though. To do this we rewrite it as
```
(
    (1 + α) * c(2α - p0)
    - 2((1 + α) * c(2α)) * c(α - p0) / c(α)
) / (1 + α)^2
```

We then consider the two remaining terms together since they mostly
cancel out.
```
T12 = (zeta(-2α - 1) - zeta(-2α + p0 - 1)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0)) +
    a0 * c(α - p0)^2 / 2 * x^(-α + 2p0 - 1) / (log(x) * (1 - x^p0))
    =
    ((zeta(-2α - 1) - zeta(-2α + p0 - 1)) +
    a0 * c(α - p0)^2 * x^(-2α + 2p0 - 2)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0))
```
We treat this in a very similar way as `T11`, by splitting it into
three factors.
```
T121 = x^(1 + α) / 2
T122 = (1 + α) / (1 - x^p0)
T123 = ((zeta(-2α - 1) - zeta(-2α + p0 - 1)) + a0 * c(α - p0)^2 * x^(-2α + 2p0 - 2)) /
    ((1 + α) * log(x))
```
The factor `T121` we can enclose directly. The factor `T122` is the
same as the factor `T112` above. We are hence left enclosing `T123`.
To do that we split it into two terms, letting
```
v(α) = zeta(-2α - 1) - zeta(-2α + p0 - 1)
w(α) = -a0 * c(α - p0)^2 = -2c(2α) * c(α - p0)^2 / c(α)^2
```
and noticing that `-2α + 2p0 - 2 = (1 + α)^2` we can write
```
T123 = (v(α) - w(α) * x^((1 + α)^2)) / ((1 + α) * log(x))
```
Adding and subtracting `w(α)` in the numerator we can split this into
two terms
```
T1231 = (v(α) - w(α)) / ((1 + α) * log(x))
T1232 = w(α) * (1 - x^((1 + α)^2)) / ((1 + α) * log(x))
```
For `T1231` it enough to handle the removable singularity of `(v(α) -
w(α)) / (1 + α)` and then multiply by an enclosure of `inv(log(x))`.

For `T1232` we can rewrite it as
```
T1232 = (1 + α) * w(α) * (1 - x^((1 + α)^2)) / ((1 + α)^2 * log(x))
```
We can compute an enclosure of
```
(1 + α) * w(α) = -2(1 + α) * c(2α) * c(α - p0)^2 / c(α)^2 =
```
by handling the removable singularity.

We can compute an enclosure of
```
(1 - x^((1 + α)^2) / ((1 + α)^2 * log(x))
```
By letting `t = (1 + α)^2 * log(x)` to write it as
```
(1 - exp(t)) / t
```
We then handle this similarly to how we do it for the tail Clausen
below.

## Handling `T2`: tail Clausen with small `j`
For small values of `j` the two terms
```
gamma(α + u0.v0.α - j * u0.v0.p0) * cospi((α + u0.v0.α - j * u0.v0.p0) / 2) *
    x^-(α + u0.v0.α - j * u0.v0.p0)
```
and
```
-zeta(-1 - α + u0.v0.α - j * u0.v0.p0) / 2 * x^2
```
in the asymptotic expansion of the Clausen functions in the tail have
very large cancellations. It is therefore beneficial to treat them
together to account for the cancellations. Which `j` we treat
separately like this is determined by the argument
`skip_singular_j_until`, it then handles `j = 1:skip_singular_j_until`
separately.


We are thus interested in bounding
```
-u0.v0.a[j] * clausenc(x, 1 - α - u0.v0.α + j * u0.v0.p0) /
    (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0))
```
We can get an enclosure of `inv(gamma(1 + α) * (1 - x^p0))` by
noticing that it is increasing in `x` and it is hence enough to
compute at the endpoints of `x`. For `x = 0` it is given by
`inv(gamma(1 + α)) = rgamma(1 + α)`. Otherwise we use the same
approach as in [`inv_u0_bound`](@ref) for enclosing it. We are then
interested in enclosing the rest.

Let `r = -u0.v0.α + j * u0.v0.p0 - 1`, then `r > 0` and for
small values of `j` it is very close to zero. We have `1 - α -
u0.v0.α + j * u0.v0.p0 = 2 - α + r`. The sum of the first two
terms in the asymptotic expansion of the Clausen is then given by
```
gamma(α - 1 - r) * cospi((α - 1 - r) / 2) * x^(1 - α + r) -
    zeta(-α + r) / 2 * x^2
```
Dividing by `log(x) * x^(1 - α)` gives us
```
(gamma(α - 1 - r) * cospi((1 - α + r) / 2) * x^r - zeta(-α + r) / 2 * x^(1 + α)) / log(x)
```
Notice that the order of the terms depends on the value of `α`, in
some cases `x^r` is leading and in some cases `x^(1 + α)`. Adding and
subtracting `zeta(-α + r) / 2 * x^r` we can rewrite this as
```
(gamma(α - 1 - r) * cospi((1 - α + r) / 2) - zeta(-α + r) / 2) * x^r / log(x) +
zeta(-α + r) / 2 * (x^r - x^(1 + α)) / log(x)
```
Let
```
T21 = (gamma(α - 1 - r) * cospi((1 - α + r) / 2) - zeta(-α + r) / 2) * x^r / log(x)
T22 = zeta(-α + r) / 2 * (x^r - x^(1 + α)) / log(x)
```

For `T21` we rewrite it using that
```
gamma(α - 1 - r) = gamma(α + 2 - r) / ((α + 1 - r) * (α - r) * (α - 1 - r))

zeta(-α + r) / 2 = zeta_deflated(-α + r) / 2 - 1 / 2(1 + α - r)
```
giving us
```
T21 = (
    gamma(α + 2 - r) / ((α - r) * (α - 1 - r)) * cospi((1 - α + r) / 2) + 1 // 2
) / (α + 1 - r) - zeta_deflated(-α + r) / 2,
```
This formulation is much more stable and can be accurately enclosed.
For the case when `α + 1 - r` overlaps zero we can handle the
removable singularity.

For `T22` we have to compute an enclosure of
```
zeta(-α + r) / 2 * (x^r - x^(1 + α)) / log(x)
```
Using [`zeta_deflated`](@ref) we have
```
zeta(-α + r) = zeta_deflated(-α + r) - 1 / (1 + α - r)
```
and can split `T22` as
```
zeta_deflated(-α + r) * (x^r - x^(1 + α)) / 2log(x) + (x^r - x^(1 + α)) / (-α + r - 1) / 2log(x)
= T221 + T222
```
The term `T221` can be enclosed directly. For `T222` we rewrite it as
```
(x^r - x^(1 + α)) / (2(r - (1 + α)) * log(x))
```
We split this into two cases, when `r >= 1 + α` and when `r < 1 + α`.
In the first case we factor out `x^(1 + α)`, giving us
```
x^(1 + α) / 2 * (x^(r - (1 + α)) - 1) / ((r - (1 + α)) * log(x))
```
If we let `t = (r - (1 + α)) * log(x)` we can write
```
(x^(r - (1 + α)) - 1) / ((r - (1 + α)) * log(x)) = (exp(t) - t) / t
```
Here we notice that since `r - (1 + α) >= 0` and `log(x) < 0` we have
`t <= 0`. Furthermore the function `(exp(t) - 1) / t` is zero at `t =
-Inf`, one at `t = 0` and increasing in `t`. It is therefore enough to
determine the endpoints of `t` and from there we can compute an
enclosure.

The second case, when `r < 1 + α`, doesn't always occur, it depends on
the value of `j` and `u0.ϵ`. We first check if `1 + α - r` contain any
positive numbers, if that is the case we proceed similar to for the
first case. We factor out `x^r`, giving us
```
x^r / 2 * (1 - x^(1 + α - r)) / ((r - (1 + α)) * log(x)) =
    x^r / 2 * (x^(1 + α - r) - 1) / ((1 + α - r) * log(x))
```
Taking `t = 1 + α - r` we in this case have `1 + α -r > 0` and `log(x)
< 0` so `t < 0`. This allows us to compute an enclosure similar to how
we did it in the first case.

## Handling `T3`: the remaining terms
Once the terms `P` and `Q` have been taken out from the expansion it
is possible to enclose the remaining ones directly. However the
expansion for the first Clausen function in the tail has very large
cancellations between the first two terms in the expansion and they
are therefore handled separately. The rest of the terms we enclose
directly by using [`eval_expansion`](@ref). If `u0.ϵ` is small this
gives good enough enclosures so that we don't have to handle it in any
more sophisticated way.

Recall that all terms are supposed to be divided by `(gamma(1 + α) *
log(x) * x^(1 - α) * (1 - x^p0))`. The `x^(1 - α)` factor will be
cancelled explicitly. For the remaining part we compute an enclosure.
We are therefore interested in computing an enclosure of
```
inv(log(x)) * inv(gamma(1 + α) * (1 - x^p0))
```
The `inv(log(x))` is easily handled using monotonicity in `x` and the
other factor is the same as in the above section.
"""
function F0(
    u0::BHKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 3,
    ϵ::Arb = Arb(0.5),
    skip_singular_j_until = 100,
    alpha_interval = :full,
)
    @assert ϵ < 1

    # Interval for α
    α = Arb((-1, -1 + u0.ϵ))
    # Interval for α + 1
    αp1 = Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))

    # This method assumes that the weight is x^(1 - u0.γ * (1 + α)) *
    # log(u0.c * inv(x)). As an extra precaution we check this.
    let x = Arb(0.5), α = Arb((-1, -1 + u0.ϵ))
        @assert Arblib.overlaps(u0.w(x), x^(1 - u0.γ * (αp1)) * log(u0.c + inv(x)))
    end

    # Function for bounding gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x)
    f1 = inv_u0_bound(u0; M, ϵ)

    # Function for enclosing log(inv(x)) / log(u0.c + inv(x))
    f2 = x -> if iszero(x)
        one(x)
    elseif Arblib.contains_zero(x)
        lower = let xᵤ = ubound(Arb, x)
            log(inv(xᵤ)) / log(u0.c + inv(xᵤ))
        end
        upper = one(x)
        Arb((lower, upper))
    else
        log(inv(x)) / log(u0.c + inv(x))
    end

    # Compute the expansion of D(u0), skipping the Clausen term in the
    # tail corresponding to j = 1 and also remove the two leading
    # term, these three terms are handled separately.
    Du0_expansion =
        D(u0, AsymptoticExpansion(); M, skip_singular_j_until, alpha_interval)(ϵ)
    delete!(Du0_expansion, (2, 0, 0, 0, 0, 0, 0))
    delete!(Du0_expansion, (0, 1, 0, 0, 0, 0, 0))

    # Divide the expansion of D(u0) by x^(1 - α)
    Du0_expansion_div_x_onemα = empty(Du0_expansion)
    for ((p, q, i, j, k, l, m), y) in Du0_expansion
        Du0_expansion_div_x_onemα[(p, q, i + 1, j, k, l, m - 1)] = y
    end

    c(a) = gamma(a) * cospi(a / 2)

    # Compute enclosures of several values depending only on α, many
    # of them with removable singularities

    # Use this for computing tighter enclosures
    extra_degree = 2

    # Enclosure of rgamma(α) / (1 + α)
    rgamma_α_div_αp1 = fx_div_x(s -> rgamma(s - 1), αp1; extra_degree)
    # Enclosure of rgamma(2α) / (1 + α)
    rgamma_2α_div_αp1 = fx_div_x(s -> rgamma(2(s - 1)), αp1; extra_degree)
    # Enclosure of rgamma(1 + α) / (1 + α)
    rgamma_1pα_div_αp1 = fx_div_x(s -> rgamma(s), αp1; extra_degree)
    # Enclosure of rgamma(α - p0) / (1 + α)^2
    rgamma_αmp0_div_αp12 = fx_div_x(s -> rgamma(-1 - s^2 / 2), αp1, 2; extra_degree)
    # Enclosure of rgamma(2α - p0) / (1 + α)
    rgamma_2αmp0_div_αp1 = fx_div_x(s -> rgamma(s - 2 - s^2 / 2), αp1; extra_degree)

    # Enclosure of cospi(α / 2) / (1 + α)
    cos_αdiv2_div_αp1 = fx_div_x(s -> cospi((s - 1) / 2), αp1; extra_degree)
    # Enclosure of cospi((α - p0) / 2) / (1 + α)^2
    cos_αmp0div2_div_αp12 = fx_div_x(s -> cospi((-1 - s^2 / 2) / 2), αp1, 2; extra_degree)

    # Enclosure of ((1 + α) * c(2α - p0) - 2(1 + α) * c(2α) * c(α - p0) / c(α)) / (1 + α)^2
    T113 = fx_div_x(αp1, 2; extra_degree, force = true) do s
        if Arblib.contains_zero(s[0])
            # Expansion of c(2α - p0) * (1 + α)
            c_2αmp0_mul_α = inv(
                fx_div_x(
                    t -> rgamma(t - 2 - t^2 / 2) / cospi((t - 2 - t^2 / 2) / 2),
                    s;
                    extra_degree,
                ),
            )

            # Enclosure of c(2α) * (1 + α)
            c_2α_mul_α =
                inv(fx_div_x(t -> rgamma(2(t - 1)) / cospi(t - 1), s; extra_degree))

            # Enclosure of c(α - p0)
            c_αmp0 =
                fx_div_x(t -> cospi((-1 - t^2 / 2) / 2), s, 2; extra_degree) /
                fx_div_x(t -> rgamma(-1 - t^2 / 2), s, 2; extra_degree)

            # Enclosure of c(α)
            c_α =
                fx_div_x(t -> cospi((t - 1) / 2), s; extra_degree) /
                fx_div_x(t -> rgamma(t - 1), s; extra_degree)

            c_2αmp0_mul_α - 2c_2α_mul_α * c_αmp0 / c_α
        else
            # Evaluate this at higher precision since it is close to
            # the removable singularity
            let α = setprecision(s - 1, 2precision(s)), p0 = (1 + α) + (1 + α)^2 / 2
                setprecision((c(2α - p0) - 2c(2α) * c(α - p0) / c(α)) / (1 + α), precision(s))
            end
        end
    end

    # Enclosure of a0 / gamma(1 + α)
    a0_div_gamma_1pα =
        rgamma_α_div_αp1^3 / (α / 2 * cospi(α) * cos_αdiv2_div_αp1^2 * rgamma_2α_div_αp1)

    # Enclosure of T1231 * log(x) = (v(α) - w(α)) / (1 + α)
    T1231_mul_logx = fx_div_x(αp1, 2; extra_degree, force = true) do s
        if Arblib.contains_zero(s[0])
            # Enclosure of s * zeta(1 - 2s)
            zeta_1m2s_mul_α = s * zeta_deflated(1 - 2s, one(Arb)) - 1 // 2

            # Enclosure of s * zeta(1 - s + s^2 / 2)
            zeta_1msps2div2_mul_α =
                s * zeta_deflated(1 - s + s^2 / 2, one(Arb)) - 1 / (1 - s / 2)

            v_mul_α = zeta_1m2s_mul_α - zeta_1msps2div2_mul_α

            # Enclosure of c(2α) * (1 + α)
            c_2α_mul_α =
                inv(fx_div_x(t -> rgamma(2(t - 1)) / cospi(t - 1), s; extra_degree))

            # Enclosure of c(α - p0)
            c_αmp0 =
                fx_div_x(t -> cospi((-1 - t^2 / 2) / 2), s, 2; extra_degree) /
                fx_div_x(t -> rgamma(-1 - t^2 / 2), s, 2; extra_degree)

            # Enclosure of c(α)
            c_α =
                fx_div_x(t -> cospi((t - 1) / 2), s; extra_degree) /
                fx_div_x(t -> rgamma(t - 1), s; extra_degree)

            w_mul_α = -2c_2α_mul_α * c_αmp0^2 / c_α^2

            v_mul_α - w_mul_α
        else
            # Evaluate this at higher precision since it is close to
            # the removable singularity
            let α = setprecision(s - 1, 2precision(s)), p0 = (1 + α) + (1 + α)^2 / 2
                v = zeta(-2α - 1) - zeta(-2α + p0 - 1)
                w = -2c(2α) * c(α - p0)^2 / c(α)^2
                setprecision(v - w, precision(s))
            end
        end
    end

    # Enclosure of w(α) * (1 + α)
    w_mul_α =
        -2(cospi(α) / rgamma_2α_div_αp1) *
        (cos_αmp0div2_div_αp12 / rgamma_αmp0_div_αp12)^2 /
        (cos_αdiv2_div_αp1 / rgamma_α_div_αp1)^2

    # α-factor of T21 for j = 1:skip_singular_j_until
    # IMPROVE: Compute tighter enclosure when α + 1 - r is close to
    # zero
    T21_α = map(1:skip_singular_j_until) do j
        let r = -u0.v0.α + j * u0.v0.p0 - 1
            ArbExtras.enclosure_series(α, degree = 4) do α
                let res = -zeta_deflated(-α + r, one(Arb)) / 2
                    if (α isa Arb && Arblib.contains_zero(α + 1 - r)) ||
                       (α isa ArbSeries && Arblib.contains_zero(α[0] + 1 - r))
                        res += fx_div_x(α + 1 - r; extra_degree) do s
                            gamma(s + 1) / ((s - 1) * (s - 2)) * cospi((2 - s) / 2) + 1 // 2
                        end
                    else
                        res +=
                            (
                                gamma(α + 2 - r) / ((α - r) * (α - 1 - r)) *
                                cospi((1 - α + r) / 2) + 1 // 2
                            ) / (α + 1 - r)
                    end
                    res
                end
            end
        end
    end

    # Enclosure of zeta_deflated(-α + r) for j = 1:skip_singular_j_until
    zeta_deflated_mαpr = map(1:skip_singular_j_until) do j
        let r = -u0.v0.α + j * u0.v0.p0 - 1
            ArbExtras.enclosure_series(α -> zeta_deflated(-α + r, one(r)), α, degree = 8)
        end
    end

    return x::Arb -> begin
        @assert x <= ϵ

        # Enclosure of inv(log(x))
        invlogx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            xᵤ = ubound(Arb, x)
            Arb((inv(log(xᵤ)), 0))
        else
            inv(log(x))
        end

        # Enclosure of inv(gamma(1 + α) * (1 - x^p0))
        invgamma1mxp0 = if iszero(x)
            rgamma(1 + α)
        elseif Arblib.contains_zero(x)
            lower = zero(x)
            upper = let xᵤ = ubound(Arb, x)
                # Enclosure of (1 - xᵤ^p0) / (1 + α)
                onemxp0_div_αp1 =
                    fx_div_x(s -> (1 - xᵤ^(s + s^2 / 2)), αp1, extra_degree = 2)
                rgamma_1pα_div_αp1 / onemxp0_div_αp1
            end
            Arb((lower, upper))
        else
            # Enclosure of (1 - x^p0) / (1 + α)
            onemxp0_div_αp1 = fx_div_x(s -> (1 - x^(s + s^2 / 2)), αp1, extra_degree = 2)
            rgamma_1pα_div_αp1 / onemxp0_div_αp1
        end

        # Enclosure for the terms P and Q
        T1 = let
            # Compute an enclosure of T11

            # Enclosure of T111 = x^(-α + p0 - 1) / log(x) using that
            # -α + p0 - 1 = (1 + α)^2 / 2
            T111 = abspow(x, Arblib.nonnegative_part!(zero(x), (αp1)^2 / 2)) * invlogx

            # Enclosure of T112 = (1 + α) / (1 - x^p0)
            T112 = if iszero(x)
                αp1
            elseif Arblib.contains_zero(x)
                lower = 1 + α
                upper = let xᵤ = ubound(Arb, x)
                    # Enclosure of inv((1 - xᵤ^p0) / (1 + α))
                    inv(fx_div_x(s -> (1 - xᵤ^(s + s^2 / 2)), αp1, extra_degree = 2))
                end
                Arb((lower, upper))
            else
                # Enclosure of inv((1 - x^p0) / (1 + α))
                inv(fx_div_x(s -> (1 - x^(s + s^2 / 2)), αp1, extra_degree = 2))
            end


            T11 = T111 * T112 * T113

            # Compute an enclosure of T12

            # Enclosure of T121 = x^(1 + α) / 2
            T121 = abspow(x, Arblib.nonnegative_part!(zero(x), 1 + α)) / 2

            # Enclosure of T122, which is the same as T112
            T122 = T112

            T123 = let
                T1231 = T1231_mul_logx * invlogx

                T1232 = let
                    # Lower and upper bounds of
                    # t = (1 + α)^2 * log(x)
                    tₗ = ubound(Arb, αp1^2) * log(abs_lbound(Arb, x))
                    tᵤ = abs_lbound(Arb, αp1^2) * log(ubound(Arb, x))
                    # Lower and upper bounds of (1 - exp(t)) / t
                    # Using that t <= 0 to handle the singular cases
                    lower = Arblib.isnegative(tᵤ) ? (1 - exp(tᵤ)) / tᵤ : -one(tᵤ)
                    upper = isfinite(tₗ) ? (1 - exp(tₗ)) / tₗ : zero(tₗ)

                    Arb((lower, upper)) * w_mul_α
                end

                T1231 + T1232
            end

            T12 = T121 * T122 * T123

            a0_div_gamma_1pα * (T11 + T12)
        end

        # Enclosure of the two singular terms in the expansion of
        # clausenc(x, 1 - α - u0.v0.α + j * u0.v0.p0) for
        # j = 1:skip_singular_j_until
        T2s = map(1:skip_singular_j_until) do j
            let r = -u0.v0.α + j * u0.v0.p0 - 1
                # Enclosure of
                # (gamma(α - 1 - r) * cospi((1 - α + r) / 2) - zeta(-α + r) / 2) * x^r / log(x)
                T21 = T21_α[j] * abspow(x, r) * invlogx

                # Enclosure of zeta_deflated(-α + r) * (x^r - x^(1 + α)) / 2log(x)
                T221 =
                    zeta_deflated_mαpr[j] * (abspow(x, r) - abspow(x, αp1)) * invlogx / 2

                # Enclosure of (x^r - x^(1 + α)) / (-α + r - 1) / 2log(x)
                T222 = let
                    # Handle the case r >= 1 + α

                    # Lower and upper bounds of
                    # t = (r - (1 + α)) * log(x)
                    tₗ = ubound(Arb, r - αp1) * log(abs_lbound(Arb, x))
                    tᵤ = abs_lbound(Arb, r - αp1) * log(ubound(Arb, x))
                    # Lower and upper bounds of (exp(t) - 1) / t
                    # Using that t <= 0 to handle the singular cases
                    lower = isfinite(tₗ) ? (exp(tₗ) - 1) / tₗ : zero(tₗ)
                    upper = Arblib.isnegative(tᵤ) ? (exp(tᵤ) - 1) / tᵤ : one(tᵤ)

                    T222 = abspow(x, αp1) / 2 * Arb((lower, upper))

                    # Handle the case r < 1 + α if it occurs
                    if Arblib.contains_positive(1 + α - r)
                        # Lower and upper bounds of
                        # t = (1 + α - r) * log(x)
                        tₗ = ubound(Arb, αp1 - r) * log(abs_lbound(Arb, x))
                        tᵤ = abs_lbound(Arb, αp1 - r) * log(ubound(Arb, x))
                        # Lower and upper bounds of (exp(t) - 1) / t
                        # Using that t <= 0 to handle the singular cases
                        lower = isfinite(tₗ) ? (exp(tₗ) - 1) / tₗ : zero(tₗ)
                        upper = Arblib.isnegative(tᵤ) ? (exp(tᵤ) - 1) / tᵤ : one(tᵤ)

                        # Add enclosure from this case to T222
                        T222 = union(T222, abspow(x, r) / 2 * Arb((lower, upper)))
                    end

                    T222
                end

                term = T21 + T221 + T222

                term *= invgamma1mxp0

                -u0.v0.a[j] * term
            end
        end

        T2 = sum(T2s, init = zero(x))

        # Enclosure of the remaining terms in the expansion
        T3 = eval_expansion(u0, Du0_expansion_div_x_onemα, x) * invlogx * invgamma1mxp0

        # (u0(x)^2 / 2 + Hu0x) / (log(x) * gamma(1 + α) * x^(1 - α) * (1 - x^p0))
        F = T1 + T2 + T3

        return f1(x) * f2(x) * F
    end
end

"""
    inv_u0_bound(u0::BHKdVAnsatz{Arb})

Return a function `f` such that `f(x)` gives an upper bound of
```
gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x)
```

It assumes that `x` is non-negative, any negative parts of `x` are
ignored.

For non-zero `x` it computes an enclosure of the value. It splits it as
```
gamma(1 + α) * (1 - x^p0) / (u0(x) / x^(-α))
```
and computes an enclosure of
```
gamma(1 + α) * (1 - x^p0) = gamma(2 + α) * (1 - x^p0) / (1 + α)
```
by handling the removable singularity at `α = -1` using
[`fx_div_x`](@ref). For `u0(x) / x^(-α)` it computes an enclosure
using the asymptotic expansion of `u0`.

We now describe how to compute an upper bound when `x` overlaps with
zero.

We start by computing the asymptotic expansion of `u0`. We then split
the expansion into the leading term and a tail. In practice both the
leading term and the tail are positive and the tail is much smaller
than the leading term. Since we are only interested in an upper bound
we can thus just skip the tail. For this to be valid we must ensure
that both the leading term and the tail are positive. For the tail
this is checked using [`expansion_ispositive`](@ref). For the leading
term this is done by checking that the end result is positive (it is
easily seen that the numerator is positive in our case)

What remains is to handle the leading term, it is given by
```
a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) * x^-α
```
Since we are interested in `gamma(1 + α) * x^-α * (1 - x^p0)` divided
by this we get
```
gamma(1 + α) * (1 - x^p0) /
    (a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0))
```
We split this further into the two factors
```
F1 = gamma(1 + α) / a0
F2 = (1 - x^p0) /
    (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0)
```

# `F1`
From the definition of `a0`
```
a0 = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```
and using `gamma(1 + α) = α * gamma(α)` we get
```
F1 = α * gamma(α)^3 * cospi(α / 2)^2 / (2gamma(2α) * cospi(α))
```
This has a removable singularity at `α = -1`. By rewriting it as
```
F1 = α / 2 * cospi(α)
    * inv(rgamma(α) / (1 + α))^3
    * (cospi(α / 2) / (1 + α))^2
    * rgamma(2α) / (1 + α)
```
, where `rgamma` is the reciprocal gamma function given by `rgamma(s)
= inv(gamma(s))`, we can use [`fx_div_x`](@ref) to compute enclosures
of the individual factors.

# `F2`
For `F2` we let `c(a) = gamma(a) * cospi(a / 2)` and then factor it
out, giving us
```
inv(c(α)) * (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
```

Similarly to `F1`, `c(α)` has a removable singularity at `α = -1`. We
can compute an enclosure using the same tools.

For the remaining part we compute an enclosure of the derivative of
`c(α - p0) / c(α)` with respect to `α` and check that this is
negative. This means that `c(α - p0) / c(α)` is decreasing in `α` and
the maximum value is hence attained at `α = -1` where it is `1`. This
means that
```
(1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
```
is upper bounded by `(1 - x^p0) / (1 - x^p0) = 1`. Since `c(α - p0) /
c(α) <= 1` we also get that
```
(1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
```
is non-increasing in `x`. A lower bound can thus be computed by
considering `xᵤ = ubound(x)`. At `xᵤ > 0` we can handle the removable
singularity in `α`
"""
function inv_u0_bound(u0::BHKdVAnsatz{Arb}; M::Integer = 3, ϵ::Arb = Arb(0.5))
    # Interval for α
    α = Arb((-1, -1 + u0.ϵ))
    # Interval for α + 1
    αp1 = Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))
    extra_degree = 2

    # Enclosure of rgamma(α) / (α + 1)
    rgamma_α_div_αp1 = fx_div_x(s -> rgamma(s - 1), αp1; extra_degree)
    # Enclosure of rgamma(2α) / (α + 1)
    rgamma_2α_div_αp1 = fx_div_x(s -> rgamma(2(s - 1)), αp1; extra_degree)
    # Enclosure of cospi(α / 2) / (α + 1)
    cos_αdiv2_div_αp1 = fx_div_x(s -> cospi((s - 1) / 2), αp1; extra_degree)

    # Enclosure of F1
    F1 =
        α / 2 * cospi(α) * inv(rgamma_α_div_αp1)^3 * cos_αdiv2_div_αp1^2 * rgamma_2α_div_αp1

    c(a) =
        if (a isa Arb && Arblib.contains_zero(a + 1)) ||
           (a isa ArbSeries && Arblib.contains_zero(a[0] + 1))
            fx_div_x(s -> cospi((s - 1) / 2), a + 1; extra_degree) /
            fx_div_x(s -> rgamma(s - 1), a + 1; extra_degree)
        else
            gamma(a) * cospi(a / 2)
        end

    inv_c_α = inv(c(α))

    # Prove that c(α - p0) / c(α) is decreasing in α
    cαmp0_div_cα = let α = ArbSeries((α, 1)), p0 = (1 + α) + (1 + α)^2 / 2
        c(α - p0) / c(α)
    end

    Arblib.isnegative(cαmp0_div_cα[1]) ||
        error("c(α - p0) / c(α) could not be proved to be decreasing")

    # Compute the expansion of u0
    u0_expansion = u0(ϵ, AsymptoticExpansion(); M)

    # Temporarily remove the leading term and check that the rest of
    # the expansion is positive, so that we can remove it from the
    # denominator and still get an upper bound.
    leading_term = u0_expansion[(1, 0, 0, 0, 0, 0, 0)]
    delete!(u0_expansion, (1, 0, 0, 0, 0, 0, 0))
    expansion_ispositive(u0, u0_expansion, ϵ) ||
        error("expansion of u0 could not prove to be positive")
    u0_expansion[(1, 0, 0, 0, 0, 0, 0)] = leading_term

    # Divide the expansion of u0 by x^-α
    u0_expansion_div_x_mα = empty(u0_expansion)
    for ((p, q, i, j, k, l, m), y) in u0_expansion
        u0_expansion_div_x_mα[(p, q, i + 1, j, k, l, m)] = y
    end

    return x::Arb -> begin
        x < 1 || throw(DomainError(x, "need 0 < x < 1"))
        x <= ϵ || throw(DomainError(x, "need x <= ϵ = $ϵ"))

        if Arblib.contains_zero(x)
            # Enclose F2
            F2 = if iszero(x)
                # We have (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0) = 1
                inv_c_α
            else
                # Compute lower and upper bounds of
                # (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
                # IMPROVE: Compute tighter enclosures for very small x
                lower = let xᵤ = ubound(Arb, x)
                    # Enclosure of (1 - x^p0) / (α + 1)
                    numerator_div_αp1 =
                        fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), αp1; extra_degree)
                    # Enclosure of (1 - c(α - p0) / c(α) * x^p0) / (α + 1)
                    denominator_div_αp1 = fx_div_x(
                        s -> 1 - c(-1 - s^2 / 2) / c(s - 1) * xᵤ^(s + s^2 / 2),
                        αp1;
                        extra_degree,
                        force = true,
                    )
                    numerator_div_αp1 / denominator_div_αp1
                end
                upper = one(Arb)

                inv_c_α * Arb((lower, upper))
            end

            F = F1 * F2

            Arblib.ispositive(F) ||
                error("leading term of u0 is not positive, this should not happen")
        else
            # Enclosure of gamma(1 + α) * (1 - x^p0) = gamma(2 + α) * (1 - x^p0) / (1 + α)
            numerator =
                gamma(2 + α) * fx_div_x(s -> (1 - x^(s + s^2 / 2)), αp1; extra_degree)

            F = numerator / eval_expansion(u0, u0_expansion_div_x_mα, x)
        end

        return F
    end
end
