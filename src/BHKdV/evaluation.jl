"""
    eval_expansion(u0::BHKdVAnsatz, expansion, x)

Evaluate the given expansion. It requires that `x < 1` and that `x` is
non-negative, any negative parts of `x` are ignored.

The terms are stored as `((p, q, i, j, k, l, m), y)`. The parameters `(i, j,
k l, m)` correspond to the term
```
y * x^(i * α + j * p0 - k*u0.v0.v0.α + l*u0.v0.v0.p0 + m)
```
where `α ∈ (-1, -1 + u0.ϵ]` and `p0 = 1 + α + (1 + α)^2 / 2`.

The parameter `p` corresponds to multiplication by the factor
```
a0 * (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0) * x^-α
```
to the power `p`, which is the part of the expansion for the main term
which is not even powers of `x`. The parameter `q` corresponds to
multiplication by the factor
```
-a0 * (
    gamma(2α) * sinpi((1 - 2α) / 2) * x^(-2α) -
    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2) * x^(-2α + p0) +
    (-zeta(1 - 2α - 2) / 2 + zeta(1 - 2α + p0 - 2) / 2) * x^2
)
```
to the power `q`, which is the part of the expansion for `H` applied
to the main term which is not on the form `x^2m` for `m >= 2`.

# Handling of `p` and `q`
The method currently does not support `q != 0` and only supports `p ==
0` and `p == 1`. These cases are handled specially in [`F0`](@ref) and
for that reason we don't bother implementing them here.

As `α` goes to `-1` we have that
```
a0 * (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0)
```
converges to
```
(1 - γ - log(x)) / π
```
- **TODO:** Prove this and compute error bounds.

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
    # k*u0.v0.v0.α + l*u0.v0.v0.p0 + m in a way that also accounts for
    # the dependence between α and p0
    exponent = zero(α)

    # Precomputed values and buffers for computing the exponent
    lower = zero(exponent)
    upper = zero(exponent)
    # Upper bounds of α and p0
    α_upper = -1 + u0.ϵ
    p0_upper = 1 + α_upper + (1 + α_upper)^2 / 2
    # Enclosure of p0 - α
    p0mα = 1 + (1 + α)^2 / 2
    _exponent!(exponent, i, j, k, l, m) = begin
        # Compute the part i * α + j * p0 + m
        # Note that i * α + j * p0 = 3j / 2 + (i + 2j) * α + α^2 / 2
        # is increasing in α if i + 2j is non-negative.
        if i + 2j >= 0
            # It is increasing in α, evaluated at endpoints

            # Lower bound at α = -1 can be done with rational numbers
            let α = -1
                Arblib.set!(lower, i * α + j * (1 + α + (1 + α)^2 // 2) + m)
            end

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

        # Add - k*u0.v0.v0.α + l*u0.v0.v0.p0
        Arblib.submul!(exponent, u0.v0.v0.α, k)
        Arblib.addmul!(exponent, u0.v0.v0.p0, l)

        return exponent
    end

    # Irrationals used
    π = Arb(Irrational{:π}())
    γ = Arb(Irrational{:γ}())

    # Precompute (1 - γ - log(x)) / π
    onemγmlogxdivπ = (1 - γ - log(x)) / π

    res = zero(x)

    for ((p, q, i, j, k, l, m), y) in expansion
        if !iszero(y)
            if iszero(p)
                _exponent!(exponent, i, j, k, l, m)

                term = abspow(x, exponent)
            elseif isone(p)
                # Add -α to the exponent coming from the p factor
                _exponent!(exponent, i - 1, j, k, l, m)

                # FIXME: Add error bounds for this term. Here we just
                # use the limiting value

                # Compute an enclosure using monotonicity
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

                    upper = let x = ubound(Arb, x)
                        (1 - γ - log(x)) / π * abspow(x, exponent)
                    end

                    term = Arb((lower, upper))
                else
                    term = onemγmlogxdivπ * abspow(x, exponent)
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
                            gamma(α) * sinpi((1 - α) / 2) -
                            gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0
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
                            gamma(2α) * sinpi((1 - 2α) / 2) * x^(-2α) -
                            gamma(2α - p0) * sinpi((1 - 2α + p0) / 2) * x^(-2α + p0) +
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
y * x^^(-u0-v0.v0.α + l*u0.v0.v0.p0)
```
Let `y₁` be the coefficient for `l = 1` and `S` be the sum of the
coefficients with `l > 1`. It checks that `S` is negative and that `y₁
+ S > 0`, ensuring us that this sum is positive for all `0 < x < 1`
and a lower bound is given by `(y₁ + S) * x^(-u0.v0.v0.α +
u0.v0.v0.p0)` for all `0 < x < 1`.

The next step is to prove that the sum of the terms with `m > 0` are
smaller than (y₁ + S) * x^(-u0.v0.v0.α + u0.v0.v0.p0)`. This is done
by noting that it's enough to check it for `x = ϵ`

**TODO:** Check that this is correct.
"""
function expansion_ispositive(
    u0::BHKdVAnsatz{Arb},
    expansion::AbstractDict{NTuple{7,Int},Arb},
    ϵ::Arb,
)
    @assert 0 < ϵ < 1
    @assert 0 < -u0.v0.v0.α + u0.v0.v0.p0 < 2

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

    a = (y₁ + S) * ϵ^(-u0.v0.v0.α + u0.v0.v0.p0)
    b = eval_expansion(u0, expansion_2, ϵ)

    a > b || return false

    return true
end


"""
    (u0::BHKdVAnsatz)(x, ::Ball)

Evaluate the ansatz `u0` at the point `x` using direct ball arithmetic
(not an asymptotic approach).

The tail term is evaluated directly.

To evaluate the main term, given by
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
we make use of the fact that this converges to
```
2 / π^2 * clausencmzeta(x, 2, 1)
```
, which is the main term for `BHAnsatz`, as `α -> -1`. We therefore
evaluate this function and add error bounds for
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) - 2 / π^2 * clausencmzeta(x, 2, 1)
```
valid for the entire range `α ∈ (-1, -1 + u0.ϵ]`.

For now we compute the value for `α = -1 + u0.ϵ` and take the union of
this result and the one computed with the limiting expression. This
works in practice since we have a monotone convergence.

This approach also works for `ArbSeries`, though it is currently less
clear if we have the same monotone convergence, probably we do.

- **TODO:** Compute rigorous error bounds. Possibly by proving the
    monotonicity of the error.
"""
function (u0::BHKdVAnsatz{Arb})(x::Union{Arb,ArbSeries}, ::Ball)
    # Main term

    # Approximation
    res = 2 / Arb(π)^2 * clausencmzeta(x, 2, 1)

    # Add error bounds
    # TODO: Implement rigorous bounds
    @warn "Non-rigorous bounds implemented for main term" maxlog = 1
    error = let α = -1 + u0.ϵ
        a0 = finda0(α)
        p0 = 1 + α + (1 + α)^2 / 2
        res2 = a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))

        if x isa Arb
            res = union(res, res2)
        elseif x isa ArbSeries
            coefficients = union.(Arblib.coeffs(res), Arblib.coeffs(res2))
            res = ArbSeries(coefficients)
        end
    end

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.v0.N0
        s = 1 - u0.v0.v0.α + j * u0.v0.v0.p0
        res += u0.v0.v0.a[j] * clausencmzeta(x, s)
    end

    # Fourier terms
    for n = 1:u0.v0.N
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
evaluation of the expansion gives an enclosure of the result.

For the tail term the expansions are easily computed exactly like for
`BHAnsatz`. For the main term we have to be a bit more careful.

For the main term the coefficients in front of `x^2m` is given by
```
a0 * (-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / factorial(2m)
```
which in the limit becomes `∞ * 0`. It converges to
```
u0.v0.a0 * (-1)^m * dzeta(2 - 2m) / factorial(2m)
```
, where `u0.v0.a0 = 2 / π^2` , which is the coefficient in front of
`x^2m` for the main term of `u0.v0`. We therefore compute the
coefficients by using this expression and then adding bounding the
error for it.

For now we bound the error by computing the value at `α = -1 + u0.ϵ`
and take the union of the result with the one computed with the
limiting expression. This works in practice since we have a monotone
convergence.
- **TODO:** Compute rigorous bounds for the coefficients. Possibly by
  proving the monotonicity of the error.
- **TODO:** Compute rigorous bounds for the error term in the
  expansion. Possibly by proving the monotonicity of it.

The only remaining part of the expansion of the main term is
```
a0 * (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0) * x^-α
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
    # TODO: Implement rigorous bounds
    @warn "Non-rigorous bounds implemented for x^2m coefficients" maxlog = 1
    let α = -1 + u0.ϵ, a0 = finda0(α), p0 = 1 + α + (1 + α)^2 / 2
        for m = 1:M-1
            coefficient = u0.v0.a0 * (-1)^m * dzeta(Arb(2 - 2m)) / factorial(2m)

            # Add error bounds
            coefficient_2 =
                a0 * (-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / factorial(2m)

            coefficient = union(coefficient, coefficient_2)

            res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
        end
    end

    # Error term for main term
    @warn "Non-rigorous error term implemented for main term" maxlog = 1
    Arblib.add_error!(
        res[(0, 0, 0, 0, 0, 0, 2M)],
        2abs(dzeta(Arb(2 - 2M)) / factorial(2M)) * u0.v0.a0,
    )

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.v0.N0
        s = 1 - u0.v0.v0.α + j * u0.v0.v0.p0
        C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)
        res[(0, 0, 0, 0, 1, j, 0)] = C * u0.v0.v0.a[j]
        for m = 1:M-1
            res[(0, 0, 0, 0, 0, 0, 2m)] += p[2m] * u0.v0.v0.a[j]
        end
        res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.v0.a[j]
    end

    # Fourier terms
    if !iszero(u0.v0.N)
        for m = 1:M-1
            res[(0, 0, 0, 0, 0, 0, 2m)] +=
                (-1)^m * sum(Arb(n)^(2m) * u0.v0.b[n] for n = 1:u0.v0.N) /
                factorial(Arb(2m))
        end
        Arblib.add_error!(
            res[(0, 0, 0, 0, 0, 0, 2M)],
            sum(Arb(n)^(2M) * abs(u0.v0.b[n]) for n = 1:u0.v0.N) / factorial(Arb(2M)),
        )
    end

    return res
end

"""
    H(u0::BHKdVAnsatz, ::Ball)

Returns a function such that `H(u0, Ball())(x)` gives an enclosure of
`H^-α[u0](x)` for all values of `α ∈ (-1, -1 + u0.ϵ]`. It uses direct
ball arithmetic (not an asymptotic approach).

The transform of the main term is given by
```
-a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
```
we make use of the fact that this converges to
```
-2 / π^2 * clausencmzeta(x, 3, 1)
```
, which is the main term for `BHAnsatz`, as `α -> -1`. We therefore
evaluate this function and bound the error.

For now we compute the value for `α = -1 + u0.ϵ` and take the union of
this result and the one computed with the limiting expression. This
works in practice since we have a monotone convergence.

This approach also works for `ArbSeries`, though it is currently less
clear if we have the same monotone convergence, probably we do.

- **TODO:** Compute rigorous error bounds. Possibly by proving the
    monotonicity of the error.

For the tail term we need to make sure that we correctly handle the
fact that the transform depends on the value of `α`.

For the Fourier terms we do this directly, the transformation takes
`cos(n * x)` to `-n^α * cos(n * x)` and in this case we just let `α`
be a ball containing `(-1, -1 + u0.ϵ]`.

For the Clausen functions we have to be a bit more careful. The
transformation takes `clausenc(x, s)` to `-clausenc(x, s - α)` but
putting `α` as a ball doesn't give good enclosures. It gives a large
overestimations, in particular when `s - α` is close to an integer,
and fails when it overlaps with an integer.

For now we take the same approach as for the main term, compute with
`α = -1 + u0.ϵ` and take the union.

- **FIXME:** The monotonicity holds when computing with `Arb` but not
    when computing with `ArbSeries` for higher order terms. At least
    this seems to be the case. So this doesn't give a rigorous bound
    but probably does give good bounds in practice.
- **TODO:** Compute rigorous error bounds. Proving the monotonicity of
    the error would be great but might be hard in practice since one
    of the coefficients of `u0.v0.v0` have a sign that differs.
"""
function H(u0::BHKdVAnsatz{Arb}, ::Ball)
    # Terms used when computing error bounds
    α = -1 + u0.ϵ
    a0 = finda0(α)
    p0 = 1 + α + (1 + α)^2 / 2

    return x::Union{Arb,ArbSeries} -> begin
        # Main term

        # Approximation
        res = -2 / Arb(π)^2 * clausencmzeta(x, 3, 1)

        # Add error bounds
        @warn "No error bounds for main term" maxlog = 1
        # TODO: Implement rigorous bounds
        res2 = -a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
        if x isa Arb
            res = union(res, res2)
        elseif x isa ArbSeries
            coefficients = union.(Arblib.coeffs(res), Arblib.coeffs(res2))
            res = ArbSeries(coefficients)
        end

        # Tail term

        # Clausen terms
        clausen_term = zero(x)

        # Approximation
        for j = 1:u0.v0.v0.N0
            term = clausencmzeta(x, 2 - u0.v0.v0.α + j * u0.v0.v0.p0)
            clausen_term -= u0.v0.v0.a[j] * term
        end

        # Add error bounds
        @warn "Non-rigorous bounds implemented for Clausen terms" maxlog = 1
        # TODO: Implement rigorous bounds
        clausen_term2 = zero(x)
        for j = 1:u0.v0.v0.N0
            term = clausencmzeta(x, 1 - α - u0.v0.v0.α + j * u0.v0.v0.p0)
            clausen_term2 -= u0.v0.v0.a[j] * term
        end
        if x isa Arb
            clausen_term = union(clausen_term, clausen_term2)
        elseif x isa ArbSeries
            coefficients = union.(Arblib.coeffs(clausen_term), Arblib.coeffs(clausen_term2))
            clausen_term = ArbSeries(coefficients)
        end

        res += clausen_term

        # Fourier terms
        let α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α
            for n = 1:u0.v0.N
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
    skip_j_one_singular = false,
    use_approx_p_and_q = false,
)
    f = H(u0, AsymptoticExpansion(); M, skip_j_one_singular)

    return x -> eval_expansion(u0, f(x), x; use_approx_p_and_q)
end

"""
    H(u0::BHKdVAnsatz, ::AsymptoticExpansion; M = 3, skip_j_one_singular = false)

Return a dictionary containing the terms in the asymptotic expansion
of `u0` which can then be evaluated with [`eval_expansion`](@ref).

The highest term, `x^2M`, is an error term is which makes sure that
evaluation of the expansion gives an enclosure of the result.

For the main term the coefficients in front of `x^2m` is given by
```
-a0 * (-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / factorial(2m)
```
which in the limit becomes `∞ * 0`. For `m >= 2` it converges to
```
-u0.v0.a0 * (-1)^m * dzeta(2 - 2m) / factorial(2m)
```
, where `u0.v0.a0 = 2 / π^2` , which is the coefficient in front of
`x^2m` for the main term of `u0.v0`. We therefore compute the
coefficients by using this expression and then bounding the error for
it. For `m = 1` we can't do this, instead we include it in the
remaining terms below.

For now we bound the error by computing the value at `α = -1 + u0.ϵ`
and take the union of the result with the one computed with the
limiting expression. This works in practice since we have a monotone
convergence.
- **TODO:** Compute rigorous bounds for the coefficients. Possibly by
  proving the monotonicity of the error.
- **TODO:** Compute rigorous bounds for the error term in the
  expansion. Possibly by proving the monotonicity of it.

The remaining part of the expansion of `H` applied to the main term is
```
-a0 * (
    gamma(2α) * sinpi((1 - 2α) / 2) * abs(x)^(-2α) -
    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2) * abs(x)^(-2α + p0)
    (zeta(1 - 2α - 2) / 2 - zeta(1 - 2α + p0 - 2) / 2) * abs(x)^2
)
```
which we don't evaluate at all yet. Instead we store it implicitly in
the expansion.

For both the Clausen terms and the Fourier terms we let `α` be a ball.
This gives good enclosures for the Fourier terms and decent enclosures
for the Clausen terms for `j >= 2`. For `j = 1` the parameter `s`
overlaps with 3 and the two leading terms in the expansion become
singular are need to be handled separately.

For now we only compute an approximation of the two singular terms by
using `α = -1 + u0.ϵ`.
- **TODO:** Compute a rigorous enclosure of the two singular terms. In
  the end we might not need this though since `F0` uses
  `skip_j_one_singular = true`.

If `skip_j_one_singular` is true then don't include the two singular
terms from the Clausen term in the tail corresponding to `j = 1`. This
is used in `F0` where these two terms are treated separately.

See [`eval_expansion`](@ref) for more details about how the
coefficients are stored.
"""
function H(
    u0::BHKdVAnsatz{Arb},
    ::AsymptoticExpansion;
    M::Integer = 3,
    skip_j_one_singular::Bool = false,
    alpha_interval = :full,
)
    @assert M >= 3

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

        # Remaining terms
        @warn "Non-rigorous bounds implemented for x^2m coefficients" maxlog = 1
        # TODO:
        for m = 2:M-1
            coefficient = -(-1)^m * dzeta(Arb(3 - 2m)) / factorial(2m) * u0.v0.a0

            # Add error bounds
            # FIXME: Implement rigorous error bounds. Currently we
            # just compute with α = -1 + u0.ϵ.
            coefficient_2 =
                let α = -1 + u0.ϵ, a0 = finda0(α), p0 = 1 + α + (1 + α)^2 / 2
                    -(-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / factorial(2m) * a0
                end

            coefficient = union(coefficient, coefficient_2)

            res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
        end

        # Error term for main term
        @warn "No error bounds for error term of main term" maxlog = 1
        # FIXME: Add error bounds for error term. Currently we don't do anything
        Arblib.add_error!(
            res[(0, 0, 0, 0, 0, 0, 2M)],
            2abs(dzeta(Arb(3 - 2M)) / factorial(2M)) * u0.v0.a0,
        )

        # Clausen terms

        # Handle the first term separately since we  since s is very close to 3 in
        # this case and it is therefore very unstable
        if u0.v0.v0.N0 >= 1
            let j = 1
                s = 1 - α - u0.v0.v0.α + j * u0.v0.v0.p0
                C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

                # Only add these terms if skip_j_one_singular is not true
                if !skip_j_one_singular
                    @warn "Using non-rigorous terms for first Clausen in tail" maxlog = 1
                    # FIXME: Compute a rigorous enclosure. In the end
                    # we might not need to use this though.
                    C2, _, p2, _ =
                        let s = 1 - (-1 + u0.ϵ) - u0.v0.v0.α + j * u0.v0.v0.p0
                            clausenc_expansion(x, s, M, skip_constant = true)
                        end

                    res[(0, 0, -1, 0, 1, j, 0)] = -C2 * u0.v0.v0.a[j]
                    res[(0, 0, 0, 0, 0, 0, 2)] -= p2[2] * u0.v0.v0.a[j]
                end

                for m = 2:M-1
                    res[(0, 0, 0, 0, 0, 0, 2m)] -= p[2m] * u0.v0.v0.a[j]
                end
                res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.v0.a[j]
            end
        end
        for j = 2:u0.v0.v0.N0
            s = 1 - α - u0.v0.v0.α + j * u0.v0.v0.p0
            C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

            res[(0, 0, -1, 0, 1, j, 0)] = -C * u0.v0.v0.a[j]
            for m = 1:M-1
                res[(0, 0, 0, 0, 0, 0, 2m)] -= p[2m] * u0.v0.v0.a[j]
            end
            res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.v0.a[j]
        end

        # Fourier terms
        if !iszero(u0.v0.N)
            for m = 1:M-1
                res[(0, 0, 0, 0, 0, 0, 2m)] -=
                    (-1)^m * sum(n^α * Arb(n)^(2m) * u0.v0.b[n] for n = 1:u0.v0.N) /
                    factorial(Arb(2m))
            end
            Arblib.add_error!(
                res[(0, 0, 0, 0, 0, 0, 2M)],
                sum(n^α * Arb(n)^(2M) * abs(u0.v0.b[n]) for n = 1:u0.v0.N) /
                factorial(Arb(2M)),
            )
        end

        return res
    end
end

function D(u0::BHKdVAnsatz, ::Asymptotic; M::Integer = 3, skip_j_one_singular = false)
    f = D(u0, AsymptoticExpansion(); M, skip_j_one_singular)

    return x -> eval_expansion(u0, f(x), x)
end

function D(
    u0::BHKdVAnsatz,
    evaltype::AsymptoticExpansion;
    M::Integer = 3,
    skip_j_one_singular = false,
    alpha_interval = :full,
)
    f = x -> u0(x, evaltype; M)
    g = H(u0, evaltype; M, skip_j_one_singular, alpha_interval)

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
    F0(u0::BHKdVAnsatz{Arb}, ::Asymptotic)

Returns a function such that an **upper bound** of `F0(u0)(x)` is
computed accurately for small values of `x`. It assumes that `x < 1`.

This method assumes that the value for `p0` is taken to be
```
p0 = 1 + α + (1 + α)^2 / 2
```

Recall that the expression we are interested in bounding is
```
abs((u0(x)^2 / 2 + H(u0)(x)) / (u0(x) * x * log(u0.c + inv(x))))
```

# Split into two factors
As a first step we split the expression into two factors which are
both bounded as `x -> 0` that we bound separately.

The first one is given by
```
F1 = abs(log(x) / log(u0.c + inv(x)) * gamma(1 + α) * x^-α * (1 - x^p0) / u0(x))
```
Notice that `log(x) / log(u0.c + inv(x))` is easily seen to be bounded
and as we will see later so is `gamma(1 + α) * x^-α * (1 - x^p0) /
u0(x)`

The second one is given by
```
F2 = abs((u0(x)^2 / 2 + H(u0)(x)) / (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0)))
```
which we will also see is bounded.

# Bounding `F1`
For `F1` we only compute an upper bound and not an enclosure, the
smaller the value of `x` the tighter it will be in general. We split
`F1` into the two factors
```
F11 = abs(log(x) / log(u0.c + inv(x)))
F12 = abs(gamma(1 + α) * x^-α * (1 - x^p0) / u0(x))
```

For `F11` we can easily compute an accurate enclosure by direct
evaluation if `x` doesn't overlap with zero. If `x` does overlap with
zero we compute an enclosure by using that `log(x) / log(u0.c + inv(x))`
converges to `-1` as `x -> 0` and is increasing in `x`.
- **PROVE:** That `log(x) / log(u0.c + inv(x))` is increasing in `x`.

For `F12` we compute the asymptotic expansion of `u0`. We then split
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
a0 * (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0) * x^-α
```
Since we are interested in `gamma(1 + α) * x^-α * (1 - x^p0)` divided
by this we get
```
gamma(1 + α) * (1 - x^p0) /
    (a0 * (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0))
```
Numerically we can see that this term converges to `π` as `α -> -1`
and `x -> 0`. We can also see that it is decreasing in both `α` and
`x`. We would therefore expect to be able to compute a rather accurate
enclosure. We don't use these observations directly though, instead we
proceed as follows. We can split this further into the two factors
```
F121 = gamma(1 + α) / a0
F122 = (1 - x^p0) /
    (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0)
```

For `F121` we can get an enclosure using that it converges to `-π^2 /
2` as `α -> -1` and is increasing in `α`.
- **PROVE:** That `gamma(1 + α) / a0` converges to `-π^2 / 2` as `α ->
  -1` and is increasing in `α`.

For `F122` we let `c(a) = gamma(a) * sinpi((1 - a) / 2)` and then
factor it out, giving us
```
inv(c(α)) * (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
```
We can enclose `c(α)` using that it converges to `-π / 2` as `α -> -1`
and is decreasing in `α`.
- **PROVE:** That `c(α)` converges to `-π / 2` as `α -> -1` and is
  decreasing in `α`.

For the remaining part we notice that `c(α - p0) / c(α)` converges to
`1` as `α -> -1` and is decreasing in `α`. An upper bound for `(1 -
x^p0) / (1 - c(α - p0) / c(α) * x^p0)` is hence given by `1` (for `α =
-1`) and a lower bound can be computed by evaluating it at `α = -1 +
u0.ϵ` and an upper bound for `x`
- **PROVE:** That `c(α - p0) / c(α)` converges to `1` as `α -> -1` and is
  decreasing in `α` and `x`.

**IMPROVE:** WE could improve the enclosure we get by incorporating
some of the terms in the tail of `u0`. For very small `x` this would
be negligible but for `x` around say `1e-10` this would slightly
improve the values. It seems to be able to give a factor of around
`0.5` for `x` close to `0.1` but only a factor `0.85` around `x =
1e-10` (checked by plotting the value for `F1` used in this method and
comparing it to the non-asymptotic version of it). Since this is not a
big improvement and for `x` values this large we can use the
non-asymptotic version anyway it is probably not worth it to implement
though.

# Bounding `F2`
Getting an accurate bound for `F2` requires some work since it's not
enough to only consider the leading term, we have to account for the
cancellation between the terms even for very small values of `x` (of
the order `1e-10000` or even smaller).

Recall that we are interested in bounding
```
F2 = abs((u0(x)^2 / 2 + H(u0)(x)) / (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0)))
```

As a first step we compute the asymptotic expansion of `u0(x)^2 / 2 +
H(u0)(x)`. We then take out the two leading terms with keys `(2, 0, 0,
0, 0, 0, 0)` and `(0, 1, 0, 0, 0, 0, 0)` which we call `P` and `Q`
respectively.

# Handling `P` and `Q`
The terms `P` and `Q` are given by
```
P = a0^2 / 2 * (c(α)^2 - 2c(α) * c(α - p0) * x^p0 + c(α - p0)^2 * x^2p0) * x^-2α
```
and
```
Q = -a0 * ((c(2α) - c(2α - p0) * x^p0) * x^-2α - (zeta(-2α - 1) / 2 - zeta(-2α + p0 - 1) / 2) * x^2)
```
where `c(a) = gamma(a) * sinpi((1 - a) / 2)` as above.

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
The factor `a0 / gamma(1 + α)` can be enclosed using that it converges
to `-2 / π^2` and is decreasing in `α`
- **PROVE:** That `a0 / gamma(1 + α)` converges to `-2 / π^2` and is decreasing.
First we focus on the term
```
F21 = (c(2α - p0) - a0 * c(α) * c(α - p0)) * x^(-α + p0 - 1) / (log(x) * (1 - x^p0))
```
This term is small and fairly stable in `α` (positive and increasing).
To compute an enclosure we split it into three factors
```
F211 = x^(-α + p0 - 1) / log(x)
F212 = (1 + α) / (1 - x^p0)
F213 = (c(2α - p0) - a0 * c(α) * c(α - p0)) / (1 + α)
```
For the first one we can directly get an enclosure using that
```
-α + p0 - 1 = -α + (1 + α + (1 + α)^2 / 2) - 1 = (1 + α)^2 / 2
```
For the second term, `F212`, we write `t = 1 + α` giving us
```
F212 = t / (1 - x^(t + t^2 / 2))
```
This is increasing in `t` so we only need to evaluate it at `t = 0`
and `t = u0.ϵ`. For `t = 0` it converges to `-inv(log(x))`. For `t =
u0.ϵ` we can evaluate it directly. To see that it converges to
`-inv(log(x))` for `t -> 0` we can rewrite it as
```
- t / (exp((t + t^2 / 2) * log(x)) - 1)
```
L'Hôpital gives us
```
- 1 / (log(x) * (1 + t) * exp((t + t^2 / 2) * log(x)))
```
and inserting `t = 0` immediately gives the limit `-inv(log(x))`.
- **PROVE:** That `t / (1 - x^(t + t^2 / 2))` is increasing in `t`.
  Possibly by multiplying with `log(x)` to get something similar to `u
  / expm1(u)`.

For the third term, `F213`, we note that `a0` can be written as
```
a0 = 2c(2α) / c(α)^2
```
which allows us to simplify it to
```
F213 = (c(2α - p0) - 2c(2α) * c(α - p0) / c(α)) / (1 + α)
```
For `α = -1` this converges to `3 // 4 - π^2 / 12`. For `α`
sufficiently close to `-1` it is also decreasing, allowing us to
compute an enclosure.
- **PROVE:** That `F213` converges to `3 // 4 - π^2 / 12` as `α ->
  -1`. Note that even though `c(α - p0) / c(α)` converges to `1` it is
  important for the asymptotics and cannot be ignored.
- **PROVE:** That `F213` is decreasing in `α` on some fixed interval
  which contains our `α`. Not that as `α` around `-0.6` it starts to
  be increasing instead, so we have to used the specified interval we
  work on.

We then consider the two remaining terms together since they mostly
cancel out.
```
F22 = (zeta(-2α - 1) - zeta(-2α + p0 - 1)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0)) +
    a0 * c(α - p0)^2 / 2 * x^(-α + 2p0 - 1) / (log(x) * (1 - x^p0))
    =
    ((zeta(-2α - 1) - zeta(-2α + p0 - 1)) +
    a0 * c(α - p0)^2 * x^(-2α + 2p0 - 2)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0))
```
We treat this in a very similar way as `F21`, by splitting it into
three factors.
```
F221 = x^(1 + α) / 2
F222 = (1 + α) / (1 - x^p0)
F223 = ((zeta(-2α - 1) - zeta(-2α + p0 - 1)) + a0 * c(α - p0)^2 * x^(-2α + 2p0 - 2)) /
    ((1 + α) * log(x))
```
The factor `F221` we can enclose directly. The factor `F222` is the
same as the factor `F212` above. We are hence left enclosing `F223`.
To do that we split it into two terms, letting
```
w1 = zeta(-2α - 1) - zeta(-2α + p0 - 1)
w2 = -a0 * c(α - p0)^2
```
and noticing that `-2α + 2p0 - 2 = (1 + α)^2` we can write
```
F223 = (w1 - w2 * x^((1 + α)^2)) / ((1 + α) * log(x))
```
Adding and subtracting `w2` in the numerator we can split this into
two terms
```
F2231 = (w1 - w2) / ((1 + α) * log(x))
F2232 = w2 * (1 - x^((1 + α)^2)) / ((1 + α) * log(x))
```
For `F2231` it enough to compute an enclosure of `(w1 - w2) / (1 + α)`
and then multiply by an enclosure of `inv(log(x))`. To get an
enclosure of `(w1 - w2) / (1 + α)` we note that it is increasing in
`α` and the limit as `α` goes to `-1` is
```
(π^2 + 8γ₁ - 4γ - 2) / 8
```
where `γᵢ` is the `i`th Stieltjes constant and `γ = γ0`.
- **PROVE:** That `(w1 - w2) / (1 + α)` converges to `(π^2 + 8γ₁ - 4γ
  - 2) / 8` and is increasing in `α`.

For `F2232` we can rewrite it as
```
F2232 = w2 * (1 + α) * (1 - x^((1 + α)^2)) / ((1 + α)^2 * log(x))
```
We can compute an enclosure of
```
(1 - x^((1 + α)^2) / ((1 + α)^2 * log(x)) =
    -(exp((1 + α)^2 * log(x)) - 1) / ((1 + α)^2 * log(x))
```
By letting `t = (1 + α)^2 * log(x)` and noticing that the function
`(exp(t) - 1) / t` is zero at `t = -Inf`, one at `t = 0` and
increasing in `t`. This is similar to how we do it for the first
Clausen term in the tail below. We are left enclosing
```
w2 * (1 + α) = -a0 * c(α - p0)^2 * (1 + α)
```
We can get an enclosure of `c(α - p0)` by using that it converges to
`-π / 2` as `α -> -1` and is increasing in `α`. We can enclose `a0 *
(1 + α)` using that it converges to `-2 / π^2` and is decreasing in
`α`
- **PROVE:** That `c(α - p0)` converges to `-π / 2` and is increasing
  in `α`.
- **PROVE:** That `a0 * (1 + α)` converges to `-2 / π^2` and is
  decreasing in `α`.

# Handling the remaining terms
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
inv(log(x) * gamma(1 + α) * (1 - x^p0))
```

## Enclosing `inv(log(x) * gamma(1 + α) * (1 - x^p0))`
We begin by noticing that `inv(log(x))` can be enclosed directly. In
the case that `x` overlaps with zero we use the monotonicity together
with that the limit is zero for `x = 0`.

We are left enclosing `inv(gamma(1 + α) * (1 - x^p0))`. For fixed `x`
this converges to `-inv(log(x))` and is increasing in `α`. To get an
enclosure it is therefore enough to compute the value of the limit as
well as the value at `α = -1 + u0.ϵ`.
- **PROVE:** That `inv(gamma(1 + α) * (1 - x^p0))` converges to
  `-inv(log(x))` and is increasing in `α`.
To get an enclosure in the case that `x` contains zero it is enough to
notice that the lower bound is zero and that it is increasing in `x`,
which is easy to see.

## Handling the first Clausen function in the tail
We are interested in bounding the first two terms in the expansion of
```
-u0.v0.v0.a[1] * clausenc(x, 1 - α - u0.v0.v0.α + u0.v0.v0.p0) /
    (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0))
```
We can get an enclosure of `inv(gamma(1 + α) * (1 - x^p0))` as
explained above. We are therefore interested in enclosing the rest.

Let `r = -u0.v0.v0.α + u0.v0.v0.p0 - 1`. Then `r > 0` and is very
small, around `1e-8` or so depending on the precise choice of
`u0.v0.v0`. We have `1 - α - u0.v0.v0.α + u0.v0.v0.p0 = 2 - α + r`.
The sum of the first two terms in the asymptotic expansion of the
Clausen is then given by
```
gamma(α - 1 - r) * sinpi((2 - α + r) / 2) * x^(1 - α + r) -
    zeta(-α + r) / 2 * x^2
```
Dividing by `log(x) * x^(1 - α)` gives us
```
(gamma(α - 1 - r) * sinpi((2 - α + r) / 2) * x^r - zeta(-α + r) / 2 * x^(1 + α)) / log(x)
```
Notice that the order of the terms depends on the value of `α`, in
some cases `x`r` is leading and in some cases `x^(1 - α)`. Adding and
subtracting `zeta(-α + r) / 2 * x^r` we can rewrite this as
```
(gamma(α - 1 - r) * sinpi((2 - α + r) / 2) - zeta(-α + r) / 2) * x^r / log(x) +
zeta(-α + r) / 2 * (x^r - x^(1 + α)) / log(x)
```
For the first term we can note that
```
gamma(α - 1 - r) * sinpi((2 - α + r) / 2) - zeta(-α + r) / 2
```
is bounded for in `α` and decreasing, so we can get an enclosures by
evaluating it at the endpoints.
- **PROVE:** That `gamma(α - 1 - r) * sinpi((2 - α + r) / 2) - zeta(-α
  + r) / 2` is decreasing in `α`.
Denoting this enclosure by `C` we get `C * x^r / log(x)` for the first
term, which is easily enclosed.

For the second term we have to compute an enclosure of
```
zeta(-α + r) / 2 * (x^r - x^(1 + α)) / log(x)
```
From the Laurent series of the zeta function we have
```
zeta(-α + r) = zeta(1 + (-α + r - 1)) = inv(-α + r - 1) +
    sum((-1)^n / factorial(n) * γₙ * (-α + r - 1)^n for n = 1:Inf)
```
where `γₙ` is the nth Stieltjes constant. Calling the sum `S` we can
thus rewrite the above as
```
(x^r - x^(1 + α)) / (-α + r - 1) / 2log(x) + S * (x^r - x^(1 + α)) / 2log(x)
```
Now `S` is bounded and increasing in `α` so we can compute an
enclosure of the second term directly. Factoring
- **PROVE:** That `S` is increasing in `α`
The remaining term to handle is, after slightly rewriting it,
```
(x^r - x^(1 + α)) / (2(r - (1 + α)) * log(x))
```
We split this into two cases, when `r >= 1 + α` and when `r < 1 + α`.
In the first case we factor out `x^(1 + α)`, giving us
```
x^(1 + α) / 2 * (x^(r - (1 + α)) - 1) / ((r - (1 + α)) * log(x)) =
    x^(1 + α) / 2 * (exp((r - (1 + α)) * log(x)) - 1) / ((r - (1 + α)) * log(x)) =
```
Here we notice that since `r - (1 + α) >= 0` and `log(x) < 0` we have
```
(r - (1 + α)) * log(x) <= 0
```
Furthermore the function `(exp(t) - 1) / t` is zero at `t = -Inf`, one
at `t = 0` and increasing in `t`. Using this we can compute an
enclosure of
```
(exp((1 + α - r) * log(x)) - 1) / ((1 + α - r) * log(x))
```
The second case, when `r < 1 + α`, doesn't always occur, it depends on
the value of `u0.ϵ`. We first check if `1 + α - r` contain any
positive numbers, if that is the case we proceed similar to for the
first case. We factor out `x^r`, giving us
```
x^r / 2 * (1 - x^(1 + α - r)) / ((r - (1 + α)) * log(x)) =
    x^r / 2 * (x^(1 + α - r) - 1) / ((1 + α - r) * log(x)) =
    x^r / 2 * (exp((1 + α - r) * log(x)) - 1) / (((1 + α) - r) * log(x))
```
And in this case `1 + α -r > 0` and `log(x) < 0` so we have
```
(1 + α - r) * log(x) < 0
```
This allows us to compute an enclosure similar to how we did it in the
first case.
"""
function F0(
    u0::BHKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 3,
    ϵ::Arb = Arb(0.5),
    alpha_interval = :full,
)
    @assert ϵ < 1

    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), x * log(u0.c + inv(x)))
    end

    # Compute the expansion of u0 and remove the leading term, which
    # is handled separately.
    u0_expansion = u0(ϵ, AsymptoticExpansion(); M)
    delete!(u0_expansion, (1, 0, 0, 0, 0, 0, 0))

    # Ensure that the tail of the expansion of u0 is positive, so that
    # we can remove it from the denominator of F1 and still get an
    # upper bound.
    #expansion_ispositive(u0, u0_expansion, ϵ) ||
    #    error("expansion of u0 not prove to be positive, this should not happen")

    # Compute the expansion of D(u0), skipping the Clausen term in the
    # tail corresponding to j = 1 and also remove the two leading
    # term, the three terms are handled separately.
    Du0_expansion =
        D(u0, AsymptoticExpansion(), skip_j_one_singular = true; M, alpha_interval)(ϵ)
    delete!(Du0_expansion, (2, 0, 0, 0, 0, 0, 0))
    delete!(Du0_expansion, (0, 1, 0, 0, 0, 0, 0))

    # Divide the expansion of D(u0) by x^(1 - α)
    Du0_expansion_div_x_onemα = empty(Du0_expansion)
    for ((p, q, i, j, k, l, m), y) in Du0_expansion
        Du0_expansion_div_x_onemα[(p, q, i + 1, j, k, l, m - 1)] = y
    end
    #return Du0_expansion_div_x_onemα
    c(a) = gamma(a) * sinpi((1 - a) / 2)

    return x::Arb -> begin
        @assert x <= ϵ

        # Compute an upper bound of F1
        F1 = let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, xᵤ = ubound(Arb, x)
            # Note that inside this statement α refers to -1 + u0.ϵ
            # and p0 to the corresponding p0 value.

            # Enclosure of abs(log(x) / log(u0.c + inv(x))), either by
            # direct evaluation or using monotonicity.
            F11 = if iszero(x)
                one(Arb)
            elseif Arblib.contains_zero(x)
                abs(Arb((-1, log(xᵤ) / log(u0.c + inv(xᵤ)))))
            else
                abs(log(x) / log(u0.c + inv(x)))
            end

            # Enclose F12
            F121_lower = -Arb(π)^2 / 2
            F122_upper = gamma(1 + α) / finda0(α)
            F121 = Arb((F121_lower, F122_upper))

            # Upper and lower bound of
            # (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
            F122_lower = (1 - xᵤ^p0) / (1 - c(α - p0) / c(α) * xᵤ^p0)
            F122_upper = one(Arb)
            # Combine upper and lower bound and multiply with
            # enclosure of inv(c(α)) to get an upper bound for F122.
            F122 = inv(Arb((c(α), -Arb(π) / 2))) * Arb((F122_lower, F122_upper))

            F12 = F121 * F122

            Arblib.ispositive(F12) ||
                error("leading term of u0 is not positive, this should not happen")

            F11 * F12
        end

        # Compute an enclosure of F2
        F2 = let α = Arb((-1, -1 + u0.ϵ))
            # Enclosure of inv(log(x))
            invlogx = if iszero(x)
                zero(x)
            elseif Arblib.contains_zero(x)
                xᵤ = ubound(Arb, x)
                Arb((inv(log(xᵤ)), 0))
            else
                inv(log(x))
            end

            # Start by handling the terms P and Q

            # Enclosure of a0 / gamma(1 + α)
            a0gamma = let α = -1 + u0.ϵ
                Arb((finda0(α) / gamma(1 + α), -2 / Arb(π)^2))
            end

            # Compute an enclosure of F21

            # Enclosure of F211 = x^(-α + p0 - 1) / log(x) using that
            # -α + p0 - 1 = (1 + α)^2 / 2
            F211 = abspow(x, Arblib.nonnegative_part!(zero(x), (1 + α)^2 / 2)) * invlogx

            # Enclosure of F212 = (1 + α) / (1 - x^p0) using t = 1 + α
            F212 = let
                lower = -invlogx
                upper = let t = u0.ϵ
                    t / (1 - abspow(x, t + t^2 / 2))
                end

                Arb((lower, upper))
            end

            # Enclosure of F213 = (c(2α - p0) - a0 * c(α) * c(α - p0)) / (1 + α)
            F213 = let
                lower = let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2
                    (c(2α - p0) - 2c(2α) * c(α - p0) / c(α)) / (1 + α)
                end
                upper = 3 // 4 - Arb(π)^2 / 12

                Arb((lower, upper))
            end

            F21 = F211 * F212 * F213

            # Compute an enclosure of F22

            # Enclosure of F221 = x^(1 + α) / 2
            F221 = abspow(x, Arblib.nonnegative_part!(zero(x), 1 + α)) / 2

            # Enclosure of F222, which is the same as F212
            F222 = F212

            F223 = let
                F2231 = let
                    # Lower and upper bound of (w1 - w2) / (1 + α)
                    lower = (Arb(π)^2 + 8stieltjes(Arb, 1) - 4stieltjes(Arb, 0) - 2) / 8
                    upper =
                        let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, a0 = finda0(α)
                            w1 = zeta(-2α - 1) - zeta(-2α + p0 - 1)
                            w2 = -a0 * c(α - p0)^2

                            (w1 - w2) / (1 + α)
                        end

                    Arb((lower, upper)) * invlogx
                end

                F2232 = let
                    # Enclosure of (1 - x^((1 + α)^2) / ((1 + α)^2 * log(x))
                    F2232 = if Arblib.contains_zero(x)
                        Arblib.unit_interval!(zero(x))
                    else
                        # We compute a lower bound of t = (1 + α)^2 * log(x)
                        t_lower = lbound(
                            Arb,
                            Arblib.nonnegative_part!(zero(x), (1 + α)^2) * log(x),
                        )

                        -Arb((expm1(t_lower) / t_lower, 1))
                    end

                    # Multiply by enclosure of c(α - p0)^2
                    F2232 *= let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2
                        Arb((-Arb(π) / 2, c(α - p0)))^2
                    end

                    # Multiply by enclosure of a0 * (1 + α)
                    F2232 *= let α = -1 + u0.ϵ
                        Arb((finda0(α) * (1 + α), -2 / Arb(π)^2))
                    end

                    -F2232
                end

                F2231 + F2232
            end

            F22 = F221 * F222 * F223

            # The enclosure of the terms coming from P + Q in the expansion
            P_plus_Q = a0gamma * (F21 + F22)

            # Enclosure of inv(gamma(1 + α) * (1 - x^p0))

            invgamma1mxp0 = let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2
                if iszero(x)
                    lower = zero(x)
                    upper = inv(gamma(1 + α))
                    Arb((lower, upper))
                elseif Arblib.contains_zero(x)
                    lower = zero(x)
                    upper = inv(gamma(1 + α) * (1 - ubound(Arb, x)^p0))
                    Arb((lower, upper))
                else
                    lower = -inv(log(x))
                    upper = inv(gamma(1 + α) * (1 - x^p0))
                    Arb((lower, upper))
                end
            end

            # Handle the two singular terms in the expansion of
            # clausenc(x, 1 - α - u0.v0.v0.α + u0.v0.v0.p0)
            if u0.v0.v0.N0 > 0
                clausen_j_one =
                    let α = Arb((-1, -1 + u0.ϵ)), r = -u0.v0.v0.α + u0.v0.v0.p0 - 1
                        # Enclosure of
                        # gamma(α - 1 - r) * sinpi((2 - α + r) / 2) - zeta(-α + r) / 2
                        C_lower = let α = -1 + u0.ϵ
                            (gamma(α - 1 - r) * sinpi((2 - α + r) / 2) - zeta(-α + r) / 2)
                        end
                        C_upper = let α = Arb(-1)
                            (gamma(α - 1 - r) * sinpi((2 - α + r) / 2) - zeta(-α + r) / 2)
                        end
                        C = Arb((C_lower, C_upper))

                        # Enclosure of C * x^r / log(x)
                        term1 = C * abspow(x, r) * invlogx

                        # Enclosure of zeta(-α + r) - inv(α + r - 1)
                        S_lower = let α = -1 + u0.ϵ
                            zeta(-α + r) - inv(-α + r - 1)
                        end
                        S_upper = let α = Arb(-1)
                            zeta(-α + r) - inv(-α + r - 1)
                        end
                        S = Arb((S_lower, S_upper))

                        # Enclosure of (x^r - x^(1 + α)) / (-α + r - 1) / 2log(x)
                        term2 = let
                            # Handle the case r > 1 + α

                            # Enclosure of (exp(t) - 1) / t with
                            # t = (r - (1 + α)) * log(x)
                            factor1 = if Arblib.contains_zero(x)
                                # We have t = (r - (1 + α)) * log(x) in interval [-Inf, 0]
                                Arblib.unit_interval!(zero(x))
                            else
                                # We compute a lower bound of t = (r - (1 + α)) * log(x)

                                # First we let s be the non-negative part
                                # of r - (1 + α), since we are assuming that
                                # r >= 1 + α.
                                s = Arblib.nonnegative_part!(zero(x), r - (1 + α))

                                # Then we take the lower bound of s * log(x)
                                t_lower = lbound(Arb, s * log(x))

                                Arb((expm1(t_lower) / t_lower, 1))
                            end

                            # Enclosure of x^(1 + α) /2 * factor1
                            term2 =
                                abspow(x, Arblib.nonnegative_part!(zero(x), 1 + α)) / 2 * factor1

                            # Handle the case r < 1 + α if it occurs
                            if Arblib.contains_positive(1 + α - r)
                                # Enclosure of (exp(t) - 1) / t with
                                # t = (1 + α - r) * log(x)
                                factor2 = if Arblib.contains_zero(x)
                                    # We have t = (1 + α - r) * log(x) in interval [-Inf, 0]
                                    Arblib.unit_interval!(zero(x))
                                else
                                    # We compute a lower bound of t = (1 + α - r) * log(x)

                                    # First we let s be the non-negative part
                                    # of 1 + α - r, since we are assuming that
                                    # r < 1 + α.
                                    s = Arblib.nonnegative_part!(zero(x), 1 + α - r)

                                    # Then we take the lower bound of s * log(x)
                                    t_lower = lbound(Arb, s * log(x))

                                    Arb((expm1(t_lower) / t_lower, 1))
                                end

                                # Set result to union of the above case
                                # and this case which is given by x^r / 2 * factor2
                                term2 = union(term2, abspow(x, r) / 2 * factor2)
                            end

                            term2
                        end

                        # Enclosure of S * (x^r - x^(1 + α)) / 2log(x)
                        term3 = let onepα = Arblib.nonnegative_part!(zero(x), 1 + α)
                            S * (abspow(x, r) - abspow(x, onepα)) * invlogx / 2
                        end

                        term = term1 + term2 + term3

                        term *= invgamma1mxp0

                        -u0.v0.v0.a[1] * term
                    end
            else
                clausen_j_one = zero(x)
            end

            # Enclosure of the remaining terms in the expansion
            remainder =
                eval_expansion(u0, Du0_expansion_div_x_onemα, x) * invlogx * invgamma1mxp0

            #(u0(x)^2 / 2 + Hu0x) / (log(x) * gamma(1 + α) * x^(1 - α) * (1 - x^p0))
            P_plus_Q + clausen_j_one + remainder
        end

        return F1 * F2
    end
end
