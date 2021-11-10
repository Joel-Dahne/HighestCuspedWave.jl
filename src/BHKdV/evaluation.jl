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
    x;
    div_a0 = false,
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

    # Compute i * α + j * p0 + m in a way that accounts for dependence
    # between α and p0 and computes more accurate enclosures in some
    # special cases. We use p0 = 1 + α + (1 + α)^2 / 2
    exponent_α_p0_m(i, j, m; verbose = false) = begin
        if i + 2j >= 0
            # It is increasing in α, evaluated at endpoints

            # α = -1 can be done with rational numbers
            lower = let α = -1
                Arb(i * α + j * (1 + α + (1 + α)^2 // 2) + m)
            end

            upper = let α = -1 + u0.ϵ
                i * α + j * (1 + α + (1 + α)^2 / 2) + m
            end

            enclosure = Arb((lower, upper))

            # If the lower bound is zero we want to avoid any spurious
            # negative parts
            iszero(lower) && Arblib.nonnegative_part!(enclosure, enclosure)

            return enclosure
        end

        return m + (i + j) * α + j * (1 + (1 + α)^2 / 2)
    end

    # Irrationals used
    π = Arb(Irrational{:π}())
    γ = Arb(Irrational{:γ}())

    res = zero(x)

    for ((p, q, i, j, k, l, m), y) in expansion
        if !iszero(y)
            if iszero(p)
                exponent = exponent_α_p0_m(i, j, m) - k * u0.v0.v0.α + l * u0.v0.v0.p0

                term = y * abspow(x, exponent)
            elseif isone(p)
                # Add -α to the exponent coming from the p factor
                exponent = exponent_α_p0_m(i - 1, j, m) - k * u0.v0.v0.α + l * u0.v0.v0.p0

                # FIXME: Add error bounds for this term. Here we just
                # use the limiting value

                # Compute an enclosure using monotonicity
                if iszero(x)
                    term = zero(x)
                elseif Arblib.contains_zero(x)
                    if Arblib.ispositive(exponent)
                        lower = zero(x)
                    else
                        lower = Arblib.indeterminate!(zero(x))
                    end

                    upper = let x = ubound(Arb, x)
                        (1 - γ - log(x)) / π * abspow(x, exponent)
                    end

                    term = y * Arb((lower, upper))
                else
                    term = y * (1 - γ - log(x)) / π * abspow(x, exponent)
                end
            else
                use_approx_p_and_q ||
                    throw(ArgumentError("only p == 0 or p == 1 supported, got p = $p"))

                exponent = exponent_α_p0_m(i, j, m) - k * u0.v0.v0.α + l * u0.v0.v0.p0

                term = y * abspow(x, exponent)

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

            res += term
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
    skip_j_one = false,
    use_approx_p_and_q = false,
)
    f = H(u0, AsymptoticExpansion(); M, skip_j_one)

    return x -> eval_expansion(u0, f(x), x; use_approx_p_and_q)
end

"""
    H(u0::BHKdVAnsatz, ::AsymptoticExpansion; M = 3, skip_j_one = false)

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

For the Clausen terms just we let `α` be a ball for `j >= 2`, this
gives good enclosures thanks to [`clausenc_expansion`](@ref). For `j =
1` this doesn't work since `1 - α - u0.v0.v0.α + u0.v0.v0.p0` contains
the integer `3` and the first two coefficients in the expansion blow
up and collapse together.
- **TODO:** Figure out how to handle the case `j = 1`. Currently we
  just take the value at `α = -1` but this is probably not good
  enough.

For the tail term the Fourier terms are handled directly by letting
`α` be a ball.

If `skip_j_one` is true then skip the Clausen term corresponding to `j
= 1`. This is used in `F0` where this term is treated separately.

See [`eval_expansion`](@ref) for more details about how the
coefficients are stored.
"""
function H(
    u0::BHKdVAnsatz{Arb},
    ::AsymptoticExpansion;
    M::Integer = 3,
    skip_j_one::Bool = false,
    alpha_interval = :full,
)
    @assert M >= 3

    # Terms used when computing error bounds
    α = -1 + u0.ϵ
    a0 = finda0(α)
    p0 = 1 + α + (1 + α)^2 / 2
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
        # TODO: Implement rigorous error bounds
        for m = 2:M-1
            coefficient = -(-1)^m * dzeta(Arb(3 - 2m)) / factorial(Arb(2m)) * u0.v0.a0

            # Add error bounds
            coefficient_2 =
                -a0 * (-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / factorial(2m)

            coefficient = union(coefficient, coefficient_2)

            res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
        end

        # Error term for main term
        @warn "No error bounds for error term of main term" maxlog = 1
        Arblib.add_error!(
            res[(0, 0, 0, 0, 0, 0, 2M)],
            2abs(dzeta(Arb(3 - 2M)) / factorial(Arb(2M))) * u0.v0.a0,
        )

        # Clausen terms

        # Handle the first term manually since s is very close to 3 in
        # this case and it is therefore very unstable
        if u0.v0.v0.N0 >= 1 && !skip_j_one
            let j = 1
                # TODO: Figure out how to handle this. Currently we just take
                # α = -1 + u0.ϵ.
                @warn "Clausen term with j = $j in tail not enclosed" maxlog = 1
                s = 1 - α - u0.v0.v0.α + j * u0.v0.v0.p0
                C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

                res[(0, 0, -1, 0, 1, j, 0)] = -C * u0.v0.v0.a[j]
                for m = 1:M-1
                    res[(0, 0, 0, 0, 0, 0, 2m)] -= p[2m] * u0.v0.v0.a[j]
                end
                res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.v0.a[j]
            end
        end
        let α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α
            for j = 2:u0.v0.v0.N0
                if alpha_interval == :full
                    s = 1 - α - u0.v0.v0.α + j * u0.v0.v0.p0
                elseif alpha_interval == :endpoint
                    s = 1 - (-1 + u0.ϵ) - u0.v0.v0.α + j * u0.v0.v0.p0
                end
                #s = 2 - u0.v0.v0.α + j * u0.v0.v0.p0
                C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

                res[(0, 0, -1, 0, 1, j, 0)] = -C * u0.v0.v0.a[j]
                for m = 1:M-1
                    res[(0, 0, 0, 0, 0, 0, 2m)] -= p[2m] * u0.v0.v0.a[j]
                end
                res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.v0.a[j]
            end
        end

        # Fourier terms
        let α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α
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
        end

        return res
    end
end

function D(u0::BHKdVAnsatz, ::Asymptotic; M::Integer = 3, skip_j_one = false)
    f = D(u0, AsymptoticExpansion(); M, skip_j_one)

    return x -> eval_expansion(u0, f(x), x)
end

function D(
    u0::BHKdVAnsatz,
    evaltype::AsymptoticExpansion;
    M::Integer = 3,
    skip_j_one = false,
    alpha_interval = :full,
)
    f = x -> u0(x, evaltype; M)
    g = H(u0, evaltype; M, skip_j_one, alpha_interval)

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

For now this method assumes that the value for `p0` is taken to be
```
p0 = 1 + α + (1 + α)^2 / 2
```
though it could be adjusted for other values.

# Split into weight part and bounded part
The asymptotic expansions for `u0` and `D(u0)` are computed. From the
former we factor out `x^-α` and from the latter `x^(-2α + p0)`. This
leaves us with the factor
```
x^(-2α + p0) / (u0.w(x) * x^-α)
```
that we have to handle. With `u0.w(x) = x * log(u0.c + inv(x))` this
simplifies to
```
x^(-1 - α + p0) / log(u0.c + inv(x))
```
Since `-1 - α + p0 > 0` for all values of `α` this can be upper
bounded by `inv(log(u0.c + inv(x))). For a lower bound we use that
```
-1 - α + p0 = (1 + α)^2 / 2 < u0.ϵ^2 / 2
```

# Simplify denominator
We can still not evaluate `u0(x) / x^-α` and `D(u0)(x) / x^(-2α + p0)`
directly since both of them diverge. As a first step we extract the
leading part of the expansion of `u0(x) / x^-α`, corresponding to the
term `(1, 0, 1, 0, 0, 0, 0)`, given by
```
a0 * (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0)
```
We then prove that the remaining part of `u0(x) / x^-α` is positive
and if we remove it we hence still get an upper bound of the value.
- **TODO:** Prove that the remaining part is positive by factoring out
    the leading power of `x` and evaluating directly.
If we let
```
w1 = gamma(α) * sinpi((1 - α) / 2)
w2 = gamma(α - p0) * sinpi((1 - α + p0) / 2)
```
it is therefore enough to consider
```
(D(u0)(x) / x^(-2α + p0)) / (a0 * (w1 - w2 * x^p0))
```

- **TODO:** Consider keeping the tail to get better bounds.

# Split numerator
The next step is to split `D(u0)(x) / x^(-2α + p0)` into three parts,
the term corresponding to `(2, 0, 2, -1, 0, 0, 0)`. which we will call
`P`, the term corresponding to `(0, 1, 2, -1, 0, 0, 0)`, which we will
call `Q`, and the remaining terms, which we will call `R`.

- **TODO:** Is this splitting good enough? Or do we need to take into
  account the terms in the tail as well? In particular the leading
  term in the tail might be useful to handle separately.

# Handle the terms `P` and `Q`
The two terms `P` and `Q` we have to treat together to account for the
cancellation that occurs between them. We have
```
P = a0^2 * / 2 * (
    (gamma(α) * sinpi((1 - α) / 2))^2 * x^-p0 -
    2gamma(α) * sinpi((1 - α) / 2) * gamma(α - p0) * sinpi((1 - α + p0) / 2) +
    (gamma(α - p0) * sinpi((1 - α + p0) / 2))^2 * x^p0
)
```
and
```
Q = -a0 * (
    gamma(2α) * sinpi((1 - 2α) / 2) * x^-p0 -
    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2) -
    (zeta(-2α - 1) / 2 - zeta(-2α + p0 - 1) / 2) * x^(2 + 2α - p0)
)
```
By construction `a0` is taken such that the terms with `x^-p0` cancel
out, leaving us with
```
P + Q = a0 * (
    a0 * gamma(α) * sinpi((1 - α) / 2) * gamma(α - p0) * sinpi((1 - α + p0) / 2) -
    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2) +
    (zeta(-2α - 1) / 2 - zeta(-2α + p0 - 1) / 2) * x^(2 + 2α - p0) +
    a0 / 2 * (gamma(α - p0) * sinpi((1 - α + p0) / 2))^2 * x^p0
)
```

We are now interested in bounding `(P + Q) / (a0 * (w1 - w2 * x^p0))`.
Cancelling the `a0` factor we have
```
((P + Q) / a0) / (w1 - w2 * x^p0)
```
However both the numerator and denominator tend to zero and we run
into problems bounding it. A more delicate approach is therefore
required to deal with it. We split `P + Q` into the part which is
constant in `x` and the rest and bound them separately.

- **TODO:** Take into account the absolute value!

## Handle the constant part of `P + Q`
For the constant part we are interested in bounding
```
(
    a0 * gamma(α) * sinpi((1 - α) / 2) * gamma(α - p0) * sinpi((1 - α + p0) / 2) -
    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2)
) / (
    (w1 - w2 * x^p0)
)
```
Since `abs(w1) > abs(w2)` and they have the same sign this is
decreasing in `x`, an upper bound is hence given by
```
(
    a0 * gamma(α) * sinpi((1 - α) / 2) * gamma(α - p0) * sinpi((1 - α + p0) / 2) -
    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2)
) / w1
```
which is zero at `α = -1` and decreasing. An enclosure can thus be
computed by evaluating it at `α = -1 + u0.ϵ`.
- **PROVE:** That `abs(w1) > abs(w2)` and they have the same sign.
- **PROVE:** That the expression indeed is decreasing and zero at `α =
    -1`.

## Handling remaining part of `P + Q`
For the other part we are interested in bounding
```
(
    a0 / 2 * (gamma(α - p0) * sinpi((1 - α + p0) / 2))^2 * x^p0 +
    (zeta(-2α - 1) / 2 - zeta(-2α + p0 - 1) / 2) * x^(2 + 2α - p0)
) / (
    (w1 - w2 * x^p0)
)
```
Let
```
v1 = a0 / 2 * (gamma(α - p0) * sinpi((1 - α + p0) / 2))^2
v2 = -(zeta(-2α - 1) / 2 - zeta(-2α + p0 - 1) / 2)
```
This gives us
```
(v1 * x^p0 - v2 * x^(2 + 2α - p0)) / (w1 - w2 * x^p0)
```
Since `abs(w1) > abs(w2)` and they have the same sign we can get an
upper bound by replacing `w2` with `w1` and then factoring it out,
giving us
```
inv(w1) * (v1 * x^p0 - v2 * x^(2 + 2α - p0)) / (1 - x^p0)
```
An enclosure of `inv(w1)` is easily computed and we are left with
bounding the remaining quotient.

Expanding `v1` and `v2` at `α = -1` gives us
```
v1 = -1 / 4(α + 1) + O(1)
v2 = -1 / 4(α + 1) + O(1)
```
and we can hence rewrite the quotient as
```
((-1 / 4(α + 1) + v1 + 1 / 4(α + 1)) * x^p0 - (-1 / 4(α + 1) + v2 + 1 / 4(α + 1)) * x^(2 + 2α - p0)) / (1 - x^p0)
```
which splits into two parts, the first one being
```
-1 // 4 * (x^p0 - x^(2 + α - p0)) / ((α + 1) * (1 - x^p0))
```
- **PROVE:** That `(x^p0 - x^(2 + α - p0)) / ((α + 1) * (1 - x^p0))`
    is contained in the interval `[-1, 0]`

The second one is
```
((v1 + 1 / 4(α + 1)) * x^p0 - (v2 + 1 / 4(α + 1)) * x^(2 + 2α - p0)) / (1 - x^p0)
```
Letting
```
a = (v1 + 1 / 4(α + 1))
b = (v2 + 1 / 4(α + 1))
```
and factoring out `-a * x^(2 + 2α - p0)` gives
```
-a * x^(2 + 2α - p0) * (b / a - x^(-2 - 2α + 2p0)) / (1 - x^p0)
```
An enclosure of `a` can be computed using that it is `-1 // 4` at `α =
-1` and increasing in `α`. The factor `x^(2 + 2α - p0)` we can upper
bound by `1`.
- **PROVE:** That `a` is `-1 // 4` at `α = 1` and increasing in `α`.
Focusing on the quotient we can rewrite it as
```
((1 - x^(-2 - 2α + 2p0)) + (b / a - 1)) / (1 - x^p0)
```
where
```
(1 - x^(-2 - 2α + 2p0)) / (1 - x^p0)
```
can be bounded by `1` since `-2 - 2α + 2p0 = (1 + α)^2 < p0`. We are
left bounding
```
A = (b / a - 1) / (1 - x^p0)
```
- **TODO:** Compute an enclosure/upper bound of `A`
Finally an upper for the whole expression is hence given by
```
-a * (1 + A)
```

# Handle rest of numerator
Finally we need to bound what we get with the terms remaining in the
expansion for `D(u0)` after we took away the first two terms. We are
interested in bounding
```
R / (a0 * (w1 - w2 * x^p0))
```
The numerator can be enclosed directly with [`eval_expansion`](@ref).
For the denominator it's enough to bound it away from zero. Since
`abs(w1) > abs(w2)` and they have the same sign (this is proved above)
we can get a value closer to zero by considering
```
a0 * (w1 - w1 * x^p0) = w1 * a0 * (1 - x^p0)
```
We can also note that the function is decreasing in `x` so it is
enough to lower bound it for `ubound(x)`.

To handle `a0 * (1 - x^p0)` we expand `a0` and `1 - x^p0` separately
in `α`. We have
```
a0 = -2 / (π^2 * (α + 1)) + O(1)
```
and
```
1 - x^p0 = -log(x) * (α + 1) + O((α + 1)^2)
```
Multiplying these two together we get
```
a0 * (1 - x^p0) = 2log(x) / π^2 + O(α + 1)
```
which we can evaluate at `ubound(x)`.
- **TODO:** Explicitly bound the error terms and incorporate that in
    the bound.
"""
function F0(u0::BHKdVAnsatz{Arb}, ::Asymptotic; M::Integer = 3, ϵ::Arb = Arb(1e-2))
    @assert ϵ < 1

    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), x * log(u0.c + inv(x)))
    end

    u0_expansion = u0(ϵ, AsymptoticExpansion(); M)
    Du0_expansion = D(u0, AsymptoticExpansion(); M)(ϵ)

    # Divide the u0 expansion by x^-α
    u0_expansion_div_xmα = empty(u0_expansion)
    for ((p, q, i, j, k, l, m), y) in u0_expansion
        u0_expansion_div_xmα[(p, q, i + 1, j, k, l, m)] = y
    end

    # Divide the Hu0 expansion by x^(-2α + p0)
    Du0_expansion_div_xm2αp0 = empty(Du0_expansion)
    for ((p, q, i, j, k, l, m), y) in Du0_expansion
        Du0_expansion_div_xm2αp0[(p, q, i + 2, j - 1, k, l, m)] = y
    end

    #return u0_expansion_div_xmα
    #return Du0_expansion_div_xm2αp0

    # Take out the term (1, 0, 1, 0, 0, 0, 0) from u0
    delete!(u0_expansion_div_xmα, (1, 0, 1, 0, 0, 0, 0))

    # Take out the two terms (2, 0, 2, -1, 0, 0, 0) and (0, 1, 2, -1,
    # 0, 0, 0) from the expansion of D(u0) and handle them separately.
    # Call the result R
    delete!(Du0_expansion_div_xm2αp0, (2, 0, 2, -1, 0, 0, 0))
    delete!(Du0_expansion_div_xm2αp0, (0, 1, 2, -1, 0, 0, 0))
    R = Du0_expansion_div_xm2αp0

    #return u0_expansion_div_xmα
    #return R

    # x^(-2α + p0) / (u0.w(x) * x^-α) = x^(-1 - α + p0) / log(u0.c + inv(x))
    weight(x) = begin
        iszero(x) && return zero(x)

        if Arblib.contains_zero(x)
            # Use that it is zero at x = 0 and monotonically
            # increasing in x

            # Lower bound is zero
            lower = zero(x)

            # Upper bound is inv(log(u0.c + inv(x))) evaluated at upper
            # bound for x
            upper = inv(log(u0.c + inv(ubound(Arb, x))))

            return Arb((lower, upper))
        end

        # (-1 - α + p0) is upper bounded by u0.ϵ^2 / 2
        lower = x^(u0.ϵ^2 / 2) / log(u0.c + inv(x))

        upper = inv(log(u0.c + inv(x)))

        return Arb((lower, upper))
    end

    return x -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && x[0] <= ϵ)

        # Compute an enclosure of w1 using the limiting value at α
        # = -1 being -π / 2 and the monotonicity in α
        # PROVE: The monotonicity in α and the limit
        w1 = let α = -1 + u0.ϵ
            lower = gamma(α) * sinpi((1 - α) / 2)
            upper = -Arb(π) / 2

            Arb((lower, upper))
        end

        # TODO: Prove that u0_expansion_div_xmα is non-negative, so
        # that we can neglect it.

        # Upper bound of the part corresponding to P + Q
        res1 = begin
            # Upper bound of the constant part, consisting of the
            # first term
            part1 = let α = -1 + u0.ϵ, w1 = gamma(α) * sinpi((1 - α) / 2)
                a0 = finda0(α)
                p0 = 1 + α + (1 + α)^2 / 2

                (
                    a0 *
                    gamma(α) *
                    sinpi((1 - α) / 2) *
                    gamma(α - p0) *
                    sinpi((1 - α + p0) / 2) -
                    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2)
                ) / w1
            end

            # Upper bound of the other part, consisting of the remaining two terms
            # FIXME: Compute a rigorous upper bound.

            part2 = begin
                main_term = -1 // 4 * Arb((-1, 0))

                # Compute an enclosure of a using the limit -1 // 4 at
                # -1 and that it is increasing.
                # PROVE: That this holds.
                a = let α = -1 + u0.ϵ
                    a0 = finda0(α)
                    p0 = 1 + α + (1 + α)^2 / 2

                    lower = -1 // 4
                    upper =
                        a0 / 2 * (gamma(α - p0) * sinpi((1 - α + p0) / 2))^2 + 1 / 4(α + 1)

                    Arb((lower, upper))
                end

                A = zero(x) # TODO: Compute an enclosure of A
                remainder_term = -a * (1 + A)

                inv(w1) * (main_term + remainder_term)
            end

            part1 + part2
        end

        # Upper bound of the remaining part
        res2 = begin
            numerator = eval_expansion(u0, R, x)

            # TODO: Lower bound of w1 * a0 * (1 / x^p0 - 1). This is an
            # approximation but we need to bound the error terms.
            denominator = w1 * 2log(ubound(Arb, x)) / Arb(π)^2

            numerator / denominator
        end

        return (res1 + res2) * weight(x)
    end
end

"""
    F0_nonzero(u0::BHKdVAnsatz{Arb}, ::Asymptotic)

Returns a function such that an **upper bound** of `F0(u0)(x)` is
computed accurately for small values of `x`. It assumes that `x < 1`.

This method is intended for the case when `x` doesn't overlap with
zero, in practice we use it for the interval `[1e-200000, 1e-5]`.

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
For `F1` we give bounds which work for `x` overlapping zero as well,
the reason for this is that the same expression comes up in other
estimates where we do need that. The bound we give is only an upper
bound and not an enclosure, the smaller the value of `x` the tighter
it will be in general. We split `F1` into the two factors
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
- **TODO:** Check expression for P and Q

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
For now we compute it by letting `α = -1 + u0.ϵ`
- **TODO:** Figure out how to bound this rigorously for the full range of `α`.
We then consider the two remaining terms together since they mostly
cancel out.
```
F22 = (zeta(-2α - 1) - zeta(-2α + p0 - 1)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0)) +
    a0 * c(α - p0)^2 / 2 * x^(-α + 2p0 - 1) / (log(x) * (1 - x^p0))
    =
    ((zeta(-2α - 1) - zeta(-2α + p0 - 1)) +
    a0 * c(α - p0)^2 * x^(-2α + 2p0 - 2)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0))
```
- **TODO:** Figure out how to bound this rigorously for the full range of `α`

# Handling the remaining terms
Once the terms `P` and `Q` have been taken out from the expansion it
is possible to enclose the remaining ones directly. However this gives
very bad enclosures, in particular for some of the terms which have
large cancellations. To handle this we extract some of the terms and
handle them explicitly as well.

The first step is to separately handle the terms in `H(u0)` coming
from `clausenc(x, 1 - α - u0.v0.v0.α + u0.v0.v0.p0)`. There
are two main reasons to handle this term separately, that the
parameter overlaps with `3` so it needs to be handled specially and
that it is the second leading term after `Q`.

The rest of the terms we enclose directly using
[`eval_expansion`](@ref). If `u0.ϵ` is small this gives good enough
enclosures so that we don't have to handle it in any more
sophisticated way.

## Handling `clausenc(x, 1 - α - u0.v0.v0.α + u0.v0.v0.p0)`
We are interested in bounding
```
-u0.v0.v0.a[1] * clausenc(x, 1 - α - u0.v0.v0.α + u0.v0.v0.p0) /
    (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0))
```
Notice the minus sign coming from the Hilbert transform. Let `r =
-u0.v0.v0.α + u0.v0.v0.p0 - 1`. Then `r > 0` and it is very small,
around `1e-8` or so depending on the precise choice of `u0.v0.v0`. We
have `1 - α - u0.v0.v0.α + u0.v0.v0.p0 = 2 - α + r` and the asymptotic
expansion of the Clausen term can then be written
```
clausenc(x, 2 - α - r) = gamma(α - 1 + r) * sinpi((2 - α + r) / 2) * x^(1 - α + r) -
    zeta(-α + r) / 2 * x^2 +
    O(x^4)
```
Ignoring the `O(x^4)` term for now and dividing by `x`(1 - α)` gives us
```
gamma(α - 1 + r) * sinpi((2 - α + r) / 2) * x^r - zeta(-α + r) / 2 * x^(1 + α)
```
Notice that the order of the terms depends on the value of `α`, in
some cases `x`r` is leading and in some cases `x^(1 - α)`.
- **TODO:** Finish computing a rigorous enclosure of this.
"""
function F0_nonzero(
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
    expansion_ispositive(u0, u0_expansion, ϵ) ||
        error("expansion of u0 not prove to be positive, this should not happen")

    # Compute the expansion of D(u0), skipping the Clausen term in the
    # tail corresponding to j = 1 and also remove the two leading
    # term, the three terms are handled separately.
    Du0_expansion = D(u0, AsymptoticExpansion(), skip_j_one = true; M, alpha_interval)(ϵ)
    delete!(Du0_expansion, (2, 0, 0, 0, 0, 0, 0))
    delete!(Du0_expansion, (0, 1, 0, 0, 0, 0, 0))

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
        F2 = let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, a0 = finda0(α)
            # Start by handling the terms P and Q

            # Enclosure of a0 / gamma(1 + α)
            a0gamma = Arb((finda0(α) / gamma(1 + α), -2 / Arb(π)^2))

            # Compute an enclosure of F21
            # FIXME: Compute a rigorous enclosure, this only computes
            # it for α = -1 + u0.ϵ. Though it seems to be quite stable.
            F21 =
                (c(2α - p0) - finda0(α) * c(α) * c(α - p0)) * x^(-α + p0 - 1) /
                (log(x) * (1 - x^p0))

            # Compute an enclosure of F22
            # FIXME: Compute a rigorous enclosure, this only computes
            # it for α = -1 + u0.ϵ. The value grows as α -> -1 so this
            # doesn't tell the full story.
            F22 =
                (
                    (zeta(-2α - 1) - zeta(-2α + p0 - 1)) +
                    finda0(α) * c(α - p0)^2 * x^(-2α + 2p0 - 2)
                ) / 2 * x^(1 + α) / (log(x) * (1 - x^p0))

            # The enclosure of the terms coming from P + Q in the expansion
            P_plus_Q = a0gamma * (F21 + F22)

            # Handle the term clausenc(x, 1 - α - u0.v0.v0.α + u0.v0.v0.p0)
            # FIXME: Compute a rigorous enclosure, this only computes
            # it for α = -1 + u0.ϵ.
            if u0.v0.v0.N0 > 0
                clausen_j_one = let r = -u0.v0.v0.α + u0.v0.v0.p0 - 1
                    s = 2 - α + r
                    term =
                        gamma(α - 1 - r) * sinpi((2 - α + r) / 2) * x^r -
                        zeta(-α + r) / 2 * x^(1 + α)
                    term /= gamma(1 + α) * log(x) * (1 - x^p0)

                    M = 3
                    _, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

                    p[2] = 0
                    remainder =
                        (p(x) + E * x^2M) / (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0))

                    -u0.v0.v0.a[1] * (term + remainder)
                end
            else
                clausen_j_one = zero(x)
            end

            # Enclosure of the remaining terms in the expansion
            # FIXME: The division by needs to be handled for the full
            # range of α
            remainder =
                eval_expansion(u0, Du0_expansion, x) /
                (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0))

            #(u0(x)^2 / 2 + Hu0x) / (log(x) * gamma(1 + α) * x^(1 - α) * (1 - x^p0))
            P_plus_Q + clausen_j_one + remainder
        end

        return F1 * F2
    end
end
