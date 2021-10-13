"""
    eval_expansion(u0::BHKdVAnsatz, expansion, x)

Evaluate the given expansion. It requires that `abs(x) < 1`.

The terms are stored as `((p, q, i, j, k, l, m), y)`. The parameters `(i, j,
k l, m)` correspond to the term
```
y * abs(x)^(i * α + j * p0 - k*u0-v0.v0.α + l*u0.v0.v0.p0 + m)
```
where `α ∈ (-1, -1 + u0.ϵ]` and `p0 = 1 + α + (1 + α)^2 / 2`.

The parameter `p` corresponds to multiplication by the factor
```
a0 * (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0) * abs(x)^-α
```
to the power `p`, which is the part of the expansion for the main term
which is not even powers of `x`. The parameter `q` corresponds to
multiplication by the factor
```
-a0 * (
    gamma(2α) * sinpi((1 - 2α) / 2) * abs(x)^(-2α) -
    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2) * abs(x)^(-2α + p0) +
    (-zeta(1 - 2α - 2) / 2 + zeta(1 - 2α + p0 - 2) / 2) * abs(x)^2
)
```
to the power `q`, which is the part of the expansion for `H` applied
to the main term which is not on the form `x^2m` for `m >= 2`.

# Handling `p`
As `α` goes to `-1` the factor corresponding to `p` converges to
```
-1 / π * abs(x) * log(abs(x)) - (γ - 1) / π * abs(x)
```
If we let
```
w1 = gamma(α) * sinpi((1 - α) / 2)
w2 = gamma(α - p0) * sinpi((1 - α + p0) / 2)
```
we can rewrite the factor as
```
a0 * (w1 - w2 * abs(x)^p0) * abs(x)^-α
```
We always have `-α > 1 - u0.ϵ` so we can factor out `abs(x)^(1 -
u0.ϵ)`. We then have to bound
```
a0 * (w1 - w2 * abs(x)^p0) * abs(x)^(-α - 1 + u0.ϵ)
```
For `abs(x)^(-α - 1 + u0.ϵ)` and upper bound is given by one and a
lower one by taking `α = -1`. For the remaining part we factor out
`w2`, giving us
```
w2 * a0 * (w1 / w2 - abs(x)^p0)
```
An enclosure of `w2` can be computed by using that it is `-π / 2` at
`α = -1` and increasing.
- **PROVE:** Prove the limit for `w2` and that it is increasing
- **TODO:** Figure out how to handle `a0 * (w1 / w2 - x^p0). It should
  be similar to what is being done in the asymptotic `F0`.
- Alternatively error bounds for the limit `-1 / π * x * log(x) - (γ -
  1) / π * x` could be computed.

# Handling `q`
As `α` goes to `-1` the factor corresponding to `q` converges to
```
u0.v0.a0 * (
    -1 // 4 * abs(x)^2 * log(abs(x))^2 + (3 // 4 - γ / 2) * abs(x)^2 * log(abs(x))
    - (3γ - γ^2 - 2γ₁ - 7 // 2 + π^2 / 12) / 4 * abs(x)^2
)
```
where `γ₁ = stieltjes(Arb, 1)`.
- **TODO:** In the end we might not need this factor since it is
    handled specially by [`F0`](@ref). But otherwise we would need to
    compute error bounds for this factor, preferably ones that depend
    on `x`
"""
function eval_expansion(
    u0::BHKdVAnsatz{Arb},
    expansion::AbstractDict{NTuple{7,Int},Arb},
    x;
    div_a0 = false,
)
    @assert abs(x) < 1

    res = zero(x)

    # Enclosure of α
    α = Arb((-1, -1 + u0.ϵ))

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
    γ₁ = stieltjes(Arb, 1)

    for ((p, q, i, j, k, l, m), y) in expansion
        if !iszero(y)
            term = y

            exponent = exponent_α_p0_m(i, j, m) - k * u0.v0.v0.α + l * u0.v0.v0.p0

            if !iszero(p)
                # Factor out x^(1 - u0.ϵ) and add to exponent
                exponent += 1 - u0.ϵ

                # This currently uses the limiting expression but
                # without adding explicit error bounds
                if Arblib.contains_zero(x)
                    # Use monotonicity
                    # TODO: Determine the interval of monotonicity
                    lower = zero(x)
                    upper = let x = ubound(Arb, x)
                        (-inv(π) * log(abs(x)) - (γ - 1) / π) * abspow(x, u0.ϵ)
                    end

                    factor = Arb((lower, upper))
                else
                    factor = (-inv(π) * log(abs(x)) - (γ - 1) / π) * abspow(x, u0.ϵ)
                end

                @warn "Error term not implemented for p = $p" maxlog = 1

                term *= factor^p
            end

            if !iszero(q)
                # This currently uses the limiting expression but
                # without adding error bounds. It is also currently
                # not needed anywhere.
                factor =
                    u0.v0.a0 * (
                        -1 // 4 * abs(x)^2 * log(abs(x))^2 +
                        (3 // 4 - γ / 2) * abs(x)^2 * log(abs(x)) +
                        (3γ - γ^2 - 2γ₁ - 7 // 2 + π^2 / 12) / 4 * abs(x)^2
                    )

                @warn "Error term not implemented for q = $q"

                term *= factor^q
            end

            term *= abspow(x, exponent)

            res += term
        end
    end

    return res
end

"""
    (u0::BHKdVAnsatz)(x, ::Ball)

Evaluate the ansatz `u0` at the point `x` using direct ball arithmetic
(not an asymptotic approach).

The tail term is evaluated directly.

To evaluate the main term, given by
```
a0 * (Ci(x, 1 - α) - Ci(x, 1 - α + p0) - (zeta(1 - α) - zeta(1 - α + p0)))
```
we make use of the fact that this converges to
```
u0.v0.a0 * (Ci(x, 2, 1) - zeta(2, d = 1))
```
, which is the main term for `BHAnsatz`, as `α -> -1`. We therefore
evaluate this function and add precomputed \$L^\\infty\$ bounds for
```
a0 * (Ci(x, 1 - α) - Ci(x, 1 - α + p0) - (zeta(1 - α) - zeta(1 - α + p0))) - u0.v0.a0 * (Ci(x, 2, 1) - zeta(2, d = 1))
```
valid for the entire range `α ∈ (-1, -1 + u0.ϵ]`.

**TODO:** Compute rigorous \$L^\\infty\$ bounds for the above
  expression. We currently use a heuristic value.
"""
function (u0::BHKdVAnsatz{Arb})(x, ::Ball)
    # Main term

    # Approximation
    res = 2 / Arb(π)^2 * (Ci(x, 2, 1) - zeta(Arb(2), d = 1))

    # Add error bounds
    @warn "L^∞ bounds not rigorously computed - using heuristic values" maxlog = 1
    if u0.ϵ <= 0.0003
        Arblib.add_error!(res, Mag(2e-4))
    else
        throw(ErrorException("no L^∞ bounds for ϵ = $ϵ"))
    end

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.v0.N0
        s = 1 - u0.v0.v0.α + j * u0.v0.v0.p0
        res += u0.v0.v0.a[j] * (Ci(x, s) - zeta(s))
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
u0.v0.a0 * (-1)^m * zeta(2 - 2m, d = 1) / factorial(2m)
```
which is the coefficient in front of `x^2m` for the main term of
`u0.v0`. We therefore compute the coefficients by using this
expression and then adding an error term for it.
- **TODO:** What error term should we add? At the moment it is hard
    coded.
- **TODO:** How should we handle the error term from the expansion?

The only remaining part of the expansion of the main term is
```
a0 * (gamma(α) * sinpi((1 - α) / 2) - gamma(α - p0) * sinpi((1 - α + p0) / 2) * x^p0) * x^-α
```
which we don't evaluate at all yet. This we don't evaluate but instead
store implicitly in the expansion.

See [`eval_expansion`](@ref) for more details.
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
    for m = 1:M-1
        coefficient = (-1)^m * zeta(Arb(2 - 2m), d = 1) / factorial(Arb(2m)) * u0.v0.a0
        if 2m == 2
            Arblib.add_error!(coefficient, Mag(2e-5))
        elseif 2m == 4
            Arblib.add_error!(coefficient, Mag(4e-8))
        else
            throw(ArgumentError("no error bound computed for m = $m"))
        end

        res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
    end

    # Error term for main term
    Arblib.add_error!(
        res[(0, 0, 0, 0, 0, 0, 2M)],
        2abs(zeta(Arb(2 - 2M), d = 1) / factorial(Arb(2M))) * u0.v0.a0,
    )

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.v0.N0
        C, _, p, E = Ci_expansion(x, 1 - u0.v0.v0.α + j * u0.v0.v0.p0, M)
        res[(0, 0, 0, 0, 1, j, 0)] = C * u0.v0.v0.a[j]
        for m = 1:M-1
            res[(0, 0, 0, 0, 0, 0, 2m)] += p[2m] * u0.v0.v0.a[j]
        end
        Arblib.add_error!(res[(0, 0, 0, 0, 0, 0, 2M)], E)
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
-a0 * (Ci(x, 1 - 2α) - Ci(x, 1 - 2α + p0) - (zeta(1 - 2α) - zeta(1 - 2α + p0)))
```
we make use of the fact that this converges to
```
-u0.v0.a0 * (Ci(x, 3, 1) - zeta(3, d = 1))
```
, which is the main term for `BHAnsatz`, as `α -> -1`. We therefore
evaluate this function and add precomputed \$L^\\infty\$ bounds for
```
-a0 * (Ci(x, 1 - 2α) - Ci(x, 1 - 2α + p0) - (zeta(1 - 2α) - zeta(1 - 2α + p0))) + u0.v0.a0 * (Ci(x, 3, 1) - zeta(3, d = 1))
```
valid for the entire range `α ∈ (-1, -1 + u0.ϵ]`.

**TODO:** Compute rigorous \$L^\\infty\$ bounds for the above
  expression. We currently use a heuristic value.

For the tail term we need to make sure that we correctly handle the
fact that the transform depends on the value of `α`.

For the Fourier terms we do this directly, the transformation takes
`cos(n * x)` to `-n^α * cos(n * x)` and in this case we just let `α`
be a ball containing `(-1, -1 + u0.ϵ]`.

For the Clausen functions we have to be a bit more careful. The
transformation takes `Ci(x, s)` to `-Ci(x, s - α)` but putting `α` as
a ball doesn't give good enclosures. It gives a large overestimations,
in particular when `s - α` is close to an integer, and fails when it
overlaps with an integer. Instead we use the approximation given by
taking the Hilbert transform and bounding the error. The Hilbert
transform of `Ci(x, s)` is given by `-Ci(x, s + 1)`. For the error we
currently use a heuristic, for `s = 2` and `ϵ = 0.0003` the error
maximum error in `x` is around `3e-4` and for larger values of `s` we
get a smaller error.

**TODO:** How to properly bound the error when using `-Ci(x, s + 1)`
instead of `-Ci(x, s - α)`?

**TODO:** We will probably have to improve on the enclosure to get a
sufficiently good defect in the end.

"""
function H(u0::BHKdVAnsatz{Arb}, ::Ball)


    @warn "L^∞ bounds for main term not rigorously computed - using heuristic values" maxlog =
        1
    @warn "L^∞ bounds for Clausen tail term not rigorously computed - using heuristic values" maxlog =
        1

    return x -> begin

        # Main term

        # Approximation
        res = -u0.v0.a0 * (Ci(x, 3, 1) - zeta(Arb(3), d = 1))

        # Add error bounds
        if u0.ϵ <= 0.0003
            Arblib.add_error!(res, Mag(1e-4))
        else
            throw(ErrorException("no L^∞ bounds for ϵ = $ϵ"))
        end

        # Tail term

        # Clausen terms
        for j = 1:u0.v0.v0.N0
            s = 2 - u0.v0.v0.α + j * u0.v0.v0.p0
            term = Ci(x, s) - zeta(s)
            if s >= 2 && u0.ϵ <= 0.0003
                Arblib.add_error!(term, Mag(3e-4))
            else
                # Error bounds only valid for s >= 2
                throw(ErrorException("no L^∞ bounds for s = $s and ϵ = $ϵ"))
            end
            res -= u0.v0.v0.a[j] * term
        end

        # Fourier terms
        let α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α
            for n = 1:u0.v0.N
                res -= u0.v0.b[n] * n^α * (cos(n * x) - 1)
            end
        end

        return res
    end
end

function H(u0::BHKdVAnsatz, ::Asymptotic; M::Integer = 3)
    f = H(u0, AsymptoticExpansion(); M)

    return x -> eval_expansion(u0, f(x), x)
end

"""
    H(u0::BHKdVAnsatz, ::AsymptoticExpansion; M = 3)

Return a dictionary containing the terms in the asymptotic expansion
of `u0` which can then be evaluated with [`eval_expansion`](@ref).

The highest term, `x^2M`, is an error term is which makes sure that
evaluation of the expansion gives an enclosure of the result.

For the main term the coefficients in front of `x^2m` is given by
```
a0 * (-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / factorial(2m)
```
which in the limit becomes `∞ * 0`. For `m >= 2` it converges to
```
u0.v0.a0 * (-1)^m * zeta(2 - 2m, d = 1) / factorial(2m)
```
which is the coefficient in front of `x^2m` for the main term of
`u0.v0`. We therefore compute the coefficients by using this
expression and then adding an error term for it. For `m = 1` we can't
do this, instead we include it in the remaining terms below.
- **TODO:** What error term should we add? At the moment it is hard
    coded.
- **TODO:** How should we handle the error term from the expansion?

The remaining part of the expansion of `H` applied to the main term is
```
-a0 * (
    gamma(2α) * sinpi((1 - 2α) / 2) * abs(x)^(-2α) -
    gamma(2α - p0) * sinpi((1 - 2α + p0) / 2) * abs(x)^(-2α + p0)
    (zeta(1 - 2α - 2) / 2 - zeta(1 - 2α + p0 - 2) / 2) * abs(x)^2
)
```
which we don't evaluate at all yet. This we don't evaluate but instead
store implicitly in the expansion.

For the tail term the Fourier terms are handled directly by letting
`α` be a ball. For the Clausen terms this gives very bad enclosures.
Instead we take `α = -1` in this case and compute error bounds.
- **TODO:** How to compute these error bounds? Is this even the right
    approach?

See [`eval_expansion`](@ref) for more details.
"""
function H(u0::BHKdVAnsatz{Arb}, ::AsymptoticExpansion; M::Integer = 3)
    @assert M >= 3

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
        @warn "No error bounds for coefficients of main term" maxlog = 1
        for m = 2:M-1
            coefficient = -(-1)^m * zeta(Arb(3 - 2m), d = 1) / factorial(Arb(2m)) * u0.v0.a0

            res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
        end

        # Error term for main term
        @warn "No error bounds for error term of main term" maxlog = 1
        Arblib.add_error!(
            res[(0, 0, 0, 0, 0, 0, 2M)],
            2abs(zeta(Arb(3 - 2M), d = 1) / factorial(Arb(2M))) * u0.v0.a0,
        )

        # Tail term
        α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α

        # Clausen terms
        @warn "No error bounds for Clausen terms in tail" maxlog = 1
        for j = 1:u0.v0.v0.N0
            s = 2 - u0.v0.v0.α + j * u0.v0.v0.p0
            C, _, p, E = Ci_expansion(x, s, M)

            res[(0, 0, -1, 0, 1, j, 0)] = -C * u0.v0.v0.a[j]
            for m = 1:M-1
                res[(0, 0, 0, 0, 0, 0, 2m)] -= p[2m] * u0.v0.v0.a[j]
            end
            Arblib.add_error!(res[(0, 0, 0, 0, 0, 0, 2M)], E)
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

function D(u0::BHKdVAnsatz, ::Asymptotic; M::Integer = 3)
    f = D(u0, AsymptoticExpansion(); M)

    return x -> eval_expansion(u0, f(x), x)
end

function D(u0::BHKdVAnsatz, evaltype::AsymptoticExpansion; M::Integer = 3)
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
that we have to handle. With `u0.w(x) = x * log(10 + inv(x))` this
simplifies to
```
x^(-1 - α + p0) / log(10 + inv(x))
```
Since `-1 - α + p0 > 0` for all values of `α` this can be upper
bounded by `inv(log(10 + inv(x))). For a lower bound we use that
```
-1 - α + p0 = (1 + α)^2 / 2 < u0.ϵ^2 / 2
```

# Simplify numerator
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

# Split numerator
The next step is to split `D(u0)(x) / x^(-2α + p0)` into three parts,
the term corresponding to `(2, 0, 2, -1, 0, 0, 0)`. which we will call
`P`, the term corresponding to `(0, 1, 2, -1, 0, 0, 0)`, which we will
call `Q`, and the remaining terms, which we will call `R`.

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
Factoring out `v2 * x^(2 + 2α - p0)` from the numerator and `w2` from
the denominator we get
```
v2 / w2 * x^(2 + 2α - p0) * (v1 / v2 - x^(-2 - 2α + 2p0)) / (w1 / w2 - x^p0)
```
**FIXME:** For now we assume that this is contained in the interval
`[0, inv(2π)]` which can be checked numerically to be the case.
However I don't know how to prove this as of yet. I started doing the
following but it **doesn't work due to the fact that `v2 / w2`
diverges**. We could notice that
1. `w1 / w2 -> 1` from above
2. `v2 / v1` -> 1` from above
3. `w1 / w2 > v2 / v1`
4. `0 < -2 - 2α + 2p0 < p0`
All of this combined allows us to conclude that
```
0 < (v2 / v1 - x^(-2 - 2α + 2p0)) / (w1 / w2 - x^p0) < 1
```
for all values of `α` and `x`. But `v2 / w2` still diverges :(

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
    gamma = SpecialFunctions.gamma

    @assert ϵ < 1

    # This uses a hard coded version of the weight so just as an extra
    # precaution we check that it seems to be the same as the one
    # used.
    let x = Arb(0.5)
        @assert isequal(u0.w(x), x * log(10 + inv(x)))
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

    # x^(-2α + p0) / (u0.w(x) * x^-α) = x^(-1 - α + p0) / log(10 + inv(x))
    weight(x) = begin
        iszero(x) && return zero(x)

        if Arblib.contains_zero(x)
            # Use that it is zero at x = 0 and monotonically
            # increasing in x

            # Lower bound is zero
            lower = zero(x)

            # Upper bound is inv(log(10 + inv(x))) evaluated at upper
            # bound for x
            upper = inv(log(10 + inv(ubound(Arb, x))))

            return Arb((lower, upper))
        end

        # (-1 - α + p0) is upper bounded by u0.ϵ^2 / 2
        lower = x^(u0.ϵ^2 / 2) / log(10 + inv(x))

        upper = inv(log(10 + inv(x)))

        return Arb((lower, upper))
    end

    return x -> begin
        @assert (x isa Arb && x <= ϵ) || (x isa ArbSeries && x[0] <= ϵ)

        # TODO: Prove that u0_expansion_div_xmα is non-negative, so
        # that we can neglect it.

        # Upper bound of the part corresponding to P + Q
        res1 = begin
            # Upper bound of the constant part, consisting of the
            # first term
            part1 = let α = -1 + u0.ϵ
                a0 = finda0(α)
                p0 = 1 + α + (1 + α)^2 / 2

                w1 = gamma(α) * sinpi((1 - α) / 2)

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
            part2 = Arb((0, inv(2Arb(π))))

            part1 + part2
        end

        # Upper bound of the remaining part
        res2 = begin
            numerator = eval_expansion(u0, R, x)

            # Compute an enclosure of w1 using the limiting value at α
            # = -1 being -π / 2 and the monotonicity in α
            # PROVE: The monotonicity in α
            w1 = let α = -1 + u0.ϵ
                lower = gamma(α) * sinpi((1 - α) / 2)
                upper = -Arb(π) / 2

                Arb((lower, upper))
            end

            # TODO: Lower bound of w1 * a0 * (1 / x^p0 - 1). This is an
            # approximation but we need to bound the error terms.
            denominator = w1 * 2log(ubound(Arb, x)) / Arb(π)^2

            numerator / denominator
        end

        return (res1 + res2) * weight(x)
    end
end
