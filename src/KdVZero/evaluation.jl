"""
    eval_expansion(u0::KdVZeroAnsatz, expansion, x)

Evaluate the given expansion. The term `((i, j, m), y)` corresponds to
```
y*abs(x)^(-i * α + j * p0 + m)
```
"""
function eval_expansion(
    u0::KdVZeroAnsatz,
    expansion::AbstractDict{NTuple{3,Int},ArbSeries},
    x::Arb;
    offset_i::Integer = 0,
    offset_m::Integer = 0,
)
    α = ArbSeries((0, 1); u0.degree)

    res = zero(α)

    for ((i, j, m), y) in expansion
        if !iszero(y)
            if iszero(i + offset_i) && iszero(j) && iszero(m + offset_m)
                # Exponent is identically equal to zero
                term = one(res)
            else
                exponent = -(i + offset_i) * α + j * u0.p0 + (m + offset_m)
                term = abspow_with_remainder(x, exponent, u0.α)
            end
            res += mul_with_remainder(y, term, u0.α)
        end
    end

    return res
end

"""
    (u0::KdVZeroAnsatz)(x::Arb, ::Ball)

Return an expansion of `u0(x)` in `α` around `α = 0`. The last term is
a remainder term which ensures that evaluating the expansion gives an
enclosure of `u0` for all values in the interval `α`.

The value is given by
```
sum(a[j] * clausencmzeta(x, 1 - α + j * p0) for j = 0:2)
```

For `j = 1, 2` we can evaluate it directly with
[`clausencmzeta`](@ref). For `j = 0` this doesn't work directly and we
have to do it manually.

We are thus interested in computing an expansion around `α = 0` of
```
a[0] * clausencmzeta(x, 1 - α) = a[0] * clausenc(x, 1 - α) - a[0] * zeta(1 - α)
```
for a fixed `x`. Using the formula
```
clausenc(x, s) = let v = 1 - s
    gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
end
```
from [`_clausenc_zeta`](@ref) we can write the Clausen function as
```
clausenc(x, 1 - α) = gamma(α) * inv(2π)^α * cospi(α / 2) * (zeta(α, x / 2π) + zeta(α, 1 - x / 2π))
```
Everything except `gamma(α)`, which has a pole of order `1` at `α =
0`, can be evaluated directly. We can also notice that the constant
term in `zeta(α, x / 2π) + zeta(α, 1 - x / 2π)` is given by
```
zeta(0, x / 2π) + zeta(0, 1 - x / 2π) = (1 / 2 - x / 2π) + (1 / 2 - (1 - x / π)) = 0
```
where we have used [Equation
25.11.13](https://dlmf.nist.gov/25.11.E13). We can thus extract a
factor `α` from this. Putting the `α` together with `gamma(α)` we can
write it as `α * gamma(α) = α / rgamma(α)` where we can handle the
removable singularity
- **TODO:** Compute remainder term in `α` for `α * gamma(α)`.

For wide values of `x` direct evaluation of `zeta(α, x / 2π) + zeta(α,
1 - x / 2π)` gives a poor enclosure. To get better enclosures we use
that all derivatives of that `zeta(α, x / 2π) + zeta(α, 1 - x / 2π)`
are monotone on the interval `(0, π)`. We don't do this for the
remainder term because it doesn't give much of an improvement
- **PROVE:** That all derivatives of `zeta(α, x / 2π) + zeta(α, 1 - x
  / 2π)` are monotone on `(0, π)`.

For computing `a[0] * zeta(1 - α)` we rewrite it as
```
(a[0] / α) * (α * zeta(1 - α))
```
and use the expansion
```
zeta(s) = 1 / (s - 1) + sum((-1)^n * stieltjes(n) / factorial(n) * (s - 1)^n for n = 0:Inf)
```
to get
```
α * zeta(1 - α) = -1 + sum(stieltjes(n) / factorial(n) * α^n for n = 0:Inf)
```
- **TODO:** Compute remainder term in `α`
"""
function (u0::KdVZeroAnsatz)(x::Arb, ::Ball)
    # The main term we handle manually
    res = let
        # Compute this to one degree higher since we divide it by α,
        # reducing the degree by one.
        clausen_term = let π = Arb(π), α = ArbSeries((0, 1), degree = u0.degree + 1)
            if iswide(x) && 0 < x < π
                z = ArbSeries(
                    union.(
                        Arblib.coeffs(
                            zeta(α, lbound(Arb, x) / 2π) + zeta(α, 1 - lbound(Arb, x) / 2π),
                        ),
                        Arblib.coeffs(
                            zeta(α, ubound(Arb, x) / 2π) + zeta(α, 1 - ubound(Arb, x) / 2π),
                        ),
                    ),
                )
            else
                z = zeta(α, x / 2π) + zeta(α, 1 - x / 2π)
            end
            inv(2π)^α * cospi(α / 2) * z
        end

        # Compute remainder term
        clausen_term_remainder =
            let π = Arb(π), α = ArbSeries((u0.α, 1), degree = u0.degree + 1)
                (inv(2π)^α*cospi(α / 2)*(zeta(α, x / 2π)+zeta(α, 1 - x / 2π)))[u0.degree+1]
            end
        clausen_term[u0.degree+1] = clausen_term_remainder

        # The constant term is exactly zero
        @assert Arblib.contains_zero(clausen_term[0])
        clausen_term[0] = 0

        # Divide clausen term by α
        clausen_term = clausen_term << 1

        # Expansion of α * gamma(α)
        gammamulα = inv(rgamma(ArbSeries((0, 1), degree = u0.degree + 1)) << 1)

        # FIXME: Properly implement this. Now we just widen the last
        # coefficient so that we get an enclosure for a lower bound of α
        if !iszero(u0.α)
            error = let α = lbound(Arb, u0.α)
                α * gamma(α) - gammamulα(α)
            end
            gammamulα[u0.degree] +=
                Arblib.add_error!(zero(error), error / lbound(Arb, u0.α)^u0.degree)
        end

        clausen_term = mul_with_remainder(clausen_term, gammamulα, u0.α)

        # Divide a[0] by α, perform the multiplication and then
        # multiply by α. This makes sure the degree after the
        # multiplication is correct.
        a0clausenterm = mul_with_remainder(u0.a[0] << 1, clausen_term, u0.α) >> 1

        zetamulα =
            ArbSeries([-1; [stieltjes(Arb, n) / factorial(n) for n = 0:u0.degree-1]])

        # FIXME: Properly implement this. Now we just widen the last
        # coefficient so that we get an enclosure for a lower bound of α
        if !iszero(u0.α)
            error = let α = lbound(Arb, u0.α)
                α * zeta(1 - α) - zetamulα(α)
            end
            zetamulα[u0.degree] +=
                Arblib.add_error!(zero(error), error / lbound(Arb, u0.α)^u0.degree)
        end

        a0zeta_term = mul_with_remainder(u0.a[0] << 1, zetamulα, u0.α)

        a0clausenterm - a0zeta_term
    end

    let α = ArbSeries((0, 1); u0.degree)
        for j = 1:2
            # term = clausencmzeta(x, 1 - α + j * u0.p0)
            term = clausencmzeta_with_remainder(x, 1 - α + j * u0.p0, u0.α)

            # res += u0.a[j] * term
            res += mul_with_remainder(u0.a[j], term, u0.α)
        end
    end

    return res
end

"""
    (u0::KdVZeroAnsatz)(x::Arb, ::AsymptoticExpansion; M::Integer)

Return an expansion of `u0(x)` in `x` around zero where the
coefficients in the expansion are themselves expansions in `α` around
zero. The last term in each of these expansions in `α` is a remainder
term which ensures that evaluating the expansion gives an enclosure of
the term for all values in the interval `α`.

It returns a dictionary `expansion` where the keys are three tuples
`(i, j, m)` and correspond to a term of the form
```
expansion[(i, j, m)] * abs(x)^(i * α + j * p0 + m)
```

The value of `M` determines the number of terms in the expansion in
`x`.

Most of the terms can be computed by evaluating them directly. The
exception is the singular term for the first Clausen function, given
by `gamma(α) * cospi(α / 2)` which has a singularity at `α = 0`.
However multiplication by `a[0]` removes the singularity and we have
the expansion
```
a[0] * gamma(α) * cospi(α / 2) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```
which is the same as the one occurring in [`expansion_p0`](@ref). The
expansion is given by
```
2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2)) =
    1 - γ * α + (γ^2 / 2 - π^2 / 8) * α^2 +
    (-4γ^3 + 3γ * π^2 + 28 * polygamma(2, 1)) / 24 * α^3 + O(α^4)
```
where `γ` is the Euler constant.
- **TODO:** Compute remainder term in `α`
"""
function (u0::KdVZeroAnsatz)(x::Arb, ::AsymptoticExpansion; M::Integer = 10)
    α = ArbSeries((0, 1); u0.degree)

    expansion = OrderedDict{NTuple{3,Int},ArbSeries}()

    # Initiate even powers of x
    for m = 1:M
        expansion[(0, 0, 2m)] = ArbSeries(; u0.degree)
    end

    # Handle main term
    # Compute the coefficient for the singular term
    a0singular_term = let γ = Arb(Irrational{:γ}()), π = Arb(π)
        ArbSeries((1, -γ, γ^2 / 2 - π^2 / 8); u0.degree)
    end
    # FIXME: Properly implement this. Now we just widen the last
    # coefficient so that we get an enclosure for a lower bound of α
    if !iszero(u0.α)
        error = let α = lbound(Arb, u0.α)
            gamma(α) * sinpi((1 - α) / 2) * finda0(α) - a0singular_term(α)
        end
        a0singular_term[u0.degree] +=
            Arblib.add_error!(zero(error), error / lbound(Arb, u0.α)^u0.degree)
    end

    expansion[(1, 0, 0)] = a0singular_term

    # Compute the coefficients for the analytic terms
    for m = 1:M-1
        term = (-1)^m * zeta(1 - α - 2m) / factorial(2m)
        term_remainder = let α = ArbSeries((u0.α, 1), degree = u0.degree)
            ((-1)^m*zeta(1 - α - 2m)/factorial(2m))[u0.degree]
        end
        term[u0.degree] = term_remainder

        expansion[(0, 0, 2m)] += mul_with_remainder(u0.a[0], term, u0.α)
    end

    # Add remainder term (in x)
    remainder_term = clausenc_expansion_remainder(x, 1 - α, M)
    expansion[(0, 0, 2M)] += mul_with_remainder(u0.a[0], remainder_term, u0.α)

    # Handle tail terms
    for j = 1:2
        s = 1 - α + j * u0.p0

        # Compute the coefficient for the singular term
        singular_term = compose_with_remainder(s -> gamma(1 - s) * sinpi(s / 2), s, u0.α)
        expansion[(1, j, 0)] = mul_with_remainder(u0.a[j], singular_term, u0.α)

        # Compute the coefficients for the analytic terms
        for m = 1:M-1
            term =
                (-1)^m * compose_with_remainder(s -> zeta(s - 2m), s, u0.α) / factorial(2m)
            expansion[(0, 0, 2m)] += mul_with_remainder(u0.a[j], term, u0.α)
        end

        # Add remainder term
        remainder_term = clausenc_expansion_remainder(x, s, M)
        expansion[(0, 0, 2M)] += mul_with_remainder(u0.a[j], remainder_term, u0.α)
    end

    return expansion
end

"""
    H(u0::KdVZeroAnsatz, ::Ball)

Return a function such that `H(u0)(x)` computes an expansion in `α`
around `α = 0`. The last term is a remainder term which ensures that
evaluating the expansion gives an enclosure of `u0` for all values in
the interval `α`.

The value is given by
```
-sum(a[j] * clausencmzeta(x, 1 - 2α + j * p0) for j = 0:2)
```

For `j = 1, 2` we can evaluate it directly with
[`clausencmzeta`](@ref). For `j = 0` this doesn't work directly and we
have to do it manually.

We are thus interested in computing an expansion around `α = 0` of
```
a[0] * clausencmzeta(x, 1 - 2α) = a[0] * clausenc(x, 1 - 2α) - a[0] * zeta(1 - 2α)
```
for a fixed `x`. Using the formula
```
clausenc(x, s) = let v = 1 - s
    gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
end
```
from [`_clausenc_zeta`](@ref) we can write the Clausen function as
```
clausenc(x, 1 - 2α) = gamma(2α) * inv(2π)^2α * cospi(α) * (zeta(2α, x / 2π) + zeta(2α, 1 - x / 2π))
```
Everything except `gamma(2α)`, which has a pole of order `1` at `α =
0`, can be evaluated directly. We can also notice that the constant
term in `zeta(2α, x / 2π) + zeta(2α, 1 - x / 2π)` is given by
```
zeta(0, x / 2π) + zeta(0, 1 - x / 2π) = (1 / 2 - x / 2π) + (1 / 2 - (1 - x / π)) = 0
```
where we have used [Equation
25.11.13](https://dlmf.nist.gov/25.11.E13). We can thus extract a
factor `α` from this. Putting the `α` together with `gamma(2α)` we can
write it as `α * gamma(2α) = α / rgamma(2α)` where we can handle the
removable singularity
- **TODO:** Compute remainder term in `α` for `α * gamma(2α)`.

remove the pole and we have the expansion
```
α * gamma(2α) = 1 / 2 - γ * α + (γ^2 + π^2 / 6) * α^2 + O(α^3)
```
- **TODO:** Compute remainder term in `α`

For wide values of `x` direct evaluation of `zeta(2α, x / 2π) +
zeta(2α, 1 - x / 2π)` gives a poor enclosure. To get better enclosures
we use that all derivatives of that `zeta(2α, x / 2π) + zeta(2α, 1 - x
/ 2π)` are monotone on the interval `(0, π)`. We don't do this for the
remainder term because it doesn't give much of an improvement
- **PROVE:** That all derivatives of `zeta(2α, x / 2π) + zeta(2α, 1 -
  x / 2π)` are monotone on `(0, π)`.

For computing `a[0] * zeta(1 - 2α)` we rewrite it as
```
(a[0] / α) * (α * zeta(1 - 2α))
```
and use the expansion
```
zeta(s) = 1 / (s - 1) + sum((-1)^n * stieltjes(n) / factorial(n) * (s - 1)^n for n = 0:Inf)
```
to get
```
α * zeta(1 - 2α) = -1 / 2 + sum(2^n * stieltjes(n) / factorial(n) * α^n for n = 0:Inf)
```
- **TODO:** Compute remainder term in `α`
"""
function H(u0::KdVZeroAnsatz, ::Ball)
    u0.degree <= 2 || throw(ArgumentError("only supports degree up to 2"))

    # Expansion of α * gamma(2α)
    gammamulα = inv(rgamma(2ArbSeries((0, 1), degree = u0.degree + 1)) << 1)

    # FIXME: Properly implement this. Now we just widen the last
    # coefficient so that we get an enclosure for a lower bound of α
    if !iszero(u0.α)
        error = let α = lbound(Arb, u0.α)
            α * gamma(2α) - gammamulα(α)
        end
        gammamulα[u0.degree] +=
            Arblib.add_error!(zero(error), error / lbound(Arb, u0.α)^u0.degree)
    end

    # Expansion of α * zeta(1 - 2α)
    zetamulα =
        ArbSeries([-1 // 2; [2^n * stieltjes(Arb, n) / factorial(n) for n = 0:u0.degree-1]])

    # FIXME: Properly implement this. Now we just widen the last
    # coefficient so that we get an enclosure for a lower bound of α
    if !iszero(u0.α)
        error = let α = lbound(Arb, u0.α)
            α * zeta(1 - 2α) - zetamulα(α)
        end
        zetamulα[u0.degree] +=
            Arblib.add_error!(zero(error), error / lbound(Arb, u0.α)^u0.degree)
    end

    return x::Arb -> begin
        # The main term we handle manually
        res = let
            clausen_term =
                let π = Arb(π), α = ArbSeries((0, 1), degree = u0.degree + 1)
                    if iswide(x) && 0 < x < π
                        z = ArbSeries(
                            union.(
                                Arblib.coeffs(
                                    zeta(2α, lbound(Arb, x) / 2π) +
                                    zeta(2α, 1 - lbound(Arb, x) / 2π),
                                ),
                                Arblib.coeffs(
                                    zeta(2α, ubound(Arb, x) / 2π) +
                                    zeta(2α, 1 - ubound(Arb, x) / 2π),
                                ),
                            ),
                        )
                    else
                        z = zeta(2α, x / 2π) + zeta(2α, 1 - x / 2π)
                    end
                    inv(2π)^2α * cospi(α) * z
                end

            # Compute remainder term
            clausen_term_remainder =
                let π = Arb(π), α = ArbSeries((u0.α, 1), degree = u0.degree + 1)
                    (inv(2π)^2α*cospi(α)*(zeta(2α, x / 2π)+zeta(2α, 1 - x / 2π)))[u0.degree+1]
                end
            clausen_term[u0.degree+1] = clausen_term_remainder

            # The constant term is exactly zero
            @assert Arblib.contains_zero(clausen_term[0])
            clausen_term[0] = 0

            # Divide clausen term by α
            clausen_term = clausen_term << 1

            clausen_term = mul_with_remainder(clausen_term, gammamulα, u0.α)

            # Divide a[0] by α, perform the multiplication and then
            # multiply by α. This makes sure the degree after the
            # multiplication is correct.
            a0clausenterm = mul_with_remainder(u0.a[0] << 1, clausen_term, u0.α) >> 1

            a0zeta_term = mul_with_remainder(u0.a[0] << 1, zetamulα, u0.α)

            -(a0clausenterm - a0zeta_term)
        end

        let α = ArbSeries((0, 1); u0.degree)
            for j = 1:2
                # term = clausencmzeta(x, 1 - α + j * u0.p0)
                term = clausencmzeta_with_remainder(x, 1 - 2α + j * u0.p0, u0.α)

                # res += u0.a[j] * term
                res -= mul_with_remainder(u0.a[j], term, u0.α)
            end
        end

        return res
    end
end

"""
    H(u0::KdVZeroAnsatz, ::AsymptoticExpansion)

Return a function such that `H(u0)(x)` computes an expansion in `x`
around zero there the coefficients are themselves expansion in `α`
around zero. The last term in each of these expansions in `α` is a
remainder term which ensures that evaluating the expansion gives an
enclosure of the term for all values in the interval `α`.

It returns a dictionary `expansion` where the keys are three tuples
`(i, j, m)` and correspond to a term of the form
```
expansion[(i, j, m)] * abs(x)^(i * α + j * p0 + m)
```

The value of `M` determines the number of terms in the expansion in
`x`.

Most of the terms can be computed by evaluating them directly. The
exception is the singular term for the first Clausen function, given
by `gamma(2α) * cospi(α)` which has a singularity at `α = 0`. However
multiplication by `a[0]` removes the singularity and we have
```
a[0] * gamma(2α) * cospi(α) = 2gamma(2α)^2 * cospi(α)^2 / (gamma(α)^2 * cospi(α / 2)^2)
```
From Mathematica the expansion is given by
```
2gamma(2α)^2 * cospi(α)^2 / (gamma(α)^2 * cospi(α / 2)^2) =
    1 / 2 - γ * α + (γ^2 - π^2 / 8) * α^2 +
    (-8γ^3 + 3γ * π^2 + 14polygamma(2, 1)) / 12 * α^3
```
where `γ` is the Euler constant.
- **TODO:** Compute remainder term in `α`
"""
function H(u0::KdVZeroAnsatz, ::AsymptoticExpansion; M::Integer = 10)
    α = ArbSeries((0, 1); u0.degree)

    return x::Arb -> begin
        expansion = OrderedDict{NTuple{3,Int},ArbSeries}()

        # Initiate even powers of x
        for m = 1:M
            expansion[(0, 0, 2m)] = ArbSeries(; u0.degree)
        end

        # Handle main term
        # Compute the coefficient for the singular term
        a0singular_term = let γ = Arb(Irrational{:γ}()), π = Arb(π)
            ArbSeries((1 // 2, -γ, γ^2 - π^2 / 8); u0.degree)
        end
        # FIXME: Properly implement this. Now we just widen the last
        # coefficient so that we get an enclosure for a lower bound of α
        if !iszero(u0.α)
            error = let α = lbound(Arb, u0.α)
                gamma(2α) * cospi(1 - α) * finda0(α) - a0singular_term(α)
            end
            a0singular_term[u0.degree] +=
                Arblib.add_error!(zero(error), error / lbound(Arb, u0.α)^u0.degree)
        end
        expansion[(2, 0, 0)] = -a0singular_term

        # Compute the coefficients for the analytic terms
        for m = 1:M-1
            term = (-1)^m * zeta(1 - 2α - 2m) / factorial(2m)
            term_remainder = let α = ArbSeries((u0.α, 1), degree = u0.degree)
                ((-1)^m*zeta(1 - 2α - 2m)/factorial(2m))[u0.degree]
            end
            expansion[(0, 0, 2m)] -= mul_with_remainder(u0.a[0], term, u0.α)
        end

        # Add remainder term
        remainder_term = clausenc_expansion_remainder(x, 1 - 2α, M)
        expansion[(0, 0, 2M)] -= mul_with_remainder(u0.a[0], remainder_term, u0.α)

        # Handle tail terms
        for j = 1:2
            s = 1 - 2α + j * u0.p0

            # Compute the coefficient for the singular term
            singular_term =
                compose_with_remainder(s -> gamma(1 - s) * sinpi(s / 2), s, u0.α)
            expansion[(2, j, 0)] = -mul_with_remainder(u0.a[j], singular_term, u0.α)

            # Compute the coefficients for the analytic terms
            for m = 1:M-1
                term =
                    (-1)^m * compose_with_remainder(s -> zeta(s - 2m), s, u0.α) /
                    factorial(2m)
                expansion[(0, 0, 2m)] -= mul_with_remainder(u0.a[j], term, u0.α)
            end

            # Add remainder term
            remainder_term = clausenc_expansion_remainder(x, s, M)
            expansion[(0, 0, 2M)] -= mul_with_remainder(u0.a[j], remainder_term, u0.α)
        end

        return expansion
    end
end

function D(u0::KdVZeroAnsatz, evaltype::Ball)
    f = H(u0, evaltype)
    return x -> begin
        p = u0(x, evaltype)
        p2 = mul_with_remainder(p, p, u0.α)
        return p2 / 2 + f(x)
    end
end

"""
    D(u0::KdVZeroAnsatz, ::AsymptoticExpansion; M::Integer)

Return an expansion of `D(u0)(x)` in `x` around zero where the
coefficients in the expansion are themselves expansions in `α` around
zero.

It returns a dictionary `expansion` where the keys are three tuples
`(i, j, m)` and correspond to a term of the form
```
expansion[(i, j, m)] * abs(x)^(i * α + j * p0 + m)
```

The value of `M` determines the number of terms in the expansion in
`x`.

The expansion is computed using `u0` and `H(u0)` but some modification
are made to it.

To begin with the leading terms in the expansion are identically equal
to zero due to the choice of `a[0]`, `a[1]`, `a[2]` and `p0`. These
are the terms with keys `(2, 0, 0), (2, 1, 0), (0, 0, 2), (1, 0, 2)`.

As is shown in [`F0`](@ref) we also know that the linear term of the
expansion in `α` is identically equal to zero. We can therefore set
all linear terms to zero since they must cancel out in the end anyway.

"""
function D(u0::KdVZeroAnsatz, evaltype::AsymptoticExpansion; M::Integer = 10)
    f = H(u0, evaltype; M)
    return x::Arb -> begin
        expansion1 = u0(x, evaltype; M)
        expansion2 = f(x)

        expansion = empty(expansion1)

        # u0^2/2 term
        let expansion1 = collect(expansion1)
            z = zero(first(expansion1)[2]) # Avoid allocating zero multiple times
            for (i, (key1, a1)) in enumerate(expansion1)
                expansion[2 .* key1] =
                    get(expansion, 2 .* key1, z) + mul_with_remainder(a1, a1, u0.α) / 2
                for j = i+1:length(expansion1)
                    (key2, a2) = expansion1[j]
                    key = key1 .+ key2
                    expansion[key] =
                        get(expansion, key, z) + mul_with_remainder(a1, a2, u0.α)
                end
            end
        end

        # H term
        merge!(+, expansion, expansion2)

        # Handle terms that are zero due to the choice of a[0], a[1],
        # a[2] and p0.
        for key in ((2, 0, 0), (2, 1, 0), (0, 0, 2), (1, 0, 2))
            y = expansion[key]
            @assert all(Arblib.contains_zero, Arblib.coeffs(y))
            expansion[key] = zero(y)
        end

        return expansion
    end
end

"""
    u0_div_xmα(u0::KdvZeroAnsatz, ::Asymptotic = Asymptotic(); ϵ)

The leading term of `u0` behaves like ` x^-α`. This method returns a
function `f` such that `f(x)` computes an enclosure of `u0(x) / x^-α`.
So an enclosure of `u0(x)` is given by
```
x^-α * f(x)
```
It does so in a way that allows evaluation at `x = 0` where `f` is
non-zero.

A value `ϵ` has to be given and this will allow evaluation for all `x
<= ϵ`.

Note that this method doesn't return an expansion in `α` but the
result is an enclosure that is valid for all values of `α`.

Recall that `u0(x)` is given by
```
u0(x) = a[0] * clausencmzeta(x, 1 - α) +
        a[1] * clausencmzeta(x, 1 - α + p0) +
        a[2] * clausencmzeta(x, 1 - α + 2p0)
```

The second and third term we can handle directly by computing
enclosures of `a[1]`, `a[2]`, `p0` and then using the expansion of the
Clausen functions around `x = 0`

For the first term direct evaluation fails computing the singular term
in `clausenc(x, 1 - α)` (it also fails computing the constant term but
that is exactly cancelled by the `zeta`). It works well for computing
the `x^2m` term for `m >= 1` and the remainder term.

The singular term is given by
```
gamma(α) * cospi(α / 2) * x^-α
```
which blows up as `α -> 0`. If we take into account the multiplication
by `a[0]` we get
```
a[0] * gamma(α) * cospi(α / 2) * x^-α =
    2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2)) * x^-α
```
The coefficient here also occurs in [`expansion_p0`](@ref) and its
expansion around `α = 0` is given by
```
2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2)) =
    1 - γ * α + (γ^2 / 2 - π^2 / 8) * α^2 +
    (-4γ^3 + 3γ * π^2 + 28 * polygamma(2, 1)) / 24 * α^3 + O(α^4)
```
which we can use to compute an enclosure

**IMPROVE:** We could attempt to compute a tighter enclosure if
needed. We could avoid some overestimations by using that the sign of
the individual terms are fixed. We could compute tighter enclosures of
`a[0]` by using the monotonicity, might be able to do similar things
for the other enclosures. We might not need to compute a tighter
enclosure in the end though.
"""
function u0_div_xmα(u0::KdVZeroAnsatz, ::Asymptotic = Asymptotic(); ϵ::Arb, M::Integer = 10)
    # We can just as well use an upper bound for ϵ
    ϵ = ubound(Arb, ϵ)

    # Compute enclosures of a[0], a[1], a[2] and p0
    a_enclosure = [a(u0.α) for a in u0.a]
    p0_enclosure = u0.p0(u0.α)

    _, _, P0, E0 = clausenc_expansion(ϵ, 1 - u0.α + p0_enclosure, M)
    C1, e1, P1, E1 = clausenc_expansion(ϵ, 1 - u0.α + p0_enclosure, M)
    C2, e2, P2, E2 = clausenc_expansion(ϵ, 1 - u0.α + 2p0_enclosure, M)

    # Enclosure of a[0] * C0 where C0 is the coefficient for the
    # singular term of the first Clausen function.
    a0C0 = let γ = Arb(Irrational{:γ}()), π = Arb(π)
        q = ArbSeries((1, -γ, γ^2 / 2 - π^2 / 8))
        q(u0.α)
    end

    return x::Union{Arb,ArbSeries} -> begin
        x isa Arb &&
            !(abs_ubound(x) <= abs_ubound(ϵ)) &&
            throw(ArgumentError("x has to be less than ϵ = $ϵ, got x = $x"))
        x isa ArbSeries &&
            !(abs_ubound(x[0]) <= abs_ubound(ϵ)) &&
            throw(ArgumentError("x has to be less than ϵ = $ϵ, got x = $x"))

        # Compute an enclosure of the first Clausen function
        clausen0_tail = zero(x)
        for m = 1:M-1
            clausen0_tail += P0[2m] * abspow(x, 2m + u0.α)
        end
        clausen0_tail += E0 * abspow(x, 2M + u0.α)

        # The exponent for the singular term is -u0.α + u0.α = 0
        res = a0C0 + a_enclosure[0] * clausen0_tail

        # Compute enclosures of the second and third Clausen functions

        # The exponent is given by e1 + u0.α = (-u0.α + p0_enclosure)
        # + u0.α = p0_enclosure
        clausen1 = C1 * abspow(x, p0_enclosure)
        for m = 1:M-1
            clausen1 += P1[2m] * abspow(x, 2m + u0.α)
        end
        clausen1 += E1 * abspow(x, 2M + u0.α)

        res += a_enclosure[1] * clausen1

        # The exponent is given by e2 + u0.α = (-u0.α + 2p0_enclosure)
        # + u0.α = 2p0_enclosure
        clausen2 = C2 * abspow(x, 2p0_enclosure)
        for m = 1:M-1
            clausen2 += P2[2m] * abspow(x, 2m + u0.α)
        end
        clausen2 += E2 * abspow(x, 2M + u0.α)

        res += a_enclosure[2] * clausen2

        return res
    end
end

"""
    F0(u0::KdVZeroAnsatz, ::Ball)

Return a function such that `F0(u0)(x)` computes an expansion in `α`
around `α = 0` of
```
(u0(x)^2 / 2 + H(u0)(x)) / (u0.w(x) * u0(x))
```

Both the constant and the linear terms in the expansion are supposed
to be zero. That the constant term is zero we get for free in the
computations, it is `1` for `u0(x)` and `-1 / 2` for `H(u0)(x)` so they
cancel exactly. For the linear term we need to prove that it is zero.

Let `p` and `q` denote the expansions in `α` from `u0(x)` and
`H(u0)(x)`. We have `p[0] = 1` and `q[0] = -1 / 2`. This means that the
linear term in `p^2 / 2 + q` is given by `p[1] + q[1]`, we want to
show that this term is exactly equal to zero.

We have
```
u0(x) = sum(a[j] * clausencmzeta(x, 1 - α + j * p0) for j = 0:2)
```
and
```
-sum(a[j] * clausencmzeta(x, 1 - 2α + j * p0) for j = 0:2)
```
It is enough to show that the linear term of the expansion in `α` of
```
a[j] * clausencmzeta(x, 1 - α + j * p0)
```
is the same as that for
```
a[j] * clausencmzeta(x, 1 - 2α + j * p0)
```
for `j = 0, 1, 2`. Since the coefficient `a[j]` is the same for both
terms we can ignore it when comparing the terms. However `a[j]` has a
zero constant coefficient in all cases so the linear term of the
product will be determined by the constant term of the Clausen
functions.

But the constant terms of the expansion in `α` are given by directly
plugging in `α = 0`, in which case the two terms are clearly
identical, which is what we wanted to show.
"""
function F0(u0::KdVZeroAnsatz, evaltype::Ball)
    f = H(u0, evaltype)

    return x -> begin
        p = u0(x, evaltype)
        q = f(x)

        @assert isone(Arblib.ref(p, 0))
        @assert Arblib.ref(q, 0) == Arb(-1 // 2)

        p2 = mul_with_remainder(p, p, u0.α)

        res = div_with_remainder(p2 / 2 + q, u0.w(x) * p, u0.α)

        @assert iszero(Arblib.ref(res, 0)) || !isfinite(Arblib.ref(res, 0))
        @assert Arblib.contains_zero(Arblib.ref(res, 1))
        res[1] = 0

        return res
    end
end

"""
    F0(u0::KdVZeroAnsatz, ::Asymptotic; ϵ)

Return a function such that `F0(u0)(x)` computes an expansion in `α`
around `α = 0` of
```
(u0(x)^2 / 2 + H(u0)(x)) / (u0.w(x) * u0(x))
```
It uses an evaluation strategy that works asymptotically in `x`.

It first computes an expansion of
```
D(u0)(x) = u0(x)^2 / 2 + H(u0)(x)
```
and `u0(x)` in `x`.

When computing
```
D(u0)(x) / (x * u0(x))
```
we then explicitly cancel `x^(1 - α)` in `D(u0)` and `x^-α` in `u0`.

Similarly to the non-asymptotic version both the constant and linear
terms in the expansion are supposed to be zero, which we enforce.
"""
function F0(u0::KdVZeroAnsatz, ::Asymptotic; ϵ::Arb = Arb(1), M::Integer = 10)
    D_expansion = D(u0, AsymptoticExpansion(); M)(ϵ)
    u0_expansion = u0(ϵ, AsymptoticExpansion(); M)

    return x::Arb -> begin
        x <= ϵ || throw(ArgumentError("need x <= ϵ, got x = $x with ϵ = $ϵ"))

        num = eval_expansion(u0, D_expansion, x, offset_i = -1, offset_m = -1)
        den = eval_expansion(u0, u0_expansion, x, offset_i = -1)

        @assert iszero(Arblib.ref(num, 0))
        @assert Arblib.contains_zero(Arblib.ref(num, 1))
        num[1] = 0

        return div_with_remainder(num, den, u0.α)
    end
end
