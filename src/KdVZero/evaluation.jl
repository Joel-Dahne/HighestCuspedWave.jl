"""
    eval_expansion(u0::KdVZeroAnsatz, expansion::AbstractDict{NTuple{3,Int},TaylorModel}, x::Arb)

Evaluate the given expansion. The term `((i, j, m), y)` corresponds to
```
y * abs(x)^(-i * α + j * p0 + m)
```
"""
function eval_expansion(
    u0::KdVZeroAnsatz,
    expansion::AbstractDict{NTuple{3,Int},TaylorModel},
    x::Arb;
    offset_i::Integer = 0,
    offset_m::Integer = 0,
)
    Mα = TaylorModel(identity, u0.α, u0.α0; u0.degree)

    res = zero(Mα)

    for ((i, j, m), y) in expansion
        if !iszero(y)
            if iszero(i + offset_i) && iszero(j) && iszero(m + offset_m)
                # Exponent is identically equal to zero
                term = one(res)
            else
                exponent = -(i + offset_i) * Mα + j * u0.p0 + (m + offset_m)
                term = abspow(x, exponent)
            end
            res += y * term
        end
    end

    return res
end

"""
    eval_expansion(u0::KdVZeroAnsatz, expansion::AbstractDict{NTuple{3,Int},Arb}, x::Union{Arb,ArbSeries})

Evaluate the given expansion. The term `((i, j, m), y)` corresponds to
```
y * abs(x)^(-i * α + j * p0 + m)
```

Note that this takes a dictionary where the values are of type `Arb`
and not `TaylorModel`. This methods doesn't compute a Taylor model in
`α` like most methods for `KdVZeroansatz`. The benefit is that it
allows evaluation with `x::ArbSeries`.
"""
function eval_expansion(
    u0::KdVZeroAnsatz,
    expansion::AbstractDict{NTuple{3,Int},Arb},
    x::Union{Arb,ArbSeries};
    offset_i::Integer = 0,
    offset_m::Integer = 0,
)
    # Enclosure of p0 on the interval u0.α
    p0 = u0.p0(u0.α)

    res = zero(x)

    for ((i, j, m), y) in expansion
        if !iszero(y)
            exponent = -(i + offset_i) * u0.α + j * p0 + (m + offset_m)
            term = abspow(x, exponent)

            res += y * term
        end
    end

    return res
end

"""
    (u0::KdVZeroAnsatz)(x::Arb, ::Ball)

Return a Taylor model of `u0(x)`

The value is given by
```
sum(a[j] * clausencmzeta(x, 1 - α + j * p0) for j = 0:2)
```

For `j = 1, 2` we can evaluate it directly with
[`clausencmzeta`](@ref). For `j = 0` we can evaluate it directly if
`u0.α0 < 0`, for `u0.α0 = 0` there is a removable singularity to
handle. We describe the procedure for handling the removable
singularity below.

# Handling the removable singularity for `α0 = 0`
We are interested in computing an expansion around `α = 0` of
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
write it as `α * gamma(α) = 1 / (rgamma(α) / α)` where we can handle
the removable singularity

For wide values of `x` direct evaluation of `zeta(α, x / 2π) + zeta(α,
1 - x / 2π)` gives a poor enclosure. To get better enclosures we use
that all derivatives of that `zeta(α, x / 2π) + zeta(α, 1 - x / 2π)`
are monotone on the interval `(0, π)`. We don't do this for the
remainder term because it doesn't give much of an improvement. To see
that all derivatives are monotone we can rewrite it as
```
zeta(α, x / 2π) + zeta(α, 1 - x / 2π) = clausenc(x, α) / (gamma(v) * inv(2π)^v * cospi(v / 2))
```
The left hand side is monotone on `(0, π)` since `clausenc(x, α)` is.

For computing `a[0] * zeta(1 - α)` we rewrite it as
```
(a[0] / α) * (α * zeta(1 - α))
```
and then use Riemann's functional equation to get
```
zeta(1 - α) = zeta(α) / (2^α * π^(α - 1) * gamma(1 - α) * sinpi(α / 2))
```
We handle the removable singularity by writing it as
```
α * zeta(1 - α) = zeta(α) / (2^α * π^(α - 1) * gamma(1 - α) * (sinpi(α / 2) / α))
```
and explicitly dealing with `sinpi(α / 2) / α`. From the expansion
```
zeta(s) = 1 / (s - 1) + sum((-1)^n * stieltjes(n) / factorial(n) * (s - 1)^n for n = 0:Inf)
```
we can also notice that the constant term in the expansion is exactly
`-1`.
"""
function (u0::KdVZeroAnsatz)(x::Arb, ::Ball)
    if iszero(u0.α0)
        # Handle the term for j = 0 manually
        res = let
            # Compute clausenc(x, 1 - α). We compute this to one
            # degree higher since we divide it by α, reducing the
            # degree by one.

            # Taylor model of
            # inv(2π)^α * cospi(α / 2) * (zeta(α, x / 2π) + zeta(α, 1 - x / 2π)) / α
            clausen_term_part1 = let π = Arb(π)
                clausen_term_part1_expansion =
                    let α = ArbSeries((u0.α0, 1), degree = u0.degree + 1)
                        if iswide(x) && 0 < x < π
                            z = ArbSeries(
                                union.(
                                    Arblib.coeffs(
                                        zeta(α, lbound(Arb, x) / 2π) +
                                        zeta(α, 1 - lbound(Arb, x) / 2π),
                                    ),
                                    Arblib.coeffs(
                                        zeta(α, ubound(Arb, x) / 2π) +
                                        zeta(α, 1 - ubound(Arb, x) / 2π),
                                    ),
                                ),
                            )
                        else
                            z = zeta(α, x / 2π) + zeta(α, 1 - x / 2π)
                        end
                        inv(2π)^α * cospi(α / 2) * z
                    end

                # The constant term is exactly zero
                @assert Arblib.contains_zero(clausen_term_part1_expansion[0])
                clausen_term_part1_expansion[0] = 0

                clausen_term_part1_remainder =
                    let α = ArbSeries((u0.α, 1), degree = u0.degree + 1)
                        (inv(2π)^α*cospi(α / 2)*(zeta(α, x / 2π)+zeta(α, 1 - x / 2π)))[u0.degree+1]
                    end

                # Extend expansion to include remainder term
                clausen_term_part1_expansion =
                    ArbSeries(clausen_term_part1_expansion, degree = u0.degree + 2)
                clausen_term_part1_expansion[u0.degree+2] = clausen_term_part1_remainder

                # Construct Taylor model and divide by α
                TaylorModel(clausen_term_part1_expansion, u0.α, u0.α0) << 1
            end

            # Taylor model of α * gamma(α) = 1 / (rgamma(α) / α)
            clausen_term_part2 =
                1 / (
                    TaylorModel(
                        rgamma,
                        u0.α,
                        u0.α0,
                        degree = u0.degree + 1,
                        enclosure_degree = 1,
                    ) << 1
                )

            # Taylor model of clausenc(x, 1 - α)
            clausen_term = clausen_term_part1 * clausen_term_part2

            # Taylor model of a[0] * clausenc(x, 1 - α)
            # Divide a[0] by α, perform the multiplication and then
            # multiply by α. This makes sure the degree after the
            # multiplication is correct.
            a0_clausen_term = ((u0.a[0] << 1) * clausen_term) >> 1

            # Taylor model of α * zeta(1 - α)
            # We first compute it to a slightly higher degree and then
            # truncate
            zetamulα = let
                # zeta(α)
                # The Arb implementation of zeta for ArbSeries doesn't
                # work well for intervals around zero which are not
                # symmetric. We therefore artificially make the
                # interval symmetric by widening it. It doesn't seem
                # to affect the computed bounds much.
                f1 = TaylorModel(
                    TaylorModel(zeta, union(u0.α, -u0.α), u0.α0, degree = u0.degree + 3).p,
                    u0.α,
                    u0.α0,
                )

                # sinpi(α / 2) / α
                f2 =
                    TaylorModel(α -> sinpi(α / 2), u0.α, u0.α0, degree = u0.degree + 4) << 1

                # 2^α * π^(α - 1) * gamma(1 - α)
                f3 = TaylorModel(
                    α -> 2^α * Arb(π)^(α - 1) * gamma(1 - α),
                    u0.α,
                    u0.α0,
                    degree = u0.degree + 3,
                )

                truncate(f1 / (f2 * f3); u0.degree)
            end

            # The constant term is exactly -1
            @assert contains(zetamulα.p[0], -1)
            zetamulα.p[0] = -1

            # Taylor model of a[0] * zeta(1 - α)
            a0zeta_term = (u0.a[0] << 1) * zetamulα

            truncate(a0_clausen_term; u0.degree) - a0zeta_term
        end
    else
        res = zero(u0.p0)
    end

    # Taylor model of α
    Mα = TaylorModel(identity, u0.α, u0.α0; u0.degree)

    # If u0.α0 is non-zero we handle the case j = 0 here
    j_start = ifelse(iszero(u0.α0), 1, 0)
    for j = j_start:2
        term = clausencmzeta(x, 1 - Mα + j * u0.p0)

        if j == 0
            # u0.a[0] has a higher degree to we truncate it first
            res += truncate(u0.a[j]; u0.degree) * term
        else
            res += u0.a[j] * term
        end
    end

    return res
end

"""
    (u0::KdVZeroAnsatz)(x::Arb, ::AsymptoticExpansion; M::Integer)

Return an expansion of `u0(x)` in `x` around zero where the
coefficients in the expansion are Taylor models.

It returns a dictionary `expansion` where the keys are three tuples
`(i, j, m)` and correspond to a term of the form
```
expansion[(i, j, m)] * abs(x)^(i * α + j * p0 + m)
```

The value of `M` determines the number of terms in the expansion in
`x`.

Most of the terms can be computed by evaluating them directly. The
exception is the singular term for the first Clausen function when
`u0.α0 = 0`. It is given by `gamma(α) * cospi(α / 2)` which has a
singularity at `α = 0`. However multiplication by `a[0]` removes the
singularity and we can compute it by rewriting it as
```
a[0] * gamma(α) * cospi(α / 2) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```
where we can handle `gamma(2α) / gamma(α) = rgamma(α) / rgamma(2α)`
similarly to how it is done in [`expansion_as`](@ref).
"""
function (u0::KdVZeroAnsatz)(x::Arb, ::AsymptoticExpansion; M::Integer = 10)
    Mα = TaylorModel(identity, u0.α, u0.α0; u0.degree)

    expansion = OrderedDict{NTuple{3,Int},TaylorModel}()

    # Initiate even powers of x
    for m = 1:M
        expansion[(0, 0, 2m)] = zero(Mα)
    end

    # Handle main term
    # Compute the coefficient for the singular term
    if iszero(u0.α0)
        # We compute a0 to a higher degree and then truncate, to
        # get a tighter enclosure
        a0singular_term = let
            # g = rgamma(α) / rgamma(2α)
            g = div_removable(
                TaylorModel(
                    rgamma,
                    u0.α,
                    u0.α0,
                    degree = u0.degree + 4,
                    enclosure_degree = 1,
                ),
                TaylorModel(
                    α -> rgamma(2α),
                    u0.α,
                    u0.α0,
                    degree = u0.degree + 4,
                    enclosure_degree = 1,
                ),
            )

            a0singular_term =
                g * TaylorModel(
                    α -> 2cospi(α) / cospi(α / 2),
                    u0.α,
                    u0.α0,
                    degree = u0.degree + 3,
                )

            truncate(a0singular_term; u0.degree)
        end
    else
        a0singular_term =
            truncate(u0.a[0]; u0.degree) * compose(α -> gamma(α) * cospi(α / 2), Mα)
    end

    expansion[(1, 0, 0)] = a0singular_term

    # Compute the coefficients for the analytic terms
    for m = 1:M-1
        term = (-1)^m * compose(α -> zeta(1 - α - 2m), Mα) / factorial(2m)

        expansion[(0, 0, 2m)] += truncate(u0.a[0]; u0.degree) * term
    end

    # Add remainder term (in x)
    remainder_term = compose(α -> clausenc_expansion_remainder(x, 1 - α, M), Mα)
    expansion[(0, 0, 2M)] += truncate(u0.a[0]; u0.degree) * remainder_term

    # Handle tail terms
    for j = 1:2
        s = 1 - Mα + j * u0.p0

        # Compute the coefficient for the singular term
        singular_term = compose(s -> gamma(1 - s) * sinpi(s / 2), s)
        expansion[(1, j, 0)] = u0.a[j] * singular_term

        # Compute the coefficients for the analytic terms
        for m = 1:M-1
            term = (-1)^m * compose(s -> zeta(s - 2m), s) / factorial(2m)
            expansion[(0, 0, 2m)] += u0.a[j] * term
        end

        # Add remainder term
        Mremainder_term = compose(s -> clausenc_expansion_remainder(x, s, M), s)
        expansion[(0, 0, 2M)] += u0.a[j] * remainder_term
    end

    return expansion
end

"""
    H(u0::KdVZeroAnsatz, ::Ball)

Return a function such that `H(u0)(x)` computes a Taylor model of
`H(u0)`.

The value is given by
```
-sum(a[j] * clausencmzeta(x, 1 - 2α + j * p0) for j = 0:2)
```

For `j = 1, 2` we can evaluate it directly with
[`clausencmzeta`](@ref). For `j = 0` we can evaluate it directly if
`u0.α0 < 0`, for `u0.α0 = 0` there is a removable singularity to
handle. We describe the procedure for handling the removable
singularity below.

# Handling the removable singularity for `α0 = 0`
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
write it as `α * gamma(2α) = 1 / (rgamma(2α) / α)` where we can handle
the removable singularity

For wide values of `x` direct evaluation of `zeta(2α, x / 2π) +
zeta(2α, 1 - x / 2π)` gives a poor enclosure. To get better enclosures
we use that all derivatives of that `zeta(2α, x / 2π) + zeta(2α, 1 - x
/ 2π)` are monotone on the interval `(0, π)`. We don't do this for the
remainder term because it doesn't give much of an improvement. For
proof of monotonicity see `u0(x)`.

For computing `a[0] * zeta(1 - 2α)` we rewrite it as
```
(a[0] / α) * (α * zeta(1 - 2α))
```
and then use Riemann's function equation to get
```
zeta(1 - 2α) = zeta(2α) / (2^2α * π^(2α - 1) * gamma(1 - 2α) * sinpi(α))
```
to get
```
α * zeta(1 - 2α) = zeta(2α) / (2^2α * π^(2α - 1) * gamma(1 - 2α) * (sinpi(α) / α))
```
and explicitly dealing with `sinpi(α) / α`. From the expansion
```
zeta(s) = 1 / (s - 1) + sum((-1)^n * stieltjes(n) / factorial(n) * (s - 1)^n for n = 0:Inf)
```
we can also notice that the constant term in the expansion is exactly
`1 / 2`.
"""
function H(u0::KdVZeroAnsatz, ::Ball)
    return x::Arb -> begin
        # The main term we handle manually
        if iszero(u0.α0)
            res = let
                # Taylor model of
                # inv(2π)^2α * cospi(α) * (zeta(2α, x / 2π) + zeta(2α, 1 - x / 2π)) / α
                clausen_term_part1 = let π = Arb(π)
                    clausen_term_part1_expansion =
                        let α = ArbSeries((u0.α0, 1), degree = u0.degree + 2)
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

                    # The constant term is exactly zero
                    @assert Arblib.contains_zero(clausen_term_part1_expansion[0])
                    clausen_term_part1_expansion[0] = 0

                    clausen_term_part1_remainder =
                        let α = ArbSeries((u0.α, 1), degree = u0.degree + 2)
                            (inv(2π)^2α*cospi(α)*(zeta(2α, x / 2π)+zeta(2α, 1 - x / 2π)))[u0.degree+2]
                        end

                    # Extend expansion to include remainder term
                    clausen_term_part1_expansion = ArbSeries(
                        clausen_term_part1_expansion;
                        degree = u0.degree + 2,
                    )
                    clausen_term_part1_expansion[u0.degree+2] =
                        clausen_term_part1_remainder

                    # Construct Taylor model and divide by α
                    TaylorModel(clausen_term_part1_expansion, u0.α, u0.α0) << 1
                end

                # # Taylor model of α * gamma(2α) =  = 1 / (rgamma(2α) / α)
                clausen_term_part2 =
                    1 / (
                        TaylorModel(
                            α -> rgamma(2α),
                            u0.α,
                            u0.α0,
                            degree = u0.degree + 1,
                            enclosure_degree = 1,
                        ) << 1
                    )

                # Taylor model of clausenc(x, 1 - α)
                clausen_term = clausen_term_part1 * clausen_term_part2

                # Taylor model of a[0] * clausenc(x, 1 - α)
                # Divide a[0] by α, perform the multiplication and then
                # multiply by α. This makes sure the degree after the
                # multiplication is correct.
                a0_clausen_term = ((u0.a[0] << 1) * clausen_term) >> 1

                # Taylor model of α * zeta(1 - 2α)
                # We first compute it to a slightly higher degree and then
                # truncate
                zetamulα = let
                    # zeta(2α)
                    # The Arb implementation of zeta for ArbSeries doesn't
                    # work well for intervals around zero which are not
                    # symmetric, we therefore artificially make the
                    # interval symmetric by widening it. It doesn't seem
                    # to affect the computed bounds much.
                    f1 = TaylorModel(
                        TaylorModel(
                            α -> zeta(2α),
                            union(u0.α, -u0.α),
                            u0.α0,
                            degree = u0.degree + 3,
                        ).p,
                        u0.α,
                        u0.α0,
                    )

                    # sinpi(α) / α
                    f2 = TaylorModel(sinpi, u0.α, u0.α0, degree = u0.degree + 4) << 1

                    # 2^2α * π^(2α - 1) * gamma(1 - 2α)
                    f3 = TaylorModel(
                        α -> 2^2α * Arb(π)^(2α - 1) * gamma(1 - 2α),
                        u0.α,
                        u0.α0,
                        degree = u0.degree + 3,
                    )

                    truncate(f1 / (f2 * f3); u0.degree)
                end

                # The constant term is exactly -1 / 2
                @assert contains(zetamulα.p[0], Arb(-1 // 2))
                zetamulα.p[0] = -1 // 2

                # Taylor model of a[0] * zeta(1 - α)
                a0zeta_term = (u0.a[0] << 1) * zetamulα

                -(truncate(a0_clausen_term; u0.degree) - a0zeta_term)
            end
        else
            res = zero(u0.p0)
        end

        Mα = TaylorModel(identity, u0.α, u0.α0; u0.degree)

        # If u0.α0 is non-zero we handle the case j = 0 here
        j_start = ifelse(iszero(u0.α0), 1, 0)
        for j = j_start:2
            term = clausencmzeta(x, 1 - 2Mα + j * u0.p0)

            if j == 0
                res -= truncate(u0.a[j]; u0.degree) * term
            else
                res -= u0.a[j] * term
            end
        end

        return res
    end
end

"""
    H(u0::KdVZeroAnsatz, ::AsymptoticExpansion)

Return a function such that `H(u0)(x)` computes an expansion in `x`
around zero where the coefficients in the expansion are Taylor models.

It returns a dictionary `expansion` where the keys are three tuples
`(i, j, m)` and correspond to a term of the form
```
expansion[(i, j, m)] * abs(x)^(i * α + j * p0 + m)
```

The value of `M` determines the number of terms in the expansion in
`x`.

Most of the terms can be computed by evaluating them directly. The
exception is the singular term for the first Clausen function when
`u0.α0 = 0`. It is given by `gamma(2α) * cospi(α)` which has a
singularity at `α = 0`. However multiplication by `a[0]` removes the
singularity and we have
```
a[0] * gamma(2α) * cospi(α) = 2gamma(2α)^2 * cospi(α)^2 / (gamma(α)^2 * cospi(α / 2)^2)
    = 2(gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2)))^2
```
where we can handle `gamma(2α) / gamma(α) = rgamma(α) / rgamma(2α)`
similarly to how it is done in [`expansion_as`](@ref).
"""
function H(u0::KdVZeroAnsatz, ::AsymptoticExpansion; M::Integer = 10)
    Mα = TaylorModel(identity, u0.α, u0.α0; u0.degree)

    return x::Arb -> begin
        expansion = OrderedDict{NTuple{3,Int},TaylorModel}()

        # Initiate even powers of x
        for m = 1:M
            expansion[(0, 0, 2m)] = zero(Mα)
        end

        # Handle main term
        # Compute the coefficient for the singular term
        if iszero(u0.α0)
            a0singular_term = let
                # We compute a0 to a higher degree and then truncate, to
                # get a tighter enclosure

                # g = rgamma(α) / rgamma(2α)
                g = div_removable(
                    TaylorModel(rgamma, u0.α, u0.α0, degree = u0.degree + 4),
                    TaylorModel(α -> rgamma(2α), u0.α, u0.α0, degree = u0.degree + 4),
                )

                # rgamma(α) * cospi(α) / (rgamma(2α) * cospi(α / 2))
                a0singular_term =
                    g * TaylorModel(
                        α -> cospi(α) / cospi(α / 2),
                        u0.α,
                        u0.α0,
                        degree = u0.degree + 3,
                    )

                # 2(rgamma(α) * cospi(α) / (rgamma(2α) * cospi(α / 2)))^2
                a0singular_term = 2a0singular_term * a0singular_term

                truncate(a0singular_term; u0.degree)
            end
        else
            a0singular_term =
                truncate(u0.a[0]; u0.degree) * compose(α -> gamma(2α) * cospi(α), Mα)
        end
        expansion[(2, 0, 0)] = -a0singular_term

        # Compute the coefficients for the analytic terms
        for m = 1:M-1
            term = (-1)^m * compose(α -> zeta(1 - 2α - 2m), Mα) / factorial(2m)

            expansion[(0, 0, 2m)] -= truncate(u0.a[0]; u0.degree) * term
        end

        # Add remainder term
        remainder_term = compose(α -> clausenc_expansion_remainder(x, 1 - 2α, M), Mα)
        expansion[(0, 0, 2M)] -= truncate(u0.a[0]; u0.degree) * remainder_term

        # Handle tail terms
        for j = 1:2
            s = 1 - 2Mα + j * u0.p0

            # Compute the coefficient for the singular term
            singular_term = compose(s -> gamma(1 - s) * sinpi(s / 2), s)
            expansion[(2, j, 0)] = -u0.a[j] * singular_term

            # Compute the coefficients for the analytic terms
            for m = 1:M-1
                term = (-1)^m * compose(s -> zeta(s - 2m), s) / factorial(2m)
                expansion[(0, 0, 2m)] -= u0.a[j] * term
            end

            # Add remainder term
            remainder_term = compose(s -> clausenc_expansion_remainder(x, s, M), s)
            expansion[(0, 0, 2M)] -= u0.a[j] * remainder_term
        end

        return expansion
    end
end

function D(u0::KdVZeroAnsatz, evaltype::Ball)
    f = H(u0, evaltype)
    return x -> begin
        p = u0(x, evaltype)
        return p * p / 2 + f(x)
    end
end

"""
    D(u0::KdVZeroAnsatz, ::AsymptoticExpansion; M::Integer)

Return a function such that `D(u0)(x)` computes an expansion in `x`
around zero of
```
u0(x)^2 / 2 + H(u0)(x)
```
where the coefficients in the expansion are Taylor models.

It returns a dictionary `expansion` where the keys are three tuples
`(i, j, m)` and correspond to a term of the form
```
expansion[(i, j, m)] * abs(x)^(i * α + j * p0 + m)
```

The value of `M` determines the number of terms in the expansion in
`x`.

The expansion is computed using `u0` and `H(u0)` but some modification
are made to it. The leading terms in the expansion are identically
equal to zero due to the choice of `a[0]`, `a[1]`, `a[2]` and `p0`.
These are the terms with keys `(2, 0, 0), (2, 1, 0), (0, 0, 2), (1, 0,
2)`.
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
                expansion[2 .* key1] = get(expansion, 2 .* key1, z) + a1 * a1 / 2
                for j = i+1:length(expansion1)
                    (key2, a2) = expansion1[j]
                    key = key1 .+ key2
                    expansion[key] = get(expansion, key, z) + a1 * a2
                end
            end
        end

        # H term
        merge!(+, expansion, expansion2)

        # Handle terms that are zero due to the choice of a[0], a[1],
        # a[2] and p0.
        for key in ((2, 0, 0), (2, 1, 0), (0, 0, 2), (1, 0, 2))
            y = expansion[key]
            @assert all(Arblib.contains_zero, Arblib.coeffs(y.p))
            expansion[key] = zero(y)
        end

        return expansion
    end
end

"""
    F0(u0::KdVZeroAnsatz, ::Ball)

Return a function such that `F0(u0)(x)` computes a Taylor model in `α`
around `u0.α0` of
```
(u0(x)^2 / 2 + H(u0)(x)) / (u0.w(x) * u0(x))
```

The computation is straight forward. However, for `u0.α0 = 0` both the
constant and linear term are supposed to be exactly zero so we have to
deal with that.

# Constant and linear term for `u0.α0 = 0`
That the constant term is zero we get for free in the computations, it
is `1` for `u0(x)` and `-1 / 2` for `H(u0)(x)` so they cancel exactly.
For the linear term we need to prove that it is zero.

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
H(u0)(x) = -sum(a[j] * clausencmzeta(x, 1 - 2α + j * p0) for j = 0:2)
```
It is enough to show that the linear term of the expansion in `α` of
```
a[j] * clausencmzeta(x, 1 - α + j * p0)
```
is the same as that for
```
a[j] * clausencmzeta(x, 1 - 2α + j * p0)
```
for `j = 0, 1, 2`. For `j = 1, 2` the Clausen term is finite and since
`a[j]` has a zero constant coefficient the linear term of the product
will be determined by the constant term of the Clausen functions,
which is trivially the same. For `j = 0` there is a removable
singularity to treat. The linear terms are the same in this case as
well, for a proof we refer to the corresponding lemma in the paper.
"""
function F0(u0::KdVZeroAnsatz, evaltype::Ball)
    f = H(u0, evaltype)

    return x -> begin
        p = u0(x, evaltype)
        q = f(x)

        res = (p * p / 2 + q) / (u0.w(x) * p)

        if iszero(u0.α0)
            # This ensures that the constant term is zero
            @assert isone(Arblib.ref(p.p, 0))
            @assert Arblib.ref(q.p, 0) == Arb(-1 // 2)
            @assert iszero(Arblib.ref(res.p, 0)) || !isfinite(Arblib.ref(res.p, 0))

            # The linear term should be zero
            @assert Arblib.contains_zero(Arblib.ref(res.p, 1))
            res.p[1] = 0
        end

        return res
    end
end

"""
    F0(u0::KdVZeroAnsatz, ::Asymptotic; ϵ)

Return a function such that `F0(u0)(x)` computes a Taylor model of
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
terms in the expansion are supposed to be zero when `u0.α0 = 0`, which
we enforce.
"""
function F0(u0::KdVZeroAnsatz, ::Asymptotic; ϵ::Arb = Arb(1), M::Integer = 10)
    D_expansion = D(u0, AsymptoticExpansion(); M)(ϵ)
    u0_expansion = u0(ϵ, AsymptoticExpansion(); M)

    return x::Arb -> begin
        x <= ϵ || throw(ArgumentError("need x <= ϵ, got x = $x with ϵ = $ϵ"))

        num = eval_expansion(u0, D_expansion, x, offset_i = -1, offset_m = -1)
        den = eval_expansion(u0, u0_expansion, x, offset_i = -1)

        if iszero(u0.α0)
            @assert iszero(Arblib.ref(num.p, 0))
            @assert Arblib.contains_zero(Arblib.ref(num.p, 1))
            num.p[1] = 0
        end

        return num / den
    end
end

"""
    F02(u0::KdVZeroAnsatz, ::Asymptotic; ϵ)

Return a function such that `F02(u0)(x)` computes an enclosure of
```
(u0(x)^2 / 2 + H(u0)(x)) / (u0.w(x) * u0(x))
```
It uses an evaluation strategy that works asymptotically in `x`.

Note that this methods doesn't compute a Taylor model in `α` like most
methods for `KdVZeroansatz`. It only computes an enclosure. For this
reason it doesn't make sense to use for `u0.α0 = 0`
"""
function F02(u0::KdVZeroAnsatz, ::Asymptotic; ϵ::Arb = Arb(1), M::Integer = 10)
    iszero(u0.α0) && @warn "F02 doesn't make sense for u0.α0 = 0"

    # Compute the expansions and evaluate their terms in α
    D_expansion = let expansion = D(u0, AsymptoticExpansion(); M)(ϵ)
        res = empty(expansion, Arb)
        for (key, value) in expansion
            res[key] = value(u0.α)
        end
        res
    end
    u0_expansion = let expansion = u0(ϵ, AsymptoticExpansion(); M)
        res = empty(expansion, Arb)
        for (key, value) in expansion
            res[key] = value(u0.α)
        end
        res
    end

    return x::Union{Arb,ArbSeries} -> begin
        if x isa Arb
            x <= ϵ || throw(ArgumentError("need x <= ϵ, got x = $x with ϵ = $ϵ"))
        elseif x isa ArbSeries
            x[0] <= ϵ || throw(ArgumentError("need x[0] <= ϵ, got x = $x with ϵ = $ϵ"))
        end

        num = eval_expansion(u0, D_expansion, x, offset_i = -1, offset_m = -1)
        den = eval_expansion(u0, u0_expansion, x, offset_i = -1)

        return num / den
    end
end
