"""
    (u0::KdVZeroAnsatz)(x::Arb, ::Ball)

Return an expansion of `u0(x)` in `α` around `α = 0`.

The value is given by
```
sum(as[j] * clausencmzeta(x, 1 - α + j * p0) for j = 0:2)
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
remove the pole and we have the expansion
```
α * gamma(α) = 1 - γ * α + (γ^2 / 2 + π^2 / 12) * α^2 + O(α^3)
```

For computing `a[0] * zeta(1 - α)` we rewrite it as
```
(a[0] / α) * (α * zeta(1 - α))
```
and use the expansion
```
α * zeta(1 - α) = -1 + stieltjes(0) * α + stieltjes(1) * α^2 + stieltjes(2) / 2 * α^3 + O(α^4)
```
"""
function (u0::KdVZeroAnsatz)(x::Arb, ::Ball)
    as = expansion_as(u0)
    p0 = expansion_p0(u0)
    α = ArbSeries((0, 1), degree = Arblib.degree(p0))

    # The main term we handle manually
    res = let
        clausen_term = let π = Arb(π)
            inv(2π)^α * cospi(α / 2) * (zeta(α, x / 2π) + zeta(α, 1 - x / 2π))
        end

        # The constant term is exactly zero
        @assert Arblib.contains_zero(clausen_term[0])
        clausen_term[0] = 0

        # Divide clausen term by α
        clausen_term = clausen_term << 1

        # Expansion of α * gamma(α)
        gammamulα = let γ = Arb(Irrational{:γ}()), π = Arb(π)
            ArbSeries((1, -γ, γ^2 / 2 + π^2 / 12))
        end

        clausen_term *= gammamulα

        # Divide as[0] by α, perform the multiplication and then
        # multiply by α. This makes sure the degree after the
        # multiplication is correct.
        a0clausenterm = ((as[0] << 1) * clausen_term) >> 1

        zetamulα =
            ArbSeries((-1, stieltjes(Arb, 0), stieltjes(Arb, 1), stieltjes(Arb, 2) / 2))

        a0zeta_term = (as[0] << 1) * zetamulα

        a0clausenterm - a0zeta_term
    end

    res += sum(as[j] * clausencmzeta(x, 1 - α + j * p0) for j = 1:2)

    return res
end

"""
    (u0::KdVZeroAnsatz)(x::Arb, ::Asymptotic; M::Integer = 3)

Return an expansion of `u0(x)` in `x` around zero where the
coefficients in the expansion are themselves expansions in `α` around
zero.

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
- **TODO:** Compute remainder term

- **TODO:** Figure out how to handle remainder terms.
"""
function (u0::KdVZeroAnsatz)(x::Arb, ::AsymptoticExpansion; M::Integer = 3)
    as = expansion_as(u0)
    p0 = expansion_p0(u0)
    α = ArbSeries((0, 1), degree = Arblib.degree(p0))

    expansion = OrderedDict{NTuple{3,Int},ArbSeries}()

    # Initiate even powers of x
    for m = 1:M
        expansion[(0, 0, 2m)] = ArbSeries(degree = 3)
    end

    # Handle main term
    s = 1 - α

    # Compute the coefficient for the singular term
    a0singular_term = let γ = Arb(Irrational{:γ}()), π = Arb(π)
        ArbSeries((1, -γ, γ^2 / 2 - π^2 / 8))
    end
    expansion[(1, 0, 0)] = a0singular_term

    # Compute the coefficients for the analytic terms
    for m = 1:M-1
        term = (-1)^m * zeta(1 - α - 2m) / factorial(2m)
        expansion[(0, 0, 2m)] += as[0] * term
    end

    # Add error term
    error_term = ArbSeries() # TODO: How to handle this?
    expansion[(0, 0, 2M)] += as[0] * error_term

    # Handle tail terms
    for j = 1:2
        s = 1 - α + j * p0

        # Compute the coefficient for the singular term
        singular_term = gamma(1 - s) * sinpi(s / 2)
        expansion[(1, j, 0)] = as[j] * singular_term

        # Compute the coefficients for the analytic terms
        for m = 1:M-1
            term = (-1)^m * zeta(s - 2m) / factorial(2m)
            expansion[(0, 0, 2m)] += as[j] * term
        end

        # Add error term
        error_term = ArbSeries() # TODO: How to handle this?
        expansion[(0, 0, 2M)] += as[j] * error_term
    end

    return expansion
end

"""
    H(u0::KdVZeroAnsatz, ::Ball)

Return a function such that `H(u0)(x)` computes an expansion in `α`
around `α = 0`.

The value is given by
```
-sum(as[j] * clausencmzeta(x, 1 - 2α + j * p0) for j = 0:2)
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
remove the pole and we have the expansion
```
α * gamma(2α) = 1 / 2 - γ * α + (γ^2 + π^2 / 6) * α^2 + O(α^3)
```

For computing `a[0] * zeta(1 - 2α)` we rewrite it as
```
(a[0] / α) * (α * zeta(1 - 2α))
```
and use the expansion
```
α * zeta(1 - 2α) = -1 / 2 + stieltjes(0) * α + 2stieltjes(1) * α^2 + 2stieltjes(2) * α^3 + O(α^4)
```
"""
function H(u0::KdVZeroAnsatz, ::Ball)
    as = expansion_as(u0)
    p0 = expansion_p0(u0)
    α = ArbSeries((0, 1), degree = Arblib.degree(p0))

    return x::Arb -> begin
        # The main term we handle manually
        res = let
            clausen_term = let π = Arb(π)
                inv(2π)^2α * cospi(α) * (zeta(2α, x / 2π) + zeta(2α, 1 - x / 2π))
            end

            # The constant term is exactly zero
            @assert Arblib.contains_zero(clausen_term[0])
            clausen_term[0] = 0

            # Divide clausen term by α
            clausen_term = clausen_term << 1

            # Expansion of α * gamma(2α)
            gammamulα = let γ = Arb(Irrational{:γ}()), π = Arb(π)
                ArbSeries((1 // 2, -γ, γ^2 + π^2 / 6))
            end

            clausen_term *= gammamulα

            # Divide as[0] by α, perform the multiplication and then
            # multiply by α. This makes sure the degree after the
            # multiplication is correct.
            a0clausenterm = ((as[0] << 1) * clausen_term) >> 1

            # Expansion of α * zeta(1 - 2α)
            zetamulα = ArbSeries((
                -1 // 2,
                stieltjes(Arb, 0),
                2stieltjes(Arb, 1),
                2stieltjes(Arb, 2),
            ))

            a0zeta_term = (as[0] << 1) * zetamulα

            -(a0clausenterm - a0zeta_term)
        end

        res -= sum(as[j] * clausencmzeta(x, 1 - 2α + j * p0) for j = 1:2)

        return res
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
function u0_div_xmα(
    u0::KdVZeroAnsatz{Arb},
    ::Asymptotic = Asymptotic();
    ϵ::Arb,
    M::Integer = 3,
)
    # We can just as well use an upper bound for ϵ
    ϵ = ubound(Arb, ϵ)

    p0 = expansion_p0(u0)
    as = expansion_as(u0)

    # Compute enclosures of a[0], a[1], a[2] and p0
    as_enclosure = [a(u0.α) for a in as]
    p0_enclosure = p0(u0.α)

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
        res = a0C0 + as_enclosure[0] * clausen0_tail

        # Compute enclosures of the second and third Clausen functions

        # The exponent is given by e1 + u0.α = (-u0.α + p0_enclosure)
        # + u0.α = p0_enclosure
        clausen1 = C1 * abspow(x, p0_enclosure)
        for m = 1:M-1
            clausen1 += P1[2m] * abspow(x, 2m + u0.α)
        end
        clausen1 += E1 * abspow(x, 2M + u0.α)

        res += as_enclosure[1] * clausen1

        # The exponent is given by e2 + u0.α = (-u0.α + 2p0_enclosure)
        # + u0.α = 2p0_enclosure
        clausen2 = C2 * abspow(x, 2p0_enclosure)
        for m = 1:M-1
            clausen2 += P2[2m] * abspow(x, 2m + u0.α)
        end
        clausen2 += E2 * abspow(x, 2M + u0.α)

        res += as_enclosure[2] * clausen2

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

        res = (p^2 / 2 + q) / (u0.w(x) * p)

        @assert iszero(Arblib.ref(res, 0))
        @assert Arblib.contains_zero(Arblib.ref(res, 1))
        res[1] = 0

        return res
    end
end

# TODO: Implement this
function F0(u0::KdVZeroAnsatz, ::Asymptotic)
    return x -> begin
        ArbSeries((0, 0, 1))
    end
end
