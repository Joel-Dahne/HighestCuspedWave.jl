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

        zetamulα = ArbSeries((-1, stieltjes(Arb, 0), stieltjes(Arb, 1), stieltjes(Arb, 2) / 2))

        a0zeta_term = (as[0] << 1) * zetamulα

        a0clausenterm - a0zeta_term
    end

    res += sum(as[j] * clausencmzeta(x, 1 - α + j * p0) for j = 1:2)

    return res
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
(a[0] / α) * (α * zeta(1 - α))
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

            # Expansion of 2α * gamma(2α)
            gammamulα = let γ = Arb(Irrational{:γ}()), π = Arb(π)
                ArbSeries((1 // 2, -γ, γ^2 + π^2 / 6))
            end

            clausen_term *= gammamulα

            # Divide as[0] by α, perform the multiplication and then
            # multiply by α. This makes sure the degree after the
            # multiplication is correct.
            a0clausenterm = ((as[0] << 1) * clausen_term) >> 1

            zetamulα = ArbSeries((-1 // 2, stieltjes(Arb, 0),24stieltjes(Arb, 1), 2stieltjes(Arb, 2)))

            a0zeta_term = (as[0] << 1) * zetamulα

            -(a0clausenterm - a0zeta_term)
        end

        res -= sum(as[j] * clausencmzeta(x, 1 - 2α + j * p0) for j = 1:2)

        return res
    end
end
