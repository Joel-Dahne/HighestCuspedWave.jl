"""
    expansion_p0(::Type{KdVZeroAnsatz}, α::Arb; degree::Integer = 2)

Compute an expansion for `p0` at `α = 0` of the given degree. The last
term is a remainder term which ensures that evaluating the expansion
gives an enclosure of `p0` for all values in the interval `α`.

We are interested in finding `p0` solving
```
gamma(2α - p0) * cospi((2α - p0) / 2) / (gamma(α - p0) * cospi((α - p0) / 2)) =
    2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```

Expanding the right hand side around `α = 0` we get from Mathematica,
with
```
Series[2 Gamma[2 a]*Cos[Pi*a]/(Gamma[a]*Cos[Pi*a/2]), {a, 0, 3}]
```
```
2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2)) =
    1 - γ * α + (γ^2 / 2 - π^2 / 8) * α^2 +
    (-4γ^3 + 3γ * π^2 + 28 * polygamma(2, 1)) / 24 * α^3 + O(α^4)
```
where `γ` is the Euler constant.

If we let `p0 = p00 + p01 * α + p02 * α^2`, insert this into the left
hand side and expand in `α` we get, using
```
Series[(Gamma[2 a - p0]* Cos[Pi*(2 a - p0)/2]/(Gamma[a - p0]*Cos[Pi (a - p0)/2])) /.
    p0 -> (p00 + p01*a + p02*a^2), {a, 0, 3}]
```
```
1 +
(digamma(-p00) + π / 2 * tan(p00 * π / 2)) * α +
(
    - 3π^2
    + 2p01 * π^2
    + 4digamma(-p00)^2
    + 12polygamma(1, -p00)
    - 8p01 * polygamma(1, -p00)
    + 4π * digamma(-p00) * tan(p00 * π / 2)
    - 2π^2 * tan(p00 * π / 2)^2
    + 2p01 * π^2 * tan(p00 * π / 2)^2
) / 8 * α^2 +
(
    + 12 * p02 * π^2
    - 18 * π^2 * digamma(-p00)
    + 12 * p01 * π^2 * digamma(-p00)
    + 8 * digamma(-p00)^3
    - 48 * p02 * polygamma(1, -p00)
    + 72 * digamma(-p00) * polygamma(1, -p00)
    - 48 * p01 * digamma(-p00) * polygamma(1, -p00)
    + 56 * polygamma(2, -p00)
    - 72 * p01 * polygamma(2, -p00)
    + 24 * p01^2 * polygamma(2, -p00)
    + 5 * π^3 * tan(p00 * π / 2)
    - 12 * p01 * π^3 * tan(p00 * π / 2)
    + 6 * p01^2 * π^3 * tan(p00 * π / 2)
    + 12 * π * digamma(-p00)^2 * tan(p00 * π / 2)
    + 36 * π * polygamma(1, -p00) * tan(p00 * π / 2)
    - 24 * p01 * π * polygamma(1, -p00) * tan(p00 * π / 2)
    + 12 * p02 * π^2 * tan(p00 * π / 2)^2
    - 12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
    + 12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
    + 6 * π^3 * tan(p00 * π / 2)^3
    - 12 * p01 * π^3 * tan(p00 * π / 2)^3
    + 6 * p01^2 * π^3 * tan(p00 * π / 2)^3
) / 48 * α^3 +
O(α^4)
```
Here we have used that `polygamma(0, x) = digamma(x)`. For `p00` we
see that we must have
```
digamma(-p00) + π / 2 * tan(p00 * π / 2) = -γ
```
For `p01` we get, after multiplying by `8`,
```
- 3π^2
+ 2p01 * π^2
+ 4digamma(-p00)^2
+ 12polygamma(1, -p00)
- 8p01 * polygamma(1, -p00)
+ 4π * digamma(-p00) * tan(p00 * π / 2)
- 2π^2 * tan(p00 * π / 2)^2
+ 2p01 * π^2 * tan(p00 * π / 2)^2
= 4γ^2 - π^2
```
And for `p02` we have, after multiplication by 48,
```
+ 12 * p02 * π^2
- 18 * π^2 * digamma(-p00)
+ 12 * p01 * π^2 * digamma(-p00)
+ 8 * digamma(-p00)^3
- 48 * p02 * polygamma(1, -p00)
+ 72 * digamma(-p00) * polygamma(1, -p00)
- 48 * p01 * digamma(-p00) * polygamma(1, -p00)
+ 56 * polygamma(2, -p00)
- 72 * p01 * polygamma(2, -p00)
+ 24 * p01^2 * polygamma(2, -p00)
+ 5 * π^3 * tan(p00 * π / 2)
- 12 * p01 * π^3 * tan(p00 * π / 2)
+ 6 * p01^2 * π^3 * tan(p00 * π / 2)
+ 12 * π * digamma(-p00)^2 * tan(p00 * π / 2)
+ 36 * π * polygamma(1, -p00) * tan(p00 * π / 2)
- 24 * p01 * π * polygamma(1, -p00) * tan(p00 * π / 2)
+ 12 * p02 * π^2 * tan(p00 * π / 2)^2
- 12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
+ 12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
+ 6 * π^3 * tan(p00 * π / 2)^3
- 12 * p01 * π^3 * tan(p00 * π / 2)^3
+ 6 * p01^2 * π^3 * tan(p00 * π / 2)^3
= -8γ^3 + 6γ * π^2 + 56 * polygamma(2, 1)
```

We can solve for `p00` using numerical methods. for `p01` we can note
that the equation is linear and we immediately get
```
p01 = -(
    4γ^2 - π^2 -
    (-3π^2 + 4digamma(-p00)^2 + 12polygamma(1, -p00) + + 4π * digamma(-p00) * tan(p00 * π / 2) - 2π^2 * tan(p00 * π / 2)^2)
    ) / (- 2π^2 + 8polygamma(1, -p00) - 2π^2 * tan(p00 * π / 2)^2)
```
similarly we get
```
p02 = (
    (-8γ^3 + 6γ * π^2 + 56 * polygamma(2, 1))
    - (
        - 18 * π^2 * digamma(-p00)
        + 12 * p01 * π^2 * digamma(-p00)
        + 8 * digamma(-p00)^3
        + 72 * digamma(-p00) * polygamma(1, -p00)
        - 48 * p01 * digamma(-p00) * polygamma(1, -p00)
        + 56 * polygamma(2, -p00)
        - 72 * p01 * polygamma(2, -p00)
        + 24 * p01^2 * polygamma(2, -p00)
        + 5 * π^3 * tan(p00 * π / 2)
        - 12 * p01 * π^3 * tan(p00 * π / 2)
        + 6 * p01^2 * π^3 * tan(p00 * π / 2)
        + 12 * π * digamma(-p00)^2 * tan(p00 * π / 2)
        + 36 * π * polygamma(1, -p00) * tan(p00 * π / 2)
        - 24 * p01 * π * polygamma(1, -p00) * tan(p00 * π / 2)
        - 12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
        + 12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
        + 6 * π^3 * tan(p00 * π / 2)^3
        - 12 * p01 * π^3 * tan(p00 * π / 2)^3
        + 6 * p01^2 * π^3 * tan(p00 * π / 2)^3
    )
) / (
    + 12 * π^2
    - 48 * polygamma(1, -p00)
    + 12 * π^2 * tan(p00 * π / 2)^2
)
```

- **IMPROVE:** Possibly compute more terms in the expansion.
- **TODO:** Compute remainder term.
"""
function expansion_p0(::Type{KdVZeroAnsatz}, α::Arb; degree::Integer = 2)
    degree <= 2 || throw(ArgumentError("only supports degree up to 2"))

    p00 = let π = Arb(π), γ = Arb(Irrational{:γ}())
        f(p00) = digamma(-p00) + π / 2 * tan(p00 * π / 2) + γ
        roots, flags = ArbExtras.isolate_roots(f, Arf(1.4), Arf(1.5), depth = 20)
        @assert only(flags)
        ArbExtras.refine_root(f, Arb(roots[1]))
    end

    p01 = let π = Arb(π), γ = Arb(Irrational{:γ}())
        -(
            4γ^2 - π^2 - (
                -3π^2 +
                4digamma(-p00)^2 +
                12real(polygamma(Acb(1), Acb(-p00))) +
                +4π * digamma(-p00) * tan(p00 * π / 2) - 2π^2 * tan(p00 * π / 2)^2
            )
        ) / (-2π^2 + 8real(polygamma(Acb(1), Acb(-p00))) - 2π^2 * tan(p00 * π / 2)^2)
    end

    p02 = let π = Arb(π), γ = Arb(Irrational{:γ}())
        (
            (-8γ^3 + 6γ * π^2 + 56 * real(polygamma(Acb(2), Acb(1)))) - (
                -18 * π^2 * digamma(-p00) +
                12 * p01 * π^2 * digamma(-p00) +
                8 * digamma(-p00)^3 +
                72 * digamma(-p00) * real(polygamma(Acb(1), Acb(-p00))) -
                48 * p01 * digamma(-p00) * real(polygamma(Acb(1), Acb(-p00))) +
                56 * real(polygamma(Acb(2), Acb(-p00))) -
                72 * p01 * real(polygamma(Acb(2), Acb(-p00))) +
                24 * p01^2 * real(polygamma(Acb(2), Acb(-p00))) +
                5 * π^3 * tan(p00 * π / 2) - 12 * p01 * π^3 * tan(p00 * π / 2) +
                6 * p01^2 * π^3 * tan(p00 * π / 2) +
                12 * π * digamma(-p00)^2 * tan(p00 * π / 2) +
                36 * π * real(polygamma(Acb(1), Acb(-p00))) * tan(p00 * π / 2) -
                24 * p01 * π * real(polygamma(Acb(1), Acb(-p00))) * tan(p00 * π / 2) -
                12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2 +
                12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2 +
                6 * π^3 * tan(p00 * π / 2)^3 - 12 * p01 * π^3 * tan(p00 * π / 2)^3 +
                6 * p01^2 * π^3 * tan(p00 * π / 2)^3
            )
        ) / (
            +12 * π^2 - 48 * real(polygamma(Acb(1), Acb(-p00))) +
            12 * π^2 * tan(p00 * π / 2)^2
        )
    end

    # Expansion without remainder term
    p0 = ArbSeries((p00, p01, p02); degree)

    # FIXME: Properly implement this. Now we just widen the last
    # coefficient so that we get an enclosure for a lower bound of α
    if !iszero(α)
        error = findp0(lbound(Arb, α)) - p0(lbound(Arb, α))
        p0[degree] += Arblib.add_error!(zero(error), error / lbound(Arb, α)^degree)
    end

    return p0
end

"""
    expansion_as(::Type{KdVZeroAnsatz}, α::Arb; degree::integer = 2)

Compute expansions for `a[i]` for `i = 1:3` at `α = 0` of the given
degree. The last term is a remainder term which ensures that
evaluating the expansion gives an enclosure of `a[i]` for all values
in the interval `α`.

The expansion of `a[0]` is computed to a degree one higher than the
other. The reason for this is that the constant term in `a[0]` is zero
and in many cases we divide `a[0]` by `α` and still want to have
sufficiently high degree.

# Computing `a[0]`
We have
```
a[0] = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```
At `α = 0` both gamma functions have a pole, we can instead rewrite it
in terms of the reciprocal gamma function as
```
a[0] = 2rgamma(α)^2 * cospi(α) / (rgamma(2α) * cospi(α / 2)^2)
```
We can compute `2cospi(α) / cospi(α / 2)^2` whereas for `rgamma(α)^2 /
rgamma(2α)` we have to handle the removable singularity.
- **TODO:** Compute remainder term. The only problematic part is
  `rgamma(α) / rgamma(2α)`. We want a very tight enclosure.

# Setting up a linear system for `a[1]` and `a[2]`
We can get the values for `a[1]` and `a[2]` in terms of a linear
system depending on `α`, `a[0]` and `p0`.

Computing the asymptotic expansion of `D(u0)` we have that the four
leading terms have the `x`-factors, in order,
- `x^(-2α)`
- `x^(-2α + p0)`
- `x^2`
- `x^(-α + 2)
The first two are zero because of the choice of `a[0]` and `p0`. We
want to find `a[1]` and `a[2]` such that the last two are also zero.

The coefficient in front of `x^2` is given by
```
-L₁⁰ = 1 / 2 * sum(a[j] * zeta(1 - 2α + j * p0 - 2) for j = 0:u0.N0) =
    1 / 2 * (a[0] * zeta(1 - 2α - 2) + a[1] * zeta(1 - 2α + p0 - 2) + a[2] * zeta(1 - 2α + 2p0 - 2))
```
and for `x^(-α + 2)`
```
a₀⁰ * K₁⁰ = gamma(α) * sinpi((1 - α) / 2) * a[0] * (-1 / 2) * sum(a[j] * zeta(1 - α + j * p0 - 2) for j = 0:u0.N0) =
    -gamma(α) * sinpi((1 - α) / 2) * a[0] / 2 * (
        a[0] * zeta(1 - α - 2) + a[1] * zeta(1 - α + p0 - 2) + a[2] * zeta(1 - α + 2p0 - 2)
    )
```
Since we are looking for zeros we cancel some of the factors and
simplify to get the two equations
```
a[0] * zeta(-1 - 2α) + a[1] * zeta(-1 - 2α + p0) + a[2] * zeta(-1 - 2α + 2p0) = 0
a[0] * zeta(-1 - α) + a[1] * zeta(-1 - α + p0) + a[2] * zeta(-1 - α + 2p0) = 0
```

If we let
```
A1 = zeta(-1 - 2α + p0)
A2 = zeta(-1 - 2α + 2p0)
C1 = a[0] * zeta(-1 - 2α)

B1 = zeta(-1 - α + p0)
B2 = zeta(-1 - α + 2p0)
C2 = a[0] * zeta(-1 - α)
```
we can write this as the linear system
```
A1 * a[1] + A2 * a[2] = -C1
B1 * a[1] + B2 * a[2] = -C2
```
Solving for `a[1]` and `a[2]` we get
```
a[1] = (-B2 * C1 + A2 * C2) / (A1 * B2 - A2 * B1)
a[2] = (B1 * C1 - A1 * C2) / (A1 * B2 - A2 * B1)
```
Or if we let
```
v1 = -B2 * C1 + A2 * C2
v2 = B1 * C1 - A1 * C2
d = A1 * B2 - A2 * B1
```
we can write it as
```
a[1] = v1 / d
a[2] = v2 / d
```

# Simplifying `d`, `v1` and `v2`
To make the linear system easier to evaluate we begin by simplifying
some of the expressions.

We have
```
d = A1 * B2 - A2 * B1 =
    zeta(-1 - 2α + p0) * zeta(-1 - α + 2p0) - zeta(-1 - 2α + 2p0) * zeta(-1 - α + p0) =
```
```
v1 = -B2 * C1 + A2 * C2 =
    a[0] * (-zeta(-1 - α + 2p0) * zeta(-1 - 2α) + zeta(-1 - 2α + 2p0) * zeta(-1 - α))
```
and
```
v2 = B1 * C1 - A1 * C2 =
    a[0] * (zeta(-1 - 2α + p0) * zeta(-1 - α) - zeta(-1 - α + p0) * zeta(-1 - 2α))
```

If we let
```
z(i, j) = zeta(-1 - i * α + j * p0)

z1 = z(2, 1) * z(1, 2) - z(2, 2) * z(1, 1)

z2 = z(2, 2) * z(1, 0) - z(1, 2) * z(2, 0)

z3 = z(2, 1) * z(1, 0) - z(1, 1) * z(2, 0)
```
we can write this as
```
d = z1
v1 = a[0] * z2
v2 = -a[0] * z3
```

# Evaluating the linear system
Unfortunately using the above formulas don't work for direct
evaluation, there are several indeterminate values we have to take
care of.

## Handling `z(i, j)`s
The constant term in the expansion of `z(i, j)` doesn't depend on the
value of `i`. This means that the constant terms of
```
z1 = z(2, 1) * z(1, 2) - z(2, 2) * z(1, 1)

z2 = z(2, 2) * z(1, 0) - z(1, 2) * z(2, 0)

z3 = z(2, 1) * z(1, 0) - z(1, 1) * z(2, 0)
```
all are zero since they exactly cancel out.

## Computing `d`, `v1` and `v2`
Recall that
```
d = z1
v1 = a[0] * z2
v2 = -a[0] * z3
```
If we naively compute these terms with `ArbSeries` the degree after
the multiplication will be the minimum degree of all factors. By
factoring out `α` from `a[0]` as well as `z1`, `z2` and `z3`,
performing the multiplication and then multiplying back the `α` we can
avoid this issue.

## Handling division by `d`
We want to compute
```
a[1] = v1 / d
a[2] = v2 / d
```
but since `v1[0] = v2[0] = d[0] = 0` we first need to cancel one `α`
from all of them. To get the higher order terms we factor out one more
`α` from `v1` and `v2` and multiply it back afterwards.
"""
function expansion_as(::Type{KdVZeroAnsatz}, α::Arb; degree::Integer = 2)
    # Expansion of a0 without remainder term
    a0 = let α = ArbSeries((0, 1), degree = degree + 1)
        # rgamma(α)^2 / rgamma(2α) handling the removable singularity
        g = let α = ArbSeries(α, degree = Arblib.degree(α) + 1)
            (rgamma(α)^2 << 1) / (rgamma(2α) << 1)
        end
        2cospi(α) / cospi(α / 2)^2 * g
    end

    # FIXME: Properly implement this. Now we just widen the last
    # coefficient so that we get an enclosure for a lower bound of α
    if !iszero(α)
        error = finda0(lbound(Arb, α)) - a0(lbound(Arb, α))
        a0[degree+1] += Arblib.add_error!(zero(error), error / lbound(Arb, α)^(degree + 1))
    end

    # Compute expansions of p0
    p0 = expansion_p0(KdVZeroAnsatz, α; degree)

    α_s = ArbSeries((0, 1); degree) # Series expansion of α
    z(i, j) = compose_with_remainder(zeta, -1 - i * α_s + j * p0, α)

    z1 = mul_with_remainder(z(2, 1), z(1, 2), α) - mul_with_remainder(z(2, 2), z(1, 1), α)
    z2 = mul_with_remainder(z(2, 2), z(1, 0), α) - mul_with_remainder(z(1, 2), z(2, 0), α)
    z3 = mul_with_remainder(z(2, 1), z(1, 0), α) - mul_with_remainder(z(1, 1), z(2, 0), α)

    # The constant coefficients for z1, z2 and z3 are all exactly
    # equal to zero.
    @assert all(Arblib.contains_zero(z[0]) for z in (z1, z2, z3))
    z1[0] = z2[0] = z3[0] = 0

    d = z1
    # Factor out α from a[0] and z2 and multiply back afterwards
    v1 = mul_with_remainder(a0 << 1, z2 << 1, α) >> 2
    # Factor out α from a[0] and z3 and multiply back afterwards
    v2 = -mul_with_remainder(a0 << 1, z3 << 1, α) >> 2

    # Factor out α^2 from v1 and v2 and α from d, multiply back one α afterwards
    a1 = div_with_remainder((v1 << 2), (d << 1), α) >> 1
    a2 = div_with_remainder((v2 << 2), (d << 1), α) >> 1

    return OffsetVector([a0, a1, a2], 0:2)
end
