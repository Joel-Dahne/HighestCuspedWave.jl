"""
    expansion_p0(::KdVZeroAnsatz)

Compute an expansion for `p0` at `α = 0`.

We are interested in finding `p0` solving
```
gamma(2α - p0) cospi((2α - p0) / 2) / (gamma(α - p0) * cospi((α - p0) / 2)) =
    2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```

Expanding the right hand side around `α = 0` we get from Mathematica
```
2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2)) =
    1 - γ * α + (γ^2 / 2 - π^2 / 8) * α^2 +
    (-4γ^3 + 3γ * π^2 + 28 * polygamma(2, 1)) / 24 * α^3 + O(α^4)
```
where `γ` is the Euler constant.
- **TODO:** Compute remainder term

If we let `p0 = p00 + p01 * α + p02 * α^2`, insert this into the left
hand side and expand in `α` we get
```
1 +
(digamma(-p00) + π / 2 * tan(p00 * π / 2)) * α +
(
    - 3π^2
    - 2p01 * π^2
    + 4digamma(-p00)^2
    + 12polygamma(1, -p00)
    + 8p01 * polygamma(1, -p00)
    + 4π * digamma(-p00) * tan(p00 * π / 2)
    - 2π^2 * tan(p00 * π / 2)^2
    - 2p01 * π^2 * tan(p00 * π / 2)^2
) / 8 * α^2 +
(
    - 12 * p02 * π^2
    - 18 * π^2 * digamma(-p00)
    - 12 * p01 * π^2 * digamma(-p00)
    + 8 * digamma(-p00)^3
    + 48 * p02 * polygamma(1, -p00)
    + 72 * digamma(-p00) * polygamma(1, -p00)
    + 48 * p01 * digamma(-p00) * polygamma(1, -p00)
    + 56 * polygamma(2, -p00)
    + 72 * p01 * polygamma(2, -p00)
    + 24 * p01^2 * polygamma(2, -p00)
    + 5 * π^3 * tan(p00 * π / 2)
    + 12 * p01 * π^3 * tan(p00 * π / 2)
    + 6 * p01^2 * π^3 * tan(p00 * π / 2)
    + 12 * π * digamma(-p00)^2 * tan(p00 * π / 2)
    + 36 * π * polygamma(1, -p00) * tan(p00 * π / 2)
    + 24 * p01 * π * polygamma(1, -p00) * tan(p00 * π / 2)
    - 12 * p02 * π^2 * tan(p00 * π / 2)^2
    - 12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
    - 12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
    + 6 * π^3 * tan(p00 * π / 2)^3
    + 12 * p01 * π^3 * tan(p00 * π / 2)^3
    + 6 * p01^2 * π^3 * tan(p00 * π / 2)^3
) / 48 * α^3 +
O(α^4)
```
Here we have used that `polygamma(0, x) = digamma(x)`. For `p00` we
see that we must have
```
digamma(-p00) + π / 2 * tan(p00 * π / 2) = -γ
```
For `p01` we get (after multiplying by `8`)
```
- 3π^2
- 2p01 * π^2
+ 4digamma(-p00)^2
+ 12polygamma(1, -p00)
+ 8p01 * polygamma(1, -p00)
+ 4π * digamma(-p00) * tan(p00 * π / 2)
- 2π^2 * tan(p00 * π / 2)^2
- 2p01 * π^2 * tan(p00 * π / 2)^2
= 4γ^2 - π^2
```
And for `p02` we have after multiplication by 48
```
- 12 * p02 * π^2
- 18 * π^2 * digamma(-p00)
- 12 * p01 * π^2 * digamma(-p00)
+ 8 * digamma(-p00)^3
+ 48 * p02 * polygamma(1, -p00)
+ 72 * digamma(-p00) * polygamma(1, -p00)
+ 48 * p01 * digamma(-p00) * polygamma(1, -p00)
+ 56 * polygamma(2, -p00)
+ 72 * p01 * polygamma(2, -p00)
+ 24 * p01^2 * polygamma(2, -p00)
+ 5 * π^3 * tan(p00 * π / 2)
+ 12 * p01 * π^3 * tan(p00 * π / 2)
+ 6 * p01^2 * π^3 * tan(p00 * π / 2)
+ 12 * π * digamma(-p00)^2 * tan(p00 * π / 2)
+ 36 * π * polygamma(1, -p00) * tan(p00 * π / 2)
+ 24 * p01 * π * polygamma(1, -p00) * tan(p00 * π / 2)
- 12 * p02 * π^2 * tan(p00 * π / 2)^2
- 12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
- 12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
+ 6 * π^3 * tan(p00 * π / 2)^3
+ 12 * p01 * π^3 * tan(p00 * π / 2)^3
+ 6 * p01^2 * π^3 * tan(p00 * π / 2)^3
= -8γ^3 + 6γ * π^2 + 56 * polygamma(2, 1)
```

We can solve for `p00` using numerical methods. For `p01` we can note
that the equation is linear and we immediately get
```
p01 = -(
    4γ^2 - π^2 -
    (-3π^2 + 4digamma(-p00)^2 + 12polygamma(1, -p00) + + 4π * digamma(-p00) * tan(p00 * π / 2) - 2π^2 * tan(p00 * π / 2)^2)
    ) / (- 2π^2 + 8polygamma(1, -p00) - 2π^2 * tan(p00 * π / 2)^2)
```
Similarly we get
```
p02 =
(
    (-8γ^3 + 6γ * π^2 + 56 * polygamma(2, 1)) -
    (
        - 18 * π^2 * digamma(-p00)
        - 12 * p01 * π^2 * digamma(-p00)
        + 8 * digamma(-p00)^3
        + 72 * digamma(-p00) * polygamma(1, -p00)
        + 48 * p01 * digamma(-p00) * polygamma(1, -p00)
        + 56 * polygamma(2, -p00)
        + 72 * p01 * polygamma(2, -p00)
        + 24 * p01^2 * polygamma(2, -p00)
        + 5 * π^3 * tan(p00 * π / 2)
        + 12 * p01 * π^3 * tan(p00 * π / 2)
        + 6 * p01^2 * π^3 * tan(p00 * π / 2)
        + 12 * π * digamma(-p00)^2 * tan(p00 * π / 2)
        + 36 * π * polygamma(1, -p00) * tan(p00 * π / 2)
        + 24 * p01 * π * polygamma(1, -p00) * tan(p00 * π / 2)
        - 12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
        - 12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2
        + 6 * π^3 * tan(p00 * π / 2)^3
        + 12 * p01 * π^3 * tan(p00 * π / 2)^3
        + 6 * p01^2 * π^3 * tan(p00 * π / 2)^3
    )
) / (
    - 12 * π^2
    + 48 * polygamma(1, -p00)
    - 12 * π^2 * tan(p00 * π / 2)^2
)


```

**TODO:** Possibly compute more terms in the expansion.
**TODO:** Compute remainder term.
"""
function expansion_p0(::KdVZeroAnsatz)
    γ = Arb(Irrational{:γ}())

    rhs = ArbSeries((1, -γ, γ^2 / 2 - Arb(π)^2 / 8))

    p00 = let π = Arb(π)
        f(p00) = digamma(-p00) + π / 2 * tan(p00 * π / 2) + γ
        roots, flags = ArbExtras.isolate_roots(f, Arf(1.4), Arf(1.5), depth = 20)
        @assert only(flags)
        ArbExtras.refine_root(f, Arb(roots[1]))
    end

    p01 = let π = Arb(π)
        -(
            4γ^2 - π^2 - (
                -3π^2 +
                4digamma(-p00)^2 +
                12real(polygamma(Acb(1), Acb(-p00))) +
                +4π * digamma(-p00) * tan(p00 * π / 2) - 2π^2 * tan(p00 * π / 2)^2
            )
        ) / (-2π^2 + 8real(polygamma(Acb(1), Acb(-p00))) - 2π^2 * tan(p00 * π / 2)^2)
    end

    p02 = let π = Arb(π)
        (
            (-8γ^3 + 6γ * π^2 + 56 * polygamma(2, 1)) - (
                -18 * π^2 * digamma(-p00) - 12 * p01 * π^2 * digamma(-p00) +
                8 * digamma(-p00)^3 +
                72 * digamma(-p00) * real(polygamma(Acb(1), Acb(-p00))) +
                48 * p01 * digamma(-p00) * real(polygamma(Acb(1), Acb(-p00))) +
                56 * real(polygamma(Acb(2), Acb(-p00))) +
                72 * p01 * real(polygamma(Acb(2), Acb(-p00))) +
                24 * p01^2 * real(polygamma(Acb(2), Acb(-p00))) +
                5 * π^3 * tan(p00 * π / 2) +
                12 * p01 * π^3 * tan(p00 * π / 2) +
                6 * p01^2 * π^3 * tan(p00 * π / 2) +
                12 * π * digamma(-p00)^2 * tan(p00 * π / 2) +
                36 * π * real(polygamma(Acb(1), Acb(-p00))) * tan(p00 * π / 2) +
                24 * p01 * π * real(polygamma(Acb(1), Acb(-p00))) * tan(p00 * π / 2) -
                12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2 -
                12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2 +
                6 * π^3 * tan(p00 * π / 2)^3 +
                12 * p01 * π^3 * tan(p00 * π / 2)^3 +
                6 * p01^2 * π^3 * tan(p00 * π / 2)^3
            )
        ) / (
            -12 * π^2 + 48 * real(polygamma(Acb(1), Acb(-p00))) -
            12 * π^2 * tan(p00 * π / 2)^2
        )
    end

    p0 = ArbSeries((p00, p01, p02))

    return p0
end

"""
    expansion_as(u0::KdVZeroAnsatz)

Compute expansions for `a[i]` for `i = 1:3` at `α = 0`.

# Computing `a[0]`
From Mathematica we have
```
a[0] = α - π^2 / 12 * a^3 + polygamma(2, 1) * a^4 - 11π^4 / 360 * a^5 + O(a^6)
```
- **TODO:** Compute remainder term.

# Setting up a linear system for `a[1]` and `a[2]`
We can get the values for `a[1]` and `a[2]` in terms of a linear
system depending on `α`, `a[0]` and `p0`.

Computing the asymptotic expansion of `D(u0)` we have that the four
leading terms have the `x`-factors are, in order,
- `x^(-2α)`
- `x^(-2α + p0)`
- `x^2`
- `x^(-α + 2)
The first two are zero because of the choice of `a[0]` and `p0`. We
want to find `a[1]` and `a[2]` such that the last two are also zero.

The coefficient in front of `x^2` is given by
```
L₁⁰ = -1 / 2 * sum(a[j] * zeta(1 - 2α + j * p0 - 2) for j = 0:u0.N0) =
    -1 / 2 * (a[0] * zeta(1 - 2α - 2) + a[1] * zeta(1 - 2α + p0 - 2) + a[2] * zeta(1 - 2α + 2p0 - 2))
```
and for `x^(-α + 2)`
```
a₀⁰ * K₁⁰ = gamma(α) * sinpi((1 - α) / 2) * a[0] * (-1 / 2) * sum(a[j] * zeta(1 - α + j * p0 - 2) for j = 0:u0.N0) =
    -gamma(α) * sinpi((1 - α) / 2) * a[0] / 2 * (
        a[0] * zeta(1 - α - 2) + a[1] * zeta(1 - α + p0 - 2) + a[2] * zeta(1 - α + 2p0 - 2)
    )
```

If we let
```
A1 = -zeta(1 - 2α + p0 - 2) / 2
A2 = -zeta(1 - 2α + 2p0 - 2) / 2
C1 = -a[0] * zeta(1 - 2α - 2) / 2

B1 = -gamma(α) * sinpi((1 - α) / 2) * a[0] * zeta(1 - α + p0 - 2) / 2
B2 = -gamma(α) * sinpi((1 - α) / 2) * a[0] * zeta(1 - α + 2p0 - 2) / 2
C2 = -gamma(α) * sinpi((1 - α) / 2) * a[0]^2 * zeta(1 - α - 2) / 2
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

# Simplifying the linear system
To make the linear system easier to evaluate we begin by simplifying
some of the expressions.

We have
```
d = A1 * B2 - A2 * B1 =
    (-zeta(1 - 2α + p0 - 2) / 2) * (-gamma(α) * sinpi((1 - α) / 2) * a[0] * zeta(1 - α + 2p0 - 2) / 2)
    - (-zeta(1 - 2α + 2p0 - 2) / 2) * (-gamma(α) * sinpi((1 - α) / 2) * a[0] * zeta(1 - α + p0 - 2) / 2) =
    zeta(1 - 2α + p0 - 2) * gamma(α) * sinpi((1 - α) / 2) * a[0] * zeta(1 - α + 2p0 - 2) / 4
    - zeta(1 - 2α + 2p0 - 2) * gamma(α) * sinpi((1 - α) / 2) * a[0] * zeta(1 - α + p0 - 2) / 4 =
    gamma(α) * sinpi((1 - α) / 2) * a[0] / 4 * (
        zeta(1 - 2α + p0 - 2) * zeta(1 - α + 2p0 - 2) -
        zeta(1 - 2α + 2p0 - 2) * zeta(1 - α + p0 - 2)
    )
```

Then we also have
```
v1 = -B2 * C1 + A2 * C2 =
    -(-gamma(α) * sinpi((1 - α) / 2) * a[0] * zeta(1 - α + 2p0 - 2) / 2) *
    (-a[0] * zeta(1 - 2α - 2) / 2) +
    (-zeta(1 - 2α + 2p0 - 2) / 2) *
    (-gamma(α) * sinpi((1 - α) / 2) * a[0]^2 * zeta(1 - α - 2) / 2) =
    gamma(α) * sinpi((1 - α) / 2) * a[0]^2 / 4 * (
        zeta(1 - 2α + 2p0 - 2) * zeta(1 - α - 2) -
        zeta(1 - α + 2p0 - 2) * zeta(1 - 2α - 2)
    )
```
and
```
v2 = B1 * C1 - A1 * C2 =
    (-gamma(α) * sinpi((1 - α) / 2) * a[0] * zeta(1 - α + p0 - 2) / 2) *
    (-a[0] * zeta(1 - 2α - 2) / 2) -
    (-zeta(1 - 2α + p0 - 2) / 2) *
    (-gamma(α) * sinpi((1 - α) / 2) * a[0]^2 * zeta(1 - α - 2) / 2) =
    -gamma(α) * sinpi((1 - α) / 2) * a[0]^2 / 4 * (
        zeta(1 - 2α + p0 - 2) * zeta(1 - α - 2) -
        zeta(1 - α + p0 - 2) * zeta(1 - 2α - 2)
    )
```

If we let
```
q = gamma(α) * sinpi((1 - α) / 2) * a[0]
```
and
```
z(i, j) = zeta(-1 - i * α + j * p0)

z1 = z(2, 1) * z(1, 2) - z(2, 2) * z(1, 1)

z2 = z(2, 2) * z(1, 0) - z(1, 2) * z(2, 0)

z3 = z(2, 1) * z(1, 0) - z(1, 1) * z(2, 0)
```
we can write this as
```
d = q / 4 * z1

v1 = q * a[0] / 4 * z2

v2 = -q * a[0] / 4 * z3
```

# Evaluating the linear system
Unfortunately using the above formulas don't work for direct
evaluation, there are several indeterminate values we have to take
care of.

## Handling `q`
To begin with we can note that
```
q = gamma(α) * sinpi((1 - α) / 2) * a[0]
```
has a removable singularity. Plugging in the expression
```
a[0] = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```
we can rewrite it as
```
q = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```
Which also occurs in [`expansion_p0`](@ref) and the expansion its given
by
```
q = 1 - γ * α + (γ^2 / 2 - π^2 / 8) * α^2 + O(α^3)
```
- **TODO:** Compute remainder term

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
d = q * z1 / 4
v1 = q * a[0] * z2 / 4
v2 = -q * a[0] * z3 / 4
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
function expansion_as(u0::KdVZeroAnsatz)
    # Compute expansion of a0
    a0 = ArbSeries((
        0,
        1,
        0,
        -Arb(π)^2 / 12,
        real(polygamma(Acb(2), Acb(1))),
        -11Arb(π)^4 / 360,
    ),)

    # Compute expansions of p0 and α
    p0 = expansion_p0(u0)
    α = ArbSeries((0, 1), degree = Arblib.degree(a0))

    # Expansion of gamma(α) * sinpi((1 - α) / 2) * a[0]
    q = let γ = Arb(Irrational{:γ}()), π = Arb(π)
        ArbSeries((1, -γ, γ^2 / 2 - π^2 / 8))
    end

    z(i, j) = zeta(-1 - i * α + j * p0)

    z1 = z(2, 1) * z(1, 2) - z(2, 2) * z(1, 1)
    z2 = z(2, 2) * z(1, 0) - z(1, 2) * z(2, 0)
    z3 = z(2, 1) * z(1, 0) - z(1, 1) * z(2, 0)

    # The constant coefficients for z1, z2 and z3 are all exactly
    # equal to zero.
    @assert all(Arblib.contains_zero(z[0]) for z in (z1, z2, z3))
    z1[0] = z2[0] = z3[0] = 0

    # Factor out α from z1 and multiply back afterwards
    d = (q * (z1 << 1) / 4) >> 1
    # Factor out α from a[0] and z2 and multiply back afterwards
    v1 = (q * (a0 << 1) * (z2 << 1) / 4) >> 2
    # Factor out α from a[0] and z3 and multiply back afterwards
    v2 = -(q * (a0 << 1) * (z3 << 1) / 4) >> 2

    # Factor out α^2 from v1 and v2 and α from d, multiply back one α afterwards
    a1 = ((v1 << 2) / (d << 1)) >> 1
    a2 = ((v2 << 2) / (d << 1)) >> 1

    return OffsetVector([a0, a1, a2], 0:2)
end
