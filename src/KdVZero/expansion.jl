"""
    expansion_p0(::Type{KdVZeroAnsatz}, α0::Arb, interval::Arb; degree::Integer = 2)

Compute an expansion for `p0` at the point `α0` of the given degree.
The last term is a remainder term which ensures that evaluating the
expansion gives an enclosure of `p0` for all `α ∈ interval`.

We are interested in finding `p0` solving
```
gamma(2α - p0) * cospi((2α - p0) / 2) / (gamma(α - p0) * cospi((α - p0) / 2)) =
    2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```
If we consider `p0` as a function of `α` with the expansion
```
p0(α) = p00 + p01 * (α - α) + p02 * (α - α0)^2 + ...
```
plugging this into the left hand side and expanding both sides at the
point `α0` we can solve for the coefficients `p00, p01, p02, ...`.

We split this into two cases, one for `α0 < 0` and one for `α0 = 0`.

# `α0 < 0`
Let
```
f(α, p0) = gamma(2α - p0) * cospi((2α - p0) / 2) / (gamma(α - p0) * cospi((α - p0) / 2))

g(α) = 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
```

We can compute `p00` as the solution to `f(α0, p00) = g(α0)`, which we
can compute directly using [`findp0`](@ref).

For `p01` we differentiate the equation ones with respect to `α`,
giving us
```
f_1(α, p0) + p0' * f_2(α, p0) = g'(α)
```
where `f_1` is `f` differentiated w.r.t to the first argument and
`f_2` w.r.t the second argument. Plugging in `α0` we get
```
f_1(α0, p00) + p01 * f_2(α0, p00) = g'(α0)
```
Solving for `p01` we have
```
p01 = (g'(α0) - f_1(α0, p00)) / f_2(α0, p00)
```
We can compute all required derivatives using `ArbSeries`.

For `p02` we differentiate the equation ones more w.r.t `α`, giving us
```
f_11(α, p0) + 2p0' * f_12(α, p0) + (p0')^2 * f_22(α, p0) + p0'' * f_2(α, p0) = g''(α)
```
Inserting `α0` we have
```
f_11(α0, p00) + 2p01 * f_12(α0, p00) + p01^2 * f_22(α0, p00) + 2 * p02 * f_2(α0, p00) = g''(α0)
```
and solving for `p02` gives
```
p02 = (g''(α0) - f_11(α0, p00) + 2p01 * f_12(α0, p00) + p01^2 * f_22(α0, p00)) / 2f_2(α0, p00)
```
We can compute most derivatives using `ArbSeries`, the exception is
`f_12` since Arb doesn't support series in two variables. Instead we
compute the derivative w.r.t. the second argument by hand and then
different this w.r.t the first argument using `ArbSeries`.
Differentiation gives
```
f_2(α, p0) =
    gamma(2α - p0) * digamma(α - p0) * cospi((2α - p0) / 2) /
    (gamma(α - p0) * cospi((α - p0) / 2)) -
    (gamma(2α - p0) * digamma(2α - p0) * cospi((2α - p0) / 2)) /
    (gamma(α - p0) * cospi((α - p0) / 2)) +
    π * (gamma(2α - p0) * sinpi((2α - p0) / 2)) /
    (2gamma(α - p0) * cospi((α - p0) / 2)) -
    π * (gamma(2α - p0) * cospi((2α - p0) / 2) * tanpi((α - p0) / 2)) /
    (2gamma(α - p0) * cospi((α - p0) / 2))
```

# `α0 = 0`
Recall that we are interested in solving
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

# Remainder term
To enclose the remainder term we want to find `R` such that
```
p00 + p01 * (α - α0) + p02 * (α - α0)^2 + R * (α - α0)^3
```
gives an enclosure of
```
p0(α) = p00 + p01 * (α - α0) + p02 * (α - α0)^2 + ...
```
for all `α ∈ interval`.

This is equivalent to saying that for all `α ∈ interval` there is `r ∈
R` such that
```
f(α, p00 + p01 * (α - α0) + p02 * (α - α0)^2 + r * (α - α0)^3) = 0
```
We denote the left hand side as a function of `α` and `R` by `g(α,
R)`. Since we expand to degree `degree` we expect the zero to be of
order `degree + 1`, we therefore solve
```
g(α, R) / (α - α0)^(degree + 1) = 0
```
If `α0 = 0` the zero is instead of degree `degree + 2`, giving us
`g(α, R) / (α - α0)^(degree + 2) = 0`.

Let `g_div_α(α, R) = g(α, R) / (α - α0)^k` where `k` is the degree of
the zero. Given a guess for `R` we can verify it by checking that
```
g_div_α(interval, lbound(R)) < 0 < g_div_α(interval, ubound(R))
```
Or possibly reversing the inequalities.

We can get a guess for `R` by computing the zeros of
`g_div_α(lbound(interval), R)` and `g_div_α(ubound(α), R)` and taking
the convex hull of them. When computing these zeros we don't have
access to derivatives w.r.t. `R` since we have to use `ArbSeries` to
compute the derivatives w.r.t. `α`, we therefore use
[`ArbExtras.refine_root_bisection`](@ref), which doesn't need access
to the derivative.

For `α0 != 0` we can evaluate `g(α, R)` directly using `ArbSeries`.
For `α0 = 0` we have to handle the removable singularity from
`gamma(2α) / gamma(α)`.
"""
function expansion_p0(::Type{KdVZeroAnsatz}, α0::Arb, interval::Arb; degree::Integer = 2)
    degree <= 2 || throw(ArgumentError("only supports degree up to 2"))

    p00 = if iszero(α0)

    else

    end

    if iszero(α0)
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
            ) /
            (-2π^2 + 8real(polygamma(Acb(1), Acb(-p00))) - 2π^2 * tan(p00 * π / 2)^2)
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
                    24 * p01 * π * real(polygamma(Acb(1), Acb(-p00))) * tan(p00 * π / 2) - 12 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2 +
                    12 * p01 * π^2 * digamma(-p00) * tan(p00 * π / 2)^2 +
                    6 * π^3 * tan(p00 * π / 2)^3 - 12 * p01 * π^3 * tan(p00 * π / 2)^3 +
                    6 * p01^2 * π^3 * tan(p00 * π / 2)^3
                )
            ) / (
                +12 * π^2 - 48 * real(polygamma(Acb(1), Acb(-p00))) +
                12 * π^2 * tan(p00 * π / 2)^2
            )
        end
    else
        p00 = findp0(α0)

        rhs = let α = ArbSeries((α0, 1); degree)
            2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
        end

        f_1, f_11 = let α = ArbSeries((α0, 1), degree = 2), p0 = p00
            res = (
                gamma(2α - p0) * cospi((2α - p0) / 2) /
                (gamma(α - p0) * cospi((α - p0) / 2))
            )
            res[1], 2res[2]
        end

        f_2, f_22 = let α = α0, p0 = ArbSeries((p00, 1), degree = 2)
            res = (
                gamma(2α - p0) * cospi((2α - p0) / 2) /
                (gamma(α - p0) * cospi((α - p0) / 2))
            )
            res[1], 2res[2]
        end

        f_12 = let p0 = p00, α = ArbSeries((α0, 1), degree = 1)
            res =
                gamma(2α - p0) * digamma(α - p0) * cospi((2α - p0) / 2) /
                (gamma(α - p0) * cospi((α - p0) / 2)) -
                (gamma(2α - p0) * digamma(2α - p0) * cospi((2α - p0) / 2)) /
                (gamma(α - p0) * cospi((α - p0) / 2)) +
                Arb(π) * (gamma(2α - p0) * sinpi((2α - p0) / 2)) /
                (2gamma(α - p0) * cospi((α - p0) / 2)) -
                Arb(π) *
                (gamma(2α - p0) * cospi((2α - p0) / 2) * tan(Arb(π) * (α - p0) / 2)) /
                (2gamma(α - p0) * cospi((α - p0) / 2))

            res[1]
        end

        p01 = (rhs[1] - f_1) / f_2

        p02 = (2rhs[2] - f_11 - 2p01 * f_12 - p01^2 * f_22) / 2f_2
    end

    # Expansion without remainder term
    p0 = ArbSeries((p00, p01, p02); degree)

    if !iszero(α0)
        # We double check that solution solves the equation, this is
        # only to catch potential bugs.
        lhs = let α = ArbSeries((α0, 1); degree)
            gamma(2α - p0) * cospi((2α - p0) / 2) / (gamma(α - p0) * cospi((α - p0) / 2))
        end
        rhs = let α = ArbSeries((α0, 1); degree)
            2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))
        end
        @assert Arblib.overlaps(lhs, rhs)
    end

    if !iszero(radius(interval))
        # Function we want to find zero of
        g(α::ArbSeries, R::Arb) =
            let p0 = p00 + p01 * (α - α0) + p02 * (α - α0)^2 + R * (α - α0)^3
                if Arblib.contains_zero(α[0])
                    # Enclosure of rgamma(α) / α
                    rgamma1_div_v = fx_div_x(α -> rgamma(α), α, extra_degree = 2)

                    # Enclosure of rgamma(2α) / α
                    rgamma2_div_v = fx_div_x(α -> rgamma(2α), α, extra_degree = 2)

                    # Enclosure of gamma(2α) / gamma(α)
                    gammadivgamma = rgamma1_div_v / rgamma2_div_v
                else
                    # Enclosure of gamma(2α) / gamma(α)
                    # Use higher precision to allow for accurate
                    # evaluation near the removable singularity
                    gammadivgamma = let α = setprecision(α, 4precision(α))
                        # Write in terms of rgamma instead of gamma
                        # since this gives better enclosures
                        setprecision(rgamma(α) / rgamma(2α), precision(Arb))
                    end
                end

                return gamma(2α - p0) * cospi((2α - p0) / 2) /
                       (gamma(α - p0) * cospi((α - p0) / 2)) -
                       2cospi(α) / cospi(α / 2) * gammadivgamma
            end

        # Degree of zero we are finding
        g_degree = ifelse(iszero(α0), degree + 2, degree + 1)

        # The function g dividing away (α - α0)^g_degree
        g_div_α(α::Union{Arb,ArbSeries}, R::Arb) =
            if (α isa Arb && Arblib.overlaps(α, α0)) ||
               (α isa ArbSeries && Arblib.overlaps(α[0], α0))
                fx_div_x(α -> g(α0 + α, R), α - α0, g_degree, force = true)
            else
                #@show g(ArbSeries(α), R)[0] (α - α0)^g_degree
                #@show g(ArbSeries(α), R)[0] / (α - α0)^g_degree
                g(ArbSeries(α), R)[0] / (α - α0)^g_degree
            end

        # Guess of lower and upper bounds for R
        R_low, R_upp = Arb(-10), Arb(10)

        # Check that the signs at the lower and upper bound differ
        sign_low = Arblib.sgn_nonzero(
            ArbExtras.enclosure_series(α -> g_div_α(α, R_low), interval, degree = 8),
        )
        sign_upp = Arblib.sgn_nonzero(
            ArbExtras.enclosure_series(α -> g_div_α(α, R_upp), interval, degree = 8),
        )
        if !(sign_low * sign_upp < 0)
            a = ArbExtras.enclosure_series(α -> g_div_α(α, R_low), interval, degree = 8)
            b = ArbExtras.enclosure_series(α -> g_div_α(α, R_upp), interval, degree = 8)
            throw(ErrorException("Sign of endpoints don't differ: $a $b $interval"))
        end

        # Compute root in R for a fixed α
        g_div_α_root(α::Arb) = Arb(
            ArbExtras.refine_root_bisection(
                R -> g_div_α(α, R),
                lbound(R_low),
                ubound(R_upp),
            ),
        )

        # Compute a guess of R by evaluating on endpoints in α
        R1 = g_div_α_root(lbound(Arb, interval))
        R2 = g_div_α_root(ubound(Arb, interval))

        # To make it easier to prove the zero we take a slightly
        # larger interval
        R_low, R_upp = let R_tmp = union(R1, R2)
            Arblib.mul!(Arblib.radref(R_tmp), Arblib.radref(R_tmp), Mag(1.01))
            getinterval(Arb, R_tmp)
        end

        # Prove that the function has a constant sign on a lower bound
        # of R for all values of α
        check_sign_low = ArbExtras.bounded_by(
            α -> -sign_low * g_div_α(α, R_low),
            getinterval(interval)...,
            Arf(0),
        )
        check_sign_low ||
            throw(ErrorException("could not determine sign on lower bound of R $interval"))
        # Prove that the function has a constant sign (opposite of the
        # above) on an upper bound of R for all values of α
        check_sign_upp = ArbExtras.bounded_by(
            α -> -sign_upp * g_div_α(α, R_upp),
            getinterval(interval)...,
            Arf(0),
        )
        check_sign_upp ||
            throw(ErrorException("could not determine sign on upper bound of R $interval"))

        R = Arb((R_low, R_upp))

        # Expansion with remainder
        p0 = ArbSeries((p00, p01, p02, R))

        # Truncate to given degree
        p0 = truncate_with_remainder(p0, interval - α0; degree)

        # Compare with old way of computing error
        # TODO: This can be removed later
        begin
            if iszero(α0)
                # Check only lower bound of α
                error = let α = lbound(Arb, interval)
                    abs(findp0(α) - p0(α - α0)) / (α - α0)^degree
                end
            else
                error1 = let α = lbound(Arb, interval)
                    abs(findp0(α) - p0(α - α0)) / abs(α - α0)^degree
                end
                error2 = let α = ubound(Arb, interval)
                    abs(findp0(α) - p0(α - α0)) / abs(α - α0)^degree
                end
                error = max(error1, error2)
            end

            @assert radius(p0[end]) < abs_ubound(error)
        end
    end

    return p0
end

"""
    expansion_as(::Type{KdVZeroAnsatz}, α0::Arb, interval::Arb; degree::integer = 2, p0)

Compute expansions for `a[i]` for `i = 1:3` at the point `α0` of the
given degree. The last term is a remainder term which ensures that
evaluating the expansion gives an enclosure of `a[i]` for all `α ∈
interval`.

The expansion of `a[0]` is computed to a degree one higher than the
other. The reason for this is that for `α0 = 0` the constant term in
`a[0]` is zero and in many cases we divide `a[0]` by `α` and still
want to have sufficiently high degree.

# Computing `a[0]`
We have
```
a[0] = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```
For `α0 < 0` we can compute the expansion directly. At `α = 0` both
gamma functions have a pole and direct evaluation fails. Instead we
proceeds as follows.

We can rewrite it in terms of the reciprocal gamma function as
```
a[0] = 2rgamma(α)^2 * cospi(α) / (rgamma(2α) * cospi(α / 2)^2)
```
We can compute `2rgamma(α) * cospi(α) / cospi(α / 2)^2` directly
whereas for `rgamma(α) / rgamma(2α)` we handle the removable
singularity by writing it as `(rgamma(α) / α) / (rgamma(2α) / α)`.

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
1 / 2 * sum(a[j] * zeta(1 - 2α + j * p0 - 2) for j = 0:u0.N0) =
    1 / 2 * (a[0] * zeta(1 - 2α - 2) + a[1] * zeta(1 - 2α + p0 - 2) + a[2] * zeta(1 - 2α + 2p0 - 2))
```
and for `x^(-α + 2)`
```
gamma(α) * sinpi((1 - α) / 2) * a[0] * (-1 / 2) * sum(a[j] * zeta(1 - α + j * p0 - 2) for j = 0:u0.N0) =
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
For `α0 < 0` we can directly evaluate `d`, `v1` and `v2`. For `α0 = 0`
we have to handle several removable singularities. We describe the
procedure for doing so below.

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
function expansion_as(
    ::Type{KdVZeroAnsatz},
    α0::Arb,
    interval::Arb;
    degree::Integer = 2,
    p0 = expansion_p0(KdVZeroAnsatz, α0, interval; degree),
)
    a0 = if iszero(α0)
        # We compute a0 to a higher degree and then truncate, to
        # get a tighter enclosure

        # rgamma(α) / rgamma(2α) = (rgamma(α) / α) / (rgamma(2α) / α)
        g = div_with_remainder(
            taylor_with_remainder(rgamma, α0, interval - α0, degree = degree + 4) << 1,
            taylor_with_remainder(
                α -> rgamma(2α),
                α0,
                interval - α0,
                degree = degree + 4,
            ) << 1,
            interval - α0,
        )

        a0 = mul_with_remainder(
            taylor_with_remainder(
                α -> rgamma(α) * 2cospi(α) / cospi(α / 2)^2,
                α0,
                interval - α0,
                degree = degree + 3,
            ),
            g,
            interval - α0,
        )

        truncate_with_remainder(a0, interval - α0, degree = degree + 1)
    else
        # We set enclosure_degree to a very high number since we don't
        # really care about performance for this method, getting a
        # good enclosure is more important.
        taylor_with_remainder(
            α -> 2gamma(2α) * cospi(α) / (gamma(α) * cospi(α / 2))^2,
            α0,
            interval,
            degree = degree + 1,
            enclosure_degree = 20,
        )
    end

    α_s = ArbSeries((α0, 1); degree) # Series expansion around α0
    z(i, j) = compose_with_remainder(zeta, -1 - i * α_s + j * p0, interval - α0)

    z1 =
        mul_with_remainder(z(2, 1), z(1, 2), interval - α0) -
        mul_with_remainder(z(2, 2), z(1, 1), interval - α0)
    z2 =
        mul_with_remainder(z(2, 2), z(1, 0), interval - α0) -
        mul_with_remainder(z(1, 2), z(2, 0), interval - α0)
    z3 =
        mul_with_remainder(z(2, 1), z(1, 0), interval - α0) -
        mul_with_remainder(z(1, 1), z(2, 0), interval - α0)

    if iszero(α0)
        # The constant coefficients for z1, z2 and z3 are all exactly
        # equal to zero.
        @assert all(Arblib.contains_zero(z[0]) for z in (z1, z2, z3))
        z1[0] = z2[0] = z3[0] = 0

        d = z1
        # Compute v1 and v2 to one degree higher since we divide it by
        # α
        v1 = mul_with_remainder(a0, z2, interval - α0, degree = degree + 1)
        v2 = -mul_with_remainder(a0, z3, interval - α0, degree = degree + 1)

        # Factor out α^2 from v1 and v2 and α from d, multiply back
        # one α afterwards
        a1 = div_with_remainder(v1 << 2, d << 1, interval - α0) >> 1
        a2 = div_with_remainder(v2 << 2, d << 1, interval - α0) >> 1
    else
        d = z1
        v1 = mul_with_remainder(a0, z2, interval - α0)
        v2 = -mul_with_remainder(a0, z3, interval - α0)

        a1 = div_with_remainder(v1, d, interval - α0)
        a2 = div_with_remainder(v2, d, interval - α0)
    end

    return OffsetVector([a0, a1, a2], 0:2)
end
