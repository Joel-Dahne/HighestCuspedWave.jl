"""
    _integrand_compute_root(u0::KdVZeroAnsatz, x::Arb, degree = 1)

Compute the unique root of
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
in `t` on the interval `[0, 1]`. It assumes that `0 <= x <= π`.

If `degree >= 0` it computes an expansion of the root in `α` around `α
= 0`. Currently it only supports `degree <= 1`. If `degree < 0` it
computes an enclosure of the root on `u0.α`, this will only work if
`u0.α` doesn't contain zero and is mostly meant for testing purposes.
In all cases it returns the root in terms of an `ArbSeries`, when
`degree < 0` only the constant term will be set.

# Enclosing the root
The first step is to compute an enclosure of the root for `α = 0`. In
the limit as `α -> 0` the function
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
tends to zero everywhere and to get a good limit we have to normalise
by dividing by `α`. This means that we are searching for the root of
```
(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) / α
```
and in the limit as `α` goes to zero this approaches
```
clausenc(x * (1 - t), 0, 1) + clausenc(x * (1 + t), 0, 1) - 2clausenc(x * t, 0, 1)
```
which is the function we compute the root of.

For wide values of `x` it uses that the root is decreasing in `x` to
only have to evaluate at the endpoints.
- **PROVE:** That the root is decreasing in `x`

If the lower bound of `x` is zero or close to zero (smaller than
`eps(x)`) it computes the root in the limiting case as `x` goes to
zero. Expanding `clausenc(x, 0, 1)` at `x = 0` gives us the expansion
```
π / 2 * abs(x)^-1 + sum(dzeta(-2m) * x^(2m) / factorial(2m) for m = 0:Inf)
```
- **PROVE:** That this is the correct expansion, we just differentiate
  the one for `clausenc(x, s)` with respect to `s`.
The limiting root can then be bound by computing the root of
```
(1 - t)^(-1) + (1 + t)^(-1) - 2t^(-1)
```
- **PROVE:** That the root of the integrand converges to the root of
  `(1 - t)^(-1) + (1 + t)^(-1) - 2t^(-1)`.

If the upper bound of `x` is close to zero, smaller than `eps(x)`, we
compute the root at `eps(x)` and use that as a lower bound. This
avoids computing with very small values of `x`.

In the case `degree < 0` it computes a root of
```
clausenc(x * (1 - t), -u0.α) + clausenc(x * (1 + t), -u0.α) - 2clausenc(x * t, -u0.α)
```
directly.

# Computing expansion in `α` of the root
Once an enclosure of the root is computed the next step is to compute
the expansion in `α`. If we see the root as a function of `α` we want
the expansion of
```
clausenc(x * (1 - root(α)), -α, 1) + clausenc(x * (1 + root(α)), -α, 1) - 2clausenc(x * root(α), -α, 1)
```
at `α = 0` to be identically equal to zero.

Differentiating with respect to `α` gives us
```
-(clausenc(x * (1 - root(α)), -α, 2) + clausenc(x * (1 + root(α)), -α, 2) - 2clausenc(x * root(α), -α, 2)) +
    -x * root'(α) * (-clausens(x * (1 - root(α)), -α - 1, 1) + clausens(x * (1 + root(α)), -α - 1, 1) - 2clausens(x * root(α), -α - 1, 1))
```
Putting this equal to zero, solving for `root'(α)` and setting `α = 0`
gives us
```
root'(0) = -(clausenc(x * (1 - root(0)), -0, 2) + clausenc(x * (1 + root(0)), -0, 2) - 2clausenc(x * root(0), -0, 2)) /
    x * (-clausens(x * (1 - root(0)), -0 - 1, 1) + clausens(x * (1 + root(0)), -0 - 1, 1) - 2clausens(x * root(0), -0 - 1, 1))
```
- **TODO:** Implement asymptotic evaluation of this. We might not need
  it in the end.
- **TODO:** Possibly compute more terms in the expansion. This might
  not be needed in the end.
"""
function _integrand_compute_root(u0::KdVZeroAnsatz, x::Arb; degree = 1)
    compute_root(x) =
        let
            # If degree < 0 compute with an enclosure of α instead of
            # for α = 0.
            f(t) =
                if degree < 0
                    clausenc(x * (1 - t), -u0.α, 1) + clausenc(x * (1 + t), -u0.α, 1) -
                    2clausenc(x * t, -u0.α, 1)
                else
                    clausenc(x * (1 - t), Arb(0), 1) + clausenc(x * (1 + t), Arb(0), 1) -
                    2clausenc(x * t, Arb(0), 1)
                end

            # The root is lower bounded by 1 / 2, take a value
            # slightly larger so that we can still isolate it even if
            # it touches 1 / 2.
            root_lower = Arf(0.5) - sqrt(eps(Arf))

            # Find a crude upper bound for the root
            δ = Arb(0.4)
            # IMPROVE: We can remove this check if we can prove that
            # the root is less than x + δ.
            Arblib.ispositive(f(root_lower + δ)) || return Arblib.indeterminate!(zero(x))
            while Arblib.ispositive(f(root_lower + δ / 2)) && δ > 1e-5
                Arblib.mul_2exp!(δ, δ, -1)
            end
            root_upper = ubound(root_lower + δ)

            # Improve the enclosure of the root
            roots, flags = ArbExtras.isolate_roots(f, root_lower, root_upper)
            if length(flags) == 1 && flags[1]
                # Refine the unique root
                root = ArbExtras.refine_root(f, Arb(only(roots)))
            else
                root = Arb((roots[1][1], roots[end][2]))
            end

            return root
        end

    compute_root_zero() =
        let
            f(t) = (1 - t)^(-1) + (1 + t)^(-1) - 2t^(-1)

            roots, flags = ArbExtras.isolate_roots(f, Arf(0.5), Arf(0.9))
            length(flags) == 1 && flags[1] || error("could not isolate root for x = 0")

            ArbExtras.refine_root(f, Arb(only(roots)))
        end

    # Compute the derivative of the root in α
    compute_derivative(t) =
        let
            num =
                clausenc(x * (1 - t), Arb(0), 2) + clausenc(x * (1 + t), Arb(0), 2) -
                2clausenc(x * t, Arb(0), 2)

            den =
                x * (
                    -clausens(x * (1 - t), Arb(-1), 1) + clausens(x * (1 + t), Arb(-1), 1) -
                    2clausens(x * t, Arb(-1), 1)
                )

            -num / den
        end

    # TODO: Implement this, if we need it in the end
    compute_derivative_zero(t) =
        let
            num =
                clausenc(x * (1 - t), Arb(0), 2) + clausenc(x * (1 + t), Arb(0), 2) -
                2clausenc(x * t, Arb(0), 2)

            den =
                x * (
                    -clausens(x * (1 - t), Arb(-1), 1) + clausens(x * (1 + t), Arb(-1), 1) -
                    2clausens(x * t, Arb(-1), 1)
                )

            -num / den
        end

    xₗ, xᵤ = getinterval(Arb, x)
    xᵤ = min(Arb(π), xᵤ) # We assume that xᵤ <= π
    ϵ = eps(x)

    if iszero(x)
        root = compute_root_zero()
    elseif !iswide(x)
        root = compute_root(x) # In this case x never overlaps zero
    elseif xᵤ < ϵ
        root = Arb((compute_root(ϵ), compute_root_zero()))
    elseif xₗ < eps(Arb)
        root = Arb((compute_root(xᵤ), compute_root_zero()))
    else
        root = Arb((compute_root(xᵤ), compute_root(xₗ)))
    end

    if degree < 0
        res = ArbSeries(root)
    elseif degree == 0
        res = ArbSeries(root)
    elseif degree == 1
        res = ArbSeries((root, compute_derivative(root)))
    else
        degree <= 1 || throw(ArgumentError("degree at most 1 is supported"))
    end

    return res
end

"""
    T0(u0::KdVZeroAnsatz, ::Ball; skip_div_u0 = false)

Compute the integral \$T_0\$.

This method is similar to [`T0_p_one`](@ref) in that it explicitly
computes the integral. It computes an expansion in `α` around `α = 0`
and treats the zero of the integrand differently, otherwise they are
very similar.

# Computing the integral
The integral is given by
```
inv(π * x * u0(x)) * ∫abs(clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * y dy
```
where we have used that the weight `u0.w(x)` is given by `x` and the
integral is taken from `0` to `π`.

The switch of variables to `t = y / x` gives us
```
x / (π * u0(x)) * ∫abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt
```
from `0` to `π / x`. Call the integral, without the factors in front,
`I`.

If we ignore the absolute value the integral can be computed
explicitly using that
```
∫(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt =
    (clausenc(x * (1 - t), 2 - α) / x^2 - t * clausens(x * (1 - t), 1 - α) / x) +
    (clausenc(x * (1 + t), 2 - α) / x^2 + t * clausens(x * (1 + t), 1 - α) / x) -
    2(clausenc(x * t, 2 - α) / x^2 + t * clausens(x * t, 1 - α) / x)
```
Call this function `primitive(t)`. The idea is to isolate the zero
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
We can then integrate explicitly using `primitive(t)` with this zero
as one of the endpoints.

On the interval `[1, π / x]` the expression inside the absolute value
is positive and we can just remove it. This gives us the integral
```
I2 = primitive(π / x) - primitive(1)
```
- **PROVE:** That the integrand is positive.

On the interval `[0, 1]` the expression inside the absolute value has
a unique root, it is negative to the left of the root and positive to
the right. If we let `root` correspond to this root then an enclosure
of the integral on `[0, 1]` is given by
```(
I1 = -(primitive(root) - primitive(0)) + (primitive(1) - primitive(root))
   = primitive(0) - 2primitive(root) + primitive(π / x)
```
Putting this together we get
```
I = I1 + I2 = primitive(0) - 2primitive(root) + primitive(π / x)
```
We also get that
```
primitive(0) = 2clausencmzeta(x, 2 - α) / x^2
```
and
```
primitive(π / x) = 2(clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x^2
```

We can notice that all terms in the result contains a division by `x`
and that we in the end multiply with `x` If we let
```
primitive_mul_x(t) =
            (
                clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) -
                2clausenc(x * t, 2 - α)
            ) / x +
            t * (
                -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) -
                2clausens(x * t, 1 - α)
            )
```
we get
```
x * I = primitive_mul_x(0) - 2primitive_mul_x(root) + primitive_mul_x(π / x)
```

# Expanding in `α`
We are now interested in computing an expansion, in `α`, of this. Most
of this can be done with direct evaluation using `ArbSeries`. However
we need to prove that the constant term in the expansion is exactly
`π` and this requires some additional work. In addition to this we
also need to handle that `root` is given as an expansion in `α`, so we
need to work a bit more to get the expansion of
`primitive_mul_x(root)`.

## Computing the constant term
To get the constant term we let `α = 0` in the above formulas. This
means that we are computing with the parameters `1` and `2` in the
Clausen functions, for which we have the explicit expressions
```
clausenc(x, 2) = π^2 / 6 - π * x / 2 + x^2 / 4
clausens(x, 1) = π / 2 - x / 2
```
when `0 <= x <= 2π`. We can use this to get the value of
`primitive_mul_x(t)` for `α = 0`, for the case when `0 <= t <= 1` (so
that the argument is always on the interval `[0, 2π]`) we get for the
`clausenc` terms
```
(clausenc(x * (1 - t), 2) + clausenc(x * (1 + t), 2) - 2clausenc(x * t, 2)) / x =

(
 + π^2 / 6 - π * x * (1 - t) / 2 + x^2 * (1 - t)^2 / 4
 + π^2 / 6 - π * x * (1 + t) / 2 + x^2 * (1 + t)^2 / 4
 - 2(π^2 / 6 - π * x * t / 2 + x^2 * t^2 / 4)
) / x =

(
 - π * (1 - t) / 2 + x * (1 - t)^2 / 4
 - π * (1 + t) / 2 + x * (1 + t)^2 / 4
 + π * t - x * t^2 / 2
) =

π / 2 * (-(1 - t) - (1 + t) + 2t) +
x / 4 * ((1 - t)^2 + (1 + t^2) - 2t^2) =

π * (t - 1) + x / 2
```
and for the `clausens` terms
```
t * (-clausens(x * (1 - t), 1) + clausens(x * (1 + t), 1) - 2clausens(x * t, 1)) =

t * (
    - (π / 2 - x * (1 - t) / 2)
    + (π / 2 - x * (1 + t) / 2)
    - 2(π / 2 - x * t / 2)
) =

t * (
    + x * (1 - t) / 2
    - x * (1 + t) / 2
    - π + x * t
) =

t * (x * ((1 - t) / 2 - (1 + t) / 2 + t) -2π) =

-π * t
```
Putting these two terms together we get
```
primitive_mul_x(t) -> π * (t - 1) + x / 2 - π * t = x / 2 - π
```
This means that
```
primitive_mul_x(0) - 2primitive_mul_x(root) =
    (x / 2 - π) - 2(x / 2 - π) =
    π - x / 2
```

For the case when `t = π / x` they above calculations do not apply
since the argument of the first `clausenc` and `clausens` functions
are negative and hence not on the interval `[0, 2π]. However we have
already established that we can write the primitive function as
```
primitive_mul_x(π / x) = 2(clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x
```
Letting `α = 0` gives us
```
2(clausenc(x + π, 2) - clausenc(π, 2)) / x =
2(
    + (π^2 / 6 - π * (x + π) / 2 + (x + π)^2 / 4)
    - (π^2 / 6 - π * π / 2 + π^2 / 4)
) / x =
2(
    - π * (x + π) / 2 + (x + π)^2 / 4
    + π^2 / 2 - π^2 / 4
) / x =
2(
    - π * (x + π) / 2 + (x + π)^2 / 4
    + π^2 / 2 - π^2 / 4
) / x =
2(x^2 / 4) / x = x / 2
```
Combining this with the above gives us that the constant function of
the integral is given by
```
primitive(0) - 2primitive(root) + primitive(π / x) =
    (π - x / 2) + (x / 2) = π
```
Which is exactly what we wanted to show. This means that after the
division by `π` the constant function should be exactly `1`.

## Computing expansion of `primitive_mul_x(root)`
The constant term in the expansion can be computed directly. For the
derivative with respect to `α` we by differentiating
`primitive_mul_x(t)` with respect to `α` and treating `t` as a
function of `α`
```
(
    - x * t' * (-clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) - 2clausens(x * t, 1 - α))
    - (-clausenc(x * (1 - t), 2 - α, 1) + clausenc(x * (1 + t), 2 - α, 1) - 2clausenc(x * t, 2 - α, 1))
) / x + t' * (
    -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) - 2clausens(x * t, 1 - α)
) + t * (
    x * t' * (clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α))
    - (-clausens(x * (1 - t), 1 - α, 1) + clausens(x * (1 + t), 1 - α, 1) - 2clausens(x * t, 1 - α, 1))
)
```
From here we can notice several simplifications, to begin with we have two copies of
```
t' * (-clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) - 2clausens(x * t, 1 - α))
```
with opposite sign that cancel out. We also have the factor
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
which is the function that we found the root of, so this will be zero
for that root. What remains is
```
- (
    clausenc(x * (1 - t), 2 - α, 1) + clausenc(x * (1 + t), 2 - α, 1) - 2clausenc(x * t, 2 - α, 1)
) / x -  t * (
    -clausens(x * (1 - t), 1 - α, 1) + clausens(x * (1 + t), 1 - α, 1) - 2clausens(x * t, 1 - α, 1)
)
```
Which is the same derivative we get if we treat `t` as a constant not
depending on `α`. The derivative of the root with respect to `α` hence
**does not** affect the derivative of the result.
- **TODO:** If we only need one derivative then we wont have to
  compute the expansion of the root at all. Saving some effort.
- **TODO:** Improve enclosure for wide values of `x`. This we will
  most likely need to do in the end.
- **TODO:** Handle remainder term in `α`.
"""
function T0(u0::KdVZeroAnsatz, ::Ball; skip_div_u0 = false)
    α = ArbSeries((0, 1), degree = 1)

    return x::Arb -> begin
        # The derivative of the result doesn't depend on the
        # derivative of the root so we only compute it to first order.
        root = _integrand_compute_root(u0, x, degree = 0)

        primitive_mul_x(t::Arb) =
            (
                clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) -
                2clausenc(x * t, 2 - α)
            ) / x +
            t * (
                -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) -
                2clausens(x * t, 1 - α)
            )

        primitive_mul_x(t::ArbSeries) =
            let
                # Compute the derivative considering the root as a constant
                res = primitive_mul_x(t[0])

                # The derivative of res doesn't depend on the derivative
                # of the root. So if the degree is at most 1 we have
                # computed the correct value
                Arblib.degree(res) <= 1 && return res

                # TODO: If we need more derivatives we have to handle them
                # manually.
                error("degree at most 1 supported")
            end

        primitive_mul_x_zero = 2clausencmzeta(x, 2 - α) / x

        # primitive(π / x)
        # If x overlaps with π this gives an indeterminate result
        # which we handle specially
        if Arblib.overlaps(x, Arb(π))
            # Use periodicity of 2π to evaluate at x - π which is
            # close to zero. Use the asymptotic expansion at x = 0 to
            # evaluate it.
            # To compute the expansion around x = 0 it uses the same
            # approach as the asymptotic version of T0 does.
            clausenc_x_plus_pi = let y = abs(x - π), s = 2 - α, M = 2
                gamma_sin = let γ = Arb(Irrational{:γ}()), π = Arb(π)
                    ArbSeries((
                        -π / 2,
                        (γ - 1) * π / 2,
                        (-24π + 24γ * π - 12γ^2 * π - π^3) / 48,
                    ))
                end

                # Singular term
                res = gamma_sin * abspow(y, 1 - α)
                # Analytic terms
                res += sum(
                    (-1)^m * zeta(s - 2m) * abspow(y, 2m) / factorial(2m) for m = 0:M-1
                )
                # Remainder term
                res += abspow(y, 2M) * clausenc_expansion_remainder(y, s, M)

                res
            end
        else
            clausenc_x_plus_pi = clausenc(x + π, 2 - α)
        end
        primitive_mul_x_pi_div_x = 2(clausenc_x_plus_pi - clausenc(Arb(π), 2 - α)) / x

        I_mul_x =
            primitive_mul_x_zero - 2primitive_mul_x(root) + primitive_mul_x_pi_div_x

        I_mul_x_div_pi = I_mul_x / π

        # The constant term should be exactly 1
        @assert Arblib.contains(Arblib.ref(I_mul_x_div_pi, 0), 1)
        I_mul_x_div_pi[0] = 1

        if skip_div_u0
            return I_mul_x_div_pi
        else
            return I_mul_x_div_pi / u0(x)
        end
    end
end

"""
    T0(u0::KdVZeroAnsatz, ::Asymptotic)

Compute the integral \$T_0\$ in a way that works for `x` close to
zero.

The method is similar to the non-asymptotic version but it evaluates
the terms that depend on `x` in an asymptotic way. From the
non-asymptotic version we get that we want to compute
```
(primitive_mul_x(0) - 2primitive_mul_x(root) + primitive_mul_x(π / x)) / (π * u0(x))
```
The leading term in `u0(x)` behaves like `x^-α` (see
[`u0_div_xmα`](@ref)) and we want to explicitly cancel this term. We
therefore rewrite the above as
```
(primitive_mul_x(0) * x^α - 2primitive_mul_x(root) * x^α + primitive_mul_x(π / x) * x^α) / (π * u0(x) * x^α)
```

For
```
primitive_mul_x(0) * x^α = 2clausencmzeta(x, 2 - α) / x^(1 - α)
```
We expand the Clausen function at `x = 0` and compute series in `α` of
the coefficients. The only problematic term is the singular one given by
```
gamma(α - 1) * sinpi(α / 2) * abspow(x, 1 - α)
```
where there is a removable singularity for `gamma(α - 1) * sinpi(α /
2)`. The expansion can be computed to be
```
gamma(α - 1) * sinpi(α / 2) =
    -π / 2 + (γ - 1) * π / 2 * α +
    (-24π + 24γ * π - 12γ^2 * π - π^3) / 48 * α^2
```
- **PROVE:** That this is the expansion, Mathematica gives this.
- **TODO:** Handle remainder term in `α`.

For
```
primitive_mul_x(π / x) * x^α = 2(clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x^(1 - α)
```
We expand `clausenc(x + π, 2 - α)` at `x = π` and cancel `clausenc(π,
2 - α)` explicitly. Since `clausenc` is analytic in `x` at `x = π` we
can compute the expansion by differentiating by hand and evaluating
with a series in `α`.

Finally for `primitive_mul_x(t)` we can expand the `clausenc` terms in
the same way as for `primitive_mul_x(0)`. For the `clausens` terms the
singular term in the expansion in `x` is given by
```
gamma(α) * sinpi(α / 2) * abspow(x, -α)
```
which also has a removable singularity. The expansion is given by
```
gamma(α) * sinpi(α / 2) =
    π / 2 - γ * π / 2 * α + (12γ^2 * π + π^3) / 48 * α^2
```
- **PROVE:** That this is the expansion, Mathematica gives this.
- **TODO:** Handle remainder term in `α`.

**TODO:** Handle remainder term in `α**.
"""
function T0(
    u0::KdVZeroAnsatz,
    ::Asymptotic;
    ϵ::Arb = one(Arb),
    M::Integer = 5,
    skip_div_u0 = false,
)
    α = ArbSeries((0, 1), degree = 1)

    u0_expansion = u0(ϵ, AsymptoticExpansion())

    return x::Arb -> begin
        x <= ϵ || throw(ArgumentError("x needs to be smaller than ϵ, got x = $x, ϵ = $ϵ"))

        # We only need an enclosure of the root, not the derivative
        root = _integrand_compute_root(u0, x, degree = 0)

        # FIXME: Figure out how to handle remainder terms in α
        # Compute primitive_mul_x(t) * x^α = primitive(t) * x^(1 + α)
        primitive_mul_x_onepα(t::Arb) =
            let γ = Arb(Irrational{:γ}()), π = Arb(π)
                part1 = let s = 2 - α
                    gamma_sin = ArbSeries((
                        -π / 2,
                        (γ - 1) * π / 2,
                        (-24π + 24γ * π - 12γ^2 * π - π^3) / 48,
                    ))

                    # Singular term
                    res =
                        gamma_sin * (
                            abspow(1 - t, 1 - α) + abspow(1 + t, 1 - α) - 2abspow(t, 1 - α)
                        )
                    # Analytic terms
                    res += sum(
                        (-1)^m *
                        zeta(s - 2m) *
                        abspow(x, 2m - 1 + α) *
                        ((1 - t)^2m + (1 + t)^2m - 2t^2m) / factorial(2m) for
                        m = 1:M-1
                    )
                    # Remainder term
                    # Here we use that max(1 - t, 1 + t, t) = 1 + t so
                    # it is enough to compute the remainder term at x * (1 + t)
                    res +=
                        abspow(x, 2M - 1 + α) *
                        ((1 - t)^2M + (1 + t)^2M - 2t^2M) *
                        clausenc_expansion_remainder(x * (1 + t), s, M)

                    res
                end

                part2 = let s = 1 - α
                    gamma_sin = ArbSeries((π / 2, -γ * π / 2, (12γ^2 * π + π^3) / 48))

                    # Singular term
                    res =
                        gamma_sin *
                        (-abspow(1 - t, -α) + abspow(1 + t, -α) - 2abspow(t, -α))
                    # Analytic terms
                    res += sum(
                        (-1)^m *
                        zeta(s - 2m - 1) *
                        abspow(x, 2m + 1 + α) *
                        (-(1 - t)^(2m + 1) + (1 + t)^(2m + 1) - 2t^(2m + 1)) /
                        factorial(2m + 1) for m = 0:M-1
                    )
                    # Remainder term
                    res +=
                        abspow(x, 2M + 1 + α) *
                        ((1 - t)^(2M + 1) + (1 + t)^(2M + 1) - 2t^(2M + 1)) *
                        clausens_expansion_remainder(x * (1 + t), s, M)

                    t * res
                end

                part1 + part2
            end

        primitive_mul_x_onepα(t::ArbSeries) =
            let
                # Compute the derivative considering the root as a constant
                res = primitive_mul_x_onepα(t[0])

                # The derivative of res doesn't depend on the derivative
                # of the root. So if the degree is at most 1 we have
                # computed the correct value
                Arblib.degree(res) <= 1 && return res

                # TODO: If we need more derivatives we have to handle them
                # manually.
                error("degree at most 1 supported")
            end

        # Compute primitive_mul_x_onepα(0) = 2clausencmzeta(x, 2 - α) / x^(1 - α)
        # FIXME: Figure out how to handle remainder terms in α
        primitive_mul_x_onepα_zero = let γ = Arb(Irrational{:γ}()), π = Arb(π)
            s = 2 - α
            gamma_sin = ArbSeries((
                -π / 2,
                (γ - 1) * π / 2,
                (-24π + 24γ * π - 12γ^2 * π - π^3) / 48,
            ))

            # Singular term
            res = gamma_sin
            # Analytic terms
            res += sum(
                (-1)^m * zeta(s - 2m) * abspow(x, 2m - 1 + α) / factorial(2m) for m = 1:M-1
            )
            # Remainder term
            res += abspow(x, 2M - 1 + α) * clausenc_expansion_remainder(x, s, M)

            2res
        end

        # Compute primitive_mul_x_onepα(π / x) =
        # 2(clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x^(1 - α)
        # FIXME: Figure out how to handle remainder terms in α
        primitive_mul_x_onepα_pi_div_x = let π = Arb(π)
            s = 2 - α
            res = zero(primitive_mul_x_onepα_zero)

            # We start from n = 1 since the first is cancelled
            for m = 1:M-1
                # m-th derivative at x = π
                deriv = (-1)^(m ÷ 2) * (iseven(m) ? clausenc(π, s - m) : clausens(π, s - m))
                # We have x^(m - 1) since we divide by x
                res += deriv / factorial(m) * abspow(x, m - 1 + α)
            end

            # Remainder term
            # Interval for the Taylor expansion
            interval = union(π, π + x)
            # Enclosure of M-th derivative on interval
            deriv =
                (-1)^(M ÷ 2) *
                (iseven(M) ? clausenc(interval, s - M) : clausens(interval, s - M))
            # Add enclosure of the remainder term
            res += deriv / factorial(M) * abspow(x, M - 1 + α)

            2res
        end

        I_mul_x =
            primitive_mul_x_onepα_zero - 2primitive_mul_x_onepα(root) +
            primitive_mul_x_onepα_pi_div_x

        I_mul_x_onepα_div_pi = I_mul_x / π

        if skip_div_u0
            res = I_mul_x_onepα_div_pi
        else
            res = I_mul_x_onepα_div_pi / eval_expansion(u0, u0_expansion, x, offset_i = -1)
        end

        # The constant term should be exactly 1. This holds not matter
        # if we divide by u0 or not since the constant term of u0 is
        # exactly 1.
        @assert Arblib.contains(Arblib.ref(I_mul_x_onepα_div_pi, 0), 1)
        res[0] = 1

        return res
    end
end
