"""
    _integrand_compute_root(::Type{KdVZeroAnsatz}, x::Arb, αₗ::Arb)

Compute enclosures of the unique root of
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
in `t` on the interval `[0, 1]`, one enclosure for `α = 0` and one
valid for `α ∈ [αₗ, 0]`. It assumes that `0 <= x <= π`.

To compute the root for `α ∈ [αₗ, 0]` it uses that the root is
decreasing in `α`. To get an enclosure it is therefore enough to
compute the root at `α = 0` and `α = αₗ`. To compute the root at `αₗ`
it uses `_integrand_compute_root(FractionalKdVAnsatz, x, αₗ)`.

# Computing the root for `α = 0`
In the limit as `α -> 0` the function
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

If the lower bound of `x` is zero or close to zero (smaller than
`eps(Arb)`) it computes the root in the limiting case as `x` goes to
zero. In that case we can use the formulation
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
and expand at `x = 0` to get that the limiting root is the root of
```
(1 - t)^(-1) + (1 + t)^(-1) - 2t^(-1)
```

If the upper bound of `x` is close to zero, smaller than `eps(Arb)`, we
compute the root at `eps(Arb)` and use that as a lower bound. This
avoids computing with very small values of `x`.
"""
function _integrand_compute_root(::Type{KdVZeroAnsatz}, x::Arb, αₗ::Arb)
    # Compute root for α = 0 and x = 0
    root0_zero = let
        f = t -> (1 - t)^(-1) + (1 + t)^(-1) - 2t^(-1)

        roots, flags = ArbExtras.isolate_roots(f, Arf(0.5), Arf(0.9))
        length(flags) == 1 && flags[1] || error("could not isolate root for x = 0")

        ArbExtras.refine_root(f, Arb(only(roots)))
    end

    # Compute root for α = 0 at a given x
    compute_root0(x::Arb) =
        let
            # If degree < 0 compute with an enclosure of α instead of
            # for α = 0.
            f =
                t ->
                    clausenc(x * (1 - t), Arb(0), 1) + clausenc(x * (1 + t), Arb(0), 1) -
                    2clausenc(x * t, Arb(0), 1)

            # The root is lower bounded by 1 / 2
            root_lower = Arf(0.5)

            # Find a crude upper bound for the root
            # root_lower + δ gives upper bound of root
            δ = root0_zero - root_lower
            while Arblib.ispositive(f(root_lower + δ / 2))
                Arblib.mul_2exp!(δ, δ, -1)
            end
            root_upper = ubound(root_lower + δ)

            # Short circuit in case the sign can't be determined on
            # the lower endpoint, this happens when x is very close to
            # π
            Arblib.contains_zero(f(Arb(root_lower))) && return Arb((root_lower, root_upper))

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

    # Compute an enclosure of the root for α = 0

    xₗ, xᵤ = getinterval(Arb, x)
    xᵤ = min(Arb(π), xᵤ) # We assume that xᵤ <= π
    ϵ = eps(Arb)

    if iszero(x)
        root0 = root0_zero
    elseif Arblib.overlaps(xᵤ, Arb(π))
        root0 = Arb((1 // 2, compute_root0(xₗ))) # Lower bound is 1 / 2
    elseif !iswide(x)
        root0 = compute_root0(x) # In this case x never overlaps zero
    elseif xᵤ < ϵ
        root0 = Arb((compute_root0(ϵ), root0_zero))
    elseif xₗ < eps(Arb)
        root0 = Arb((compute_root0(xᵤ), root0_zero))
    else
        root0 = Arb((compute_root0(xᵤ), compute_root0(xₗ)))
    end

    if iszero(αₗ)
        rootₗ = root0
    elseif -1 < αₗ < 0
        # Compute root at α = αₗ
        rootₗ = _integrand_compute_root(FractionalKdVAnsatz, x, αₗ)
    else
        throw(ArgumentError("requires that -1 < αₗ < 0"))
    end

    return root0, Arb((root0, rootₗ))
end

"""
    T0(u0::KdVZeroAnsatz, ::Ball; skip_div_u0 = false)

Return a Taylor model of the integral ``T_0`` in `α` around `α = 0` of
degree `1`.

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
Call this function `primitive(t)`. The idea is to isolate the zero of
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
We can then integrate explicitly using `primitive(t)` with this zero
as one of the endpoints.

On the interval `[1, π / x]` the expression inside the absolute value
is positive, due to lemma [`lemma_integrand_2`](@ref) and we can just
remove it. This gives us the integral ```
I2 = primitive(π / x) - primitive(1)
```

On the interval `[0, 1]` the expression inside the absolute value has
a unique root, it is negative to the left of the root and positive to
the right. If we let `root` correspond to this root then the integral
on `[0, 1]` is given by
```(
I1 = -(primitive(root) - primitive(0)) + (primitive(1) - primitive(root))
   = primitive(0) - 2primitive(root) + primitive(1)
```
Putting this together we get
```
I = I1 + I2 = primitive(0) - 2primitive(root) + primitive(π / x)
```
Furthermore we get
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
of this can be done with direct evaluation using Taylor models.
However we need to prove that the constant term in the expansion is
exactly `π` and this requires some additional work. In addition to
this we also need to handle that `root` is given as an expansion in
`α`, so we need to work a bit more to get the expansion of
`primitive_mul_x(root)`.

## Computing the constant term
To get the constant term we let `α = 0` in the above formulas. This
means that we are computing with the parameters `s = 1` and `s = 2` in
the Clausen functions, for which we have the explicit expressions
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

For the case when `t = π / x` the above calculations do not apply
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

## Computing Taylor model of `primitive_mul_x(root)`
We want to compute a Taylor model of degree `0`, we do it by manually
computing the expansion and the remainder. The constant term in the
expansion can be computed directly, to get a remainder term we bound
the derivative in `α`. Differentiating `primitive_mul_x(t)` with
respect to `α` while treating `t` as a function of `α` gives us
```
(
    - x * dt(α) * (-clausens(x * (1 - t(α)), 1 - α) + clausens(x * (1 + t(α)), 1 - α) - 2clausens(x * t(α), 1 - α))
    - (-clausenc(x * (1 - t(α)), 2 - α, 1) + clausenc(x * (1 + t(α)), 2 - α, 1) - 2clausenc(x * t(α), 2 - α, 1))
) / x + dt(α) * (
    -clausens(x * (1 - t(α)), 1 - α) + clausens(x * (1 + t(α)), 1 - α) - 2clausens(x * t(α), 1 - α)
) + t(α) * (
    x * dt(α) * (clausenc(x * (1 - t(α)), -α) + clausenc(x * (1 + t(α)), -α) - 2clausenc(x * t(α), -α))
    - (-clausens(x * (1 - t(α)), 1 - α, 1) + clausens(x * (1 + t(α)), 1 - α, 1) - 2clausens(x * t(α), 1 - α, 1))
)
```
Where `dt(α)` is the derivative of `t(α)` with respect to `α`. From
here we can notice several simplifications, to begin with we have two
copies of
```
dt(α) * (-clausens(x * (1 - t(α)), 1 - α) + clausens(x * (1 + t(α)), 1 - α) - 2clausens(x * t(α), 1 - α))
```
with opposite sign that cancel out. We also have the factor
```
clausenc(x * (1 - t(α)), -α) + clausenc(x * (1 + t(α)), -α) - 2clausenc(x * t(α), -α)
```
which is the function that we found the root of, so this will be zero
for that root. What remains is
```
- (
    clausenc(x * (1 - t(α)), 2 - α, 1) + clausenc(x * (1 + t(α)), 2 - α, 1) - 2clausenc(x * t(α), 2 - α, 1)
) / x -  t(α) * (
    -clausens(x * (1 - t(α)), 1 - α, 1) + clausens(x * (1 + t(α)), 1 - α, 1) - 2clausens(x * t(α), 1 - α, 1)
)
```
We get an enclosure of the derivative by enclosing this in `α`. An
important observation is that the derivative of the root with respect
to `α` **does not** affect the derivative of the result. It is
therefore enough to compute only an enclosure of the root and we do
not need to compute its derivative in `α`.
"""
function T0(u0::KdVZeroAnsatz, ::Ball; skip_div_u0 = false)
    iszero(u0.α0) || throw(ArgumentError("only works for u0.α0 = 0, got u0.α0 = $(u0.α0)"))

    Mα = TaylorModel(identity, u0.α, u0.α0, degree = 0)

    return x::Arb -> begin
        # primitive_mul_x(0)
        primitive_mul_x_zero = 2clausencmzeta(x, 2 - Mα) / x

        # primitive_mul_x(root)
        primitive_mul_x_root = let
            r0, r = _integrand_compute_root(typeof(u0), x, lbound(Arb, u0.α))

            # Constant term in expansion
            c =
                (
                    clausenc(x * (1 - r0), 2) + clausenc(x * (1 + r0), 2) -
                    2clausenc(x * r0, 2)
                ) / x +
                r0 * (
                    -clausens(x * (1 - r0), 1) + clausens(x * (1 + r0), 1) -
                    2clausens(x * r0, 1)
                )

            # Remainder term, enclosure of derivative
            Δ =
                -(
                    clausenc(x * (1 - r), 2 - u0.α, 1) +
                    clausenc(x * (1 + r), 2 - u0.α, 1) -
                    2clausenc(x * r, 2 - u0.α, 1)
                ) / x -
                r * (
                    -clausens(x * (1 - r), 1 - u0.α, 1) +
                    clausens(x * (1 + r), 1 - u0.α, 1) -
                    2clausens(x * r, 1 - u0.α, 1)
                )

            TaylorModel(ArbSeries((c, Δ), degree = 1), u0.α, u0.α0)
        end

        # primitive_mul_x(π / x)
        # If x overlaps with π this gives an indeterminate result
        # which we handle specially
        if Arblib.overlaps(x, Arb(π))
            # Use periodicity of 2π to evaluate at x - π which is
            # close to zero. Use the asymptotic expansion at x = 0 to
            # evaluate it.
            # To compute the expansion around x = 0 it uses the same
            # approach as the asymptotic version of T0 does.
            clausenc_x_plus_pi = let y = abs(x - π), Ms = 2 - Mα, M = 2
                # Expansion of gamma(α - 1) * sinpi(α / 2) = (sinpi(α
                # / 2) / α) / (rgamma(α - 1) / α).
                gamma_sin = truncate(
                    div_removable(
                        TaylorModel(α -> sinpi(α / 2), u0.α, u0.α0, degree = 3),
                        TaylorModel(α -> rgamma(α - 1), u0.α, u0.α0, degree = 3),
                    ),
                    degree = 0,
                )

                # Singular term
                res = gamma_sin * compose(e -> abspow(y, e), 1 - Mα)

                # Analytic terms
                res += sum(
                    (-1)^m * compose(zeta, Ms - 2m) * abspow(y, 2m) / factorial(2m) for m = 0:M-1
                )
                # Remainder term
                res +=
                    abspow(y, 2M) * compose(s -> clausenc_expansion_remainder(y, s, M), Ms)

                res
            end
        else
            clausenc_x_plus_pi = clausenc(x + π, 2 - Mα)
        end

        primitive_mul_x_pi_div_x = 2(clausenc_x_plus_pi - clausenc(Arb(π), 2 - Mα)) / x

        I_mul_x =
            primitive_mul_x_zero - 2primitive_mul_x_root + primitive_mul_x_pi_div_x

        I_mul_x_div_pi = I_mul_x / Arb(π)

        # The constant term should be exactly 1
        @assert Arblib.contains(I_mul_x_div_pi.p[0], 1) || !isfinite(I_mul_x_div_pi.p[0])
        I_mul_x_div_pi.p[0] = 1

        if skip_div_u0
            return I_mul_x_div_pi
        else
            return I_mul_x_div_pi / truncate(u0(x), degree = 0)
        end
    end
end

"""
    T0(u0::KdVZeroAnsatz, ::Asymptotic)

Return a Taylor model of the integral ``T_0`` in `α` around `α = 0` of
degree `1`. Computed in a way that works for `x` close to zero.

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
2)`. We can handle this by rewriting it as `(sinpi(α / 2) / α) /
(rgamma(α - 1) / α)` and explicitly dealing with it.

For
```
primitive_mul_x(π / x) * x^α = 2(clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x^(1 - α)
```
We expand `clausenc(x + π, 2 - α)` at `x = π` and cancel `clausenc(π,
2 - α)` explicitly. Since `clausenc` is analytic in `x` at `x = π` we
can compute the expansion by differentiating by hand and evaluating
with a series in `α`.

# Computing Taylor model of `primitive_mul_x(root) * x^α`
Similarly to the non-asymptotic version we manually compute the
expansion and the remainder. We have
```
primitive_mul_x(t) * x^α =
            (
                clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) -
                2clausenc(x * t, 2 - α)
            ) / x^(1 - α) +
            t * (
                -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) -
                2clausens(x * t, 1 - α)
            ) * x^α
```
The derivative with respect to `α` at `t = root` is
```
primitive_mul_x'(t) * x^α + primitive_mul_x(t) * log(x) * x^α
```
With `primitive_mul_x'(t)` given by
```
- (
    clausenc(x * (1 - t), 2 - α, 1) + clausenc(x * (1 + t), 2 - α, 1) - 2clausenc(x * t, 2 - α, 1)
) / x -  t(α) * (
    -clausens(x * (1 - t), 1 - α, 1) + clausens(x * (1 + t), 1 - α, 1) - 2clausens(x * t, 1 - α, 1)
)
```
We can split the derivative into four parts
```
primitive_mul_x'(t) * x^α + primitive_mul_x(t) * log(x) * x^α =
    -(clausenc(x * (1 - t), 2 - α, 1) + clausenc(x * (1 + t), 2 - α, 1) - 2clausenc(x * t, 2 - α, 1)) * x^(α - 1) +

    -t * (-clausens(x * (1 - t), 1 - α, 1) + clausens(x * (1 + t), 1 - α, 1) - 2clausens(x * t, 1 - α, 1)) * x^α +

    (clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) - 2clausenc(x * t, 2 - α)) * log(x) * x^(α - 1) +

    t * (-clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) - 2clausens(x * t, 1 - α)) * log(x) * x^α
```
Expanding the Clausen functions we get for each part
```
(clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) - 2clausenc(x * t, 2 - α)) * log(x) * x^(α - 1) =
    gamma(α - 1) * sinpi((2 - α) / 2) * log(x) * ((1 - t)^(1 - α) + (1 + t)^(1 - α) - 2t^(1 - α)) +
    sum((-1)^m * zeta(2 - α - 2m) * log(x) * x^(2m - 1 + α) * ((1 - t)^2m + (1 + t)^2m - 2t^2m) / factorial(2m) for m = 1:Inf)

t * (-clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) - 2clausens(x * t, 1 - α)) * log(x) * x^α =
    t * (
        gamma(α) * cospi((1 - α) / 2) * log(x) * ((1 - t)^(-α) + (1 + t)^(-α) - 2t^(-α)) +
        sum((-1)^m * zeta(-α - 2m) * log(x) * x^(2m + 1 + α) * ((1 - t)^(2m + 1) + (1 + t)^(2m + 1) - 2t^(2m + 1)) / factorial(2m + 1) for m = 0:Inf)
    )

-(clausenc(x * (1 - t), 2 - α, 1) + clausenc(x * (1 + t), 2 - α, 1) - 2clausenc(x * t, 2 - α, 1)) * x^(α - 1) =
    -ds₁(gamma(α - 1) * sinpi((2 - α) / 2) * x^-α * ((1 - t)^(1 - α) + (1 + t)^(1 - α) - 2t^(1 - α))) * x^α +
    -sum((-1)^m * dzeta(2 - α - 2m) * x^(2m - 1 + α) * ((1 - t)^2m + (1 + t)^2m - 2t^2m) / factorial(2m) for m = 1:Inf)

-t * (-clausens(x * (1 - t), 1 - α, 1) + clausens(x * (1 + t), 1 - α, 1) - 2clausens(x * t, 1 - α, 1)) * x^α =
    -t * (
        ds₂(gamma(α) * cospi((1 - α) / 2) * x^-α * ((1 - t)^(-α) + (1 + t)^(-α) - 2t^(-α))) * x^α +
        sum((-1)^m * zeta(-α - 2m) * x^(2m + 1 + α) * ((1 - t)^(2m + 1) + (1 + t)^(2m + 1) - 2t^(2m + 1)) / factorial(2m + 1) for m = 0:Inf)
    )
```
Here we use `ds₁` and `ds₂` to represent differentiation with respect
to `2 - α` and `1 - α` respectively. To evaluate this close to `x = 0`
we need to handle the cancellations between the first and third term
as well as between the second and fourth term. We can write the
singular term in the third term as
```
ds₁(gamma(α - 1) * sinpi((2 - α) / 2) * ((1 - t)^(1 - α) + (1 + t)^(1 - α) - 2t^(1 - α))) +
    gamma(α - 1) * sinpi((2 - α) / 2) * log(x) * ((1 - t)^(1 - α) + (1 + t)^(1 - α) - 2t^(1 - α))
```
and note that the term containing `log(x)` exactly cancels the
corresponding term in the expansion of
```
(clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) - 2clausenc(x * t, 2 - α)) * log(x) * x^(α - 1)
```
Similarly we get for the singular part of the fourth term
```
ds₂(gamma(α) * cospi((1 - α) / 2) * ((1 - t)^(-α) + (1 + t)^(-α) - 2t^(-α))) +
    gamma(α) * cospi((1 - α) / 2) * log(x) * ((1 - t)^(-α) + (1 + t)^(-α) - 2t^(-α))
```
which cancels the singular term in the expansion of
```
(clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) - 2clausens(x * t, 1 - α)) * log(x) * x^α
```

We hence get the expansion
```
primitive_mul_x'(t) * x^α + primitive_mul_x(t) * log(x) * x^α =

-(
    ds₁(gamma(α - 1) * sinpi((2 - α) / 2) * ((1 - t)^(1 - α) + (1 + t)^(1 - α) - 2t^(1 - α))) -
    sum(
        (-1)^m * zeta(2 - α - 2m) * log(x) * x^(2m - 1 + α) *
            ((1 - t)^2m + (1 + t)^2m - 2t^2m) / factorial(2m)
        for m = 1:Inf
    ) +
    sum(
        (-1)^m * dzeta(2 - α - 2m) * x^(2m - 1 + α) *
            ((1 - t)^2m + (1 + t)^2m - 2t^2m) / factorial(2m)
        for m = 1:Inf
) -

t * (
    ds₂(gamma(α) * cospi((1 - α) / 2) * (-(1 - t)^(-α) + (1 + t)^(-α) - 2t^(-α))) -
    sum(
        (-1)^m * zeta(-α - 2m) * log(x) * x^(2m + 1 + α) *
            (-(1 - t)^(2m + 1) + (1 + t)^(2m + 1) - 2t^(2m + 1)) / factorial(2m + 1)
        for m = 0:Inf
    ) +
    sum(
        (-1)^m * zeta(-α - 2m) * x^(2m + 1 + α) *
            (-(1 - t)^(2m + 1) + (1 + t)^(2m + 1) - 2t^(2m + 1)) / factorial(2m + 1)
         for m = 0:Inf
    )
)
```
"""
function T0(
    u0::KdVZeroAnsatz,
    ::Asymptotic;
    ϵ::Arb = one(Arb),
    M::Integer = 10,
    skip_div_u0 = false,
)
    iszero(u0.α0) || throw(ArgumentError("only works for u0.α0 = 0, got u0.α0 = $(u0.α0)"))

    Mα = TaylorModel(identity, u0.α, u0.α0, degree = 0)

    u0_expansion = u0(ϵ, AsymptoticExpansion())

    # Expansion of gamma(α - 1) * sinpi(α / 2) = (sinpi(α / 2) / α) /
    # (rgamma(α - 1) / α). We compute it to higher degree and then
    # truncate to get a better enclosure.
    gamma_sin = truncate(
        div_removable(
            TaylorModel(α -> sinpi(α / 2), u0.α, u0.α0, degree = 3),
            TaylorModel(α -> rgamma(α - 1), u0.α, u0.α0, degree = 3),
        ),
        degree = 0,
    )

    return x::Arb -> begin
        x <= ϵ || throw(ArgumentError("x needs to be smaller than ϵ, got x = $x, ϵ = $ϵ"))

        # Compute primitive_mul_x_onepα(0) = 2clausencmzeta(x, 2 - α) / x^(1 - α)
        primitive_mul_x_onepα_zero = let Ms = 2 - Mα
            # Singular term
            res = gamma_sin

            # Analytic terms
            res += sum(
                (-1)^m * compose(zeta, Ms - 2m) * abspow(x, 2m - 1 + Mα) / factorial(2m) for m = 1:M-1
            )

            # Remainder term
            res +=
                abspow(x, 2M - 1 + Mα) *
                compose(s -> clausenc_expansion_remainder(x, s, M), Ms)

            2res
        end

        # primitive_mul_x(root) * x^α
        primitive_mul_x_onepα_root = let
            r0, r = _integrand_compute_root(typeof(u0), x, lbound(Arb, u0.α))

            # Constant term at α = 0
            c = let t = r0
                # Enclosure of
                # (
                #   clausenc(x * (1 - t), 2) +
                #   clausenc(x * (1 + t), 2) -
                #   2clausenc(x * t, 2)
                # ) / x
                part1 = let s = Arb(2)
                    # Singular term
                    res =
                        gamma_sin.p[0] * (
                            abspow(1 - t, s - 1) + abspow(1 + t, s - 1) - 2abspow(t, s - 1)
                        )
                    # Analytic terms
                    res += sum(
                        (-1)^m *
                        zeta(s - 2m) *
                        abspow(x, 2m - 1) *
                        ((1 - t)^2m + (1 + t)^2m - 2t^2m) /
                        factorial(2m) for m = 1:M-1
                    )
                    # Remainder term
                    # Here we use that max(1 - r0, 1 + r0, r0) = 1 + r0 so
                    # it is enough to compute the remainder term at x * (1 + r0)
                    res +=
                        ((1 - t)^2M + (1 + t)^2M - 2t^2M) *
                        abspow(x, 2M - 1) *
                        clausenc_expansion_remainder(x * (1 + t), s, M)

                    res
                end

                # Enclosure of
                # t * (
                #   -clausens(x * (1 - t), 1) +
                #   clausens(x * (1 + t), 1) -
                #   2clausenc(x * t, 1)
                # )
                part2 = let s = Arb(1)
                    # Singular term
                    res =
                        -gamma_sin.p[0] * (
                            -abspow(1 - t, s - 1) + abspow(1 + t, s - 1) -
                            2abspow(t, s - 1)
                        )
                    # Analytic terms
                    res += sum(
                        (-1)^m *
                        zeta(s - 2m - 1) *
                        abspow(x, 2m + 1) *
                        (-(1 - t)^(2m + 1) + (1 + t)^(2m + 1) - 2t^(2m + 1)) / factorial(2m + 1) for m = 0:M-1
                    )
                    # Remainder term
                    res +=
                        ((1 - t)^(2M + 1) + (1 + t)^(2M + 1) - 2t^(2M + 1)) *
                        abspow(x, 2M + 1) *
                        clausens_expansion_remainder(x * (1 + t), s, M)

                    t * res
                end

                part1 + part2
            end

            # Enclosure of derivative over interval
            Δ = let t = r
                # Enclosure of
                # (
                #   clausenc(x * (1 - t), 2 - α) +
                #   clausenc(x * (1 + t), 2 - α) -
                #   2clausenc(x * t, 2 - α)
                # ) * log(x) * x^(α - 1) -
                # (
                #    clausenc(x * (1 - t), 2 - α, 1) +
                #    clausenc(x * (1 + t), 2 - α, 1) -
                #    2clausenc(x * t, 2 - α, 1)
                #  ) * x^(α - 1)
                part1 = let s = 2 - u0.α
                    # Singular term
                    # Derivative of
                    # gamma(1 - s) * sinpi(s / 2) * (
                    #   abspow(1 - t, s - 1) + abspow(1 + t, s - 1) - 2abspow(t, s - 1)
                    # )
                    # with respect to s
                    res = let s = ArbSeries((s, 1))
                        gamma_sin_tmp =
                            fx_div_x(y -> sinpi((y + 2) / 2), s - 2) /
                            fx_div_x(y -> rgamma(1 - (y + 2)), s - 2)

                        expansion =
                            gamma_sin_tmp * (
                                abspow(1 - t, s - 1) + abspow(1 + t, s - 1) -
                                2abspow(t, s - 1)
                            )

                        expansion[1]
                    end

                    # Analytic terms
                    # First sum
                    res -= sum(1:M-1) do m
                        (-1)^m *
                        zeta(s - 2m) *
                        logabspow(x, 1, 2m - 1 + u0.α) *
                        ((1 - t)^2m + (1 + t)^2m - 2t^2m) /
                        factorial(2m)
                    end
                    # Second sum
                    res += sum(1:M-1) do m
                        (-1)^m *
                        dzeta(s - 2m) *
                        abspow(x, 2m - 1 + u0.α) *
                        ((1 - t)^2m + (1 + t)^2m - 2t^2m) /
                        factorial(2m)
                    end

                    # Remainder term
                    # Here we use that max(1 - t, 1 + t, t) = 1 + t so
                    # it is enough to compute the remainder term at x * (1 + t)
                    res -=
                        ((1 - t)^2M + (1 + t)^2M - 2t^2M) *
                        logabspow(x, 1, 2M - 1 + u0.α) *
                        clausenc_expansion_remainder(x * (1 + t), s, M)
                    res +=
                        ((1 - t)^2M + (1 + t)^2M - 2t^2M) *
                        abspow(x, 2M - 1 + u0.α) *
                        clausenc_expansion_remainder(x * (1 + t), s, 1, M)

                    res
                end

                # Enclosure of
                # t * (
                #   (
                #     -clausens(x * (1 - t), 1 - α) +
                #     clausens(x * (1 + t), 1 - α) -
                #     2clausens(x * t, 1 - α)
                #   ) * log(x) * x^α -
                #   (
                #     -clausens(x * (1 - t), 1 - α, 1) +
                #     clausens(x * (1 + t), 1 - α, 1) -
                #     2clausens(x * t, 1 - α, 1)
                #   ) * x^α
                # )
                part2 = let s = 1 - u0.α
                    # Singular term
                    # Derivative of
                    # gamma(1 - s) * cospi(s / 2) * (
                    #   -abspow(1 - t, s - 1) + abspow(1 + t, s - 1) - 2abspow(t, s - 1)
                    # )
                    # with respect to s
                    res = let s = ArbSeries((s, 1))
                        gamma_sin_tmp =
                            fx_div_x(y -> cospi((y + 1) / 2), s - 1) /
                            fx_div_x(y -> rgamma(1 - (y + 1)), s - 1)

                        expansion =
                            gamma_sin_tmp * (
                                -abspow(1 - t, s - 1) + abspow(1 + t, s - 1) -
                                2abspow(t, s - 1)
                            )

                        expansion[1]
                    end

                    # Analytic terms
                    # First sum
                    res -= sum(
                        (-1)^m *
                        zeta(s - 2m - 1) *
                        logabspow(x, 1, 2m + 1 + u0.α) *
                        (-(1 - t)^(2m + 1) + (1 + t)^(2m + 1) - 2t^(2m + 1)) / factorial(2m + 1) for m = 0:M-1
                    )
                    # Second sum
                    res += sum(
                        (-1)^m *
                        dzeta(s - 2m - 1) *
                        abspow(x, 2m + 1 + u0.α) *
                        (-(1 - t)^(2m + 1) + (1 + t)^(2m + 1) - 2t^(2m + 1)) / factorial(2m + 1) for m = 0:M-1
                    )

                    # Remainder term
                    # Here we use that max(1 - t, 1 + t, t) = 1 + t so
                    # it is enough to compute the remainder term at x * (1 + t)
                    res -=
                        (-(1 - t)^(2M + 1) + (1 + t)^(2M + 1) - 2t^(2M + 1)) *
                        logabspow(x, 1, 2M + 1 + u0.α) *
                        clausens_expansion_remainder(x * (1 + t), s, M)

                    res +=
                        (-(1 - t)^(2M + 1) + (1 + t)^(2M + 1) - 2t^(2M + 1)) *
                        abspow(x, 2M + 1 + u0.α) *
                        clausens_expansion_remainder(x * (1 + t), s, 1, M)

                    t * res
                end

                -part1 - part2
            end

            TaylorModel(ArbSeries((c, Δ), degree = 1), u0.α, u0.α0)
        end

        # Compute primitive_mul_x_onepα(π / x) =
        # 2(clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x^(1 - α)
        primitive_mul_x_onepα_pi_div_x = let Ms = 2 - Mα
            res = zero(primitive_mul_x_onepα_zero)

            # The function is even around π so it is beneficial to
            # always take an even number of terms
            N = iseven(M) ? M : M + 1

            # We skip m = 0 since it is cancelled. We take only even
            # values of m since the function is even around x = π
            for m = 2:2:N-1
                # m-th derivative at x = π, note that m is always even
                deriv = (-1)^(m ÷ 2) * clausenc(Arb(π), Ms - m)

                # We have x^(m - 1 + α) since we divide by x^(1 - α)
                res += deriv * abspow(x, m - 1 + Mα) / factorial(m)
            end

            # Remainder term
            # Interval for the Taylor expansion
            J = union(Arb(π), Arb(π) + x)
            # Enclosure of N-th derivative on interval
            deriv = (-1)^(N ÷ 2) * clausenc(J, Ms - N)

            # Add enclosure of the remainder term
            res += deriv * abspow(x, N - 1 + Mα) / factorial(N)

            2res
        end

        I_mul_x_onepα =
            primitive_mul_x_onepα_zero - 2primitive_mul_x_onepα_root +
            primitive_mul_x_onepα_pi_div_x

        I_mul_x_onepα_div_pi = I_mul_x_onepα / Arb(π)

        if skip_div_u0
            res = I_mul_x_onepα_div_pi * abspow(x, -Mα)
        else
            res =
                I_mul_x_onepα_div_pi /
                truncate(eval_expansion(u0, u0_expansion, x, offset_i = -1), degree = 0)
        end

        # The constant term should be exactly 1. This holds no matter
        # if we divide by u0 or not since the constant term of u0 is
        # exactly 1.
        @assert Arblib.contains(I_mul_x_onepα_div_pi.p[0], 1) ||
                !isfinite(I_mul_x_onepα_div_pi.p[0])
        res.p[0] = 1

        return res
    end
end
