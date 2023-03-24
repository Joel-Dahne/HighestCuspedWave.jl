"""
    _integrand_compute_root(::Type{KdVZeroAnsatz}, x::Arb, αₗ::Arb)

Compute enclosures of the unique root of
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
in `t` on the interval `[0, 1]`, one enclosure for `α = 0` and one
valid for `α ∈ [αₗ, 0]`. It assumes that `0 <= x <= π`.

This is the root occurring in [`lemma_I_hat_root`](@ref). However the
lemma only discusses `α < 0` and to consider the case `α = 0` we have
to normalise to get a well defined limit, more precisely we consider
the root of the function
```
(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) / α
```
See the proof of [`lemma_kdvzero_U0_constant_term`](@ref) where this
limiting root is discussed. Except for that it has to handle the case
`α = 0` the implementation is similar to the version for
[`FractionalKdVAnsatz`](@ref).

To compute the root for `α ∈ [αₗ, 0]` it uses that the root in
practice is decreasing in `α`. It computes the root at `α = 0` and `α
= αₗ` and does an a posteriori check for the root being decreasing in
`α`. The a posteriori check is based on
[`lemma_I_hat_root_alpha`](@ref) and done in the same was as in the
corresponding `FractionalKdVAnstaz` version of this method, by
checking that the derivative w.r.t. `α` is positive. To be able to
compute the derivative near `α` we have to normalise the function by
dividing by `α` as mentioned above.

To compute the root at `αₗ` it uses
`_integrand_compute_root(FractionalKdVAnsatz, x, αₗ)`.

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
    # For x = 0 the version for FractionalKdVAnsatz is able to compute
    # the root without further work.

    # Compute root for α = 0 and x = 0
    root0_zero = _integrand_compute_root(FractionalKdVAnsatz, Arb(0), Arb(0))

    if iszero(x)
        # Compute root for α = [αₗ, 0] and x = 0
        root_zero = _integrand_compute_root(FractionalKdVAnsatz, Arb(0), Arb((αₗ, 0)))

        return root0_zero, root_zero
    end

    # If x is wide we don't need to compute the root to very high
    # precision. We therefore take the absolute tolerance to depend on
    # its radius.
    atol = Arblib.mul_2exp!(Mag(), radius(x), -20)

    # Compute root for α = 0 at a given x
    compute_root0(x::Arb) =
        let
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
                root = ArbExtras.refine_root(f, Arb(only(roots)); atol)
            else
                root = Arb((roots[1][1], roots[end][2]))
            end

            return root
        end

    root_is_decreasing(x::Arb, t::Arb, α::Arb) =
        let
            # Use when α overlaps zero
            f1 = ArbExtras.derivative_function() do α
                fx_div_x(
                    α -> _integrand_I_hat(x, t, α),
                    α,
                    force = true,
                    enclosure_degree = -1,
                    extra_degree = 2,
                )
            end

            # Use when α doesn't overlap zero
            f2 = ArbExtras.derivative_function(α -> _integrand_I_hat(x, t, α) / α)

            dfdα = ArbExtras.enclosure_series(α, degree = 4) do α
                if α isa Arb && abs(α) < 1e-10
                    α = union(α, zero(α)) # Avoid using very small non-zero α
                end

                if α isa ArbSeries && Arblib.contains_zero(α[0]) ||
                   α isa Arb && Arblib.contains_zero(α)
                    f1(α)
                else
                    f2(α)
                end
            end

            Arblib.isnegative(dfdα)
        end

    # Compute an enclosure of the root for α = 0 as well as a lower
    # bound of the root on the interval [αₗ, 0].

    xₗ, xᵤ = getinterval(Arb, x)
    xᵤ = min(Arb(π), xᵤ) # We assume that xᵤ <= π
    ϵ = eps(Arb)

    if iszero(x)
        rootₗ = root0_zero
        root0 = rootₗ
    elseif Arblib.overlaps(xᵤ, Arb(π))
        rootₗ = Arb(1 // 2)
        root0 = Arb((rootₗ, compute_root0(xₗ))) # Lower bound is 1 / 2
    elseif !iswide(x)
        rootₗ = compute_root0(x) # In this case x never overlaps zero
        root0 = rootₗ
    elseif xᵤ < ϵ
        rootₗ = compute_root0(ϵ)
        root0 = Arb((rootₗ, root0_zero))
    elseif xₗ < ϵ
        rootₗ = compute_root0(xᵤ)
        root0 = Arb((rootₗ, root0_zero))
    else
        rootₗ = compute_root0(xᵤ)
        root0 = Arb((rootₗ, compute_root0(xₗ)))
    end

    # Compute an upper bound of root on the interval [αₗ, 0] by
    # computing the root at αₗ
    if iszero(αₗ)
        rootᵤ = root0
    elseif -1 < αₗ < 0
        # Compute root at αₗ
        rootᵤ = _integrand_compute_root(FractionalKdVAnsatz, xₗ, αₗ)
    else
        throw(ArgumentError("requires that -1 < αₗ < 0"))
    end

    if !Arblib.overlaps(xᵤ, Arb(π)) && !root_is_decreasing(xᵤ, rootₗ, Arb((αₗ, 0)))
        # Root is always lower bounded by 1 / 2. In practice we expect
        # this to fail only for x very close to π, in which case the
        # lower bound is fairly tight.
        rootₗ = Arb(1 // 2)
    end
    if !root_is_decreasing(xₗ, rootᵤ, Arb((αₗ, 0)))
        # Root is always upper bounded by the value at x = 0. In
        # practice we don't expect this case to occur.
        rootᵤ = _integrand_compute_root(FractionalKdVAnsatz, Arb(0), Arb((αₗ, 0)))
    end

    return root0, Arb((rootₗ, rootᵤ))
end

"""
    T0(u0::KdVZeroAnsatz, ::Ball)

Return a function such that `T0(u0)(x)` computes a Taylor model of
```
inv(π * x * u0(x)) * ∫abs(clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * y dy
```
, where the integration is taken from `0` to `π`, in `α` around `α =
0` of degree `0`. Note that this only handles the case `u0.α0 = 0`.

This method is similar to [`T0_p_one`](@ref) in that it explicitly
computes the integral. It computes an expansion in `α` around `α = 0`
and treats the zero of the integrand differently, otherwise they are
very similar.

# Formula for the integral
The switch of variables to `t = y / x` gives us
```
x / (π * u0(x)) * ∫abs(clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)) * t dt
```
from `0` to `π / x`. Call the integral, without the factors in front,
`I`.

As in [`T0_p_one`](@ref), see also
[`lemma_U0_primitive_weight_x`](@ref), we have
```
x * I = primitive_mul_x(0) - 2primitive_mul_x(root) + primitive_mul_x(π / x)
```
where
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

Similarly to in [`T0_p_one`](@ref) we can simplify
`primitive_mul_x(t)` when `t = 0` and `t = π / x`, we get
```
primitive_mul_x(0) = 2clausencmzeta(x, 2 - α) / x
```
and
```
primitive_mul_x(π / x) = 2(clausenc(x + π, 2 - α) - clausenc(Arb(π), 2 - α)) / x
```

# Computing Taylor model of `x * I`
For `primitive_mul_x(0)` we can directly compute a Taylor model. This
is the case also for `primitive_mul_x(π / x)`, as long as `x` doesn't
overlap `π`. If `x` overlaps with `π` we use that [`clausenc`](@ref)
is `2π` periodic and even to instead evaluate
```
primitive_mul_x(π / x) = 2(clausenc(π - x, 2 - α) - clausenc(Arb(π), 2 - α)) / x
```
We then handle `clausenc(π - x, 2 - α)` by expanding at zero and using
the approach from the asymptotic version of this method.

Handling `primitive_mul_x(root)` requires slightly more work. The main
reason is that the root also depends on `α`. In this case we compute
the constant term and the remainder term in the Taylor model
separately.

For the constant term of `primitive_mul_x(root)` we simply compute it
for `α = 0`, giving us
```
(
    clausenc(x * (1 - r0), 2) + clausenc(x * (1 + r0), 2) - 2clausenc(x * r0, 2)
) / x + r0 * (
    -clausens(x * (1 - r0), 1) + clausens(x * (1 + r0), 1) - 2clausens(x * r0, 1)
)
```
where `r0` is the root for `α = 0`.

For the remainder term we want to enclose the derivative w.r.t. `α` on
the interval `u0.α`. This is similar to
[`equation_kdvzero_K_dα`](@ref), except that we here have an extra
division by `x`. We also explain the procedure below.

Differentiating `primitive_mul_x(t)` with respect
to `α` while treating `t` as a function of `α` gives us
```
(
    - x * dt(α) * (
        -clausens(x * (1 - t(α)), 1 - α) +
        clausens(x * (1 + t(α)), 1 - α) -
        2clausens(x * t(α), 1 - α))
    - (
        -clausenc(x * (1 - t(α)), 2 - α, 1) +
        clausenc(x * (1 + t(α)), 2 - α, 1) -
        2clausenc(x * t(α), 2 - α, 1)
    )
) / x + dt(α) * (
    -clausens(x * (1 - t(α)), 1 - α) + clausens(x * (1 + t(α)), 1 - α) - 2clausens(x * t(α), 1 - α)
) + t(α) * (
    x * dt(α) * (clausenc(x * (1 - t(α)), -α) + clausenc(x * (1 + t(α)), -α) - 2clausenc(x * t(α), -α))
    - (
        -clausens(x * (1 - t(α)), 1 - α, 1) +
        clausens(x * (1 + t(α)), 1 - α, 1) -
        2clausens(x * t(α), 1 - α, 1)
    )
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
for that root. If we let `t(α) = root` what remains is therefore
```
- (
    clausenc(x * (1 - root), 2 - α, 1) +
    clausenc(x * (1 + root), 2 - α, 1) -
    2clausenc(x * root, 2 - α, 1)
) / x -  t(α) * (
    -clausens(x * (1 - root), 1 - α, 1) +
    clausens(x * (1 + root), 1 - α, 1) -
    2clausens(x * root, 1 - α, 1)
)
```
We get an enclosure of the derivative by enclosing this in `α`. An
important observation is that the derivative of the root with respect
to `α` **does not** affect the derivative of the result. It is
therefore enough to compute only an enclosure of the root and we do
not need to compute its derivative in `α`.

# Constant term of Taylor model
From [`lemma_kdvzero_T0_constant_term`](@ref) we have that the
constant term in the final Taylor model should be exactly `1`, we
hence explicitly set it to this in the end.

In principle we could avoid computing the constant term since we know
that it should be exactly `1`. However we choose to compute it anyway
and verify that it encloses `1`.
"""
function T0(u0::KdVZeroAnsatz, ::Ball)
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
        primitive_mul_x_pi_div_x = let
            # If x overlaps with π this gives an indeterminate result
            # which we handle specially
            if Arblib.overlaps(x, Arb(π))
                # Use periodicity and evenness of 2π to instead
                # evaluate at π - x which is close to zero. Use the
                # asymptotic expansion at x = 0 to evaluate it. To
                # compute the expansion around x = 0 it uses the same
                # approach as the asymptotic version of T0 does.
                clausenc_x_plus_pi = let y = π - x, Ms = 2 - Mα, M = 2
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
                    clausenc_x_plus_pi = gamma_sin * compose(e -> abspow(y, e), 1 - Mα)

                    # Analytic terms
                    clausenc_x_plus_pi += sum(
                        (-1)^m * compose(zeta, Ms - 2m) * abspow(y, 2m) / factorial(2m) for m = 0:M-1
                    )
                    # Remainder term
                    clausenc_x_plus_pi +=
                        abspow(y, 2M) *
                        compose(s -> clausenc_expansion_remainder(y, s, M), Ms)

                    clausenc_x_plus_pi
                end
            else
                clausenc_x_plus_pi = clausenc(x + π, 2 - Mα)
            end

            2(clausenc_x_plus_pi - clausenc(Arb(π), 2 - Mα)) / x
        end

        I_mul_x =
            primitive_mul_x_zero - 2primitive_mul_x_root + primitive_mul_x_pi_div_x

        res = I_mul_x / (Arb(π) * truncate(u0(x), degree = 0))

        # The constant term should be exactly 1
        @assert Arblib.contains(res.p[0], 1) || !isfinite(res.p[0])
        res.p[0] = 1

        return res
    end
end

"""
    T0(u0::KdVZeroAnsatz, ::Asymptotic)

Return a function such that `T0(u0)(x)` computes a Taylor model of
```
inv(π * x * u0(x)) * ∫abs(clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * y dy
```
, where the integration is taken from `0` to `π`, in `α` around `α =
0` of degree `0`. It is computed in a way that works for `x` close to
zero.

The method is similar to the non-asymptotic version but it evaluates
the terms that depend on `x` in an asymptotic way. From the
non-asymptotic version we get that we want to compute
```
(primitive_mul_x(0) - 2primitive_mul_x(root) + primitive_mul_x(π / x)) / (π * u0(x))
```
The leading term in `u0(x)` behaves like `x^-α`and we want to
explicitly cancel this term. We therefore rewrite the above as
```
(primitive_mul_x(0) * x^α - 2primitive_mul_x(root) * x^α + primitive_mul_x(π / x) * x^α) / (π * u0(x) * x^α)
```

The division by `π * u0(x) * x^α` is straight forward to handle. We
now describe how to compute the numerator.

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

What remains is computing a Taylor model of `primitive_mul_x(root) *
x^α`, this requires slightly more work.

# Computing Taylor model of `primitive_mul_x(root) * x^α`
Similarly to the non-asymptotic version we manually compute the
constant term and the remainder. To begin with we note that
```
primitive_mul_x(root) * x^α =
            (
                clausenc(x * (1 - root), 2 - α) + clausenc(x * (1 + root), 2 - α) -
                2clausenc(x * root, 2 - α)
            ) / x^(1 - α) +
            root * (
                -clausens(x * (1 - root), 1 - α) + clausens(x * (1 + root), 1 - α) -
                2clausens(x * root, 1 - α)
            ) * x^α
```
is exactly equal to `K(x, α) * x^(α - 1)` with `K(x, α)` as in
[`equation_kdvzero_K`](@ref).

For the constant term in the Taylor model `primitive_mul_x(root) *
x^α` we simply expand the Clausen functions and evaluate.

For the remainder term we need to enclose the derivative w.r.t. `α`.
An expansion of the derivative is given in
[`lemma_kdvzero_U0_primitive_K_expansion`](@ref). It consists of six
terms. The first two are from the leading terms of the `clausenc` and
`clausens` terms respectively, the removable singularities of `gamma(s
- 1) * sinpi(s / 2)` and `gamma(s - 1) * cospi(s / 2)` are handled as
described above. The remaining four terms are all infinite sums, where
we explicitly sum the first `M - 1` terms. To enclose the remainder we
use [`clausenc_expansion_remainder`](@ref) and
[`clausens_expansion_remainder`](@ref).


# Constant term of Taylor model
Similar to in the non-asymptotic version the constant term in the
final Taylor model should be exactly `1`, we hence explicitly set it
to this in the end.

In principle we could avoid computing the constant term since we know
that it should be exactly `1`. However we choose to compute it anyway
and verify that it encloses `1`.
"""
function T0(u0::KdVZeroAnsatz, ::Asymptotic; ϵ::Arb = one(Arb), M::Integer = 10)
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

    # Precompute clausenc(Arb(π), 2 - α - m) for m = 2:2:M
    clausenc_pi = OrderedDict(m => clausenc(Arb(π), 2 - Mα - m) for m = 2:2:M)

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

        # Compute primitive_mul_x(root) * x^α
        # This requires a significant amount of work
        primitive_mul_x_onepα_root = let
            r0, r = _integrand_compute_root(typeof(u0), x, lbound(Arb, u0.α))

            # Constant term at α = 0
            c = let
                # Enclosure of
                # (
                #   clausenc(x * (1 - r0), 2) +
                #   clausenc(x * (1 + r0), 2) -
                #   2clausenc(x * r0, 2)
                # ) / x
                part1 = let s = Arb(2)
                    # Singular term
                    res =
                        gamma_sin.p[0] * (
                            abspow(1 - r0, s - 1) + abspow(1 + r0, s - 1) -
                            2abspow(r0, s - 1)
                        )
                    # Analytic terms
                    res += sum(
                        (-1)^m *
                        zeta(s - 2m) *
                        abspow(x, 2m - 1) *
                        ((1 - r0)^2m + (1 + r0)^2m - 2r0^2m) /
                        factorial(2m) for m = 1:M-1
                    )
                    # Remainder term
                    # Here we use that max(1 - r0, 1 + r0, r0) = 1 + r0 so
                    # it is enough to compute the remainder term at x * (1 + r0)
                    res +=
                        ((1 - r0)^2M + (1 + r0)^2M - 2r0^2M) *
                        abspow(x, 2M - 1) *
                        clausenc_expansion_remainder(x * (1 + r0), s, M)

                    res
                end

                # Enclosure of
                # r0 * (
                #   -clausens(x * (1 - r0), 1) +
                #   clausens(x * (1 + r0), 1) -
                #   2clausenc(x * r0, 1)
                # )
                part2 = let s = Arb(1)
                    # Singular term
                    res =
                        -gamma_sin.p[0] * (
                            -abspow(1 - r0, s - 1) + abspow(1 + r0, s - 1) -
                            2abspow(r0, s - 1)
                        )
                    # Analytic terms
                    res += sum(
                        (-1)^m *
                        zeta(s - 2m - 1) *
                        abspow(x, 2m + 1) *
                        (-(1 - r0)^(2m + 1) + (1 + r0)^(2m + 1) - 2r0^(2m + 1)) / factorial(2m + 1) for m = 0:M-1
                    )
                    # Remainder term
                    res +=
                        ((1 - r0)^(2M + 1) + (1 + r0)^(2M + 1) - 2r0^(2M + 1)) *
                        abspow(x, 2M + 1) *
                        clausens_expansion_remainder(x * (1 + r0), s, M)

                    r0 * res
                end

                part1 + part2
            end

            # Enclosure of derivative over interval
            Δ = let
                # Enclosure of leading terms from clausenc functions
                leading_1 = let s = ArbSeries((2 - u0.α, 1))
                    gamma_sin_tmp =
                        fx_div_x(sm2 -> sinpi((sm2 + 2) / 2), s - 2) /
                        fx_div_x(sm2 -> rgamma(1 - (sm2 + 2)), s - 2)

                    expansion =
                        gamma_sin_tmp * (
                            abspow(1 - r, s - 1) + abspow(1 + r, s - 1) - 2abspow(r, s - 1)
                        )

                    expansion[1]
                end

                # Enclosure of leading terms from clausens functions
                leading_2 = let s = ArbSeries((1 - u0.α, 1))
                    gamma_sin_tmp =
                        fx_div_x(sm1 -> cospi((sm1 + 1) / 2), s - 1) /
                        fx_div_x(sm1 -> rgamma(1 - (sm1 + 1)), s - 1)

                    expansion =
                        gamma_sin_tmp * (
                            -abspow(1 - r, s - 1) + abspow(1 + r, s - 1) -
                            2abspow(r, s - 1)
                        )

                    expansion[1]
                end

                # For the remainder terms we use that max(1 - t, 1 +
                # t, t) = 1 + t so it is enough to compute the
                # remainder term at x * (1 + t)

                # Enclosure of first sum times log(x)
                S1_log = sum(1:M-1) do m
                    (-1)^m *
                    zeta(2 - u0.α - 2m) *
                    logabspow(x, 1, 2m - 1 + u0.α) *
                    ((1 - r)^2m + (1 + r)^2m - 2r^2m) / factorial(2m)
                end
                S1_log +=
                    ((1 - r)^2M + (1 + r)^2M - 2r^2M) *
                    logabspow(x, 1, 2M - 1 + u0.α) *
                    clausenc_expansion_remainder(x * (1 + r), 2 - u0.α, M)

                # Enclosure of second sum times log(x)
                S2_log = sum(0:M-1) do m
                    (-1)^m *
                    zeta(1 - u0.α - 2m - 1) *
                    logabspow(x, 1, 2m + 1 + u0.α) *
                    (-(1 - r)^(2m + 1) + (1 + r)^(2m + 1) - 2r^(2m + 1)) /
                    factorial(2m + 1)
                end
                S2_log +=
                    (-(1 - r)^(2M + 1) + (1 + r)^(2M + 1) - 2r^(2M + 1)) *
                    logabspow(x, 1, 2M + 1 + u0.α) *
                    clausens_expansion_remainder(x * (1 + r), 1 - u0.α, M)

                # Third sum
                S3 = sum(1:M-1) do m
                    (-1)^m *
                    dzeta(2 - u0.α - 2m) *
                    abspow(x, 2m - 1 + u0.α) *
                    ((1 - r)^2m + (1 + r)^2m - 2r^2m) / factorial(2m)
                end
                S3 +=
                    ((1 - r)^2M + (1 + r)^2M - 2r^2M) *
                    abspow(x, 2M - 1 + u0.α) *
                    clausenc_expansion_remainder(x * (1 + r), 2 - u0.α, 1, M)

                # Fourth sum
                S4 = sum(0:M-1) do m
                    (-1)^m *
                    dzeta(1 - u0.α - 2m - 1) *
                    abspow(x, 2m + 1 + u0.α) *
                    (-(1 - r)^(2m + 1) + (1 + r)^(2m + 1) - 2r^(2m + 1)) /
                    factorial(2m + 1)
                end
                S4 +=
                    (-(1 - r)^(2M + 1) + (1 + r)^(2M + 1) - 2r^(2M + 1)) *
                    abspow(x, 2M + 1 + u0.α) *
                    clausens_expansion_remainder(x * (1 + r), 1 - u0.α, 1, M)

                # Combine the terms
                -leading_1 - r * leading_2 + S1_log + r * S2_log - S3 - r * S4
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
                # Note that clausenc_pi[m] = clausenc(Arb(π), Ms - m)
                deriv = (-1)^(m ÷ 2) * clausenc_pi[m]

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

        # res = I_mul_x_onepα / (π * u0(x) * x^α)
        res =
            I_mul_x_onepα / (
                Arb(π) * truncate(
                    eval_expansion(u0, u0_expansion, x, offset_i = -1),
                    degree = 0,
                )
            )

        # The constant term should be exactly 1. This holds no matter
        # if we divide by u0 or not since the constant term of u0 is
        # exactly 1.
        @assert Arblib.contains(res.p[0], 1) || !isfinite(res.p[0])
        res.p[0] = 1

        return res
    end
end
