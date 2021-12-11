"""
    T0(u0::KdVZeroAnsatz, ::Ball; skip_div_u0 = false)

Compute the integral \$T_0\$.

This method is similar to [`T0_p_one`](@ref) in that it explicitly
computes the integral. It computes an expansion in `α` around `α = 0`
and treats the zero of the integrand differently, otherwise they are
very similar.

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
the right. If we let `root` be an enclosure of the root then an
enclosure of the integral on `[0, 1]` is given by
```(
I1 = -(primitive(root) - primitive(0)) + (primitive(1) - primitive(root))
   = primitive(0) - 2primitive(root) + primitive(π / x)
```
We also get that
```
primitive(0) = 2clausencmzeta(x, 2 - α) / x^2
```
and
```
primitive(π / x) = 2(clausenc(x + π, 2 - α) - clausenc(π, 2 - α)) / x^2
```
- **TODO:** Handle computation of `clausenc(x + π, 2 - α)` when `x`
  overlaps with `π`.

We can notice that all terms in the result, except the `I12` term,
contains a division by `x` and that we in the end multiply with `x` If
we let
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

We are now interested in computing an expansion, in `α`, of this. This
can be done with direct evaluation using `ArbSeries`. However we need
to prove that the constant term in the expansion is exactly `π` and
this requires some additional work.

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

- **TODO:** Figure out how to handle the enclosure of the root. It
  seems like we don't need to expand it in `α` but it is enough to use
  an enclosure of the root. However the equation for the root becomes
  singular at `α = 0` so we need to handle that.
- **TODO:** Improve enclosure for wide values of `x`. This we will
  most likely need to do in the end.
- **TODO:** Compute the remainder terms.
"""
function T0(u0::KdVZeroAnsatz, ::Ball; skip_div_u0 = false)
    α = ArbSeries((0, 1), degree = 2)

    return x::Arb -> begin
        # Isolate the root of clausenc(x * (1 - t), -α) + clausenc(x *
        # (1 + t), -α) - 2clausenc(x * t, -α)
        # FIXME: For now we use α = lbound(u0.α), this should be
        # updated to work for the whole enclosure in α
        f(t) =
            let α = lbound(u0.α)
                clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
            end

        # The root is lower bounded by 1 / 2
        root_lower = Arf(0.5)

        # Find a crude upper bound for the root
        δ = Arb(0.4)
        Arblib.ispositive(f(root_lower + δ)) || return ArbSeries((1, NaN))
        while Arblib.ispositive(f(root_lower + δ / 2)) && δ > 1e-5
            Arblib.mul_2exp!(δ, δ, -1)
        end
        root_upper = ubound(root_lower + δ)

        # Improve the enclosure of the root
        roots, flags = ArbExtras.isolate_roots(f, root_lower, root_upper, depth = 5)
        if length(flags) == 1 && flags[1]
            # Refine the unique root
            root = ArbExtras.refine_root(f, Arb(only(roots)))
        else
            root = Arb((roots[1][1], roots[end][2]))
        end

        primitive_mul_x(t) =
            (
                clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) -
                2clausenc(x * t, 2 - α)
            ) / x +
            t * (
                -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) -
                2clausens(x * t, 1 - α)
            )

        primitive_mul_x_zero = 2clausencmzeta(x, 2 - α) / x

        # primitive(π / x)
        # If x overlaps with π this gives an indeterminate result
        # which we handle specially
        if Arblib.overlaps(x, Arb(π))
            # TODO: Handle this case. We probably have to use the
            # asymptotic expansion and expand that in α
            @warn "x overlaps π, not implemented yet"
            clausenc_x_plus_pi = clausenc(x + π, 2 - α)
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

**TODO:** Implement this.
"""
function T0(u0::KdVZeroAnsatz, ::Asymptotic)
    return x::Arb -> begin
        # TODO: Implement this
        return ArbSeries((1, 0))
    end
end
