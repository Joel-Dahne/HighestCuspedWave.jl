"""
    T0(u0::FractionalKdVAnsatz{Arb}, ::Ball)

**TODO:** Tune the values for `δ0, δ1, δ2` depending on `u0.α` and
`u0.p`.
"""
function T0(
    u0::FractionalKdVAnsatz{Arb},
    evaltype::Ball;
    δ0::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    δ1::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    δ2::Arf = Arf(1e-2),
    ϵ::Arb = 1 + u0.α,
    skip_div_u0 = false,
)
    # Use specialised implementation in the case the weight is x
    isone(u0.p) && return T0_p_one(u0, evaltype; skip_div_u0)

    f = T01(u0, evaltype, skip_div_u0 = true; δ0, δ1)
    g = T02(u0, evaltype, skip_div_u0 = true; δ2, ϵ)

    return x -> begin
        ## Integral on [0, x]
        part1 = f(x)

        # Short circuit on a non-finite result
        isfinite(part1) || return part1

        ## Integral on [x, π]
        part2 = g(x)

        isfinite(part2) || return part2

        if skip_div_u0
            return part1 + part2
        else
            return (part1 + part2) / u0(x)
        end
    end
end

"""
    T0_p_one(u0, Ball())

Compute the integral \$T_0\$ for `u0` with `u0.p == 1`.

In this case the integral is given by
```
inv(π * x * u0(x)) * ∫abs(clausenc(x - y, -α) + clausenc(x + y, -α) - 2clausenc(y, -α)) * y dy
```
with the integral taken from `0` to `π`. The switch of coordinates to
`t = y / x` gives us
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
Call this function `primitive(t)`. The idea is to isolate the parts of
the interval where
```
clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
is zero. For the other parts of the interval we can remove the
absolute value (taking the appropriate sign) and integrate explicitly
using `primitive(t)`.

On the interval `[1, π / x]` the expression inside the absolute value
is positive and we can just remove it. This gives us the integral
```
I2 = primitive(π / x) - primitive(1)
```
- **PROVE:** That the integrand is positive.

On the interval `[0, 1]` the expression inside the absolute value has
a unique root which we can isolate with
[`_integrand_compute_root`](@ref), it is negative to the left of the
root and positive to the right. If we let `root_lower` and
`root_upper` be lower and upper bounds of the roots respectively we
have that the integral is given by
```
I1 = -I11 + I12 + I13
```
where
```
I11 = primitive(root_lower) - primitive(0)
```
is the integral on `[0, root_lower]`.
```
I12
```
is the integral on `[root_lower, root_upper]` which we enclose
directly.
```
I13 = primitive(1) - primitive(root_upper)
```

Putting all of this together we get
```
I1 + I2 = -(primitive(root_lower) - primitive(0)) + I12 +
    primitive(1) - primitive(root_upper) + primitive(π / x) - primitive(1)
```
which we can simplify to
```
I = I1 + I2 = primitive(0) - primitive(root_lower) + I12 -
    primitive(root_upper) + primitive(π / x)
```
We also get that
```
primitive(0) = 2clausencmzeta(x, 2 - α) / x^2
```
and
```
primitive(π / x) = 2(clausenc(x + π, 2 - α) - clausenc(Arb(π), 2 - α)) / x^2
```
where `clausenc(Arb(π), 2 - α)` can also be given as the negated
alternating zeta function, `-eta(2 - α)`.

We can notice that all terms in the result, except the `I12` term,
contains a division by `x` and that we in the end multiply with `x` If
we let
```
primitive_mul_x(t) = (clausenc(x * (1 - t), 2 - α) / x - t * clausens(x * (1 - t), 1 - α)) +
    (clausenc(x * (1 + t), 2 - α) / x + t * clausens(x * (1 + t), 1 - α)) -
    2(clausenc(x * t, 2 - α) / x + t * clausens(x * t, 1 - α))
```
we get
```
x * I = primitive_mul_x(0) - primitive_mul_x(root_lower) + x * I12 -
    primitive_mul_x(root_upper) + primitive_mul_x(π / x)
```
- **IMPROVE:** We could get much better enclosures for wide values of
  `x` by expanding with `ArbSeries` and enclosing that. This would
  require us to rewrite everything as a function explicitly depending
  on `x` but should otherwise be straight forward.
"""
function T0_p_one(u0::FractionalKdVAnsatz, evaltype::Ball = Ball(); skip_div_u0 = false)
    @assert isone(u0.p)

    α = u0.α

    return x::Arb -> begin
        root = _integrand_compute_root(u0, x)
        root_lower, root_upper = getinterval(Arb, root)

        primitive_mul_x(x, t) =
            (
                clausenc(x * (1 - t), 2 - α) + clausenc(x * (1 + t), 2 - α) -
                2clausenc(x * t, 2 - α)
            ) / x +
            t * (
                -clausens(x * (1 - t), 1 - α) + clausens(x * (1 + t), 1 - α) -
                2clausens(x * t, 1 - α)
            )

        # primitive(0)
        primitive_mul_x_zero(x) = 2clausencmzeta(x, 2 - α) / x

        # primitive(π / x)
        # If x overlaps with π this gives an indeterminate result
        # which we handle specially
        if Arblib.overlaps(x, Arb(π))
            # Use periodicity of 2π to evaluate at x - π which is
            # close to zero. Use the asymptotic expansion to evaluate
            # it.
            y = x - π
            M = 3
            C, e, P, E = clausenc_expansion(y, 2 - α, M)
            # Enclosure of clausenc(x + π, 2 - α)
            # Note that it doesn't use the argument, it just returns
            # an enclosure valid for the global x
            clausenc_x_plus_pi = _ -> C * abspow(y, e) + P(y) + E * abspow(y, 2M)
        else
            clausenc_x_plus_pi = x -> clausenc(x + π, 2 - α)
        end
        # eta is only implemented for Acb in Arblib
        primitive_mul_x_pi_div_x(x) =
            2(clausenc_x_plus_pi(x) + real(eta(Acb(2 - α)))) / x

        # Compute a tighter enclosure by expanding in x
        I_mul_x(x) =
            primitive_mul_x_zero(x) - 2primitive_mul_x(x, root) +
            primitive_mul_x_pi_div_x(x)

        res = ArbExtras.enclosure_series(I_mul_x, x)

        if skip_div_u0
            return res / π
        else
            return res / (π * u0(x))
        end
    end
end
