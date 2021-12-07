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
a unique root, it is negative to the left of the root and positive to
the right. If we let `root_lower` and `root_upper` be lower and upper
bounds of the roots respectively we have that the integral is given by
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

# Isolating the root
Let
```
f(t) = clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)
```
be the expression inside the absolute value. For isolating the root on
`[0, 1]` we make use of the fact that `f` is increasing on the interval
- **PROVE:** That `f` is increasing on the interval `[0, 1]`.

Further we notice that for `t = 1 / 2` we have
```
f(1 / 2) = clausenc(x / 2, -α) + clausenc(3x / 2, -α) - 2clausenc(x / 2, -α) =
    clausenc(3x / 2, -α) - clausenc(x / 2, -α)
```
which is negative for all `0 < x < π`.
- **PROVE:** Prove the above using the monotonicity of `clausenc`, it
    is straight forward.
It follows that the root is lower bounded by `1 / 2`. In practice the
root is upper bounded by something around `0.8` or even lower. To find
the root we therefore start at some point above this upper bound and
check that the value is positive, we then check smaller and smaller
points until we find a place where it becomes negative. As long as the
value at the first point was positive this gives a crude enclosure of
the root. The enclosure is then be refined using
[`ArbExtras.isolate_roots`](@ref) and [`ArbExtras.refine_root`](@ref).
- **IMPROVE:** We could possibly get better enclosures for the root
  for wide values of `x`.
"""
function T0_p_one(u0::FractionalKdVAnsatz, evaltype::Ball = Ball(); skip_div_u0 = false)
    @assert isone(u0.p)

    α = u0.α

    return x::Arb -> begin
        # Attempt to isolate the root of clausenc(x * (1 - t), mα) +
        # clausenc(x * (1 + t), mα) - 2clausenc(x * t, mα)
        f(t) = clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) - 2clausenc(x * t, -α)

        # The root is lower bounded by 1 / 2
        root_lower = Arf(0.5)

        # Find a crude upper bound for the root
        δ = Arb(0.4)
        Arblib.ispositive(f(root_lower + δ)) || return Arblib.indeterminate!(zero(x))
        while Arblib.ispositive(f(root_lower + δ / 2)) && δ > 1e-5
            Arblib.mul_2exp!(δ, δ, -1)
        end
        root_upper = ubound(root_lower + δ)

        # Improve the enclosure of the root
        roots, flags = ArbExtras.isolate_roots(f, root_lower, root_upper, depth = 5)
        if length(flags) == 1 && flags[1]
            # Refine the unique root
            root = ArbExtras.refine_root(f, Arb(only(roots)))
            # Get lower and upper bounds for the root
            root_lower, root_upper = getinterval(root)
        else
            # Get lower and upper bounds for possible roots
            root_lower, root_upper = roots[1][1], roots[end][2]
        end

        root_lower, root_upper = Arb(root_lower), Arb(root_upper)

        integrand(t; analytic) = begin
            # Check that the real part of t is strictly between 0 and
            # 1 or return an indeterminate result
            Arblib.ispositive(Arblib.realref(t)) && Arblib.realref(t) < 1 ||
            return Arblib.indeterminate!(zero(t))

            if isreal(t)
                rt = Arblib.realref(t)

                res = Acb(
                    clausenc(x * (1 - rt), -α) + clausenc(x * (1 + rt), -α) -
                    2clausenc(x * rt, -α),
                )
            else
                res =
                    clausenc(x * (1 - t), -α) + clausenc(x * (1 + t), -α) -
                    2clausenc(x * t, -α)
            end

            return Arblib.real_abs!(res, res, analytic) * t
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

        #@assert Arblib.overlaps(x * primitive(Arb(0.5)), primitive_mul_x(Arb(0.5)))

        # primitive(0)
        primitive_mul_x_zero = 2clausencmzeta(x, 2 - α) / x

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
            clausenc_x_plus_pi = C * abspow(y, e) + P(y) + E * abspow(y, 2M)
        else
            clausenc_x_plus_pi = clausenc(x + π, 2 - α)
        end
        # eta is only implemented for Acb in Arblib
        primitive_mul_x_pi_div_x = 2(clausenc_x_plus_pi + real(eta(Acb(2 - α)))) / x

        I12 = real(
            Arblib.integrate(
                integrand,
                root_lower,
                root_upper,
                check_analytic = true,
                rtol = 1e-5,
                atol = 1e-5,
                opts = Arblib.calc_integrate_opt_struct(0, 1_000, 0, 0, 0),
            ),
        )

        I_mul_x = (
            primitive_mul_x_zero - primitive_mul_x(root_lower) + x * I12 -
            primitive_mul_x(root_upper) + primitive_mul_x_pi_div_x
        )

        if skip_div_u0
            return I_mul_x / π
        else
            return I_mul_x / (π * u0(x))
        end
    end
end
