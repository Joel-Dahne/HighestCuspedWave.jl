"""
    _T0_asymptotic_main_1(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G1(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫ abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
where the integration is taken from `0` to `1`, defined in
[`T0_asymptotic`](@ref).

Using that `1 - t >= 0` the problem reduces to computing an upper
bound of
```
G1(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```

As `α -> -1` and `x -> 0` we have that `G1(x)` tends to zero and we
therefore don't need to compute a very accurate upper bound.

As `α` goes to `-1` the factor `(1 - t)^(-α - 1) + (t + 1)^(-α - 1) -
2t^(-α - 1)` of the integrand goes to zero, by factoring out `(1 + α)`
we make it instead go to a non-zero function. This means we want to
analyse
```
G1(x) = (1 + α) / (1 - x^p0) * 1 / log(inv(x))
            ∫ abs((1 - t)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
For the factor `(1 + α) / (1 - x^p0)` we note that it is increasing in
`x` and given by `1 + α` at `x = 0`, for non-zero `x` we can handle
the removable singularity. What remains to handle is to handle
```
inv(log(inv(x))) * ∫ abs((1 - t)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) *
    t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt =
        inv(log(inv(x))) * J
```
Splitting the log-term in the weight as
```
log(c + inv(x * t)) = log((c * x * t + 1) / (x * t)) = log(1 + c * x * t) - log(x) - log(t)
```
we can split `J` into three integrals
```
J = ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) * log(1 + c * x * t) dt
    - log(x) * ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) dt
    - ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) * log(t) dt
```
By using that `log(1 + c * x * t) < log(1 + c * x)` we get
```
J <= log(c + inv(x)) * ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) dt
    - ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) * log(t) dt
```
Let
```
J1 = ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) dt

J2 = ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) * log(t) dt
```
This gives us
```
J <= log(c + inv(x)) * J1 - J2
```
and
```
G1(x) <= (1 + α) / (1 - x^p0) * (log(c + inv(x)) / log(inv(x)) * J1 - inv(log(inv(x))) * J2)
```

# Enclosing `J1` and `J2`
For computing `J1` and `J2` there are three things to take care of
1. The removable singularity of `((1 - t)^(-α - 1) + (1 + t)^(-α - 1)
  - 2t^(-α - 1)) / (1 + α)`
2. The endpoint `t = 0` where the integrands are zero but not analytic
3. The endpoint `t = 1` where the integrands are singular but
  integrable.

The removable singularity can be handled using [`fx_div_x`](@ref). For
`t` overlapping zero this doesn't quite work since the expression
blows up. Instead we rewrite it as (for `J1`)
```
abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α))
= abs(t^a * (1 - t)^(-α - 1) + t^a * (1 + t)^(-α - 1) - 2t^(a - α - 1)) / (1 + α) * t^(1 - γ * (1 + α) - a)
```
with `a = 1 / 2`. In this case both factors are finite and can be
enclosed directly. For `J2` we do the same thing but also have to
handle
```
t^(1 - γ * (1 + α) - a) * log(t)
```
Differentiating we get the derivative
```
t^(-γ * (1 + α) - a) * ((1 - γ * (1 + α) - a) * log(t) + 1)
```
which is negative for `t < exp(inv(1 - γ * (1 + α) - a))`. We hence
check so that `t` satisfies this inequality and in that case we use
monotonicity to get the upper bound zero and lower bound at
`ubound(t)`.

To handle the right endpoint we note that `t^(-γ * (1 + α))` as well
as `log(t)` are bounded for `t` close to `-1`. We can thus bound them
and factor them out of the integrals. This leaves us with the integral
```
∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t dt
```
If we remove the absolute value and let `s = -(α + 1)` a primitive
function is given by
```
(
    t^2 * (2t^s - (1 - t)^s - (1 + t)^s) / s
    + ((1 - t)^s + (1 + t)^s - 2) / s
    + t * ((1 - t)^(1 + s) - (1 + t)^(1 + s) + 2t^(1 + s))
) / ((1 + s) * (2 + s))
```
where the constant is chosen to make it finite as `s` goes to zero.
For `t = 1` this simplifies to
```
-2 * (1 - 2^-(α + 1)) / (α * (1 - α))
```
This allows us to compute the integral as long as we can remove the
absolute value. Since
```
(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)
```
is increasing in `t` it is enough to check that it is positive at the
left endpoint of the interval of integration.
"""
function _T0_asymptotic_main_1(α::Arb, γ::Arb, c::Arb)
    αp1 = Arblib.nonnegative_part!(zero(α), α + 1)

    # Primitive function of
    # ((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t
    primitive = let s = -αp1
        t -> begin
            t1 = fx_div_x(v -> 2t^v - (1 - t)^v - (1 + t)^v, s)
            t2 = fx_div_x(v -> (1 - t)^v + (1 + t)^v - 2, s)
            (t^2 * t1 + t2 + t * ((1 - t)^(1 + s) - (1 + t)^(1 + s) + 2t^(1 + s))) /
            ((1 + s) * (2 + s))
        end
    end

    # primitive(1)
    primitive_one = -2 * (1 - 2^-αp1) / (α * (1 - α))

    # Enclosure of
    # ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) dt
    J1 = begin
        integrand_J1(t; analytic = false) = begin
            if isreal(t)
                t = real(t)
                @assert !analytic
                if Arblib.contains_zero(t)
                    # Note that t is always positive in this case
                    let a = Arb(1 // 2)
                        return abs(
                            fx_div_x(
                                s ->
                                    abspow(t, a) * (1 - t)^-s + abspow(t, a) * (1 + t)^-s - 2abspow(t, a - s),
                                αp1,
                                force = true,
                            ),
                        ) * abspow(t, (1 - γ * αp1 - a))
                    end
                else
                    return abs(fx_div_x(s -> (1 - t)^-s + (1 + t)^-s - 2t^-s, αp1)) * t^(1 - γ * αp1)
                end
            else
                res = fx_div_x(s -> (1 - t)^-s + (1 + t)^-s - 2t^-s, Acb(αp1))
                Arblib.real_abs!(res, res, analytic)

                return res * t^(1 - γ * αp1)
            end
        end

        # We integrate numerically on [0, b] and explicitly on [b, 1]
        b = Arb(0.99)

        # Check that expression inside absolute value is positive on
        # [a, 1]
        @assert Arblib.ispositive(fx_div_x(s -> (1 - b)^-s + (1 + b)^-s - 2b^-s, -αp1))

        # Compute for the interval [0, b]
        part1 = real(
            Arblib.integrate(integrand_J1, 0, b, check_analytic = true, rtol = 1e-5),
        )

        # Compute for the interval [b, 1]. Note that we factor out
        # t^(-γ * (1 + α))
        part2 = let T = Arb((b, 1))
            T^(-γ * αp1) * (primitive_one - primitive(b))
        end

        part1 + part2
    end

    # Enclosure of
    # ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) * log(t) dt
    J2 = begin
        integrand_J2(t; analytic = false) = begin
            if isreal(t)
                t = real(t)
                @assert !analytic
                if Arblib.contains_zero(t)
                    # Note that t is always positive in this case
                    let a = Arb(1 // 2), tᵤ = ubound(Arb, t)
                        if t < exp(inv(1 - γ * αp1 - a))
                            return abs(
                                fx_div_x(
                                    s ->
                                        abspow(t, a) * (1 - t)^-s +
                                        abspow(t, a) * (1 + t)^-s -
                                        2abspow(t, a - s),
                                    αp1,
                                    force = true,
                                ),
                            ) * Arb((abspow(tᵤ, (1 - γ * αp1 - a)) * log(tᵤ), 0))
                        else
                            return Arblib.indeterminate!(t)
                        end
                    end
                else
                    return abs(fx_div_x(s -> (1 - t)^-s + (1 + t)^-s - 2t^-s, αp1)) *
                           t^(1 - γ * αp1) *
                           log(t)
                end
            else
                res = fx_div_x(s -> (1 - t)^-s + (1 + t)^-s - 2t^-s, Acb(αp1))
                Arblib.real_abs!(res, res, analytic)

                return res * t^(1 - γ * αp1) * log(t)
            end
        end

        # We integrate numerically on [0, b] and explicitly on [b, 1]
        b = Arb(0.99)

        # Check that expression inside absolute value is positive on
        # [a, 1]
        @assert Arblib.ispositive(fx_div_x(s -> (1 - b)^-s + (1 + b)^-s - 2b^-s, -αp1))

        # Compute for the interval [0, b]
        part1 = real(
            Arblib.integrate(integrand_J2, 0, b, check_analytic = true, rtol = 1e-5),
        )

        # Compute for the interval [b, 1]. Note that we factor out
        # t^(-γ * (1 + α)) * log(t)
        part2 = let T = Arb((b, 1))
            T^(-γ * αp1) * log(T) * (primitive_one - primitive(b))
        end

        part1 + part2
    end

    return x::Arb -> begin
        x < 1 || throw(DomainError(x, "must have x < 1"))

        # Enclosure of log(c + inv(x)) / log(inv(x))
        J1_factor = let xₗ = abs_lbound(Arb, x), xᵤ = ubound(Arb, x)
            lower = iszero(xₗ) ? one(x) : log(c + inv(xₗ)) / log(inv(xₗ))
            upper = iszero(xᵤ) ? one(x) : log(c + inv(xᵤ)) / log(inv(xᵤ))
            Arb((lower, upper))
        end

        # Enclosure of inv(log(inv(x)))
        J2_factor = let xₗ = abs_lbound(Arb, x), xᵤ = ubound(Arb, x)
            lower = iszero(xₗ) ? zero(x) : inv(log(inv(xₗ)))
            upper = iszero(xᵤ) ? zero(x) : inv(log(inv(xᵤ)))
            Arb((lower, upper))
        end

        # Enclosure of (1 + α) / (1 - x^p0)
        αp1_div_1mxp0 = if iszero(x)
            αp1
        elseif Arblib.contains_zero(x)
            lower = αp1
            upper = let xᵤ = ubound(Arb, x)
                inv(fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), αp1, extra_degree = 2))
            end
            Arb((lower, upper))
        else
            inv(fx_div_x(s -> 1 - x^(s + s^2 / 2), αp1, extra_degree = 2))
        end

        return αp1_div_1mxp0 * (J1_factor * J1 - J2_factor * J2)
    end
end

"""
    _T0_asymptotic_main_2(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G2(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫ abs(abs(1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
where the integration is taken from `1` to `π / x`, defined in
[`T0_asymptotic`](@ref).

Using that `1 - t <= 0` and that
```
(t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)
```
is positive on the interval ``[1, ∞]`` (see
[`lemma_integrand_2`](@ref)), the problem reduces to computing an upper
bound of
```
G2(x) = 1 / ((1 - x^p0) * log(inv(x))) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```

If `x` is not too small it uses [`_T0_asymptotic_main_2_1`](@ref) and
[`_T0_asymptotic_main_2_2`](@ref) to compute the integral by splitting
the interval into ``[1, 2]`` and ``[2, π / x]``. Otherwise it uses the
approach described below.
- **TODO:** What does "too small" mean? For now we always use this
  approach since the approach described below is not finished. It
  seems like it might be possible to use the below approach to handle
  the interval `[0, 1e-10]` and the above approach for larger `x`.

## Expanding integrand
We want to get rid of the fractional exponents inside the integral. To
do this we focus on the part of the integrand given by
```
((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-γ * (1 + α))
```
We can rewrite this as
```
((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-γ * (1 + α)) =
    (exp(-log(t - 1) * (1 + α)) + exp(-log(t + 1) * (1 + α)) - 2exp(-log(t) * (1 + α))) *
        exp(-γ * log(t) * (1 + α)) =
    exp(-(log(t - 1) + γ * log(t)) * (1 + α)) +
        exp(-(log(t + 1) + γ * log(t)) * (1 + α)) -
        2exp(-(1 + γ) * log(t) * (1 + α))
```
Expanding the exponentials and joining the three sums we arrive at
```
sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n)
    for n = 1:Inf
)
```
Inserting this into the integral and switching the order of
integration and summation gives us
```
G2(x) = -1 / ((1 - x^p0) * log(x)) * sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
        t * log(c + inv(x * t)) dt
    for n = 1:Inf
)
```
Where the integration is done from `1` to `π / x`. We are hence
interested in computing the integrals
```
I(n, x) = ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
    t * log(c + inv(x * t)) dt
```
for `n = 1, 2, ...`. From which we can recover `G2(x)` as
```
G2(x) = -1 / ((1 - x^p0) * log(x)) * sum((-1)^n * (1 + α)^n / factorial(n) * I(n, x) for n = 1:Inf)
```

## Computing `I(n, x)`

### Split `I(n, x)` into three parts
As a first step we split `I(n, x)` into three parts by using that
```
log(c + inv(x * t)) = log((c * x * t + 1) / (x * t)) = log(1 + c * x * t) - log(x) - log(t)
```
Letting
```
I1(n, x) = ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
    t dt

I2(n, x) = ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
    t * log(t) dt

I3(n, x) = ∫ ((log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n) *
    t * log(1 + c * x * t) dt
```
we then get
```
I(n, x) = -log(x) * I1(n, x) - I2(n, x) + I3(n, x)
```

### Bound `I3(n, x)` in terms of `I1(n, x)`
For `I3(n, x)` the factor `log(1 + c * x * t)` in the integrand is
bounded and an upper bound is given by `log(1 + c * π)`. We thus have
```
I3(n, x) <= log(1 + c * π) * I1(n, x)
```
giving us
```
I(n, x) <= (log(1 + c * π) - log(x)) * I1(n, x) - I2(n, x)
```
**FIXME:** This upper bound only holds if `I(n, x)` is positive. This
is not always the case. Most likely it is also not always the case
that `I1(n, x)` and `I3(n, x)` have the same sign. We don't even have
an enclosure if we instead multiply with the enclosure of `log(1 + c *
x * t)` on the interval since the sign changes.

### Simplifying the integrand for `I1(n, x)` and `I2(n, x)`
To begin with we want to study the part of the integrand given by
```
f(n, t) = (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n - 2(1 + γ)^n * log(t)^n
```
Using the binomial theorem we can write this as
```
f(n, t) = sum(binomial(n, k) γ^k * log(t)^k * log(t - 1)^(n - k) for k = 0:n) +
    sum(binomial(n, k) γ^k * log(t)^k * log(t + 1)^(n - k) for k = 0:n) -
    2(1 + γ)^n * log(t)^n
```
By joining the sums we can write this as
```
f(n, t) = sum(binomial(n, k) * γ^k * log(t)^k * (log(t - 1)^(n - k) + log(t + 1)^(n - k)) for k = 0:n) -
    2(1 + γ)^n * log(t)^n
```

For large values of `t` the factor
```
log(t - 1)^(n - k) + log(t + 1)^(n - k)
```
behaves almost like `2log(t)^(n - k)`. We can therefore split it into
this term plus a remainder and handle them separately. If we let
```
R(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l
```
we can split the sum as
```
sum(binomial(n, k) * γ^k * log(t)^k * (log(t - 1)^(n - k) + log(t + 1)^(n - k)) for k = 0:n) =
    sum(binomial(n, k) * γ^k * log(t)^k * 2log(t)^(n - k) for k = 0:n) +
    sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n)
```

We can now notice that the first of these sums
```
sum(binomial(n, k) * γ^k * log(t)^k * 2log(t)^(n - k) for k = 0:n)
```
is exactly equal to
```
2(1 + γ)^n * log(t)^n
```
and hence cancels the corresponding term in `f(n, t)`. This leaves us
with
```
f(n, t) = sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n)
```
Noticing that `R(0, t) = 0` we can simplify this to
```
f(n, t) = sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1)
```

Inserting this back into `I1(n, x)` and `I2(n, x)` we get
```
I1(n, x) = ∫ sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1) * t dt

I2(n, x) = ∫ sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1) * t * log(t) dt
```
where we recall that
```
R(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l
```

### Computing `I1(n, x)` and `I2(n, x)`
Switching the summation and integration in `I1(n, x)` and `I2(n, x)` we
arrive at
```
I1(n, x) = sum(binomial(n, k) * γ^k * ∫ log(t)^k * R(n - k, t) * t dt for k = 0:n-1)

I2(n, x) = sum(binomial(n, k) * γ^k * ∫ log(t)^(k + 1) * R(n - k, t) * t dt for k = 0:n-1)
```
Using that `R(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l` we can
write these integrals as
```
∫ log(t)^k * (log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)) * t dt
```
and
```
∫ log(t)^(k + 1) * (log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)) * t dt
```

#### Split integrand for `I1(n, x)` and `I2(n, x)` into main part and remainder
Asymptotically as `t` goes to infinity we have that
```
log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)
```
behaves like
```
(n - k) * log(t)^(n - k - 2) * ((n - k - 1) - log(t)) / t^2 + O(1 / t^4)
```
If we let
```
h(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l
    - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2
```
Then we have
```
log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k) =
    (n - k) * log(t)^(n - k - 2) * ((n - k - 1) - log(t)) / t^2
    + h(n - k, t)
```

Inserting the main term into the integrals for `I1(n, x)` and `I2(n, x)`
gives us
```
∫ log(t)^k * (n - k) * log(t)^(n - k - 2) * ((n - k - 1) - log(t)) / t^2 * t dt
= (n - k) * ∫ log(t)^(n - 2) * ((n - k - 1) - log(t)) / t dt
= (n - k) * (n - k - 1) * ∫ log(t)^(n - 2) / t dt -
    (n - k) * ∫ log(t)^(n - 1) / t dt

∫ log(t)^(k + 1) * (n - k) * log(t)^(n - k - 2) * ((n - k - 1) - log(t)) / t^2 * t dt
= (n - k) * ∫ log(t)^(n - 1) * ((n - k - 1) - log(t)) / t dt
= (n - k) * (n - k - 1) * ∫ log(t)^(n - 1) / t dt -
    (n - k) * ∫ log(t)^n / t dt
```
The integrals can be computed exactly to be
```
(n - k) * (n - k - 1) * ∫ log(t)^(n - 2) / t dt =
    (n - k) * (n - k - 1) / (n - 1) * log(π / x)^(n - 1)

(n - k) * ∫ log(t)^(n - 1) / t dt =
    (n - k) / n * log(π / x)^n

(n - k) * (n - k - 1) * ∫ log(t)^(n - 1) / t dt =
    (n - k) * (n - k - 1) / n * log(π / x)^n

(n - k) * ∫ log(t)^n / t dt =
    (n - k) / (n + 1) * log(π / x)^(n + 1)
```
**FIXME:** The above integrals are only valid when the power of the
logarithm is not `-1`. So the first one fails for `n = 1`, but in that
case `k = 0` and the factor in front is zero. So it doesn't actually
matter.

For the `h(n - k, t)` term inserted into the integral we let
```
H1(n, k, x) = ∫ log(t)^k * h(n - k, t) * t dt
H2(n, k, x) = ∫ log(t)^(k + 1) * h(n - k, t) * t dt
```

Inserting this into `I1(n, x)` and `I2(n, x)` gives us
```
I1(n, x) = sum(
    binomial(n, k) * γ^k * (n - k) * log(π / x)^(n - 1) * ((n - k - 1) / (n - 1) - log(π / x) / n)
    for k = 0:n-1
) + sum(binomial(n, k) * γ^k * H1(n, k, x) for k = 0:n-1)
= log(π / x)^(n - 1) * (
    1 / (n - 1) * sum(binomial(n, k) * γ^k * (n - k) * (n - k - 1)  for k = 0:n-1) -
    log(π / x) / n * sum(binomial(n, k) * γ^k * (n - k) for k = 0:n-1)
) + sum(binomial(n, k) * γ^k * H1(n, k, x) for k = 0:n-1)

I2(n, x) = sum(
    binomial(n, k) * γ^k * (n - k) * log(π / x)^n * ((n - k - 1) / n - log(π / x) / (n + 1))
    for k = 0:n-1
) + sum(binomial(n, k) * γ^k * H2(n, k, x) for k = 0:n-1)
= log(π / x)^n * (
    1 / n * sum(binomial(n, k) * γ^k * (n - k) * (n - k - 1)  for k = 0:n-1) -
    log(π / x) / (n + 1) * sum(binomial(n, k) * γ^k * (n - k) for k = 0:n-1)
) + sum(binomial(n, k) * γ^k * H2(n, k, x) for k = 0:n-1)
```

The sums in `I1(n, x)` and `I2(n, x)` without `H1` and `H2` above are
the same and they can be explicitly computed to be
```
sum(binomial(n, k) * γ^k * (n - k) * (n - k - 1)  for k = 0:n-1) =
    (1 + γ)^(n - 2) * (n - 1) * n

sum(binomial(n, k) * γ^k * (n - k) for k = 0:n-1) =
    (1 + γ)^(n - 1) * n
```
Inserting this backs gives us
```
I1(n, x) = log(π / x)^(n - 1) * (1 + γ)^(n - 2) * (n - log(π / x) * (1 + γ))
     + sum(binomial(n, k) * γ^k * H1(n, k, x) for k = 0:n-1)

I2(n, x) = log(π / x)^n * (1 + γ)^(n - 2) * ((n - 1) - n / (n + 1) * log(π / x) * (1 + γ))
     + sum(binomial(n, k) * γ^k * H2(n, k, x) for k = 0:n-1)
```

## Inserting `I(n, x)` back into the sum
Recall that we are interested in computing
```
G2(x) = -1 / ((1 - x^p0) * log(x)) * sum((-1)^n * (1 + α)^n / factorial(n) * I(n, x) for n = 1:Inf)
```
with
```
I(n, x) = (log(1 + c * π) - log(x)) * I1(n, x) - I2(n, x)
```
We can split `G2` into two sums as
```
G2(x) <= -1 / ((1 - x^p0) * log(x)) * (
    + (log(1 + c * π) - log(x)) * sum((-1)^n * (1 + α)^n / factorial(n) * I1(n, x) for n = 1:Inf)
    - sum((-1)^n * (1 + α)^n / factorial(n) * I2(n, x) for n = 1:Inf)
)
```
We are now interested in computing these two sums.

Using the above expressions of `I1(n, x)` and `I2(n, x)` we get for the
sum with `I1(n, x)`
```
sum((-1)^n * (1 + α)^n / factorial(n) * I1(n, x) for n = 1:Inf)

= sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^(n - 1) * (1 + γ)^(n - 2) * (n - log(π / x) * (1 + γ))
    for n = 1:Inf
) + sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    sum(binomial(n, k) * γ^k * H1(n, k, x) for k = 0:n-1)
    for n = 1:Inf
)
```
Focusing on the first of these two sums we can write it as
```
sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^(n - 1) * (1 + γ)^(n - 2) * (n - log(π / x) * (1 + γ))
    for n = 1:Inf
)

= sum(
    (-1)^n * (1 + α)^n / factorial(n - 1) * log(π / x)^(n - 1) * (1 + γ)^(n - 2)
    for n = 1:Inf
)
- sum(
    (-1)^n * (1 + α)^n / factorial(n) * log(π / x)^(n) * (1 + γ)^(n - 1)
    for n = 1:Inf
)

= -(1 + α) / (1 + γ) * sum(
    (-1)^(n - 1) * (1 + α)^(n - 1) / factorial(n - 1) * log(π / x)^(n - 1) * (1 + γ)^(n - 1)
    for n = 1:Inf
)
- 1 / (1 + γ) * sum(
    (-1)^n * (1 + α)^n / factorial(n) * log(π / x)^(n) * (1 + γ)^n
    for n = 1:Inf
)

= -(1 + α) / (1 + γ) * (π / x)^(-(1 + α) * (1 + γ))
- 1 / (1 + γ) * ((π / x)^(-(1 + α) * (1 + γ)) - 1)

= (1 - (2 + α) * (π / x)^(-(1 + α) * (1 + γ))) / (1 + γ)
```

For the sum with `I2(n, x)` we get
```
sum((-1)^n * (1 + α)^n / factorial(n) * I2(n, x) for n = 1:Inf)

= sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^n * (1 + γ)^(n - 2) * ((n - 1) - n / (n + 1) * log(π / x) * (1 + γ))
    for n = 1:Inf
) + sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    sum(binomial(n, k) * γ^k * H2(n, k, x) for k = 0:n-1)
    for n = 1:Inf
)
```
Focusing on the first of these two sums we can write it as
```
sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^n * (1 + γ)^(n - 2) * ((n - 1) - n / (n + 1) * log(π / x) * (1 + γ))
    for n = 1:Inf
)

= sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^n * (1 + γ)^(n - 2) * (n - 1)
    for n = 1:Inf
)
- sum(
    (-1)^n * (1 + α)^n / factorial(n) *
    log(π / x)^n * (1 + γ)^(n - 2) * n / (n + 1) * log(π / x) * (1 + γ)
    for n = 1:Inf
)

= 1 / (1 + γ)^2 * sum(
    (-1)^n * (1 + α)^n / factorial(n) * log(π / x)^n * (1 + γ)^n * (n - 1) for n = 1:Inf
)
+ 1 / ((1 + α) * (1 + γ)^2) * sum(
    (-1)^(n + 1) * (1 + α)^(n + 1) / factorial(n + 1) * log(π / x)^(n + 1) * (1 + γ)^(n + 1) * n
    for n = 1:Inf
)

= 1 / (1 + γ)^2 * (1 - (π / x)^(-(1 + α) * (1 + γ)) * ((1 + α) * (1 + γ) * log(π / x) + 1))
+ 1 / ((1 + α) * (1 + γ)^2) * (1 - (π / x)^(-(1 + α) * (1 + γ)) * ((1 + α) * (1 + γ) * log(π / x) + 1))

= (2 + α) / ((1 + α) * (1 + γ)^2) *
    (1 - (π / x)^(-(1 + α) * (1 + γ)) * ((1 + α) * (1 + γ) * log(π / x) + 1))

= ((2 + α) * (1 - (π / x)^(-(1 + α) * (1 + γ)))) / ((1 + α) * (1 + γ)^2)
- (2 + α) / (1 + γ) * (π / x)^(-(1 + α) * (1 + γ)) * log(π / x)
```

If we now let `s = (π / x)^(-(1 + α) * (1 + γ))` and use the
above we can write `G2` as
```
G2(x) <= -1 / ((1 - x^p0) * log(x)) * (
    (log(1 + c * π) - log(x)) * (1 - (2 + α) * s) / (1 + γ)
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ)^2)
    + (2 + α) / (1 + γ) * s * (1 + α) * (1 + γ) * log(π / x)
    + (log(1 + c * π) - log(x)) * sum(
    (-1)^n * (1 + α)^n / factorial(n)
    * sum(binomial(n, k) * γ^k * H1(n, k, x) for k = 0:n-1)
    for n = 1:Inf
    )
    - sum(
    (-1)^n * (1 + α)^n / factorial(n)
    * sum(binomial(n, k) * γ^k * H2(n, k, x) for k = 0:n-1)
    for n = 1:Inf
    )
)
```
We can split it into one part with the explicit terms and one part
with the sums as
```
G21(x) = -1 / ((1 - x^p0) * log(x)) * (
    (log(1 + c * π) - log(x)) * (1 - (2 + α) * s) / (1 + γ)
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ)^2)
    + (2 + α) / (1 + γ) * s * (1 + α) * (1 + γ) * log(π / x)
)

G22(x) = -1 / ((1 - x^p0) * log(x)) * (
    (log(1 + c * π) - log(x)) * sum(
    (-1)^n * (1 + α)^n / factorial(n)
    * sum(binomial(n, k) * γ^k * H1(n, k, x) for k = 0:n-1)
    for n = 1:Inf
    )
    - sum(
    (-1)^n * (1 + α)^n / factorial(n)
    * sum(binomial(n, k) * γ^k * H2(n, k, x) for k = 0:n-1)
    for n = 1:Inf
    )
)
```
with `G2(x) <= G21(x) + G22(x)`.

## Bounding `G21(x)`
We now compute a bound `G21(x)`.

### Simplifying `G21(x)`
Writing `log(π / x) = log(π) - log(x)` and putting the log-terms
together, as well as factoring out `1 / (1 + γ)`, we get
```
G21(x) = -1 / ((1 - x^p0) * log(x) * (1 + γ)) * (
    (log(1 + c * π) - log(x)) * (1 - (2 + α) * s)
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ))
    + (2 + α) * (1 + α) * (1 + γ) * s * (log(π) - log(x))
)

= -1 / ((1 - x^p0) * log(x) * (1 + γ)) * (
    log(1 + c * π) * (1 - (2 + α) * s)
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ))
    + (2 + α) * (1 + α) * (1 + γ) * s * log(π)
    - log(x) * (1 - (2 + α) * s)
    - log(x) * (2 + α) * (1 + α) * (1 + γ) * s
)

= -1 / ((1 - x^p0) * log(x) * (1 + γ)) * (
    log(1 + c * π) * (1 - (2 + α) * s)
    - ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ))
    + (2 + α) * (1 + α) * (1 + γ) * s * log(π)
    - log(x) * (1 + (2 + α) * s * ((1 + α) * (1 + γ) - 1))
)
```

Next we cancel the `log(x)` explicitly and reorder the signs a bit
```
G21(x) = 1 / ((1 - x^p0) * (1 + γ)) * (
    - log(1 + c * π) * (1 - (2 + α) * s) / log(x)
    + ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ) * log(x))
    - (2 + α) * (1 + α) * (1 + γ) * s * log(π) / log(x)
    + (1 + (2 + α) * s * ((1 + α) * (1 + γ) - 1))
)
```
Factoring out `(2 + α) / log(x)` from the first three terms gives
```
G21(x) = 1 / ((1 - x^p0) * (1 + γ)) * (
    (2 + α) / log(x) * (
        - log(1 + c * π) * (1 / (2 + α) - s)
        + (1 - s) / ((1 + α) * (1 + γ))
        - (1 + α) * (1 + γ) * s * log(π)
    )
    + (1 + (2 + α) * s * ((1 + α) * (1 + γ) - 1))
)
```

If we let `q0 = (1 + α) * (1 + γ)` we have `s = (x / π)^q0` and we can
write this as
```
G21(x) = 1 / ((1 - x^p0) * (1 + γ)) * (
    (2 + α) / log(x) * (
        - log(1 + c * π) * (1 / (2 + α) - (x / π)^q0)
        + (1 - (x / π)^q0) / q0
        - q0 * (x / π)^q0 * log(π)
    )
    + (1 + (2 + α) * (x / π)^q0 * (q0 - 1))
)
```

### Uniform bound for `G21(x)`
We now show that `G21(x)` is bounded by a factor not depending on `x`.
More precisely we want to show that
```
G21(x) <= 1 / (1 + γ) + ϵ
```
For some `ϵ` yet to be determined.

Using the expression for `G21(x)` from the above section we can write
the inequality as
```
1 / ((1 - x^p0) * (1 + γ)) * (
    (2 + α) / log(x) * (
        - log(1 + c * π) * (1 / (2 + α) - (x / π)^q0)
        + (1 - (x / π)^q0) / q0
        - q0 * (x / π)^q0 * log(π)
    )
    + (1 + (2 + α) * (x / π)^q0 * (q0 - 1))
) <= 1 / (1 + γ) + ϵ
```

Multiplying with `(1 + γ) * (1 - x^p0) * log(x)` and moving everything
to one side we get the inequality
```
(2 + α) * (
    - log(1 + c * π) * (1 / (2 + α) - (x / π)^q0)
    + (1 - (x / π)^q0) / q0
    - q0 * (x / π)^q0 * log(π)
) + (1 + (2 + α) * (x / π)^q0 * (q0 - 1)) * log(x)
- (1 + ϵ) * (1 - x^p0) * log(x) >= 0
```
Which we can rewrite slightly as
```
(2 + α) * (
    - log(1 + c * π) * (1 / (2 + α) - x^q0 / π^q0)
    + (1 - x^q0 / π^q0) / q0
    - q0 / π^q0 * x^q0 * log(π)
) + (2 + α) / π^q0 * (q0 - 1) * x^q0 * log(x)
+ (1 + ϵ) * x^p0 * log(x) - ϵ * log(x)
>= 0
```
Our goal is to show that the left hand side, `lhs(x)` is decreasing in
`x`, so that it is enough to check the inequality for an upper bound
of `x`.

Differentiating the left hand side we get
```
lhs'(x) = (2 + α) * (
    log(1 + c * π) * q0 * x^(q0 - 1) / π^q0
    - x^(q0 - 1) / π^q0
    - q0^2 / π^q0 * x^(q0 - 1) * log(π)
)
+ (2 + α) / π^q0 * (q0 - 1) * x^(q0 - 1) * (q0 * log(x) + 1)
+ (1 + ϵ) * x^(p0 - 1) * (p0 * log(x) + 1)
- ϵ / x
```
Which we can rewrite as
```
lhs'(x) = (2 + α) * x^(q0 - 1) / π^q0 * (
        log(1 + c * π) * q0 - 1 - q0^2 * log(π)
        + (q0 - 1) * (q0 * log(x) + 1)
    )
    + (1 + ϵ) * x^(p0 - 1) * (p0 * log(x) + 1)
    - ϵ / x
```
and furthermore as
```
lhs'(x) = (2 + α) * x^(q0 - 1) / π^q0 * (
        -2 + q0 * (1 + log((1 + c * π) / x)) - q0^2 * log(π / x)
    )
    + (1 + ϵ) * x^(p0 - 1) * (p0 * log(x) + 1)
    - ϵ / x
```

Since we want to prove that `lhs(x)` is decreasing we have to prove
that `lhs'(x)` is negative. Since the term `-ϵ / x` is always negative
it is enough to prove that the remaining part also is. We are hence
interested in looking at
```
(2 + α) * x^(q0 - 1) / π^q0 * (
    -2 + q0 * (1 + log((1 + c * π) / x)) - q0^2 * log(π / x)
)
+ (1 + ϵ) * x^(p0 - 1) * (p0 * log(x) + 1)
```
Which we can rewrite as
```
x^(p0 - 1) * (
    (2 + α) * x^(q0 - p0) / π^q0 * (
        -2 + q0 * (1 + log((1 + c * π) / x)) - q0^2 * log(π / x)
    )
    + (1 + ϵ) * (p0 * log(x) + 1)
)
```
Since `x^(p0 - 1)` is positive it is enough to check that the other
factor is negative, we let
```
g(x) = (2 + α) * x^(q0 - p0) / π^q0 * (
        -2 + q0 * (1 + log((1 + c * π) / x)) - q0^2 * log(π / x)
    )
    + (1 + ϵ) * (p0 * log(x) + 1)
```
denote the other factor.

**FIXME:** The below is not correct. It is not increasing in `x` so we
  can't just check the value at `x = 0`. Our goal is to prove that
  `g(x)` is negative.

As a first step we check that `g(1) < 0`, we have
```
g(1) = 1 + ϵ + (2 + α) / π^q0 * (
        -2 + q0 * (1 + log(1 + c * π)) - q0^2 * log(π)
    )
```
Which can be checked to be negative by a simple evaluation.

Next we check that `g(x)` is increasing in `x` so that `g(1)` is an
upper bound. We have
```
g'(x) = (2 + α) * (q0 - p0) * x^(q0 - p0 - 1) / π^q0 * (
        -2 + q0 * (1 + log((1 + c * π) / x)) - q0^2 * log(π / x)
    )
    + (2 + α) * x^(q0 - p0) / π^q0 * (-q0 / x + q0^2 / x)
    + (1 + ϵ) * p0 / x

    = (2 + α) * x^(q0 - p0 - 1) / π^q0 * (
        (q0 - p0) * (
            -2 + q0 * (1 + log((1 + c * π) / x)) - q0^2 * log(π / x)
        )
        -q0 + q0^2
    )
    + (1 + ϵ) * p0 / x
```
We want to prove that this is positive

**FIXME:** Recall that the above up until the other FIXME is wrong!

## Bounding `G22(x)`
We now compute a bound for `G22(x)`. Recall that it is given by
```
G22(x) = -1 / ((1 - x^p0) * log(x)) * (
    (log(1 + c * π) - log(x)) * sum(
    (-1)^n * (1 + α)^n / factorial(n)
    * sum(binomial(n, k) * γ^k * H1(n, k, x) for k = 0:n-1)
    for n = 1:Inf
    )
    - sum(
    (-1)^n * (1 + α)^n / factorial(n)
    * sum(binomial(n, k) * γ^k * H2(n, k, x) for k = 0:n-1)
    for n = 1:Inf
    )
)
```
with
```
H1(n, k, x) = ∫ log(t)^k * h(n - k, t) * t dt
H2(n, k, x) = ∫ log(t)^(k + 1) * h(n - k, t) * t dt
```
and
```
h(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l
    - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2
```

Factoring out `1 + α` from the two sums we can rewrite this as
```
G22(x) = (1 + α) / (1 - x^p0) * (
    (1 - log(1 + c * π) / log(x)) * S221(x)
    + 1 / log(x) * S222(x)
)
```
with
```
S221(x) = sum(
    (-1)^n * (1 + α)^(n - 1) / factorial(n)
    * sum(binomial(n, k) * γ^k * H1(n, k, x) for k = 0:n-1)
    for n = 1:Inf
)

S222(x) = sum(
    (-1)^n * (1 + α)^(n - 1) / factorial(n)
    * sum(binomial(n, k) * γ^k * H2(n, k, x) for k = 0:n-1)
    for n = 1:Inf
)
```

In particular, since we are only interested in an upper bound we can
consider
```
G22(x) <= (1 + α) / (1 - x^p0) * (
    (1 - log(1 + c * π) / log(x)) * abs(S221(x))
    - 1 / log(x) * abs(S222(x))
)
```
and
```
abs(S221(x)) <= sum(
    (1 + α)^(n - 1) / factorial(n)
    * sum(binomial(n, k) * γ^k * abs(H1(n, k, x)) for k = 0:n-1)
    for n = 1:Inf
)

abs(S222(x)) <= sum(
    (1 + α)^(n - 1) / factorial(n)
    * sum(binomial(n, k) * γ^k * abs(H2(n, k, x)) for k = 0:n-1)
    for n = 1:Inf
)
```
Our goal is to compute bounds for `abs(S221(x))` and `abs(S222(x))`.

### Bounding `abs(S221(x))` and `abs(S222(x))`
As a first step we are interested in bounding `abs(H1(n, k, x))` and
`abs(H2(n, k, x))`. Recall that they are given by
```
H1(n, k, x) = ∫ log(t)^k * h(n - k, t) * t dt
H2(n, k, x) = ∫ log(t)^(k + 1) * h(n - k, t) * t dt
```
integrated from `1` to `π / x`, where
```
h(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l
    - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2
```
An upper bound is given by
```
abs(H1(n, k, x)) = ∫ log(t)^k * abs(h(n - k, t)) * t dt
abs(H2(n, k, x)) = ∫ log(t)^(k + 1) * abs(h(n - k, t)) * t dt
```

Since `h(l, t)` has a singularity at `t = 1` we split the integration
at `t = 2` to avoid having to deal with both the singularity and
integrating to infinity. We then give upper bounds for `abs(h(l, t))`
on the interval `1 < t < 2` and `2 < t < π / x` separately.

#### Computing `H1(1, 0, x)` and `H2(1, 0, x)`
Since `H1(1, 0, x)` and `H2(1, 0, x)` are the main contributors to
`S221` and `S222` it makes sense to try and compute slightly more
accurate bounds for them. We have
```
H1(1, 0, x) = ∫ h(1, t) * t dt
H2(1, 0, x) = ∫ log(t) * h(1, t) * t dt
```
with
```
h(1, t) = log(t - 1) + log(t + 1) - 2log(t) + 1 / t^2
```

Note that `h(1, t)` is negative for `t > 1`. To see this write it as
```
h(1, t) = log(t - 1) + log(t + 1) - 2log(t) + 1 / t^2
    = (log(t) + log(1 - 1 / t)) + (log(t) + log(1 + 1 / t)) - 2log(t) + 1 / t^2
    = log(1 - 1 / t) + log(1 + 1 / t) + 1 / t^2
    = (-1 / t - 1 / 2t^2 - 1 / 3t^3 - 1 / 4t^4 - ...)
      + (1 / t - 1 / 2t^2 + 1 / 3t^3 - 1 / 4t^4 - ...)
      - 1 / t^2
    = -1 / 2t^4 - ...
```
Here all terms in the expansion that don't cancel are negative so the
result is clearly negative.

Since `h(1, t)` is negative an upper bound of `abs(H1(1, 0, x))` is
given by integrating from `1` to `∞`. Which can be done explicitly to
give us
```
abs(H1(1, 0, x)) = -∫ log(t - 1) * t + log(t + 1) * t - 2log(t) * t + 1 / t dt
    <= 1 // 2
```

For `H2` we get, again integrating from `1` to `∞` to get an upper
bound,
```
abs(H2(1, 0, x)) = -∫ log(t) * (log(t - 1) * t + log(t + 1) * t - 2log(t) * t + 1 / t) dt
    <= (6 - π^2) / 24
```

**TODO:** Both above integrals were computed using Mathematica. Might
  want to add some more details about it.

#### Bounding `abs(h(l, t))` for `1 < t < 2`
The only unbounded term is `log(t - 1)^l`, which we keep as it is.
What remains is
```
log(t + 1)^l - 2log(t)^l - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2

= log(t + 1)^l - 2log(t)^l - l * (l - 1) * log(t)^(l - 2) / t^2 + l * log(t)^(l - 1) / t^2
```
For which a bound for the absolute value is, using that `l >= 0`,
```
log(3)^l + 2log(2)^l + l * (l - 1) * log(2)^(l - 2) + l * log(2)^(l - 1)
```

Since `log(2) < 1` and `log(3) < 2` this is further bounded by
```
2^l + 2 + l * (l - 1) + l = 2^l + l^2 + 2
```

Summarizing we have the bound
```
abs(h(l, t)) <= abs(log(t - 1)^l) + 2^l + l^2 + 2
```
valid for `1 < t < 2`.

#### Bounding `H1(n, k, x)` and `H2(n, k, x)` integrated from `1` to `2`
For `H1(n, k, x)` we have
```
∫_1^2 log(t)^k * abs(h(n - k, t)) * t dt <=
    2log(2)^k * ∫_1^2 abs(h(n - k, t)) dt <=
    2log(2)^k * (2^(n - k) + (n - k)^2 + 2 + ∫_1^2 abs(log(t - 1)^(n - k)) dt) =
    2log(2)^k * (2^(n - k) + (n - k)^2 + 2 + factorial(n - k))
```
and for `H2(n, k, x)` we similarly get
```
∫_1^2 log(t)^(k + 1) * abs(h(n - k, t)) * t dt <=
    2log(2)^(k + 1) * (2^(n - k) + (n - k)^2 + 2 + factorial(n - k))
```

**IMPROVE:** We can get better bounds by keeping the `log(3)` and
`log(2)` from `h(l, t)`. We'll see if this is needed or not.

#### Bounding `abs(h(l, t))` for `2 < t`
We are interested in bounding
```
h(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l
    - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2
```
Focusing on `log(t - 1)^l + log(t + 1)^l - 2log(t)^l` and using that
```
log(t - 1) = log(t) + log(1 - 1 / t)
log(t + 1) = log(t) + log(1 + 1 / t)
```
we can rewrite it as
```
log(t - 1)^l + log(t + 1)^l - 2log(t)^l
= (log(t) + log(1 - 1 / t))^l + (log(t) + log(1 + 1 / t))^l - log(t)^l
= sum(0:l-1) do j
    binomial(l, j) * log(t)^j * (log(1 - 1 / t)^(l - j) + log(1 + 1 / t)^(l - j))
  end
```

Expanding the logarithms we have
```
log(1 - 1 / t) = -1 / t - 1 / 2t^2 - 1 / 3t^3 + C1(t) / t^4
log(1 + 1 / t) = 1 / t - 1 / 2t^2 + 1 / 3t^3 + C2(t) / t^4
```
where both `C1(t)` and `C2(t)` are bounded for `2 > t`.
**TODO:** Give bounds for `C1(t)` and `C2(t)`

**IN PROGRESS**

Using the above we can write
```
log(1 - 1 / t)^(l - j) + log(1 + 1 / t)^(l - j) =
```
We have that `c(l, 1) = c(l, 3) = 0` and `c(l, 2) = `


**NEXT:** Insert this into the above sum. Use the multinomial theorem
  to expand. Reorder the sums to have `1 / t^k` as the common factor.
  Ensure that the `1 / t` and `1 / t^3` terms both are zero and that
  the `1 / t^2` term is cancelled by the subtracted term. Give a bound
  in the form of `2^l * log(t)^l / t^4` or similar.

#### Bounding `H1(n, k, x)` and `H2(n, k, x)` integrated from `2` to `π / x`
**TODO**



## Fixed `x`
**FIXME: Not update**
To begin with we consider the case when `x` is fixed and non-zero,
then the only asymptotics we have to do is in `α`.

By multiplying and dividing by `1 + α` we can write `G2(x)` as
```
G2(x) = (1 + α) / ((1 - x^p0) * (1 + γ)) * (
    + (2 + α) / log(x) * (
        (1 - s) / ((1 + α) * (1 + γ))
        - log(π) * s
    ) + 1
    ) / (1 + α)
```
and then split it into the two factors
```
G21(x) = (1 + α) / ((1 - x^p0) * (1 + γ))

G22(x) = ((2 + α) / log(x) * ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s) + 1) / (1 + α)
```
Our goal will be to enclose these two separately.

### Enclosing `G21(x)`
Expanding `1 - x^p0` as
```
1 - x^p0 = 1 - exp(p0 * log(x))
= -sum(p0^n * log(x)^n / factorial(n) for n = 1:Inf)
= -p0 * sum(p0^(n - 1) * log(x)^n / factorial(n) for n = 1:Inf)
```
we can write `G21(x)` as
```
G21(x) = -1 / (1 + γ) * (1 + α) / p0 * 1 / sum(p0^(n - 1) * log(x)^n / factorial(n) for n = 1:Inf)
```
Since `p0 = (1 + α) * (1 + (1 + α) / 2)` we get
```
G21(x) = -1 / ((1 + γ) * (1 + (1 + α) / 2)) * 1 / sum(p0^(n - 1) * log(x)^n / factorial(n) for n = 1:Inf)
```
and what remains is to enclose the sum. Since `log(x) < 0` the sum is
alternating. For large enough `n` it is also decreasing, more
precisely we have
```
abs(p0^(n - 1) * log(x)^n / factorial(n)) > abs(p0^n * log(x)^(n + 1) / factorial(n + 1))
⟺
1 > abs(p0 * log(x) / (n + 1))
⟺
n > p0 * abs(log(x)) - 1
```
So for every `n` greater than `p0 * abs(log(x)) - 1` the terms in the
sum are decreasing. We can thus explicitly sum the first `n` and get
an error which is bounded by the magnitude of the `n + 1`-th term. To
get a slightly better enclosure we can explicitly sum a few extra
terms after the `n`-th.

### Enclosing `G22(x)`
To be able to explicitly cancel the division by `1 + α` we rewrite
`G22(x)` as.
```
G22(x) = (
    (1 + α) / log(x) * ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s)
    + 1 / log(x) * ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s)
    + 1
) / (1 + α)

= ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s) / log(x)
+ (1 / log(x) * ((1 - s) / ((1 + α) * (1 + γ)) - log(π) * s) + 1) / (1 + α)
```
Now, using that `s = (x / π)^((1 + α) * (1 + γ))` and letting `q0 = (1
+ α) * (1 + γ)` we can write this as
```
G22(x) = ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) / log(x)
+ (1 / log(x) * ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) + 1) / (1 + α)
```
The second term we can rewrite as
```
(1 / log(x) * ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) + 1) / (1 + α)

= ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0 + log(x)) / ((1 + α) * log(x))

= ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0 + log(x / π) + log(π)) / ((1 + α) * log(x))

= ((1 + q0 * log(x / π) - (x / π)^q0) / q0 + log(π) * (1 - (x / π)^q0)) / ((1 + α) * log(x))
```
Giving us
```
G22(x) = ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) / log(x)
+ ((1 + q0 * log(x / π) - (x / π)^q0) / q0 + log(π) * (1 - (x / π)^q0)) / ((1 + α) * log(x))

= inv(log(x)) * (
    (1 - (x / π)^q0) / q0
    - log(π) * (x / π)^q0
    + (1 + q0 * log(x / π) - (x / π)^q0) / (q0 * (1 + α))
    + log(π) * (1 - (x / π)^q0) / (1 + α)
)
```
We will now bound each of these four terms separately. The term
`log(π) * (x / π)^q0` can be enclosed directly. For the other three
terms we make use of the expansion
```
(x / π)^q0 = sum(q0^n * log(x / π)^n / factorial(n) for n = 0:Inf)
```
For the three terms this gives us
```
(1 - (x / π)^q0) / q0 = -sum(q0^(n - 1) * log(x / π)^n / factorial(n) for n = 1:Inf)

(1 + q0 * log(x / π) - (x / π)^q0) / (q0 * (1 + α)) =
    -sum(q0^(n - 1) * log(x / π)^n / factorial(n) for n = 2:Inf) / (1 + α) =
    -(1 + γ) * sum(q0^(n - 2) * log(x / π)^n / factorial(n) for n = 2:Inf) =

log(π) * (1 - (x / π)^q0) / (1 + α) =
    -log(π) * sum(q0^n * log(x / π)^n / factorial(n) for n = 1:Inf) / (1 + α) =
    -(1 + γ) * log(π) * sum(q0^(n - 1) * log(x / π)^n / factorial(n) for n = 1:Inf)
```
Note that the first and last sum are identical. Similar to the sum in
`G21(x)` all of these sums are alternating and to compute enclosures
it is enough to find when the terms start to decrease. In this case we
get
```
abs(q0^(n - 1) * log(x / π)^n / factorial(n)) > abs(q0^n * log(x / π)^(n + 1) / factorial(n + 1))
⟺
1 > abs(q0 * log(x / π) / (n + 1))
⟺
n > q0 * abs(log(x / π)) - 1
```
for the first sum and a similar calculation gives the same for the
second sum. So for every `n` greater than `q0 * abs(log(x / π)) - 1`
the terms in both sums are decreasing. We can thus explicitly sum the
first `n` and get an error which is bounded by the magnitude of the `n
+ 1`-th term. To get a slightly better enclosure we can explicitly sum
a few extra terms after the `n`-th.

## Monotonicity in `x`
Above we showed how to compute an enclosure of `G2(x)` for a fixed
non-zero `x`. We will now show a certain type of monotonicity in `x`,
so that it is enough to work with an upper bound of `x`. We want to
show that
```
G2(x) < 1 / (1 + γ)
```
We can rewrite this as
```
(2 + α) * ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) + log(x) >= (1 - x^p0) * log(x)
```
where we have multiplied both sides by `(1 + γ) * (1 - x^p0) *
log(x)`, since this is negative we also reversed the inequality. If we
move everything to one side we get
```
(2 + α) * ((1 - (x / π)^q0) / q0 - log(π) * (x / π)^q0) + x^p0 * log(x) >= 0
```
We now want to show that the left hand side is decreasing in `x`, so
that it is enough to check the inequality for an upper bound of `x`.

Differentiating the left hand side we get
```
(2 + α) * (-x^(q0 - 1) / π^q0 - q0 * log(π) * x^(q0 - 1) / π^q0) + x^(p0 - 1) * (p0 * log(x) + 1)
```
Grouping terms we have
```
x^(p0 - 1) * (
    - (2 + α) * x^(q0 - p0) / π^q0
    - (2 + α) * q0 * log(π) * x^(q0 - p0) / π^q0
    + p0 * log(x)
    + 1
)
```
It is therefore enough to show that
```
(
    - (2 + α) * x^(q0 - p0) / π^q0
    - (2 + α) * q0 * log(π) * x^(q0 - p0) / π^q0
    + p0 * log(x)
    + 1
)
```
is negative. We will prove that this is increasing in `x` and
non-positive at `x = 1`. At `x = 1` we get
```
(
    - (2 + α) / π^q0
    - (2 + α) * q0 * log(π) / π^q0
    + 1
) =
1 - (2 + α) * (1 + q0 * log(π)) / π^q0
```
which we want to prove is upper bounded by `0`. To prove move one of
the terms to the other side and multiply both sides by `π^q0` and use
the expansion `π^q0 = sum(q0^n * log(π)^n / factorial(n) for n =
0:Inf)` to get
```
sum(q0^n * log(π)^n / factorial(n) for n = 0:Inf) <= (2 + α) * (1 + q0 * log(π))
⟺
sum(q0^n * log(π)^n / factorial(n) for n = 2:Inf) <= (1 + α) * (1 + q0 * log(π))
⟺
(1 + α) * (1 + γ) * log(π) * sum(q0^(n - 1) * log(π)^n / factorial(n) for n = 2:Inf) <=
    (1 + α) * (1 + q0 * log(π))
⟺
(1 + γ) * log(π) * sum(q0^(n - 1) * log(π)^n / factorial(n) for n = 2:Inf) <=
    1 + q0 * log(π)
```
It therefore suffices to show that
```
sum(q0^(n - 1) * log(π)^n / factorial(n) for n = 2:Inf) <= 1 / (1 + γ) * log(π)
```
Since the sum is increasing in `q0` we only have to show this for an
upper bound of `q0`. The sum can be computed explicitly to be
```
sum(q0^(n - 1) * log(π)^n / factorial(n) for n = 2:Inf) = (π^q0 - q0 * log(π) - 1) / q0
```
and since we now know that this is increasing in `q0` we check that
```
(π^q0 - q0 * log(π) - 1) / q0 <= 1 / (1 + γ) * log(π)
```
holds for an upper bound of `q0`.

For the derivative we get
```
(
    - (q0 - p0) * (2 + α) * x^(q0 - p0 - 1) / π^q0
    - (q0 - p0) * (2 + α) * q0 * log(π) * x^(q0 - p0 - 1) / π^q0
    + p0 / x
)
```
Multiplying by `x`, which doesn't change the sign, and reordering we
get
```
p0 - (q0 - p0) * (2 + α) * (1 + q0 * log(π)) / π^q0 * x^(q0 - p0)
```
Since `q0 - p0 > 0` and hence `0 <= x^(q0 - p0) <= 1`, it is enough to
show
```
p0 - (q0 - p0) * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
Factoring out `1 + α` and cancelling it we get
```
(1 + (1 + α) / 2) - ((1 + γ) - (1 + (1 + α) / 2)) * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
⟺
1 + (1 + α) / 2 - (γ - (1 + α) / 2) * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
Since `(1 + α) / 2` is positive we can remove it, giving us
```
1 - (γ - (1 + α) / 2) * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
Similarly we can replace `(γ - (1 + α) / 2)` with `γ` without losing
anything.
```
1 - γ * (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
Finally `γ < 1` so it is enough to show
```
1 - (2 + α) * (1 + q0 * log(π)) / π^q0 >= 0
```
which is exactly what we already had above for the `x = 1` case.
"""
function _T0_asymptotic_main_2(α::Arb, γ::Arb, c::Arb)

    f1 = _T0_asymptotic_main_2_1(α, γ, c)
    f2 = _T0_asymptotic_main_2_2(α, γ, c)

    return x::Arb -> begin
        x < 1 || throw(DomainError(x, "must have x < 1"))

        if true#Arblib.contains_zero(x)
            # Tests for integrand
            integrand1(t) =
                let α = ubound(Arb, α)
                    (((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) * t^(-γ * (1 + α))) *
                    t *
                    log(c + inv(x * t))
                end

            integrand2(t) =
                let α = ubound(Arb, α)
                    t *
                    log(c + inv(x * t)) *
                    sum(1:10) do n
                        (-1)^n * (1 + α)^n / factorial(n) * (
                            (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                            2(1 + γ)^n * log(t)^n
                        )
                    end
                end

            #@show integrand1(Arb(3)) integrand2(Arb(3))
            @assert integrand1(Arb(3)) ≈ integrand2(Arb(3))

            # Integration limits for testing
            a, b = Arb(2), Arb(100)

            integral1 = Arblib.integrate(integrand1, a, b)

            I_part(n, x) =
                Arblib.integrate(a, b) do t
                    (
                        (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                        2(1 + γ)^n * log(t)^n
                    ) *
                    t *
                    log(c + inv(x * t))
                end |> real

            integral2 = let α = ubound(Arb, α)
                sum(1:10) do n
                    (-1)^n * (1 + α)^n / factorial(n) * I_part(n, x)
                end
            end

            @show integral1 integral2
            @assert integral1 ≈ integral2

            I1_part(n, x) =
                Arblib.integrate(a, b) do t
                    (
                        (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                        2(1 + γ)^n * log(t)^n
                    ) * t
                end |> real

            I2_part(n, x) =
                Arblib.integrate(a, b) do t
                    (
                        (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                        2(1 + γ)^n * log(t)^n
                    ) *
                    t *
                    log(t)
                end |> real

            I3_part(n, x) =
                Arblib.integrate(a, b) do t
                    (
                        (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                        2(1 + γ)^n * log(t)^n
                    ) *
                    t *
                    log(1 + c * x * t)
                end |> real

            for n = 1:10
                q1 = I_part(n, x)
                q2 = -log(x) * I1_part(n, x) - I2_part(n, x) + I3_part(n, x)

                #@show n q1 q2 I1_part(n, x) I2_part(n, x) I3_part(n, x)
                @assert Arblib.overlaps(q1, q2)
            end

            for n = 1:10
                q1 = I3_part(n, x)#-log(x) * I1_part(n, x) - I2_part(n, x) + I3_part(n, x)
                q2 = I1_part(n, x)#(Arb((log(1 + c * x), log(1 + c * π))) - log(x)) * I1_part(n, x) - I2_part(n, x)
                #@show q1 q2
                #@assert Arblib.overlaps(q1, q2)
            end

            R(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l

            I1_part_v2(n, x) =
                Arblib.integrate(a, b) do t
                    t * sum(0:n-1) do k
                        binomial(n, k) * γ^k * log(t)^k * R(n - k, t)
                    end
                end |> real

            I2_part_v2(n, x) =
                Arblib.integrate(a, b) do t
                    t * log(t) * sum(0:n-1) do k
                        binomial(n, k) * γ^k * log(t)^k * R(n - k, t)
                    end
                end |> real

            for n = 1:10
                q1 = I1_part(n, x)
                q2 = I1_part_v2(n, x)

                @assert Arblib.overlaps(q1, q2)

                q1 = I2_part(n, x)
                q2 = I2_part_v2(n, x)

                @assert Arblib.overlaps(q1, q2)
            end

            I1_part_v3(n, x) =
                sum(0:n-1) do k
                    binomial(n, k) * γ^k * Arblib.integrate(a, b) do t
                        log(t)^k * R(n - k, t) * t
                    end
                end |> real

            I2_part_v3(n, x) =
                sum(0:n-1) do k
                    binomial(n, k) * γ^k * Arblib.integrate(a, b) do t
                        log(t)^k * R(n - k, t) * t * log(t)
                    end
                end |> real

            for n = 1:10
                q1 = I1_part_v2(n, x)
                q2 = I1_part_v3(n, x)

                @assert Arblib.overlaps(q1, q2)

                q1 = I2_part_v2(n, x)
                q2 = I2_part_v3(n, x)

                @assert Arblib.overlaps(q1, q2)
            end

            return zero(α)
            # New version

            # Bound G21(x)
            # TODO: Work in progress!

            G21 =
                let α = ubound(Arb, α),
                    p0 = 1 + α + (1 + α)^2 / 2,
                    s = (π / x)^(-(1 + α) * (1 + γ))

                    -1 / ((1 - x^p0) * log(x)) * (
                        (log(1 + c * π) - log(x)) * (1 - (2 + α) * s) / (1 + γ) -
                        ((2 + α) * (1 - s)) / ((1 + α) * (1 + γ)^2) +
                        (2 + α) / (1 + γ) * s * (1 + α) * (1 + γ) * log(π / x)
                    )
                end

            T21 =
                let α = ubound(Arb, α), p0 = 1 + α + (1 + α)^2 / 2, q0 = (1 + α) * (1 + γ)
                    1 / ((1 - x^p0) * (1 + γ)) * (
                        (2 + α) / log(x) * (
                            -log(1 + c * π) * (1 / (2 + α) - (x / π)^q0) +
                            (1 - (x / π)^q0) / q0 - q0 * (x / π)^q0 * log(Arb(π))
                        ) + (1 + (2 + α) * (x / π)^q0 * (q0 - 1))
                    )
                end

            @show G21 T21

            # Compute g(1)
            # FIXME: Choose value to use for ϵ
            g_1 = let q0 = (1 + α) * (1 + γ), ϵ = Arb(0.5)
                1 + ϵ + (2 + α) / π^q0 * (-2 + q0 * (1 + log(1 + c * π)) - q0^2 * log(π))
            end
            #@show g_1
            # Check that g(1) indeed is negative
            Arblib.isnegative(g_1) || error("expected g(1) to be negative")

            g(x) =
                let p0 = (1 + α) * (1 + (1 + α) / 2), q0 = (1 + α) * (1 + γ), ϵ = Arb(0.5)
                    (
                        (2 + α) * x^(q0 - p0) / π^q0 *
                        (-2 + q0 * (1 + log((1 + c * π) / x)) - q0^2 * log(π / x)) +
                        (1 + ϵ) * (p0 * log(x) + 1)
                    )
                end

            #@show g(Arb(0.1)) g(Arb(0.5)) g(Arb(0.9)) getball(Arb, g(Arb("1e-100")))

            g_diff =
                let p0 = (1 + α) * (1 + (1 + α) / 2),
                    q0 = (1 + α) * (1 + γ),
                    ϵ = Arb(0.5),
                    x = Arb(0.1)

                    (
                        (2 + α) * x^(q0 - p0) / π^q0 * (
                            (q0 - p0) *
                            (-2 + q0 * (1 + log((1 + c * π) / x)) - q0^2 * log(π / x)) -
                            q0 + q0^2
                        ) + (1 + ϵ) * p0
                    )
                end

            #@show getball(Arb, g_diff)

            # Bound G22(x)
            # TODO: Work in progress!

            # Enclosure of (1 + α) / (1 - x^p0)
            # PROVE: Add comment about where this is already proved.
            factor_G22 = if iszero(x)
                α + 1
            elseif Arblib.contains_zero(x)
                lower = α + 1
                upper = let xᵤ = ubound(Arb, x)
                    inv(fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), α + 1, extra_degree = 2))
                end
                Arb((lower, upper))
            else
                inv(fx_div_x(s -> 1 - x^(s + s^2 / 2), α + 1, extra_degree = 2))
            end

            # Enclosure of inv(log(x))
            invlogx = if iszero(x)
                zero(x)
            elseif Arblib.contains_zero(x)
                Arb((inv(log(ubound(Arb, x))), 0))
            else
                inv(log(x))
            end

            # Enclosure of 1 - log(1 + c * π) / log(x)
            factor_S221 = 1 - log(1 + c * π) * invlogx

            # Enclosure of 1 / log(x)
            factor_S222 = invlogx

            # Bound for integral from 1 to 2
            H11(n, k, x) =
                2log(Arb(2))^k * let l = n - k
                    # This is not the simplest for of the upper bound
                    # but the one which gives slightly better bounds
                    log(Arb(3))^l +
                    2log(Arb(2))^l +
                    l * (l - 1) * log(Arb(2))^(l - 2) +
                    l * log(Arb(2))^(l - 1) +
                    factorial(l)
                end

            # Bound for integral from 2 to π / x
            H12(n, k, x) = 0

            H1(n, k, x) =
                if n == 1 && k == 0
                    # Use explicitly computed upper bound
                    abs(Arb(-1 // 2))
                else
                    H11(n, k, x) + H12(n, k, x)
                end

            # Bound for integral from 1 to 2
            H21(n, k, x) =
                2log(Arb(2))^(k + 1) * let l = n - k
                    # This is not the simplest for of the upper bound
                    # but the one which gives slightly better bounds
                    log(Arb(3))^l +
                    2log(Arb(2))^l +
                    l * (l - 1) * log(Arb(2))^(l - 2) +
                    l * log(Arb(2))^(l - 1) +
                    factorial(l)
                end

            # Bound for integral from 2 to π / x
            H22(n, k, x) = 0

            H2(n, k, x) =
                if n == 1 && k == 0
                    # Use explicitly computed upper bound
                    abs((6 - Arb(π)^2) / 24)
                else
                    H21(n, k, x) + H22(n, k, x)
                end

            # Bound of abs(S221(x))
            S221 = sum(1:10) do n
                term =
                    (1 + α)^(n - 1) / factorial(n) * sum(0:n-1) do k
                        binomial(n, k) * γ^k * H1(n, k, x)
                    end
                q = inv(factorial(n)) * sum(0:n-1) do k
                    binomial(n, k) * γ^k * H1(n, k, x)
                end
                @show q
                term
            end

            # Bound of abs(S222(x))
            S222 = sum(1:10) do n
                term =
                    (1 + α)^(n - 1) / factorial(n) * sum(0:n-1) do k
                        binomial(n, k) * γ^k * H2(n, k, x)
                    end
                q = inv(factorial(n)) * sum(0:n-1) do k
                    binomial(n, k) * γ^k * H2(n, k, x)
                end
                @show q
                term
            end

            G22 = factor_G22 * (factor_S221 * S221 + factor_S222 * S222)

            @show S221 S222 G22

            G2 = G21 + G22

            T2 = f1(x) + f2(x)

            @show G2 T2

            # In this case we prove that G2 is bounded by 1 / (1 +
            # γ). By above we only have to check two things
            # 1. 1 - (2 + α) * (1 + q0 * log(π)) / π^q0 <= 0
            # 2. 0 <= p0 - (q0 - p0) * (2 + α) * (1 + q0 * log(π)) / π^q0

            # 1. This was reduced to checking that
            # (π^q0 - q0 * log(π) - 1) / q0 <= 1 / (1 + γ) * log(π)
            # holds for an upper bound of q0
            lemma1 = let q0 = ubound(Arb, (1 + α) * (1 + γ)), π = Arb(π)
                (π^q0 - q0 * log(π) - 1) / q0 <= 1 / (1 + γ) * log(π)
            end

            # 2. This was reduced to checking the same inequality
            # as for 1. and using that γ <= 1
            lemma2 = γ <= 1 && lemma1

            if lemma1 && lemma2
                return 1 / (1 + γ)
            else
                return indeterminate(α)
            end
        else
            # TODO: Only use f1(x) + f2(x) for x not too small
            return f1(x) + f2(x)

            # FIXME: Handle the remaining remainder term in the below approach
            extra_degree = 2

            # Enclosure of (1 + α) / (1 - x^p0)
            factor1 = inv(fx_div_x(ϵ -> 1 - x^(ϵ + ϵ^2 / 2), 1 + α; extra_degree))

            factor2 = let
                term1 =
                    inv(1 + γ) * fx_div_x(1 + α, 2, force = true; extra_degree) do ϵ
                        (x / π)^(ϵ * (1 + γ)) - ϵ * (1 + γ) * log(x / π) - 1
                    end

                term2 = log(Arb(π)) * fx_div_x(1 + α; extra_degree) do ϵ
                    (x / π)^(ϵ * (1 + γ)) - 1
                end

                term3 = inv(1 + γ) * fx_div_x(1 + α; extra_degree) do ϵ
                    (x / π)^(ϵ * (1 + γ)) - 1
                end

                term4 = log(Arb(π)) * (x / π)^((1 + α) * (1 + γ))

                term1 + term2 + term3 + term4
            end

            return factor1 * factor2 / (log(inv(x)) * (1 + γ))
        end
    end
end

"""
T0_asymptotic_main_2_testing(x::Arb = Arb(1e-10), α::Arb = Arb(-0.9999), γ::Arb = Arb(0.5), c::Arb = 2Arb(ℯ))

Method for testing the development of [`_T0_asymptotic_main_2`](@ref).
The background for that method is very long and it is easy to make
mistakes. This method tests a lot of the rewrites and simplifications
that is done for [`_T0_asymptotic_main_2`](@ref) to catch any
mistakes.

Instead of integrating from `1` to `π / x` it integrates from `2` to
`10`. This avoids any singularities and allows for computing most
integrals numerically.
"""
function _T0_asymptotic_main_2_testing(
    x::Arb = Arb(1e-10),
    α::Arb = Arb(-0.9999),
    γ::Arb = Arb(0.5),
    c::Arb = 2Arb(ℯ),
)
    # Integration limits
    a, b = Arb(2), Arb(10)

    # Maximum value of n to include
    N = 5

    p0 = 1 + α + (1 + α)^2

    # Enclosure of G2 when integrating from a to b
    G2_v1 =
        1 / ((1 - x^p0) * log(inv(x))) * Arblib.integrate(a, b) do t
            ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
            t^(1 - γ * (1 + α)) *
            log(c + inv(x * t))
        end |> real

    # Only integral part from G2
    F2_v1 =
        Arblib.integrate(a, b) do t
            ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
            t^(1 - γ * (1 + α)) *
            log(c + inv(x * t))
        end |> real
    @show F2_v1

    #####
    # Rewrite integrand of F2 as infinite sum and switch integration
    # and sum
    #####

    # Integrals in sum
    I(n) =
        Arblib.integrate(a, b) do t
            (
                (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                2(1 + γ)^n * log(t)^n
            ) *
            t *
            log(c + inv(x * t))
        end |> real

    F2_v2 = sum(1:N) do n
        (-1)^n * (1 + α)^n / factorial(n) * I(n)
    end

    @show F2_v2
    @assert isapprox(F2_v1, F2_v2, rtol = Arb(1e-15))

    #####
    # Split I(n) into three parts
    #####

    I1(n) =
        Arblib.integrate(a, b) do t
            (
                (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                2(1 + γ)^n * log(t)^n
            ) * t
        end |> real

    I2(n) =
        Arblib.integrate(a, b) do t
            (
                (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                2(1 + γ)^n * log(t)^n
            ) *
            t *
            log(t)
        end |> real

    I3(n) =
        Arblib.integrate(a, b) do t
            (
                (log(t - 1) + γ * log(t))^n + (log(t + 1) + γ * log(t))^n -
                2(1 + γ)^n * log(t)^n
            ) *
            t *
            log(1 + c * x * t)
        end |> real

    for n = 1:N
        r1 = I(n)
        r2 = -log(x) * I1(n) - I2(n) + I3(n)
        @assert isfinite(r2)
        @assert Arblib.overlaps(r1, r2)
    end

    #####
    # Simplify integrands for I1, I2 and I3
    #####

    R(l, t) = log(t - 1)^l + log(t + 1)^l - 2log(t)^l

    I1_v2(n) =
        Arblib.integrate(a, b) do t
            sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1) * t
        end |> real

    I2_v2(n) =
        Arblib.integrate(a, b) do t
            sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1) * t * log(t)
        end |> real

    I3_v2(n) =
        Arblib.integrate(a, b) do t
            sum(binomial(n, k) * γ^k * log(t)^k * R(n - k, t) for k = 0:n-1) *
            t *
            log(1 + c * x * t)
        end |> real

    for n = 1:N
        r11, r12, r13 = I1(n), I2(n), I3(n)
        r21, r22, r23 = I1_v2(n), I2_v2(n), I3_v2(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert isfinite(r23)
        @assert Arblib.overlaps(r11, r21)
        @assert Arblib.overlaps(r12, r22)
        @assert Arblib.overlaps(r13, r23)
    end

    #####
    # Switch summation and integration for I1, I2 and I3
    #####

    I1_v3(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            Arblib.integrate(a, b) do t
                log(t)^k * (log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)) * t
            end
        end |> real

    I2_v3(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            Arblib.integrate(a, b) do t
                log(t)^(k + 1) *
                (log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)) *
                t
            end
        end |> real

    I3_v3(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            Arblib.integrate(a, b) do t
                log(t)^k *
                (log(t - 1)^(n - k) + log(t + 1)^(n - k) - 2log(t)^(n - k)) *
                t *
                log(1 + c * x * t)
            end
        end |> real

    for n = 1:N
        r11, r12, r13 = I1_v2(n), I2_v2(n), I3_v2(n)
        r21, r22, r23 = I1_v3(n), I2_v3(n), I3_v3(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert isfinite(r23)
        @assert Arblib.overlaps(r11, r21)
        @assert Arblib.overlaps(r12, r22)
        @assert Arblib.overlaps(r13, r23)
    end

    #####
    # Split integrands for I1, I2 and I3 into main term and remainder
    #####

    h(l, t) =
        log(t - 1)^l + log(t + 1)^l - 2log(t)^l -
        l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2

    I1_v4_main(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            (n - k) *
            ((n - k - 1) * Arblib.integrate(a, b) do t
                log(t)^(n - 2) / t
            end - Arblib.integrate(a, b) do t
                log(t)^(n - 1) / t
            end)

        end |> real

    H1(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * h(n - k, t) * t
    end |> real

    I1_v4_remainder(n) = sum(0:n-1) do k
        binomial(n, k) * γ^k * H1(n, k)
    end |> real

    I2_v4_main(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            (n - k) *
            ((n - k - 1) * Arblib.integrate(a, b) do t
                log(t)^(n - 1) / t
            end - Arblib.integrate(a, b) do t
                log(t)^n / t
            end)
        end |> real

    H2(n, k) = Arblib.integrate(a, b) do t
        log(t)^(k + 1) * h(n - k, t) * t
    end |> real

    I2_v4_remainder(n) = sum(0:n-1) do k
        binomial(n, k) * γ^k * H2(n, k)
    end |> real

    I3_v4_main(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            (n - k) *
            (
                (n - k - 1) * Arblib.integrate(a, b) do t
                    log(t)^(n - 2) / t * log(1 + c * x * t)
                end - Arblib.integrate(a, b) do t
                    log(t)^(n - 1) / t * log(1 + c * x * t)
                end
            )
        end |> real

    H3(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * h(n - k, t) * t * log(1 + c * x * t)
    end |> real

    I3_v4_remainder(n) = sum(0:n-1) do k
        binomial(n, k) * γ^k * H3(n, k)
    end |> real

    I1_v4(n) = I1_v4_main(n) + I1_v4_remainder(n)
    I2_v4(n) = I2_v4_main(n) + I2_v4_remainder(n)
    I3_v4(n) = I3_v4_main(n) + I3_v4_remainder(n)

    for n = 1:N
        r11, r12, r13 = I1_v3(n), I2_v3(n), I3_v3(n)
        r21, r22, r23 = I1_v4(n), I2_v4(n), I3_v4(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert isfinite(r23)
        @assert Arblib.overlaps(r11, r21)
        @assert Arblib.overlaps(r12, r22)
        @assert Arblib.overlaps(r13, r23)
    end

    #####
    # Primitive functions for integrals in main terms for I1 and I2
    #####

    # Primitive function of log(t)^l / t
    primitive(l, t) =
        if l == -1
            log(log(t)) # This case does occur but is always multiplied by zero
        else
            log(t)^(l + 1) / (l + 1)
        end

    I1_v5_main(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            (n - k) *
            (
                (n - k - 1) * (primitive(n - 2, b) - primitive(n - 2, a)) -
                (primitive(n - 1, b) - primitive(n - 1, a))
            )
        end

    I2_v5_main(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            (n - k) *
            (
                (n - k - 1) * (primitive(n - 1, b) - primitive(n - 1, a)) -
                (primitive(n, b) - primitive(n, a))
            )
        end

    for n = 1:N
        r11, r12 = I1_v4_main(n), I2_v4_main(n)
        r21, r22 = I1_v5_main(n), I2_v5_main(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert Arblib.overlaps(r11, r21)
        @assert Arblib.overlaps(r12, r22)
    end

    #####
    # Simplify sums for main terms for I1 and I2
    #####

    I1_v6_main(n) =
        (primitive(n - 2, b) - primitive(n - 2, a)) * sum(0:n-1) do k
            binomial(n, k) * γ^k * (n - k) * (n - k - 1)
        end -
        (primitive(n - 1, b) - primitive(n - 1, a)) * sum(0:n-1) do k
            binomial(n, k) * γ^k * (n - k)
        end

    I2_v6_main(n) =
        (primitive(n - 1, b) - primitive(n - 1, a)) * sum(0:n-1) do k
            binomial(n, k) * γ^k * (n - k) * (n - k - 1)
        end - (primitive(n, b) - primitive(n, a)) * sum(0:n-1) do k
            binomial(n, k) * γ^k * (n - k)
        end


    for n = 1:N
        r11, r12 = I1_v5_main(n), I2_v5_main(n)
        r21, r22 = I1_v6_main(n), I2_v6_main(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert Arblib.overlaps(r11, r21)
        @assert Arblib.overlaps(r12, r22)
    end

    #####
    # Explicitly compute sums in main terms for I1 and I2
    #####

    # Note that if a = 1 then primitive(l, a) = 0 and the expressions
    # simplify. Since we don't have a = 1 we can't make this
    # simplification here unfortunately.

    I1_v7_main(n) =
        (primitive(n - 2, b) - primitive(n - 2, a)) * (1 + γ)^(n - 2) * (n - 1) * n -
        (primitive(n - 1, b) - primitive(n - 1, a)) * (1 + γ)^(n - 1) * n

    I2_v7_main(n) =
        (primitive(n - 1, b) - primitive(n - 1, a)) * (1 + γ)^(n - 2) * (n - 1) * n -
        (primitive(n, b) - primitive(n, a)) * (1 + γ)^(n - 1) * n

    for n = 1:N
        r11, r12 = I1_v6_main(n), I2_v6_main(n)
        r21, r22 = I1_v7_main(n), I2_v7_main(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert Arblib.overlaps(r11, r21)
        @assert Arblib.overlaps(r12, r22)
    end

    #####
    # Insert back to get one expression for I1 and I2
    #####

    I1_v8(n) =
        (1 + γ)^(n - 2) *
        n *
        (
            (primitive(n - 2, b) - primitive(n - 2, a)) * (n - 1) -
            (primitive(n - 1, b) - primitive(n - 1, a)) * (1 + γ)
        ) + sum(binomial(n, k) * γ^k * H1(n, k) for k = 0:n-1)

    I2_v8(n) =
        (1 + γ)^(n - 2) *
        n *
        (
            (primitive(n - 1, b) - primitive(n - 1, a)) * (n - 1) -
            (primitive(n, b) - primitive(n, a)) * (1 + γ)
        ) + sum(binomial(n, k) * γ^k * H2(n, k) for k = 0:n-1)

    for n = 1:N
        r11, r12 = I1(n), I2(n)
        r21, r22 = I1_v8(n), I2_v8(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert Arblib.overlaps(r11, r21)
        @assert Arblib.overlaps(r12, r22)
    end

    #####
    # Insert back into sum for F2 and split it into two sums
    #####

    I1_v9_main(n) =
        (1 + γ)^(n - 2) *
        n *
        (
            (primitive(n - 2, b) - primitive(n - 2, a)) * (n - 1) -
            (primitive(n - 1, b) - primitive(n - 1, a)) * (1 + γ)
        )

    I1_v9_remainder(n) = sum(binomial(n, k) * γ^k * H1(n, k) for k = 0:n-1)

    I2_v9_main(n) =
        (1 + γ)^(n - 2) *
        n *
        (
            (primitive(n - 1, b) - primitive(n - 1, a)) * (n - 1) -
            (primitive(n, b) - primitive(n, a)) * (1 + γ)
        )

    I2_v9_remainder(n) = sum(binomial(n, k) * γ^k * H2(n, k) for k = 0:n-1)

    # Note that this has not been simplified very much
    I3_v9_main(n) =
        sum(0:n-1) do k
            binomial(n, k) *
            γ^k *
            (n - k) *
            (
                (n - k - 1) * Arblib.integrate(a, b) do t
                    log(t)^(n - 2) / t * log(1 + c * x * t)
                end - Arblib.integrate(a, b) do t
                    log(t)^(n - 1) / t * log(1 + c * x * t)
                end
            )
        end |> real

    I3_v9_remainder(n) = sum(binomial(n, k) * γ^k * H3(n, k) for k = 0:n-1)

    F21_v1 = sum(1:N) do n
        (-1)^n * (1 + α)^n / factorial(n) *
        (-log(x) * I1_v9_main(n) - I2_v9_main(n) + I3_v9_main(n))
    end

    F22_v1 = sum(1:N) do n
        (-1)^n * (1 + α)^n / factorial(n) *
        (-log(x) * I1_v9_remainder(n) - I2_v9_remainder(n) + I3_v9_remainder(n))
    end

    F2_v3 = F21_v1 + F22_v1

    @show F21_v1
    @show F2_v3
    @assert isapprox(F2_v1, F2_v3, rtol = Arb(1e-15))
    @assert Arblib.overlaps(F2_v2, F2_v3)

    #####
    # Split F21 into several sums
    #####

    F21_v2_1 = -log(x) * sum(1:N) do n
        (-1)^n * (1 + α)^n / factorial(n) * I1_v9_main(n)
    end

    F21_v2_2 = -sum(1:N) do n
        (-1)^n * (1 + α)^n / factorial(n) * I2_v9_main(n)
    end

    F21_v2_3 = sum(1:N) do n
        (-1)^n * (1 + α)^n / factorial(n) * I3_v9_main(n)
    end

    F21_v2 = F21_v2_1 + F21_v2_2 + F21_v2_3

    @show F21_v2
    @assert isfinite(F21_v2)
    @assert Arblib.overlaps(F21_v1, F21_v2)

    #####
    # Splits sums for F21 further
    #####

    F21_v3_1 =
        -log(x) * (
            sum(1:N) do n
                (-1)^n * (1 + α)^n / factorial(n) *
                (1 + γ)^(n - 2) *
                n *
                (n - 1) *
                (primitive(n - 2, b) - primitive(n - 2, a))
            end - sum(1:N) do n
                (-1)^n * (1 + α)^n / factorial(n) *
                (1 + γ)^(n - 1) *
                n *
                (primitive(n - 1, b) - primitive(n - 1, a))
            end
        )

    F21_v3_2 =
        -sum(1:N) do n
            (-1)^n * (1 + α)^n / factorial(n) *
            (1 + γ)^(n - 2) *
            n *
            (n - 1) *
            (primitive(n - 1, b) - primitive(n - 1, a))
        end + sum(1:N) do n
            (-1)^n * (1 + α)^n / factorial(n) *
            (1 + γ)^(n - 1) *
            n *
            (primitive(n, b) - primitive(n, a))
        end

    # TODO: We can't simplify F21_v2_3 further for now

    @assert isfinite(F21_v3_1)
    @assert isfinite(F21_v3_2)
    @assert Arblib.overlaps(F21_v2_1, F21_v3_1)
    @assert Arblib.overlaps(F21_v2_2, F21_v3_2)

    #####
    # Simplify the split sums by inserting the primitive function
    #####

    # Note that primitive(l, t) = log(t)^(l + 1) / (l + 1) unless l =
    # -1. But we only have l = 1 in one case which is zero anyway.

    F21_v4_1 =
        -log(x) * (
            -(1 + α) / (1 + γ) * sum(2:N) do n
                (-1)^(n - 1) * (1 + α)^(n - 1) / factorial(n - 1) *
                (1 + γ)^(n - 1) *
                (log(b)^(n - 1) - log(a)^(n - 1))
            end -
            inv(1 + γ) * sum(1:N) do n
                (-1)^n * (1 + α)^n / factorial(n) * (1 + γ)^n * (log(b)^n - log(a)^n)
            end
        )

    F21_v4_2 =
        -inv(1 + γ)^2 * sum(1:N) do n
            (-1)^n * (1 + α)^n / factorial(n) *
            (1 + γ)^(n) *
            (n - 1) *
            (log(b)^n - log(a)^n)
        end +
        inv((1 + α) * (1 + γ)^2) * sum(1:N) do n
            (-1)^n * (1 + α)^(n + 1) / factorial(n + 1) *
            (1 + γ)^(n + 1) *
            n *
            (log(b)^(n + 1) - log(a)^(n + 1))
        end

    @assert isfinite(F21_v4_1)
    @assert isfinite(F21_v4_2)
    @assert Arblib.overlaps(F21_v3_1, F21_v4_1)
    @assert Arblib.overlaps(F21_v3_2, F21_v4_2)

    #####
    # Compute the sums explicitly.
    #####

    F21_v5_1 =
        -log(x) * (
            (1 + α) / (1 + γ) * (a^(-(1 + α) * (1 + γ)) - b^(-(1 + α) * (1 + γ))) -
            inv(1 + γ) * (b^(-(1 + α) * (1 + γ)) - a^(-(1 + α) * (1 + γ)))
        )

    F21_v5_2 =
        -inv(1 + γ)^2 * (
            (1 + (1 + α) * (1 + γ) * log(a)) * a^(-(1 + α) * (1 + γ)) -
            (1 + (1 + α) * (1 + γ) * log(b)) * b^(-(1 + α) * (1 + γ))
        ) +
        inv((1 + α) * (1 + γ)^2) * (
            b^(-(1 + α) * (1 + γ)) * (1 + (1 + α) * (1 + γ) * log(b)) -
            a^(-(1 + α) * (1 + γ)) * (1 + (1 + α) * (1 + γ) * log(a))
        )

    @show F21_v4_1 F21_v5_1
    @show F21_v4_2 F21_v5_2

    @assert isfinite(F21_v5_1)
    @assert isfinite(F21_v5_2)
    # We now include the infinite sum so the results are note exactly
    # the same. The new one is the correct one.
    @assert isapprox(F21_v4_1, F21_v5_1, rtol = Arb(1e-15))
    @assert isapprox(F21_v4_2, F21_v5_2, rtol = Arb(1e-15))

    #####
    # Simplify after explicitly having computed the sums
    # TODO: This expressions are possibly not the same as in the
    # documentation above, this should be checked.
    #####

    F21_v6_1 =
        -log(x) * ((2 + α) / (1 + γ) * (a^(-(1 + α) * (1 + γ)) - b^(-(1 + α) * (1 + γ))))

    F21_v6_2 =
        (2 + α) / ((1 + α) * (1 + γ)^2) *
        (b^(-(1 + α) * (1 + γ)) - a^(-(1 + α) * (1 + γ))) +
        (2 + α) / (1 + γ) *
        (log(b) * b^(-(1 + α) * (1 + γ)) - log(a) * a^(-(1 + α) * (1 + γ)))

    @assert isfinite(F21_v6_1)
    @assert isfinite(F21_v6_2)
    @assert Arblib.overlaps(F21_v5_1, F21_v6_1)
    @assert Arblib.overlaps(F21_v5_2, F21_v6_2)

    #####
    # Insert back into F2 and G2
    #####

    F21_v3 = F21_v6_1 + F21_v6_2 + F21_v2_3

    F2_v4 = F21_v3 + F22_v1

    G2_v2 = 1 / ((1 - x^p0) * log(inv(x))) * F2_v4

    @assert isapprox(F21_v2, F21_v3, rtol = Arb(1e-15))
    @assert isapprox(F2_v1, F2_v4, rtol = Arb(1e-15))
    @assert isapprox(G2_v1, G2_v2, rtol = Arb(1e-15))

    #####
    # Enclose I3 in terms of I1
    # Note that this gives very poor enclosures if b is large.
    #####

    I3_v10_main(n) = Arb((-log(1 + c * x * b), log(1 + c * x * b))) * I1_v9_main(n)

    I3_v10_remainder(n) =
        Arb((-log(1 + c * x * b), log(1 + c * x * b))) * I1_v9_remainder(n)

    for n = 1:N
        r11, r12 = I3_v9_main(n), I3_v9_remainder(n)
        r21, r22 = I3_v10_main(n), I3_v10_remainder(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert Arblib.overlaps(r11, r21)
        @assert Arblib.overlaps(r12, r22)
    end

    #####
    # Give closed form expression for F21
    # Note that this gives very poor enclosures if b is large.
    #####

    # Same as before
    F21_v7_1 =
        -log(x) * ((2 + α) / (1 + γ) * (a^(-(1 + α) * (1 + γ)) - b^(-(1 + α) * (1 + γ))))

    # Same as before
    F21_v7_2 =
        (2 + α) / ((1 + α) * (1 + γ)^2) *
        (b^(-(1 + α) * (1 + γ)) - a^(-(1 + α) * (1 + γ))) +
        (2 + α) / (1 + γ) *
        (log(b) * b^(-(1 + α) * (1 + γ)) - log(a) * a^(-(1 + α) * (1 + γ)))

    F21_v7_3 =
        Arb((-log(1 + c * x * b), log(1 + c * x * b))) *
        ((2 + α) / (1 + γ) * (a^(-(1 + α) * (1 + γ)) - b^(-(1 + α) * (1 + γ))))

    F21_v4 = F21_v7_1 + F21_v7_2 + F21_v7_3

    @assert Arblib.overlaps(F21_v3, F21_v4)

    #####
    # Write G2 = G21 + G22 in terms of F21 and F22
    #####

    G21_v1 = 1 / ((1 - x^p0) * log(inv(x))) * F21_v4

    G22_v1 = 1 / ((1 - x^p0) * log(inv(x))) * F22_v1
    @show G22_v1
    G2_v3 = G21_v1 + G22_v1

    @assert isfinite(G2_v3)
    # In theory they might not overlap due to truncating infinite sum.
    # However G2_v3 has large overestimations so they should in
    # practice overlap.
    @assert Arblib.overlaps(G2_v2, G2_v3)

    #####
    # Simplify closed form expression for G21
    #####

    G21_v2 = let D = Arb((-log(1 + c * x * b), log(1 + c * x * b))), q0 = (1 + α) * (1 + γ)
        (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) * (
            (D - log(x) - inv(q0)) * (a^(-q0) - b^(-q0)) +
            (log(b) * b^(-q0) - log(a) * a^(-q0))
        )
    end

    @show G21_v1 G21_v2
    @assert isfinite(G21_v2)
    @assert Arblib.overlaps(G21_v1, G21_v2)

    #####
    # Simplify G21 with a = 1
    #####

    # Putting it in the let makes formatting in Emacs weird!
    D = Arb((-log(1 + c * x * b), log(1 + c * x * b)))

    G21_a_eq_1_v1 = let a = Arb(1), q0 = (1 + α) * (1 + γ)
        (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) * (
            (D - log(x) - inv(q0)) * (a^(-q0) - b^(-q0)) +
            (log(b) * b^(-q0) - log(a) * a^(-q0))
        )
    end

    G21_a_eq_1_v2 = let q0 = (1 + α) * (1 + γ)
        (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) *
        ((D - log(x) - 1 / q0) * (1 - b^(-q0)) + log(b) * b^(-q0))
    end

    G21_a_eq_1_v3 = let q0 = (1 + α) * (1 + γ)
        (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) *
        ((D - 1 / q0) * (1 - b^(-q0)) + log(inv(x)) + log(b * x) * b^(-q0))
    end

    @show G21_a_eq_1_v1 G21_a_eq_1_v2 G21_a_eq_1_v3
    @assert Arblib.overlaps(G21_a_eq_1_v1, G21_a_eq_1_v2)
    @assert Arblib.overlaps(G21_a_eq_1_v1, G21_a_eq_1_v3)

    #####
    # Goal: Prove that G21 is bounded by (2 + α) / (1 + γ)
    # Meaning inv((1 - x^p0) * log(inv(x))) *
    # ((D - 1 / q0) * (1 - b^(-q0)) + log(inv(x)) + log(b * x) * b^(-q0))
    # is bounded by 1.
    #####

    # Function for computing G21_a_eq_1_v2 for correct b and given x and α
    return (x, α) ->
        let b = π / x,
            q0 = (1 + α) * (1 + γ),
            p0 = 1 + α + (1 + α)^2 / 2,
            D = Arb((-log(1 + c * x * b), log(1 + c * x * b)))

            (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) *
            ((D - log(x) - 1 / q0) * (1 - b^(-q0)) + log(b) * b^(-q0))
        end

    # NEXT: Start looking at computing bounds for F22
end

"""
T0_asymptotic_main_2_testing_remainder(x::Arb = Arb(1e-10), α::Arb = Arb(-0.9999), γ::Arb = Arb(0.5), c::Arb = 2Arb(ℯ))

Method for testing the development of [`_T0_asymptotic_main_2`](@ref).
The background for that method is very long and it is easy to make
mistakes. This method tests a lot of the rewrites and simplifications
that is done for [`_T0_asymptotic_main_2`](@ref) to catch any
mistakes.

Instead of integrating from `1` to `π / x` it integrates from `2` to
`10`. This avoids any singularities and allows for computing most
integrals numerically.

This method looks at how to bound `G22`.
"""
function _T0_asymptotic_main_2_testing_remainder(
    x::Arb = Arb(1e-10),
    α::Arb = Arb(-0.9999),
    γ::Arb = Arb(0.5),
    c::Arb = 2Arb(ℯ),
)
    # Integration limits
    a, b = Arb(2), Arb(10)

    # Maximum value of n to include
    N = 10

    p0 = 1 + α + (1 + α)^2

    #####
    # Starting point is same bounds as in
    # _T0_asymptotic_main_2_testing
    # In particular I1_v4_remainder, I2_v4_remainder and I3_v4_remainder
    #####

    h_v1(l, t) =
        log(t - 1)^l + log(t + 1)^l - 2log(t)^l -
        l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2

    H1_v1(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * h_v1(n - k, t) * t
    end |> real

    H2_v1(n, k) = Arblib.integrate(a, b) do t
        log(t)^(k + 1) * h_v1(n - k, t) * t
    end |> real

    H3_v1(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * h_v1(n - k, t) * t * log(1 + c * x * t)
    end |> real

    I1_remainder_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H1_v1(n, k)
        end

    I2_remainder_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H2_v1(n, k)
        end

    I3_remainder_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H3_v1(n, k)
        end

    F22_v1 = sum(1:N) do n
        (-1)^n * (1 + α)^n / factorial(n) *
        (-log(x) * I1_remainder_v1(n) - I2_remainder_v1(n) + I3_remainder_v1(n))
    end

    G22_v1 = 1 / ((1 - x^p0) * log(inv(x))) * F22_v1

    #####
    # Split G22 into two individually finite factors
    #####

    G22_1_v1 = (1 + α) / ((1 - x^p0))

    G22_2_v1 = sum(1:N) do n
        (-1)^n * (1 + α)^(n - 1) / factorial(n) * (
            I1_remainder_v1(n) - inv(log(inv(x))) * I2_remainder_v1(n) +
            inv(log(inv(x))) * I3_remainder_v1(n)
        )
    end

    G22_v2 = G22_1_v1 * G22_2_v1
    @show G22_v1 G22_v2
    @assert isfinite(G22_v2)
    @assert Arblib.overlaps(G22_v1, G22_v2)

    #####
    # Ensure that G22_1 is finite for α overlapping -1 and x overlapping zero
    #####

    G22_1_v2 = let α = Arb((-1, α)), x = Arb((0, x))
        lower = α + 1
        upper = let xᵤ = ubound(Arb, x)
            inv(fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), α + 1, extra_degree = 2))
        end
        Arb((lower, upper))
    end

    @assert isfinite(G22_1_v2)
    @assert Arblib.overlaps(G22_1_v1, G22_1_v2)

    #####
    # Simplify I1_remainder, I2_remainder and I3_remainder for n = 1
    #####

    I1_remainder_n1_v1 = Arblib.integrate(a, b) do t
        (log(t - 1) + log(t + 1) - 2log(t) + 1 / t^2) * t
    end |> real

    I2_remainder_n1_v1 =
        Arblib.integrate(a, b) do t
            log(t) * (log(t - 1) + log(t + 1) - 2log(t) + 1 / t^2) * t
        end |> real

    I3_remainder_n1_v1 =
        Arblib.integrate(a, b) do t
            (log(t - 1) + log(t + 1) - 2log(t) + 1 / t^2) * t * log(1 + c * x * t)
        end |> real

    r11, r12, r13 = I1_remainder_v1(1), I2_remainder_v1(1), I3_remainder_v1(1)
    r21, r22, r23 = I1_remainder_n1_v1, I2_remainder_n1_v1, I3_remainder_n1_v1

    @assert isfinite(r21)
    @assert isfinite(r22)
    @assert isfinite(r23)
    @assert Arblib.overlaps(r11, r21)
    @assert Arblib.overlaps(r12, r22)
    @assert Arblib.overlaps(r13, r23)

    #####
    # Explicitly compute I1_remainder, I2_remainder and I3_remainder for n = 1
    #####

    I1_remainder_n1_primitive_v1(t) = (t^2 - 1) * (log(t - 1) + log(t + 1) - 2log(t)) / 2
    I1_remainder_n1_v2 = I1_remainder_n1_primitive_v1(b) - I1_remainder_n1_primitive_v1(a)

    # This is the primitive given by Mathematica. In this form some of
    # the terms are complex, but cancel out.
    I2_remainder_n1_primitive_v1(t) =
        (
            -3 + 4t^2 * log(t) - 4log(Acb(1 - t)) * log(t) + 4log(t)^2 - 8t^2 * log(t)^2 +
            log(t - 1) * (2 - 2t^2 + 4t^2 * log(t)) +
            2(t^2 - 1) * (2log(t) - 1) * log(t + 1) - 2polylog(2, Acb(t^2))
        ) / 8 |> real
    I2_remainder_n1_v2 = I2_remainder_n1_primitive_v1(b) - I2_remainder_n1_primitive_v1(a)

    # We have that log(t - 1) + log(t + 1) - 2log(t) + 1 / t^2 is
    # negative for t > 1. This can be seen from the Taylor expansion.
    # The integrand thus has a constant sign and we can factor out an
    # enclosure of log(1 + c * x * t)
    I3_remainder_n1_v2 = Arb((log(1 + c * x * a), log(1 + c * x * b))) * I1_remainder_n1_v2

    r11, r12, r13 = I1_remainder_n1_v1, I2_remainder_n1_v1, I3_remainder_n1_v1
    r21, r22, r23 = I1_remainder_n1_v2, I2_remainder_n1_v2, I3_remainder_n1_v2

    @assert isfinite(r21)
    @assert isfinite(r22)
    @assert isfinite(r23)
    @assert Arblib.overlaps(r11, r21)
    @assert Arblib.overlaps(r12, r22)
    @assert Arblib.overlaps(r13, r23)

    #####
    # Rewrite I2_remainder(1) to not use Acb
    #####

    # Uses https://dlmf.nist.gov/25.12.E4 to rewrite the polylog
    # polylog still only support evaluation on Acb, but the result is
    # real
    I2_remainder_n1_primitive_v2(t) =
        (
            -3 - Arb(π)^2 / 3 +
            2(log(t - 1) + log(t + 1) + 2log(t)^2) +
            2t^2 * (
                -log(t - 1) - log(t + 1) +
                2log(t) +
                log(t) * (+2log(t - 1) + 2log(t + 1) - 4log(t))
            ) +
            2real(polylog(2, Acb(1 - t^2)))
        ) / 8 |> real
    I2_remainder_n1_v3 = I2_remainder_n1_primitive_v2(b) - I2_remainder_n1_primitive_v2(a)

    r1 = I2_remainder_n1_v2
    r2 = I2_remainder_n1_v3
    @assert isfinite(r1)
    @assert Arblib.overlaps(r1, r2)

    #####
    # Write G22_2 as sum of first term plus remaining sum
    #####

    G22_2_1_v1 = -(
        I1_remainder_n1_v2 - inv(log(inv(x))) * I2_remainder_n1_v3 +
        inv(log(inv(x))) * I3_remainder_n1_v2
    )

    G22_2_2_v1 = sum(2:N) do n
        (-1)^n * (1 + α)^(n - 1) / factorial(n) * (
            I1_remainder_v1(n) - inv(log(inv(x))) * I2_remainder_v1(n) +
            inv(log(inv(x))) * I3_remainder_v1(n)
        )
    end

    G22_2_v2 = G22_2_1_v1 + G22_2_2_v1

    @assert isfinite(G22_2_v2)
    @assert Arblib.overlaps(G22_2_v1, G22_2_v2)

    #####
    # The goal is now to Compute an upper bound of G22_2_2
    #####

    #####
    # Upper bound G22_2_2 by taking the absolute value term wise
    #####

    G22_2_2_upper_v1 = sum(2:N) do n
        (1 + α)^(n - 1) / factorial(n) * (
            abs(I1_remainder_v1(n)) +
            inv(log(inv(x))) * abs(I2_remainder_v1(n)) +
            inv(log(inv(x))) * abs(I3_remainder_v1(n))
        )
    end

    @show G22_2_2_v1

    @assert isfinite(G22_2_2_upper_v1)
    @assert G22_2_2_v1 < G22_2_2_upper_v1

    #####
    # Upper bound I1_remainder, I2_remainder and I3_remainder by
    # taking the absolute value term wise
    #####

    I1_remainder_upper_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * abs(H1_v1(n, k))
        end

    I2_remainder_upper_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * abs(H2_v1(n, k))
        end

    I3_remainder_upper_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * abs(H3_v1(n, k))
        end

    for n = 1:N
        r11, r12, r13 = I1_remainder_v1(n), I2_remainder_v1(n), I3_remainder_v1(n)
        r21, r22, r23 =
            I1_remainder_upper_v1(n), I2_remainder_upper_v1(n), I3_remainder_upper_v1(n)
        @assert isfinite(r21)
        @assert isfinite(r22)
        @assert isfinite(r23)
        @assert !(abs(r11) > r21)
        @assert !(abs(r12) > r22)
        @assert !(abs(r13) > r23)
        @assert abs_ubound(r11) <= ubound(r21)
        @assert abs_ubound(r12) <= ubound(r22)
        @assert abs_ubound(r13) <= ubound(r23)
    end

    #####
    # Upper bound I3_remainder using H1 instead of H3
    #####

    I3_remainder_upper_v2(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * log(1 + c * x * b) * abs(H1_v1(n, k))
        end

    for n = 1:N
        r1 = I3_remainder_upper_v1(n)
        r2 = I3_remainder_upper_v2(n)
        @assert isfinite(r2)
        @assert r1 < r2
    end

    #####
    # Use upper bounds for I1_remainder, I2_remainder and I3_remainder
    # by to get upper bound for G22_2_2
    #####

    G22_2_2_upper_v2 = sum(2:N) do n
        (1 + α)^(n - 1) / factorial(n) * (
            I1_remainder_upper_v1(n) +
            inv(log(inv(x))) * I2_remainder_upper_v1(n) +
            inv(log(inv(x))) * I3_remainder_upper_v2(n)
        )
    end

    @show G22_2_2_upper_v1
    @show G22_2_2_upper_v2

    @assert isfinite(G22_2_2_upper_v2)
    @assert G22_2_2_upper_v1 < G22_2_2_upper_v2

    #####
    # Upper bound abs(H1) and abs(H2) by moving the absolute value
    # inside the integral
    #####

    # In principle we want to compute this
    H1_upper_wrong_v1(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * abs(h_v1(n - k, t)) * t
    end |> real
    H2_upper_wrong_v1(n, k) = Arblib.integrate(a, b) do t
        log(t)^(k + 1) * abs(h_v1(n - k, t)) * t
    end |> real
    # But since the integrand is not analytic we need to add
    # check_analytic. However this makes it converge extremely slowly
    # and we can't really get any meaningful information from it.
    H1_upper_slow_v1(n, k) =
        let f(t; analytic) =
                log(t)^k * Arblib.real_abs!(zero(t), h_v1(n - k, t), analytic) * t
            Arblib.integrate(f, a, b, check_analytic = true) |> real
        end

    H2_upper_slow_v1(n, k) =
        let f(t; analytic) =
                log(t)^(k + 1) * Arblib.real_abs!(zero(t), h_v1(n - k, t), analytic) * t
            Arblib.integrate(f, a, b, check_analytic = true) |> real
        end

    for n = 1:N
        for k = 0:n-1
            # Don't run the tests since the integration doesn't
            # converge

            #r11, r12 = H1_v1(n, k), H2_v1(n, k)
            #r21, r22 = H1_upper_slow_v1(n, k), H2_upper_slow_v1(n, k)
            #@assert isfinite(r21)
            #@assert isfinite(r22)
            #@assert !(abs(r11) > r21)
            #@assert !(abs(r12) > r22)
            #@assert abs_ubound(r11) <= ubound(r21)
            #@assert abs_ubound(r12) <= ubound(r22)
        end
    end

    #####
    # Bound h(l, t) for 1 < t <= t₀
    #####

    # Split into one part with log(t - 1) and one constant part
    # Note that the sign depends on t, negative for t < 2 and positive
    # for t > 2.
    h_upper_1_v1(l, t) = log(t - 1)^l

    h_upper_2_v1(l, t₀) = begin
        log(t₀ + 1)^l + 2log(t₀)^l + l * (l - 1) * log(t₀)^(l - 2) + l * log(t₀)^(l - 1)
    end

    h_upper_v1(l, t, t₀ = Arb(2)) = begin
        @assert t <= t₀
        abs(h_upper_1_v1(l, t)) + h_upper_2_v1(l, t₀)
    end

    for l = 1:N
        for t in range(Arb(1), 2, 100)[2:end]
            @assert abs(h_v1(l, t)) < h_upper_v1(l, t)
        end
    end

    #####
    # Bound abs(H1), abs(H2) and abs(H3) using bound of h(l, t) valid
    # on 1 < t <= b
    #####

    # Everything in the integrals that is constant we factor out out.
    # The only non-trivial part to integrate is log(t - 1)^(n - k).
    # We want to integrate abs(log(t - 1)^(n - k)) but absolute value
    # makes the Arblib integrator sad. Here we just assume that a >= 2
    # so that it is positive.
    a >= 2 || @warn "Integration in H1_upper_v1 and similar assumes that a >= 2"
    H1_upper_v1(n, k) =
        (
            log(b)^k * b * Arblib.integrate(a, b) do t
                h_upper_1_v1(n - k, t)
            end + h_upper_2_v1(n - k, b) * log(b)^k * b * (b - a)
        ) |> real

    H2_upper_v1(n, k) =
        (
            log(b)^(k + 1) * b * Arblib.integrate(a, b) do t
                h_upper_1_v1(n - k, t)
            end + h_upper_2_v1(n - k, b) * log(b)^(k + 1) * b * (b - a)
        ) |> real

    H3_upper_v1(n, k) =
        log(1 + c * x * b) * (
            log(b)^k * b * Arblib.integrate(a, b) do t
                h_upper_1_v1(n - k, t)
            end + h_upper_2_v1(n - k, b) * log(b)^k * b * (b - a)
        ) |> real

    for n = 1:N
        for k = 0:n-1
            r11, r12, r13 = H1_v1(n, k), H2_v1(n, k), H3_v1(n, k)
            r21, r22, r23 = H1_upper_v1(n, k), H2_upper_v1(n, k), H3_upper_v1(n, k)
            @assert isfinite(r21)
            @assert isfinite(r22)
            @assert isfinite(r23)
            @assert Arblib.ispositive(r21)
            @assert Arblib.ispositive(r22)
            @assert Arblib.ispositive(r23)
            @assert !(abs(r11) > r21)
            @assert !(abs(r12) > r22)
            @assert !(abs(r13) > r23)
            @assert abs_ubound(r11) <= ubound(r21)
            @assert abs_ubound(r12) <= ubound(r22)
            @assert abs_ubound(r13) <= ubound(r23)
        end
    end

    #####
    # Give upper bounds of H1, H2 and H3 when integrated from 1 to 2
    #####

    let a = Arb(1), b = Arb(2)
        # The integral ∫ abs(log(t - 1)^(n - k)) dt can in this case be
        # explicitly computed to be factorial(n - k)

        H1_upper_1to2_v1(n, k) =
            log(b)^k * b * factorial(n - k) +
            h_upper_2_v1(n - k, b) * log(b)^k * b * (b - a)

        H2_upper_1to2_v1(n, k) =
            log(b)^(k + 1) * b * factorial(n - k) +
            h_upper_2_v1(n - k, b) * log(b)^(k + 1) * b * (b - a)

        H3_upper_1to2_v1(n, k) =
            log(1 + c * x * b) * (
                log(b)^k * b * factorial(n - k) +
                h_upper_2_v1(n - k, b) * log(b)^k * b * (b - a)
            )

        for n = 1:N
            for k = 0:n-1
                r1, r2, r3 =
                    H1_upper_1to2_v1(n, k), H2_upper_1to2_v1(n, k), H3_upper_1to2_v1(n, k)
                @assert Arblib.ispositive(r1)
                @assert Arblib.ispositive(r2)
                @assert Arblib.ispositive(r3)
            end
        end

        H1_upper_1to2_v2(n, k) =
            log(b)^k * b * (factorial(n - k) + h_upper_2_v1(n - k, b) * (b - a))

        H2_upper_1to2_v2(n, k) =
            log(b)^(k + 1) * b * (factorial(n - k) + h_upper_2_v1(n - k, b) * (b - a))

        H3_upper_1to2_v2(n, k) =
            log(1 + c * x * b) *
            log(b)^k *
            b *
            (factorial(n - k) + h_upper_2_v1(n - k, b) * (b - a))

        for n = 1:N
            for k = 0:n-1
                r1, r2, r3 =
                    H1_upper_1to2_v2(n, k), H2_upper_1to2_v2(n, k), H3_upper_1to2_v2(n, k)
                @assert Arblib.ispositive(r1)
                @assert Arblib.ispositive(r2)
                @assert Arblib.ispositive(r3)
            end
        end
    end

    # NEXT: Give upper bounds for I1_remainder, I2_remainder and
    # I3_remainder when integrated from 1 to 2 using the above bounds
    # for H1, H2 and H3.

    #####
    # NEXT: Bound h(l, t) for t > t₀ > 1
    #####
end

"""
    _T0_asymptotic_main_2_1(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G21(x) = 1 / ((1 - x^p0) * log(inv(x))) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
where the integration is taken from `1` to `2`.

For ``t ∈ [1, 2]`` we have
```
log(c + inv(2x)) <= log(c + inv(x * t)) <= log(c + inv(x))
```
allowing us to factor out `log(c + inv(x * t))`. If we also factor out
`1 + α` we get
```
G21(x) <= (1 + α) / ((1 - x^p0) * log(inv(x))) *
         log(c + inv(x)) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) *
                t^(1 - γ * (1 + α)) dt
```
The integral no longer depends on `x` and converges to a non-zero
number in `α`, which we can compute.

# Computing the integral
Let
```
I = ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) dt
```
Since `1 + α` is small the factor `t^(-γ * (1 + α))` in the integrand
will be close to one. If we bound it on the interval and factor it out
we are left with
```
∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t dt
```
If we let `s = -(α + 1)` a primitive function is given by
```
(
    - t^2 * ((t - 1)^s + (t + 1)^s - 2t^s) / s
    + ((t - 1)^s + (t + 1)^s - 2) / s
    - t * ((t - 1)^(1 + s) + (t + 1)^(1 + s) - 2t^(1 + s))
) / ((1 + s) * (2 + s))
```
where the constant is chosen to make it finite as `s` goes to zero.
For `t = 1` this simplifies to
```
-(2^-α - 2) / (α * (1 - α))
```
This allows us to compute the integral
"""
function _T0_asymptotic_main_2_1(α::Arb, γ::Arb, c::Arb)
    αp1 = Arblib.nonnegative_part!(zero(α), α + 1)

    # Primitive function of
    # ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t
    primitive = let s = -αp1
        t -> begin
            t1 = fx_div_x(v -> (t - 1)^v + (t + 1)^v - 2t^v, s)
            t2 = fx_div_x(v -> (t - 1)^v + (t + 1)^v - 2, s)
            #t1 = ((t - 1)^s + (t + 1)^s - 2t^s) / s
            #t2 = ((t - 1)^s + (t + 1)^s - 2) / s
            (-t^2 * t1 + t2 - t * ((t - 1)^(1 + s) + (t + 1)^(1 + s) - 2t^(1 + s))) / ((1 + s) * (2 + s))
        end
    end

    # primitive(1)
    primitive_one = (2^-α - 2) / (α * (1 - α))

    # Enclosure of
    # ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t^(1 - γ * (1 + α)) dt
    I = Arb((1, 2))^(-γ * (1 + α)) * (primitive(Arb(2)) - primitive_one)

    return x::Arb -> begin
        # This function assumes that x is less than or equal to 1 // 2
        @assert x <= 1 // 2

        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        # Enclosure of log(c + inv(x * Arb((1, 2)))) / log(inv(x))
        # = 1 + log(cx + inv(Arb((1, 2)))) / log(inv(x))
        logfactor_lower = 1 + log(c * x + 1 // 2) * invloginvx
        logfactor_upper = 1 + log1p(c * x) * invloginvx
        logfactor = Arb((logfactor_lower, logfactor_upper))

        # Enclosure of (1 + α) / (1 - x^p0)
        αp1_div_1mxp0 = if iszero(x)
            αp1
        elseif Arblib.contains_zero(x)
            lower = αp1
            upper = let xᵤ = ubound(Arb, x)
                inv(fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), αp1, extra_degree = 2))
            end
            Arb((lower, upper))
        else
            inv(fx_div_x(s -> 1 - x^(s + s^2 / 2), αp1, extra_degree = 2))
        end

        return logfactor * αp1_div_1mxp0 * I
    end
end

"""
    _T0_asymptotic_main_2_2(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G22(x) = 1 / ((1 - x^p0) * log(inv(x))) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
where the integration is taken from `2` to `π / x`.

This method is meant for non-zero `x`.

We factor out `1 + α` from the integral and write it as
```
G22(x) = (1 + α) / ((1 - x^p0) * log(inv(x))) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
and compute the integral and the factor outside of the integral
separately. The integral is computed by integrating numerically.

To compute better enclosures in the numerical integration we compute a
lower bound by integrating to the lower bound of `π / x` and an upper
bound by integrating to the upper bound of it. We also use that
```
log(c + inv(ubound(x) * t)) <= log(c + inv(x * t)) <= log(c + inv(lbound(x) * t))
```
"""
function _T0_asymptotic_main_2_2(α::Arb, γ::Arb, c::Arb)
    return x::Arb -> begin
        # This function assumes that x is less than or equal to 1
        @assert x <= 1 // 2

        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        # Enclosure of (1 + α) / (1 - x^p0)
        αp1_div_1mxp0 = if iszero(x)
            α + 1
        elseif Arblib.contains_zero(x)
            lower = α + 1
            upper = let xᵤ = ubound(Arb, x)
                inv(fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), 1 + α, extra_degree = 2))
            end
            Arb((lower, upper))
        else
            inv(fx_div_x(s -> 1 - x^(s + s^2 / 2), 1 + α, extra_degree = 2))
        end

        # Enclosure of the integral
        extra_degree = 2
        xₗ, xᵤ = getinterval(Arb, x)

        I_lower =
            Arblib.integrate(2, lbound(Arb, π / x), warn_on_no_convergence = false) do t
                if isreal(t)
                    t = real(t)

                    return ArbExtras.enclosure_series(α, degree = 4) do α
                        if α isa ArbSeries && contains(α[0], -1)
                            fx_div_x(α + 1; extra_degree) do αp1
                                (t - 1)^-αp1 + (t + 1)^-αp1 - 2t^-αp1
                            end * t^(1 - γ * (α + 1))
                        else
                            # Large cancellations so do the
                            # computations at a higher precision
                            let αp1 = setprecision(α + 1, 2precision(α))
                                ((t - 1)^-αp1 + (t + 1)^-αp1 - 2t^-αp1) / αp1 *
                                t^(1 - γ * αp1)
                            end
                        end
                    end * log(c + inv(xᵤ * t))
                else
                    return fx_div_x(s -> (t - 1)^-s + (t + 1)^-s - 2t^-s, Acb(1 + α)) *
                           t^(1 - γ * (α + 1)) *
                           log(c + inv(xᵤ * t))
                end
            end |> real

        I_upper =
            Arblib.integrate(2, ubound(π / x), warn_on_no_convergence = false) do t
                if isreal(t)
                    t = real(t)

                    return ArbExtras.enclosure_series(α, degree = 4) do α
                        if α isa ArbSeries && contains(α[0], -1)
                            fx_div_x(α + 1; extra_degree) do αp1
                                (t - 1)^-αp1 + (t + 1)^-αp1 - 2t^-αp1
                            end * t^(1 - γ * (α + 1))
                        else
                            # Large cancellations so do the
                            # computations at a higher precision
                            let αp1 = setprecision(α + 1, 2precision(α))
                                ((t - 1)^-αp1 + (t + 1)^-αp1 - 2t^-αp1) / αp1 *
                                t^(1 - γ * αp1)
                            end
                        end
                    end * log(c + inv(xₗ * t))
                else
                    return fx_div_x(s -> (t - 1)^-s + (t + 1)^-s - 2t^-s, Acb(1 + α)) *
                           t^(1 - γ * (α + 1)) *
                           log(c + inv(xₗ * t))
                end
            end |> real

        I = Arb((I_lower, I_upper))

        return invloginvx * αp1_div_1mxp0 * I
    end
end

"""
    _T0_asymptotic_remainder(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
R(x) = ∫ abs(R(x * abs(1 - t)) + R(x * (1 + t)) - 2R(x * t)) *
            t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
defined in [`T0_asymptotic`](@ref).

From the expansion of the Clausen functions we have
```
R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t) = sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * ((1 - t)^2m + (1 + t)^2m - 2t^2m) for m = 1:Inf
)
```
where we have removed the absolute value around `1 - t` since all
powers are even. We have that
```
(1 - t)^2m + (1 + t)^2m - 2t^2m =
    sum(binomial(2m, k) * (1 + (-1)^k) * t^k for k = 0:2m) - 2t^2m =
    2sum(binomial(2m, 2k) * t^2k for k = 0:m-1)
```
Giving us
```
R(x * (t - 1)) + R(x * (1 + t)) - 2R(x * t) = 2sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m * sum(binomial(2m, 2k) * t^2k for k = 0:m-1) for m = 1:Inf
)
```
The factor `zeta(-α - 2m) * (-1)^m` is positive since `0 < -α < 1` and
`m >= 1`. All terms are therefore positive and adding the integrating
we can write it as
```
R(x) = 2sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * x^2m *
        sum(binomial(2m, 2k) * ∫ t^(2k + 1 - γ * (1 + α)) * log(c + inv(x * t)) dt for k = 0:m-1) for m = 1:Inf
)
```
where the integration is taken from `0` to `π / x`. For the integral
we first split the log-factor as
```
log(1 + inv(x * t)) = log(1 + c * x * t) - log(x * t)
```
giving us
```
∫ t^(2k + 1 - γ * (1 + α)) * log(c + inv(x * t)) dt =
    ∫ t^(2k + 1 - γ * (1 + α)) * log(1 + c * x * t) dt
    - ∫ t^(2k + 1 - γ * (1 + α)) * log(x * t) dt
```
Using that `log(1 + c * x * t)` is bounded by `log(1 + c * π)` and
integrating explicitly gives us
```
∫ t^(2k + 1 - γ * (1 + α)) * log(1 + c * x * t) dt <=
    log(1 + c * π) * ∫ t^(2k + 1 - γ * (1 + α)) dt =
    log(1 + c * π) * (π / x)^(2(k + 1) - γ * (1 + α)) / (2(k + 1) - γ * (1 + α)) <=
    log(1 + c * π) * (π / x)^(2(k + 1)) / (2(k + 1) - γ * (1 + α))
```
```
∫ t^(2k + 1 - γ * (1 + α)) * log(x * t) dt =
    ((2(k + 1) - γ * (1 + α)) * log(π) - 1) * (π / x)^(2(k + 1) - γ * (1 + α)) / (2(k + 1) - γ * (1 + α))^2
    <= log(π) * (π / x)^(2(k + 1)) / (2(k + 1) - γ * (1 + α))
```
And hence
```
∫ t^(2k + 1 - γ * (1 + α)) * log(c + inv(x * t)) dt <=
    log(c + inv(π)) * (π / x)^(2(k + 1)) / (2(k + 1) - γ * (1 + α))
```
Inserting this back into `R(x)` and factoring out `(π / x)^2m` we get
```
R(x) <= 2log(c + inv(π)) * sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * π^2m *
        sum(binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α)) for k = 0:m-1) for m = 1:Inf
)
```
We split this sum into `m = 1:M-1` and `m = M:Inf`. For the finite
part we sum directly. For the infinite tail we note that
```
sum(binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α)) for k = 0:m-1) <=
    sum(binomial(2m, 2k) * (1 / π)^(2(m - 1 - k)) / (2k + 1) for k = 0:m-1) <=
    sum(binomial(2m, 2k) * (1 / π)^(2(m - 1 - k)) / (2k + 1) for k = 0:m) =
    π^(2 - 2m) * ((π - 1)^2m + (π + 1)^2m) / 2 <=
    π^(2 - 2m) * (π + 1)^2m <=
    = π^2 * (1 + inv(π))^2m
```
Giving us
```
sum(
    zeta(-α - 2m) * (-1)^m / factorial(2m) * π^2m *
        sum(binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α)) for k = 0:m-1) for m = M:Inf
) <=
π^2 * sum(zeta(-α - 2m) * (-1)^m / factorial(2m) * (π * (1 + inv(π)))^2m for m = M:Inf)
```
This sum is exactly the same as in the remainder for `clausenc` in
[`clausenc_expansion_remainder`](@ref), substituting `x` for `π * (1 +
inv(π))`, and we can get a bound from that function.
"""
function _T0_asymptotic_remainder(α::Arb, γ::Arb, c::Arb; M = 20)
    # Precompute the factor
    # (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m))
    # in the finite sum
    factors = [(-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) for m = 1:M-1]

    # Compute tail of sum, which is the same as for clausenc_expansion
    # It is valid for x < 1
    tail = let y = π * (1 + inv(Arb(π)))
        Arb(π)^2 * clausenc_expansion_remainder(y, -α, M) * y^2M
    end

    return x::Arb -> begin
        x < 1 || throw(DomainError(x, "x must be less than 1"))

        # Sum directly for m = 1:M-1
        main = sum(1:M-1, init = zero(x)) do m
            factors[m] * sum(
                binomial(2m, 2k) * abspow(x / π, 2(m - 1 - k)) / (2(k + 1) - γ * (1 + α)) for k = 0:m-1
            )
        end

        return 2log(c + inv(Arb(π))) * (main + tail)
    end
end
