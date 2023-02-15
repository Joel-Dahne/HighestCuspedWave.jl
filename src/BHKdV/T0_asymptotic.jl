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
        # [b, 1]
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
        # [b, 1]
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

If `x` does not overlap zero it uses [`_T0_asymptotic_main_2_1`](@ref)
and [`_T0_asymptotic_main_2_2`](@ref) to compute the integral.
Factoring out `1 + α` from the integral we get
```
G2(x) = (1 + α) / (1 - x^p0) * 1 / log(inv(x)) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α)*
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
The factor `(1 + α) / (1 - x^p0)` can be enclosed by handling the
removable singularity. For the remaining part we split the interval of
integration into ``[1, 2]`` and ``[2, π / x]``, which are handled in
[`_T0_asymptotic_main_2_1`](@ref) and
[`_T0_asymptotic_main_2_2`](@ref). This only works well as long as `x`
is not too small, in practice it is only used for `x > 1e-10`.

When `x` overlaps zero significantly more work is required to get a
bound. The details are given in the corresponding section of the
paper. We here give an overview of the approach used.

# Split into `G2_M`, `G2_R_1` and `G2_R_2`
We have
```
G2(x) <= G2_M(x) + G2_R_factor(x) * (G2_R_n1(x) + G2_R_1(x) + G2_R_2(x))
```
Here
```
G2_M(x) = inv(1 - x^(1 + α + (1 + α)^2 / 2)) * inv(log(inv(x))) * (4 + 2α) / 3 *
    (
        (D - log(x) - 2 / 3 * inv(1 + α)) * (1 - (x / π)^(3 / 2 * (1 + α)))
        + log(π / x) * (x / π)^(3 / 2 * (1 + α))
    )
```
with `D = Arb((-log(1 + cπ), log(1 + cπ)))`.
```
G2_R_factor(x) = (1 + α) / (1 - x^(1 + α + (1 + α)^2 / 2)) * (1 + log(1 + c * x) / log(1 / x))
```
```
G2_R_n1 = (1 + α)^(n - 1) / factorial(n) * sum(0:n-1) do k
        binomial(n, k) / 2^k * ∫_1^(π / x) log(t)^k * abs(h(n - k, t)) * t dt
    end
```
```
G2_R_1(x) = sum(2:Inf) do n
    (1 + α)^(n - 1) / factorial(n) * sum(0:n-1) do k
        binomial(n, k) / 2^k * ∫_1^2 log(t)^k * abs(h(n - k, t)) * t dt
    end
end
```
and
```
G2_R_2(x) = sum(2:Inf) do n
    (1 + α)^(n - 1) / factorial(n) * sum(0:n-1) do k
        binomial(n, k) / 2^k * ∫_2^(π / x) log(t)^k * abs(h(n - k, t)) * t dt
    end
end
```
With
```
h(k, t) = log(t - 1)^k + log(t + 1)^k - 2log(t)^k - k * (k - 1 - log(t)) * log(t)^(k - 2) / t
```

Bounds for the functions `G2_M(x)`, `G2_R_factor(x)`, `G2_R_1(x)` and
`G2_R_2(x)` are computed individually.

Note that compared to the version in the paper we use the variable `c`
instead of the explicit value `2ℯ`. In practice `c` will be equal to
`2ℯ` but the code supports using different values as well. We avoid
specifying the exact requirements on `c` for these bounds to hold
though, since in the actual proof we only use `c = 2ℯ`.

# Bounding `G2_M(x)`
We have the following bound for `G2_M(x)`, from a lemma in the paper.
```
G2_M(x) <= (4 + 2α) / 3 * (2log(1 + cπ) / log(inv(x)) + 1)
```
Which holds for `x < inv(π^3)`.

# Enclosing `G2_R_factor(x)`
We can enclose `(1 + α) / (1 - x^(1 + α + (1 + α)^2 / 2))` using that
it is increasing in `x`. It is equal to `1 + α` for `x = 0` and for `x
> 0` we can handle the removable singularity.

# Bounding `G2_R_n1(x)`
From a lemma in the paper we have that this is bounded by `1 / 2` for
`x < 1`.

# Bounding `G2_R_1(x)`
We have the following bound for `G2_R_1(x)`, from a lemma in the
paper.
```
G2_R_1(x) <= 2(
    sqrt(ℯ) * (1 + α) / (-α)
    + 4(exp(3(1 + α) - 3(1 + α) - 1)) / (3(1 + α))
    + (2 + α) * exp(3 / 2 * (1 + α))
    - 1
)
```

# Bounding `G2_R_2(x)`
We have the following bound for `G2_R_2(x)`, from a lemma in the
paper.
```
G2_R_2(x) <= 192log(2) * (1 + α) / (-α)
```
"""
function _T0_asymptotic_main_2(α::Arb, γ::Arb, c::Arb)

    f1 = _T0_asymptotic_main_2_1(α, γ, c)
    f2 = _T0_asymptotic_main_2_2(α, γ, c)

    return x::Arb -> begin
        x < 1 || throw(DomainError(x, "must have x < 1"))

        if Arblib.contains_zero(x)
            @assert γ == 1 // 2 # The paper assumes this

            # Enclosure of inv(log(inv(x)))
            invloginvx = if iszero(x)
                zero(x)
            elseif Arblib.contains_zero(x)
                -Arb((inv(log(ubound(Arb, x))), 0))
            else
                -inv(log(x))
            end

            G2_M = let D = Arb((-log(1 + c * π), log(1 + c * π)))
                G21_1 = if x < inv(Arb(π)^3) # Required for the bound we use
                    2invloginvx
                else
                    indeterminate(x)
                end
                G21_2 = one(G21_1)

                (2 + α) / (1 + γ) * (D * G21_1 + G21_2)
            end

            # First factor of remainder term.
            # Given by (1 + α) / ((1 - x^p0)) * (1 + log(1 + c * x) / log(1 / x))
            G2_R_factor = begin
                lower = α + 1
                upper = let xᵤ = ubound(Arb, x)
                    inv(fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), α + 1, extra_degree = 2))
                end
                Arb((lower, upper)) * (1 + log(1 + c * x) * invloginvx)
            end

            # First term in sum of second factor of remainder term, n = 1
            G2_R_n1 = Arb(1 // 2)

            # Remaining terms in second factor of remainder terms
            # integrated from 1 to 2
            G2_R_1 = begin
                s = fx_div_x(3(1 + α)) do t
                    exp(t) - t - 1
                end

                2(sqrt(Arb(ℯ)) * (1 + α) / (-α) + 4s + (2 + α) * exp(3(1 + α) / 2) - 1)
            end

            # Remaining terms in second factor of remainder terms
            # integrated from 2 to π / x
            G2_R_2 = 192log(Arb(2)) * (1 + α) / (-α)

            # Combine factors and all parts of sum
            G2_R = G2_R_factor * (G2_R_n1 + G2_R_1 + G2_R_2)

            G2 = G2_M + G2_R

            return G2
        else
            # Enclosure of (1 + α) / (1 - x^p0)
            αp1_div_1mxp0 = inv(fx_div_x(s -> 1 - x^(s + s^2 / 2), 1 + α, extra_degree = 2))

            return αp1_div_1mxp0 * (f1(x) + f2(x))
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

    h_upper_2_v1(l, t₀) =
        log(t₀ + 1)^l + 2log(t₀)^l + l * (l - 1) * log(t₀)^(l - 2) + l * log(t₀)^(l - 1)

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
                r11, r12, r13 =
                    H1_upper_1to2_v1(n, k), H2_upper_1to2_v1(n, k), H3_upper_1to2_v1(n, k)
                r21, r22, r23 =
                    H1_upper_1to2_v2(n, k), H2_upper_1to2_v2(n, k), H3_upper_1to2_v2(n, k)
                @assert Arblib.ispositive(r21)
                @assert Arblib.ispositive(r22)
                @assert Arblib.ispositive(r23)
                @assert Arblib.overlaps(r11, r21)
                @assert Arblib.overlaps(r12, r22)
                @assert Arblib.overlaps(r13, r23)
            end
        end

        H1_upper_1to2_v3(n, k) =
            log(b)^k * b * (factorial(n - k) + h_upper_2_v1(n - k, b) * (b - a))

        H2_upper_1to2_v3(n, k) =
            log(b)^(k + 1) * b * (factorial(n - k) + h_upper_2_v1(n - k, b) * (b - a))

        H3_upper_1to2_v3(n, k) =
            log(1 + c * x * b) *
            log(b)^k *
            b *
            (factorial(n - k) + h_upper_2_v1(n - k, b) * (b - a))

        for n = 1:N
            for k = 0:n-1
                r11, r12, r13 =
                    H1_upper_1to2_v2(n, k), H2_upper_1to2_v2(n, k), H3_upper_1to2_v2(n, k)
                r21, r22, r23 =
                    H1_upper_1to2_v3(n, k), H2_upper_1to2_v3(n, k), H3_upper_1to2_v3(n, k)
                @assert Arblib.ispositive(r21)
                @assert Arblib.ispositive(r22)
                @assert Arblib.ispositive(r23)
                @assert Arblib.overlaps(r11, r21)
                @assert Arblib.overlaps(r12, r22)
                @assert Arblib.overlaps(r13, r23)
            end
        end

        #####
        # Give upper bounds of I1_remainder, I2_remainder and I3_remainder
        # when integrated from 1 to 2.
        #####

        I1_remainder_upper_1to2_v1(n) =
            sum(0:n-1) do k
                binomial(n, k) * γ^k * H1_upper_1to2_v2(n, k)
            end

        I2_remainder_upper_1to2_v1(n) =
            sum(0:n-1) do k
                binomial(n, k) * γ^k * H2_upper_1to2_v2(n, k)
            end

        I3_remainder_upper_1to2_v1(n) =
            sum(0:n-1) do k
                binomial(n, k) * γ^k * H3_upper_1to2_v2(n, k)
            end

        for n = 1:N
            r1 = I1_remainder_upper_1to2_v1(n)
            r2 = I2_remainder_upper_1to2_v1(n)
            r3 = I3_remainder_upper_1to2_v1(n)

            @assert Arblib.ispositive(r1)
            @assert Arblib.ispositive(r2)
            @assert Arblib.ispositive(r3)
        end

        # Insert definitions of H1, H2 and H3 and simplify the sums

        I1_remainder_upper_1to2_v2(n) =
            b * (
                sum(0:n-1) do k
                    binomial(n, k) * (γ * log(b))^k * factorial(n - k)
                end +
                (b - a) * (
                    sum(
                        binomial(n, k) * γ^k * log(b)^k * log(b + 1)^(n - k) for k = 0:n-1
                    ) +
                    2log(b)^n * sum(binomial(n, k) * γ^k for k = 0:n-1) +
                    log(b)^(n - 2) *
                    sum(binomial(n, k) * γ^k * (n - k) * (n - k - 1) for k = 0:n-1) +
                    log(b)^(n - 1) * sum(binomial(n, k) * γ^k * (n - k) for k = 0:n-1)
                )
            )

        # Both of these can be written in terms of I1_remainder
        I2_remainder_upper_1to2_v2(n) = log(b) * I1_remainder_upper_1to2_v2(n)
        I3_remainder_upper_1to2_v2(n) = log(1 + c * x * b) * I1_remainder_upper_1to2_v2(n)

        for n = 1:N
            r11 = I1_remainder_upper_1to2_v1(n)
            r12 = I2_remainder_upper_1to2_v1(n)
            r13 = I3_remainder_upper_1to2_v1(n)
            r21 = I1_remainder_upper_1to2_v2(n)
            r22 = I2_remainder_upper_1to2_v2(n)
            r23 = I3_remainder_upper_1to2_v2(n)

            @assert Arblib.ispositive(r21)
            @assert Arblib.ispositive(r22)
            @assert Arblib.ispositive(r23)
            @assert Arblib.overlaps(r11, r21)
            @assert Arblib.overlaps(r12, r22)
            @assert Arblib.overlaps(r13, r23)
        end

        # For the sums we have
        # sum(binomial(n, k) * (γ * log(b))^k * factorial(n - k) for k in 0:n-1) = b^γ * n * gamma(n, γ * log(b))
        # sum(binomial(n, k) * γ^k for k = 0:n-1) = (1 + γ)^n - γ^n
        # sum(binomial(n, k) * γ^k * (n - k) for k = 0:n-1) = n * (1 + γ)^(n - 1)
        # sum(binomial(n, k) * γ^k * (n - k) * (n - k - 1) for k = 0:n-1) = n * (n - 1) * (1 + γ)^(n - 2)
        # sum(binomial(n, k) * γ^k * log(b)^k * log(b + 1)^(n - k) for k = 0:n-1)
        # = log(b + 1)^n * (1 + γ * log(b) / log(b + 1))^n - (γ * log(b))^n

        I1_remainder_upper_1to2_v3(n) =
            b * (
                b^γ * n * gamma(Arb(n), γ * log(b)) +
                (b - a) * (
                    (log(b + 1) + γ * log(b))^n - 3(γ * log(b))^n +
                    2(log(b) * (1 + γ))^n +
                    (log(b) * (1 + γ))^(n - 1) * n +
                    (log(b) * (1 + γ))^(n - 2) * n * (n - 1)
                )
            )

        I2_remainder_upper_1to2_v3(n) = log(b) * I1_remainder_upper_1to2_v3(n)
        I3_remainder_upper_1to2_v3(n) = log(1 + c * x * b) * I1_remainder_upper_1to2_v3(n)

        for n = 1:N
            r11 = I1_remainder_upper_1to2_v2(n)
            r12 = I2_remainder_upper_1to2_v2(n)
            r13 = I3_remainder_upper_1to2_v2(n)
            r21 = I1_remainder_upper_1to2_v3(n)
            r22 = I2_remainder_upper_1to2_v3(n)
            r23 = I3_remainder_upper_1to2_v3(n)
            @assert Arblib.ispositive(r21)
            @assert Arblib.ispositive(r22)
            @assert Arblib.ispositive(r23)
            @assert Arblib.overlaps(r11, r21)
            @assert Arblib.overlaps(r12, r22)
            @assert Arblib.overlaps(r13, r23)
        end

        #####
        # Insert I1_remainder, I2_remainder and I3_remainder back in G22_2_2
        #####

        G22_2_2_upper_1to2_v1 = sum(2:N) do n
            (1 + α)^(n - 1) / factorial(n) * (
                I1_remainder_upper_1to2_v3(n) +
                inv(log(inv(x))) * I2_remainder_upper_1to2_v3(n) +
                inv(log(inv(x))) * I3_remainder_upper_1to2_v3(n)
            )
        end

        @show G22_2_2_upper_1to2_v1
        @assert Arblib.ispositive(G22_2_2_upper_1to2_v1)

        #####
        # Simplify G22_2_2
        #####

        G22_2_2_upper_1to2_v2 =
            (1 + inv(log(inv(x))) * (log(b) + log(1 + c * x * b))) * sum(2:N) do n
                (1 + α)^(n - 1) / factorial(n) * I1_remainder_upper_1to2_v3(n)
            end

        @assert Arblib.ispositive(G22_2_2_upper_1to2_v2)
        @assert Arblib.overlaps(G22_2_2_upper_1to2_v1, G22_2_2_upper_1to2_v2)

        #####
        # Insert I1_remainder_upper_1to2 directly and split into
        # several sums
        #####

        G22_2_2_upper_1to2_v3 =
            (1 + inv(log(inv(x))) * (log(b) + log(1 + c * x * b))) *
            b *
            (
                b^γ * sum(
                    (1 + α)^(n - 1) * n * gamma(Arb(n), γ * log(b)) / factorial(n) for
                    n = 2:N
                ) +
                (b - a) * (
                    sum(
                        (1 + α)^(n - 1) / factorial(n) * (log(b + 1) + γ * log(b))^n for
                        n = 2:N
                    ) - 3sum((1 + α)^(n - 1) / factorial(n) * (γ * log(b))^n for n = 2:N) +
                    2sum(
                        (1 + α)^(n - 1) / factorial(n) * (log(b) * (1 + γ))^n for n = 2:N
                    ) +
                    sum(
                        (1 + α)^(n - 1) / factorial(n) * (log(b) * (1 + γ))^(n - 1) * n for
                        n = 2:N
                    ) +
                    sum(
                        (1 + α)^(n - 1) / factorial(n) *
                        (log(b) * (1 + γ))^(n - 2) *
                        n *
                        (n - 1) for n = 2:N
                    )
                )
            )

        @assert Arblib.ispositive(G22_2_2_upper_1to2_v3)
        @assert Arblib.overlaps(G22_2_2_upper_1to2_v2, G22_2_2_upper_1to2_v3)

        #####
        # Simplify the sums
        #####

        G22_2_2_upper_1to2_v4 =
            (1 + inv(log(inv(x))) * (log(b) + log(1 + c * x * b))) *
            b *
            (
                b^γ * sum(
                    (1 + α)^(n - 1) * n * gamma(Arb(n), γ * log(b)) / factorial(n) for
                    n = 2:N
                ) +
                (b - a) * (
                    (log(b + 1) + γ * log(b)) * sum(
                        ((1 + α) * (log(b + 1) + γ * log(b)))^(n - 1) / factorial(n) for
                        n = 2:N
                    ) -
                    3γ *
                    log(b) *
                    sum(((1 + α) * γ * log(b))^(n - 1) / factorial(n) for n = 2:N) +
                    2log(b) *
                    (1 + γ) *
                    sum(((1 + α) * log(b) * (1 + γ))^(n - 1) / factorial(n) for n = 2:N) +
                    sum(
                        ((1 + α) * log(b) * (1 + γ))^(n - 1) / factorial(n - 1) for n = 2:N
                    ) +
                    (1 + α) *
                    sum(((1 + α) * log(b) * (1 + γ))^(n - 2) / factorial(n - 2) for n = 2:N)
                )
            )

        @assert Arblib.ispositive(G22_2_2_upper_1to2_v4)
        @assert Arblib.overlaps(G22_2_2_upper_1to2_v3, G22_2_2_upper_1to2_v4)

        #####
        # Explicitly compute the sums
        #####

        # Use that
        # sum(c^(n - 1) / factorial(n) for n = 2:Inf) = (exp(c) - c - 1) / c
        # sum(c^(n - 1) / factorial(n - 1) for n = 2:Inf) = exp(c) - 1
        # sum(c^(n - 2) / factorial(n - 2) for n = 2:Inf) = exp(c)

        G22_2_2_upper_1to2_v5 =
            (1 + inv(log(inv(x))) * (log(b) + log(1 + c * x * b))) *
            b *
            (
                b^γ * sum(
                    (1 + α)^(n - 1) * n * gamma(Arb(n), γ * log(b)) / factorial(n) for
                    n = 2:N
                ) +
                (b - a) * (
                    (
                        exp((1 + α) * (log(b + 1) + γ * log(b))) -
                        (1 + α) * (log(b + 1) + γ * log(b)) - 1
                    ) / (1 + α) -
                    3(exp((1 + α) * γ * log(b)) - (1 + α) * γ * log(b) - 1) / (1 + α) +
                    2(exp((1 + α) * log(b) * (1 + γ)) - (1 + α) * log(b) * (1 + γ) - 1) /
                    (1 + α) +
                    (exp((1 + α) * log(b) * (1 + γ)) - 1) +
                    (1 + α) * exp((1 + α) * log(b) * (1 + γ))
                )
            )

        @show G22_2_2_upper_1to2_v5
        @assert Arblib.ispositive(G22_2_2_upper_1to2_v5)
        @assert isapprox(G22_2_2_upper_1to2_v4, G22_2_2_upper_1to2_v5, rtol = Arb(1e-20))

        #####
        # Handle gamma(Arb(n), γ * log(b)) and compute the last sum
        #####

        # We have gamma(Arb(n), γ * log(b)) < gamma(n) = factorial(n - 1)
        # Leaving us with the sum sum((1 + α)^(n - 1) for n = 2:N) = (1 + α) / (-α)

        G22_2_2_upper_1to2_v6 =
            (1 + inv(log(inv(x))) * (log(b) + log(1 + c * x * b))) *
            b *
            (
                b^γ * (1 + α) / (-α) +
                (b - a) * (
                    (
                        exp((1 + α) * (log(b + 1) + γ * log(b))) -
                        (1 + α) * (log(b + 1) + γ * log(b)) - 1
                    ) / (1 + α) -
                    3(exp((1 + α) * γ * log(b)) - (1 + α) * γ * log(b) - 1) / (1 + α) +
                    2(exp((1 + α) * log(b) * (1 + γ)) - (1 + α) * log(b) * (1 + γ) - 1) /
                    (1 + α) +
                    (exp((1 + α) * log(b) * (1 + γ)) - 1) +
                    (1 + α) * exp((1 + α) * log(b) * (1 + γ))
                )
            )

        @show G22_2_2_upper_1to2_v6
        @assert Arblib.ispositive(G22_2_2_upper_1to2_v6)
        @assert G22_2_2_upper_1to2_v5 < G22_2_2_upper_1to2_v6


        #####
        # Isolate terms which have removable singularities
        #####

        G22_2_2_upper_1to2_v7 = let
            t1 = expm1((1 + α) * (log(b + 1) + γ * log(b))) / (1 + α)
            t2 = expm1((1 + α) * γ * log(b)) / (1 + α)
            t3 = expm1((1 + α) * log(b) * (1 + γ)) / (1 + α)

            (1 + inv(log(inv(x))) * (log(b) + log(1 + c * x * b))) *
            b *
            (
                b^γ * (1 + α) / (-α) +
                (b - a) * (
                    t1 - (log(b + 1) + γ * log(b)) - 3t2 + 3γ * log(b) + 2t3 -
                    2log(b) * (1 + γ) +
                    expm1((1 + α) * log(b) * (1 + γ)) +
                    (1 + α) * exp((1 + α) * log(b) * (1 + γ))
                )
            )
        end

        @show G22_2_2_upper_1to2_v7
        @assert Arblib.ispositive(G22_2_2_upper_1to2_v7)
        @assert Arblib.overlaps(G22_2_2_upper_1to2_v6, G22_2_2_upper_1to2_v7)

        #####
        # Ensure that it works for α overlapping -1 and x overlapping 0
        #####

        G22_2_2_upper_1to2_v8 = let α = union(Arb(-1), α), x = union(Arb(0), x)
            t1 = fx_div_x(1 + α) do s
                exp(s * (log(b + 1) + γ * log(b))) - 1
            end
            t2 = fx_div_x(1 + α) do s
                exp(s * γ * log(b)) - 1
            end
            t3 = fx_div_x(1 + α) do s
                exp(s * log(b) * (1 + γ)) - 1
            end

            invloginvx = if iszero(x)
                zero(x)
            elseif Arblib.contains_zero(x)
                -Arb((inv(log(ubound(Arb, x))), 0))
            else
                -inv(log(x))
            end

            (1 + invloginvx * (log(b) + log(1 + c * x * b))) *
            b *
            (
                b^γ * (1 + α) / (-α) +
                (b - a) * (
                    t1 - (log(b + 1) + γ * log(b)) - 3t2 + 3γ * log(b) + 2t3 -
                    2log(b) * (1 + γ) +
                    expm1((1 + α) * log(b) * (1 + γ)) +
                    (1 + α) * exp((1 + α) * log(b) * (1 + γ))
                )
            )
        end

        @show G22_2_2_upper_1to2_v8
        @assert Arblib.overlaps(G22_2_2_upper_1to2_v7, G22_2_2_upper_1to2_v8)
    end
end

"""
T0_asymptotic_main_2_testing_remainder_tail(x::Arb = Arb(1e-10), α::Arb = Arb(-0.9999), γ::Arb = Arb(0.5), c::Arb = 2Arb(ℯ))

Method for testing the development of [`_T0_asymptotic_main_2`](@ref).
The background for that method is very long and it is easy to make
mistakes. This method tests a lot of the rewrites and simplifications
that is done for [`_T0_asymptotic_main_2`](@ref) to catch any
mistakes.

This tests the procedure for bounding what in
[`_T0_asymptotic_main_2_testing_remainder`](@ref) is called `G22_2_2`
when the integration is done for the tail of the interval. In the real
code this will be from `2` to `π / x`. For the testing it is from `2`
to `10`.
"""
function _T0_asymptotic_main_2_testing_remainder_tail(
    x::Arb = Arb(1e-10),
    α::Arb = Arb(-0.9999),
    γ::Arb = Arb(0.5),
    c::Arb = 2Arb(ℯ),
)
    # Integration limits
    a, b = Arb(2), Arb(10)

    # Maximum value of n to include
    N = 10

    #####
    # Starting point is G22_2_2_v1 from _T0_asymptotic_main_2_testing_remainder
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

    G22_2_2_v1 = sum(2:N) do n
        (-1)^n * (1 + α)^(n - 1) / factorial(n) * (
            I1_remainder_v1(n) - inv(log(inv(x))) * I2_remainder_v1(n) +
            inv(log(inv(x))) * I3_remainder_v1(n)
        )
    end
    @show G22_2_2_v1

    #####
    # Expand h(l, t) as binomial sum
    #####

    h_v2(l, t) =
        sum(0:l-1) do j
            binomial(l, j) * log(t)^j * (log(1 - 1 / t)^(l - j) + log(1 + 1 / t)^(l - j))
        end - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2

    for l = 1:N
        for t in range(a, b, 10)
            r1 = h_v1(l, t)
            r2 = h_v2(l, t)
            @assert Arblib.overlaps(r1, r2)
        end
    end

    #####
    # Expand logarithms into truncated sums
    #####

    K = 40

    h_v3(l, t) =
        sum(0:l-1) do j
            binomial(l, j) *
            log(t)^j *
            (
                sum(-inv(k * t^k) for k = 1:K)^(l - j) +
                sum(-(-1)^k * inv(k * t^k) for k = 1:K)^(l - j)
            )
        end - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2

    for l = 1:N
        for t in range(a, b, 10)
            r1 = h_v2(l, t)
            r2 = h_v3(l, t)
            @assert isapprox(r1, r2, atol = Arb(1e-10))
        end
    end

    #####
    # Factor out inv(t)^(l - j)
    #####

    K = 40

    h_v4(l, t) =
        sum(0:l-1) do j
            binomial(l, j) *
            log(t)^j *
            inv(t)^(l - j) *
            (
                sum(-inv((k + 1) * t^k) for k = 0:K-1)^(l - j) +
                sum((-1)^k * inv((k + 1) * t^k) for k = 0:K-1)^(l - j)
            )
        end - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2

    for l = 1:N
        for t in range(a, b, 10)
            r1 = h_v3(l, t)
            r2 = h_v4(l, t)
            @assert Arblib.overlaps(r1, r2)
        end
    end

    #####
    # Simplify the sums
    #####

    S1_v1(p, t) = sum(-inv((k + 1) * t^k) for k = 0:K-1)^p
    S2_v1(p, t) = sum((-1)^k * inv((k + 1) * t^k) for k = 0:K-1)^p
    S_v1(n, t) = S1_v1(n, t) + S2_v1(n, t)

    S1_cs_v1(n) =
        let cs = Vector{Arb}(undef, K)
            cs[1] = (-1)^n
            for m = 1:K-1
                cs[m+1] =
                    inv(-Arb(m)) *
                    sum((k * n - m + k) * (-inv(Arb(k + 1))) * cs[m-k+1] for k = 1:m)
            end
            cs
        end

    S2_cs_v1(n) =
        let cs = Vector{Arb}(undef, K)
            cs[1] = 1
            for m = 1:K-1
                cs[m+1] =
                    inv(Arb(m)) *
                    sum((k * n - m + k) * inv((-1)^k * Arb(k + 1)) * cs[m-k+1] for k = 1:m)
            end
            cs
        end


    S1_v2(n, t) =
        let cs = S1_cs_v1(n)
            sum(cs[q+1] / t^q for q = 0:K-1)
        end
    S2_v2(n, t) =
        let cs = S2_cs_v1(n)
            sum(cs[q+1] / t^q for q = 0:K-1)
        end
    S_v2(n, t) = S1_v2(n, t) + S2_v2(n, t)

    for p = 0:N
        for t in range(a, b, 4)
            # The degree is truncated so they will not be exactly
            # equal
            tol = a < 3 ? Arb(1e-3) : Arb(1e-5)
            @assert isapprox(S1_v1(p, t), S1_v2(p, t), rtol = tol, atol = tol)
            @assert isapprox(S2_v1(p, t), S2_v2(p, t), rtol = tol, atol = tol)
        end
    end

    #####
    # Join the sums
    #####

    S_v3(n, t) =
        let cs1 = S1_cs_v1(n), cs2 = S2_cs_v1(n)
            sum((cs1[q+1] + cs2[q+1]) / t^q for q = 0:K-1)
        end

    for n = 0:N
        for t in range(a, b, 4)
            # The degree is truncated so they will not be exactly
            # equal
            @assert Arblib.overlaps(S_v2(n, t), S_v3(n, t))
        end
    end

    #####
    # Check that coefficients cancel for even q and are equal for odd
    # q if n is odd and the opposite if n is even
    #####

    for n = 0:N
        cs1 = S1_cs_v1(n)
        cs2 = S2_cs_v1(n)
        for q = 0:K-1
            if isodd(n) && iseven(q) || iseven(n) && isodd(q)
                @assert Arblib.overlaps(cs1[q+1], -cs2[q+1])
            else
                @assert Arblib.overlaps(cs1[q+1], cs2[q+1])
            end
        end
    end

    #####
    # Use the above to simplify the coefficients
    #####

    S_cs_v1(n) =
        let cs1 = S1_cs_v1(n)
            cs = 2cs1
            if isodd(n)
                for q = 0:2:K-1
                    cs[q+1] = 0
                end
            else
                for q = 1:2:K-1
                    cs[q+1] = 0
                end
            end

            cs
        end

    for n = 0:N
        cs1 = S1_cs_v1(n)
        cs2 = S2_cs_v1(n)
        cs = S_cs_v1(n)
        for q = 0:K-1
            @assert Arblib.overlaps(cs1[q+1] + cs2[q+1], cs[q+1])
        end
    end

    S_v4(n, t) =
        let cs = S_cs_v1(n)
            sum(cs[q+1] / t^q for q = 0:K-1)
        end

    for n = 0:N
        for t in range(a, b, 4)
            @assert Arblib.overlaps(S_v3(n, t), S_v4(n, t))
        end
    end

    for n = 0:N
        for t in range(a, b, 4)
            tol = a < 3 ? Arb(1e-3) : Arb(1e-5)
            @assert isapprox(S_v1(n, t), S_v4(n, t), atol = tol)
        end
    end

    #####
    # Insert sum back into h(l, t)
    #####

    h_v5(l, t) =
        sum(0:l-1) do j
            binomial(l, j) * log(t)^j * inv(t)^(l - j) * S_v4(l - j, t)
        end - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2

    for l = 1:N
        for t in range(a, b, 4)
            r1 = h_v1(l, t)
            r2 = h_v4(l, t)
            r3 = h_v5(l, t)
            @assert isapprox(r1, r3, atol = Arb(1e-5))
            @assert isapprox(r2, r3, atol = Arb(1e-5))
        end
    end

    #####
    # Switch order of sums
    #####

    h_v6(l, t) =
        let css = [S_cs_v1(n) for n = 1:l]
            sum(1:K) do m
                inv(t)^m * sum(0:l-1) do j
                    q = m - l + j
                    if 0 <= q
                        binomial(l, j) * log(t)^j * css[l-j][q+1]
                    else
                        Arb(0)
                    end
                end
            end - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2
        end

    for l = 1:N
        for t in range(a, b, 4)[2:end]
            r1 = h_v1(l, t)
            r2 = h_v5(l, t)
            r3 = h_v6(l, t)
            @assert isapprox(r1, r3, atol = Arb(1e-5))
            @assert isapprox(r2, r3, atol = Arb(1e-5))
        end
    end

    #####
    # Write coefficients for sum in different form
    #####

    # As matrix of coefficients
    # For performance reasons we precompute this for l = 1:N
    h_cs_v1 = map(1:N) do l
        let css = [S_cs_v1(n) for n = 1:l], cs = Matrix{Arb}(undef, l, K)
            for n = 1:l
                for m = 0:K-1
                    cs[n, m+1] = css[n][m+1]
                end
            end
            cs
        end
    end
    h_v7(l, t) =
        let cs = h_cs_v1[l]
            sum(1:K) do m
                coefficient = sum(max(0, l - m):l-1, init = zero(t)) do j
                    binomial(l, j) * log(t)^j * cs[l-j, m-l+j+1]
                end
                inv(t)^m * coefficient
            end - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2
        end

    for l = 1:N
        for t in range(a, b, 4)[2:end]
            r1 = h_v6(l, t)
            r2 = h_v7(l, t)
            @assert Arblib.overlaps(r1, r2)
        end
    end

    # As vector of coefficients
    h_cs_v2(l, t) =
        let css = h_cs_v1[l], cs = Vector{typeof(t)}(undef, K)
            for m = 1:K
                cs[m] = sum(max(0, l - m):l-1, init = zero(t)) do j
                    binomial(l, j) * log(t)^j * css[l-j, m-l+j+1]
                end
            end
            cs
        end

    h_v8(l, t) =
        let cs = h_cs_v2(l, t)
            sum(1:K) do m
                inv(t)^m * cs[m]
            end - l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2
        end

    for l = 1:N
        for t in range(a, b, 4)
            r1 = h_v7(l, t)
            r2 = h_v8(l, t)
            @assert Arblib.overlaps(r1, r2)
        end
    end

    # Assert that odd terms in the expansion are zero and that second
    # term matches l * log(t)^(l - 2) * (l - 1 - log(t)) / t^2
    for l = 1:N
        cs = h_cs_v2(l, b)
        @assert iszero(cs[1:2:end])
        @assert Arblib.overlaps(cs[2], l * log(b)^(l - 2) * (l - 1 - log(b)))
    end

    #####
    # Sum only even terms and start from 4
    #####

    h_v9(l, t) =
        let cs = h_cs_v2(l, t)
            sum(4:2:K) do m
                inv(t)^m * cs[m]
            end
        end

    for l = 1:N
        for t in range(a, b, 4)
            r1 = h_v8(l, t)
            r2 = h_v9(l, t)
            @assert Arblib.overlaps(r1, r2)
        end
    end

    #####
    # Insert back into G22_2_2
    #####

    H1_v2(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * h_v9(n - k, Acb(t)) * t
    end |> real

    H2_v2(n, k) = Arblib.integrate(a, b) do t
        log(t)^(k + 1) * h_v9(n - k, Acb(t)) * t
    end |> real

    H3_v2(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * h_v9(n - k, Acb(t)) * t * log(1 + c * x * t)
    end |> real

    I1_remainder_v2(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H1_v2(n, k)
        end

    I2_remainder_v2(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H2_v2(n, k)
        end

    I3_remainder_v2(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H3_v2(n, k)
        end

    if false
        # This takes some time to compute so we skip it in general
        G22_2_2_v2 = sum(2:N) do n
            (-1)^n * (1 + α)^(n - 1) / factorial(n) * (
                I1_remainder_v2(n) - inv(log(inv(x))) * I2_remainder_v2(n) +
                inv(log(inv(x))) * I3_remainder_v2(n)
            )
        end

        @show G22_2_2_v2
        @assert isapprox(G22_2_2_v1, G22_2_2_v2, rtol = Arb(1e-10))
    end

    #####
    # Start working on upper bound
    #####

    h_upper_v1(l, t) =
        let cs = h_cs_v2(l, t)
            sum(4:2:K) do m
                inv(t)^m * abs(cs[m])
            end
        end

    for l = 1:N
        for t in range(a, b, 4)
            r1 = h_v9(l, t)
            r2 = h_upper_v1(l, t)
            @assert !(r1 > r2)
        end
    end

    #####
    # Put the absolute value inside h_cs
    #####

    # As vector of coefficients
    h_cs_upper_v1(l, t) =
        let css = h_cs_v1[l], cs = Vector{typeof(t)}(undef, K)
            for m = 1:K
                cs[m] = sum(max(0, l - m):l-1, init = zero(t)) do j
                    binomial(l, j) * log(t)^j * abs(css[l-j, m-l+j+1])
                end
            end
            cs
        end

    h_upper_v2(l, t) =
        let cs = h_cs_upper_v1(l, t)
            sum(4:2:K) do m
                inv(t)^m * cs[m]
            end
        end

    for l = 1:N
        for t in range(a, b, 4)
            r1 = h_upper_v1(l, t)
            r2 = h_upper_v2(l, t)
            @assert !(r1 > r2)
        end
    end

    #####
    # Insert back into G22_2_2
    #####

    # We have to be careful here because depending on where we
    # put the absolute values we might or might not get an analytic
    # function.
    H1_upper_v1(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * h_upper_v2(n - k, Acb(t)) * t
    end |> real

    H2_upper_v1(n, k) = Arblib.integrate(a, b) do t
        log(t)^(k + 1) * h_upper_v2(n - k, Acb(t)) * t
    end |> real

    H3_upper_v1(n, k) =
        Arblib.integrate(a, b) do t
            log(t)^k * h_upper_v2(n - k, Acb(t)) * t * log(1 + c * x * t)
        end |> real

    I1_remainder_upper_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H1_upper_v1(n, k)
        end

    I2_remainder_upper_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H2_upper_v1(n, k)
        end

    I3_remainder_upper_v1(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H3_upper_v1(n, k)
        end

    if false
        # This takes some time to compute so we skip it in general
        G22_2_2_upper_v1 = sum(2:N) do n
            (1 + α)^(n - 1) / factorial(n) * (
                I1_remainder_upper_v1(n) - inv(log(inv(x))) * I2_remainder_upper_v1(n) +
                inv(log(inv(x))) * I3_remainder_upper_v1(n)
            )
        end

        @show G22_2_2_upper_v1
        @assert G22_2_2_v1 < G22_2_2_upper_v1
    else
        # Set to previous version so that further checks still work
        G22_2_2_upper_v1 = G22_2_2_v1
    end

    #####
    # We want to use that log(t)^j <= log(t)^(l - 1) to factor it out
    # all the way to h(l, t). But that only holds if t >= ℯ. We might
    # therefore want to split the integration further. For now we use
    # that log(t)^j <= 1 + log(t)^(l - 1) for t >= 1.
    #####

    h_cs_upper_v2(l) =
        let css = h_cs_v1[l], cs = Vector{Arb}(undef, K)
            for m = 1:K
                cs[m] = sum(max(0, l - m):l-1, init = zero(Arb)) do j
                    binomial(l, j) * abs(css[l-j, m-l+j+1])
                end
            end
            cs
        end

    h_upper_v3(l, t) =
        let cs = h_cs_upper_v2(l)
            #p1 = inv(t)^4 * cs[4]
            #p2 = sum(6:2:K) do m
            #    inv(t)^m * cs[m]
            #end
            #p1, p2
            (1 + log(t)^(l - 1)) * sum(4:2:K) do m
                inv(t)^m * cs[m]
            end
        end

    for l = 1:N
        for t in range(a, b, 4)
            r1 = h_upper_v2(l, t)
            r2 = h_upper_v3(l, t)
            @assert !(r1 > r2)
        end
    end

    #####
    # Factor out inv(t)^4 and use that inv(t) < inv(a) to bound the
    # rest
    #####

    h_cs_upper_v3(l) =
        let css = h_cs_v1[l], cs = Vector{Arb}(undef, K)
            for m = 1:K
                cs[m] = sum(max(0, l - m):l-1, init = zero(Arb)) do j
                    binomial(l, j) * abs(css[l-j, m-l+j+1])
                end
            end
            cs
        end

    h_upper_v4(l, t) =
        let cs = h_cs_upper_v2(l)
            (1 + log(t)^(l - 1)) * inv(t)^4 * a^4 * sum(4:2:K) do m
                cs[m] / a^m
            end
        end

    for l = 1:N
        for t in range(a, b, 4)
            r1 = h_upper_v3(l, t)
            r2 = h_upper_v4(l, t)
            @assert !(r1 > r2)
        end
    end

    #####
    # Take out part of h_upper(l, t) that doesn't depend on t
    #####

    # Precompute for performance reasons
    h_upper_constant_v1 = map(1:N) do l
        let cs = h_cs_upper_v2(l)
            a^4 * sum(4:2:K) do m
                cs[m] / a^m
            end
        end
    end

    h_upper_v5(l, t) = (1 + log(t)^(l - 1)) * inv(t)^4 * h_upper_constant_v1[l]

    for l = 1:N
        for t in range(a, b, 4)
            r1 = h_upper_v4(l, t)
            r2 = h_upper_v5(l, t)
            @assert Arblib.overlaps(r1, r2)
        end
    end

    #####
    # Insert back into G22_2_2
    #####

    H1_upper_v2(n, k) = Arblib.integrate(a, b) do t
        log(t)^k * h_upper_v5(n - k, Acb(t)) * t
    end |> real

    H2_upper_v2(n, k) = Arblib.integrate(a, b) do t
        log(t)^(k + 1) * h_upper_v5(n - k, Acb(t)) * t
    end |> real

    H3_upper_v2(n, k) =
        Arblib.integrate(a, b) do t
            log(t)^k * h_upper_v5(n - k, Acb(t)) * t * log(1 + c * x * t)
        end |> real

    for n = 1:N
        for k = 0:n-1
            r11 = H1_upper_v1(n, k)
            r12 = H2_upper_v1(n, k)
            r13 = H3_upper_v1(n, k)
            r21 = H1_upper_v2(n, k)
            r22 = H2_upper_v2(n, k)
            r23 = H3_upper_v2(n, k)
            @assert Arblib.ispositive(r21)
            @assert Arblib.ispositive(r22)
            @assert Arblib.ispositive(r23)
            @assert !(r11 > r21)
            @assert !(r12 > r22)
            @assert !(r13 > r23)
        end
    end

    I1_remainder_upper_v2(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H1_upper_v2(n, k)
        end

    I2_remainder_upper_v2(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H2_upper_v2(n, k)
        end

    I3_remainder_upper_v2(n) =
        sum(0:n-1) do k
            binomial(n, k) * γ^k * H3_upper_v2(n, k)
        end

    if false
        # This takes some time to compute so we skip it in general
        G22_2_2_upper_v2 = sum(2:N) do n
            (1 + α)^(n - 1) / factorial(n) * (
                I1_remainder_upper_v2(n) - inv(log(inv(x))) * I2_remainder_upper_v2(n) +
                inv(log(inv(x))) * I3_remainder_upper_v2(n)
            )
        end

        @show G22_2_2_upper_v2
        @assert G22_2_2_v1 < G22_2_2_upper_v2
    else
        # Set to previous version so that further checks still work
        G22_2_2_upper_v2 = G22_2_2_upper_v1
    end

    #####
    # Insert h(l, t) into H1, H2 and H3 to be able to simplify
    #####

    # Recall that if a > ℯ we don't need the 1 + in front of the inner log
    H1_upper_v3(n, k) =
        h_upper_constant_v1[n-k] * Arblib.integrate(a, b) do t
            log(t)^k * (1 + log(t)^(n - k - 1)) / t^3
        end |> real

    H2_upper_v3(n, k) =
        h_upper_constant_v1[n-k] * Arblib.integrate(a, b) do t
            log(t)^(k + 1) * (1 + log(t)^(n - k - 1)) / t^3
        end |> real

    H3_upper_v3(n, k) =
        h_upper_constant_v1[n-k] * Arblib.integrate(a, b) do t
            log(t)^k * (1 + log(t)^(n - k - 1)) / t^3 * log(1 + c * x * t)
        end |> real

    for n = 1:N
        for k = 0:n-1
            r11 = H1_upper_v2(n, k)
            r12 = H2_upper_v2(n, k)
            r13 = H3_upper_v2(n, k)
            r21 = H1_upper_v3(n, k)
            r22 = H2_upper_v3(n, k)
            r23 = H3_upper_v3(n, k)
            @assert Arblib.ispositive(r21)
            @assert Arblib.ispositive(r22)
            @assert Arblib.ispositive(r23)
            @assert Arblib.overlaps(r11, r21)
            @assert Arblib.overlaps(r12, r22)
            @assert Arblib.overlaps(r13, r23)
        end
    end

    #####
    # Explicitly compute the integrals
    #####

    H1_upper_v4(n, k) =
        h_upper_constant_v1[n-k] * (
            (gamma(Arb(k + 1), 2log(a)) - gamma(Arb(k + 1), 2log(b))) / Arb(2)^(k + 1) +
            (gamma(Arb(n), 2log(a)) - gamma(Arb(n), 2log(b))) / Arb(2)^n
        )

    H2_upper_v4(n, k) =
        h_upper_constant_v1[n-k] * (
            (gamma(Arb(k + 2), 2log(a)) - gamma(Arb(k + 2), 2log(b))) / Arb(2)^(k + 2) +
            (gamma(Arb(n + 1), 2log(a)) - gamma(Arb(n + 1), 2log(b))) / Arb(2)^(n + 1)
        )

    # In this case we factor out the bound log(1 + c * x * b) for
    # log(1 + c * x * t)
    H3_upper_v4(n, k) = log(1 + c * x * b) * H1_upper_v4(n, k)

    for n = 1:N
        for k = 0:n-1
            r11 = H1_upper_v3(n, k)
            r12 = H2_upper_v3(n, k)
            r13 = H3_upper_v3(n, k)
            r21 = H1_upper_v4(n, k)
            r22 = H2_upper_v4(n, k)
            r23 = H3_upper_v4(n, k)
            @assert Arblib.ispositive(r21)
            @assert Arblib.ispositive(r22)
            @assert Arblib.ispositive(r23)
            @assert Arblib.overlaps(r11, r21)
            @assert Arblib.overlaps(r12, r22)
            @assert !(r13 > r23)
        end
    end

    #####
    # Upper bound the explicitly computed integrals. Uses that
    # gamma(n, c) - gamma(n, d) <= gamma(n) for 0 < c < d and that k <
    # n. FIXME: We don't have factorial(k) / 2^(k + 1) < factorial(n -
    # 1) / 2^n in all cases though! So the below is not fully correct.
    # The paper does this correctly.
    #####

    H1_upper_v5(n, k) = h_upper_constant_v1[n-k] * factorial(n - 1) / Arb(2)^(n - 1)

    H2_upper_v5(n, k) = h_upper_constant_v1[n-k] * factorial(n) / Arb(2)^n

    # In this case we factor out the bound log(1 + c * x * b) for
    # log(1 + c * x * t)
    H3_upper_v5(n, k) = log(1 + c * x * b) * H1_upper_v5(n, k)

    for n = 1:N
        for k = 0:n-1
            r11 = H1_upper_v4(n, k)
            r12 = H2_upper_v4(n, k)
            r13 = H3_upper_v4(n, k)
            r21 = H1_upper_v5(n, k)
            r22 = H2_upper_v5(n, k)
            r23 = H3_upper_v5(n, k)
            @assert Arblib.ispositive(r21)
            @assert Arblib.ispositive(r22)
            @assert Arblib.ispositive(r23)
            @assert !(r11 > r21)
            @assert !(r12 > r22)
            @assert !(r13 > r23)
        end
    end

    #####
    # Insert H1, H2 and H3 into I1_remainder, I2_remainder and I3_remainder
    #####

    I1_remainder_upper_v3(n) =
        factorial(n - 1) / Arb(2)^(n - 1) * sum(0:n-1) do k
            binomial(n, k) * γ^k * h_upper_constant_v1[n-k]
        end

    I2_remainder_upper_v3(n) =
        factorial(n) / Arb(2)^n * sum(0:n-1) do k
            binomial(n, k) * γ^k * h_upper_constant_v1[n-k]
        end

    I3_remainder_upper_v3(n) = log(1 + c * x * b) * I1_remainder_upper_v3(n)

    G22_2_2_upper_v3 = sum(2:N) do n
        (1 + α)^(n - 1) / factorial(n) * (
            I1_remainder_upper_v3(n) - inv(log(inv(x))) * I2_remainder_upper_v3(n) +
            inv(log(inv(x))) * I3_remainder_upper_v3(n)
        )
    end

    @show G22_2_2_upper_v3
    @assert G22_2_2_upper_v2 < G22_2_2_upper_v3

    #####
    # Insert I1_remainder, I2_remainder and I3_remainder into H22_2_2
    # to simplify further
    #####

    G22_2_2_upper_v4 = sum(2:N) do n
        (1 + α)^(n - 1) / Arb(2)^n *
        (2(1 + inv(log(inv(x))) * log(1 + c * x * b)) / n - inv(log(inv(x)))) *
        sum(0:n-1) do k
            binomial(n, k) * γ^k * h_upper_constant_v1[n-k]
        end
    end

    @show G22_2_2_upper_v4
    @assert Arblib.overlaps(G22_2_2_upper_v3, G22_2_2_upper_v4)

    #####
    # Simplify sum further
    #####

    G22_2_2_upper_v5 =
        ((1 + inv(log(inv(x))) * log(1 + c * x * b)) - inv(log(inv(x)))) * sum(2:N) do n
            (1 + α)^(n - 1) / Arb(2)^n * sum(0:n-1) do k
                binomial(n, k) * γ^k * h_upper_constant_v1[n-k]
            end
        end

    @show G22_2_2_upper_v5
    @assert !(G22_2_2_upper_v4 > G22_2_2_upper_v5)

    #####
    # Next: Bound the sum in G22_2_2 which no longer depends on a, b or x
    #####

    G22_2_2_upper_constant_v1 = sum(2:N) do n
        (1 + α)^(n - 1) / Arb(2)^n * sum(0:n-1) do k
            binomial(n, k) * γ^k * h_upper_constant_v1[n-k]
        end
    end
end

"""
T0_asymptotic_main_2_testing_remainder_tail_sum(α::Arb = Arb(-0.9999), γ::Arb = Arb(0.5))

Method for testing the development of [`_T0_asymptotic_main_2`](@ref).
The background for that method is very long and it is easy to make
mistakes. This method tests a lot of the rewrites and simplifications
that is done for [`_T0_asymptotic_main_2`](@ref) to catch any
mistakes.

This tests the procedure for bounding what in
[`_T0_asymptotic_main_2_testing_remainder_tail`](@ref) is called
`G22_2_2_upper_constant`. We fix `a = 2`.
"""
function _T0_asymptotic_main_2_testing_remainder_tail_sum(
    α::Arb = Arb(-0.9999),
    γ::Arb = Arb(0.5),
)
    N = 50
    K = 50

    #####
    # Setup the sum we want to bound in the way it is done in
    # _T0_asymptotic_main_2_testing_remainder_tail
    #####

    S1_cs_v1 = map(1:N) do n
        let cs = Vector{Arb}(undef, K)
            cs[1] = (-1)^n
            for m = 1:K-1
                cs[m+1] =
                    inv(-Arb(m)) *
                    sum((k * n - m + k) * (-inv(Arb(k + 1))) * cs[m-k+1] for k = 1:m)
            end
            cs
        end
    end

    S_cs_v1 = map(1:N) do n
        let cs1 = S1_cs_v1[n]
            cs = 2cs1
            if isodd(n)
                for q = 0:2:K-1
                    cs[q+1] = 0
                end
            else
                for q = 1:2:K-1
                    cs[q+1] = 0
                end
            end

            cs
        end
    end

    h_cs_v1 = map(1:N) do l
        let css = [S_cs_v1[n] for n = 1:l], cs = Matrix{Arb}(undef, l, K)
            for n = 1:l
                for m = 0:K-1
                    cs[n, m+1] = css[n][m+1]
                end
            end
            cs
        end
    end

    h_cs_upper_v1 = map(1:N) do l
        let css = h_cs_v1[l], cs = Vector{Arb}(undef, K)
            for m = 1:K
                cs[m] = sum(max(0, l - m):l-1, init = zero(Arb)) do j
                    binomial(l, j) * abs(css[l-j, m-l+j+1])
                end
            end
            cs
        end
    end

    h_upper_constant_v1 = map(1:N) do l
        let cs = h_cs_upper_v1[l]
            2^4 * sum(4:2:K) do m
                cs[m] / Arb(2)^m
            end
        end
    end

    G22_2_2_upper_constant_v1 = sum(2:N) do n
        (1 + α)^(n - 1) / Arb(2)^n * sum(0:n-1) do k
            binomial(n, k) * γ^k * h_upper_constant_v1[n-k]
        end
    end

    @show G22_2_2_upper_constant_v1

    #####
    # Rewrite S1_cs and S_cs as matrices
    #####

    S1_cs_v2 = let cs = Matrix{Arb}(undef, N, K)
        for n = 1:N
            cs[n, 1] = (-1)^n
            for m = 1:K-1
                cs[n, m+1] =
                    inv(-Arb(m)) *
                    sum((k * n - m + k) * (-inv(Arb(k + 1))) * cs[n, m-k+1] for k = 1:m)
            end
        end
        cs
    end

    S_cs_v2 = let cs = 2S1_cs_v2
        for n = 1:N
            for m = 0:K-1
                if isodd(n + m)
                    cs[n, m+1] = 0
                end
            end
        end
        cs
    end

    for n = 1:N
        for m = 0:K-1
            @assert isequal(S1_cs_v1[n][m+1], S1_cs_v2[n, m+1])
            @assert isequal(S_cs_v1[n][m+1], S_cs_v2[n, m+1])
        end
    end


    #####
    # Simplify S1_cs and S_cs
    #####

    S1_cs_v3 = let cs = Matrix{Arb}(undef, N, K)
        for n = 1:N
            cs[n, 1] = (-1)^n
            for m = 2:K
                cs[n, m] =
                    sum((k * n - m + k + 1) // (k + 1) * cs[n, m-k] for k = 1:m-1) / (m - 1)
            end
        end
        cs
    end

    S_cs_v3 = let cs = 2S1_cs_v3
        for n = 1:N
            for m = 1:K
                if iseven(n + m)
                    cs[n, m] = 0
                end
            end
        end
        cs
    end

    @assert all(Arblib.overlaps.(S1_cs_v2, S1_cs_v3))
    @assert all(Arblib.overlaps.(S_cs_v2, S_cs_v3))

    #####
    # Rewrite h_cs_upper as a matrix
    #####

    h_cs_upper_v2 = let cs = similar(S_cs_v3)
        for n = 1:N
            for m = 1:K
                cs[n, m] = sum(max(0, n - m):n-1, init = zero(Arb)) do j
                    binomial(n, j) * abs(S_cs_v3[n-j, m-n+j+1])
                end
            end
        end
        cs
    end

    for n = 1:N
        for m = 0:K-1
            @assert Arblib.overlaps(h_cs_upper_v1[n][m+1], h_cs_upper_v2[n, m+1])
        end
    end

    #####
    # Rewrite h_upper_constant, factoring out 2^4 and including terms
    # in the sum which are zero.
    #####

    h_upper_constant_v2 = map(1:N) do n
        sum(1:K) do m
            if m == 2
                Arb(0) # This is the term we explicitly subtract
            else
                h_cs_upper_v2[n, m] / Arb(2)^m
            end
        end
    end

    @assert all(Arblib.overlaps.(h_upper_constant_v1, 2^4 * h_upper_constant_v2))

    #####
    # Insert back into G22_2_2_upper_constant
    #####


    G22_2_2_upper_constant_v2 =
        2^3 * sum(2:N) do n
            ((1 + α) / 2)^(n - 1) * sum(0:n-1) do k
                binomial(n, k) * γ^k * h_upper_constant_v2[n-k]
            end
        end

    @show G22_2_2_upper_constant_v2
    @assert Arblib.overlaps(G22_2_2_upper_constant_v1, G22_2_2_upper_constant_v2)

    #####
    # Check what would happen if we had h_upper_constant_v2[n] < 2^n
    #####

    all(n -> h_upper_constant_v2[n] < 2^n, 1:N) ||
        @warn "h_upper_constant_v2[n] < 2^n not satisfied"

    G22_2_2_upper_constant_test_v1 =
        2^3 * sum(2:N) do n
            ((1 + α) / 2)^(n - 1) * sum(0:n-1) do k
                binomial(n, k) * γ^k * 2^(n - k)
            end
        end

    @show G22_2_2_upper_constant_test_v1

    G22_2_2_upper_constant_test_v2 = 2^3 * sum(2:N) do n
        ((1 + α) / 2)^(n - 1) * ((2 + γ)^n - γ^n)
    end

    @show G22_2_2_upper_constant_test_v2

    G22_2_2_upper_constant_test_v3 =
        2^3 * (1 + α) / 2 * ((2 + γ)^2 * sum(0:N) do n
            ((1 + α) * (2 + γ) / 2)^n
        end - γ^2 * sum(0:N) do n
            ((1 + α) * γ / 2)^n
        end)

    @show G22_2_2_upper_constant_test_v3

    G22_2_2_upper_constant_test_v4 =
        2^3 * (1 + α) / 2 *
        ((2 + γ)^2 / (1 - (1 + α) * (2 + γ) / 2) - γ^2 / (1 - (1 + α) * γ / 2))

    @show G22_2_2_upper_constant_test_v4

    # Check that it works for α overlapping -1
    G22_2_2_upper_constant_test_v5 = let α = union(Arb(-1), α)
        2^3 * (1 + α) / 2 *
        ((2 + γ)^2 / (1 - (1 + α) * (2 + γ) / 2) - γ^2 / (1 - (1 + α) * γ / 2))
    end

    @show G22_2_2_upper_constant_test_v5

    #####
    # The above seems to be good enough! So the goal is to prove
    # h_upper_constant_v2[n] < 2^n
    #####

    #####
    # Switch notation
    #####

    D_v1 = let cs = Matrix{Arb}(undef, N, K)
        for n = 1:N
            cs[n, 1] = (-1)^n
            for m = 2:K
                cs[n, m] =
                    sum((k * n - m + k + 1) // (k + 1) * cs[n, m-k] for k = 1:m-1) / (m - 1)
            end
        end
        cs
    end

    C_v1 = let cs = 2D_v1
        for n = 1:N
            for m = 1:K
                if iseven(n + m)
                    cs[n, m] = 0
                end
            end
        end
        cs
    end

    B_v1 = let cs = similar(C_v1)
        for n = 1:N
            for m = 1:K
                cs[n, m] = sum(max(0, n - m):n-1, init = zero(Arb)) do j
                    binomial(n, j) * abs(C_v1[n-j, m-n+j+1])
                end
            end
        end
        cs
    end

    A_v1 = map(1:N) do n
        # Terms for m = 1:3 are all zero
        sum(4:K) do m
            B_v1[n, m] / Arb(2)^m
        end
    end

    @assert isequal(D_v1, S1_cs_v3)
    @assert isequal(C_v1, S_cs_v3)
    @assert isequal(B_v1, h_cs_upper_v2)
    @assert isequal(A_v1, h_upper_constant_v2)

    #return A_v1, B_v1, C_v1, D_v1

    #####
    # Simplify initial value for D. Setting it to 1 instead of (-1)^n
    # only changes the sign, we take the absolute value in the end
    # anyway.
    #####

    D_v2 = let cs = Matrix{Arb}(undef, N, K)
        for n = 1:N
            cs[n, 1] = 1
            for m = 2:K
                cs[n, m] =
                    sum((k * n - m + k + 1) // (k + 1) * cs[n, m-k] for k = 1:m-1) / (m - 1)
            end
        end
        cs
    end

    @assert isequal(abs.(D_v1), D_v2)

    C_v2 = let cs = 2D_v2
        for n = 1:N
            for m = 1:K
                if iseven(n + m)
                    cs[n, m] = 0
                end
            end
        end
        cs
    end

    B_v2 = let cs = similar(C_v1)
        for n = 1:N
            for m = 1:K
                cs[n, m] = sum(max(0, n - m):n-1, init = zero(Arb)) do j
                    binomial(n, j) * abs(C_v2[n-j, m-n+j+1])
                end
            end
        end
        cs
    end

    A_v2 = map(1:N) do n
        # Terms for m = 1:3 are all zero
        sum(4:K) do m
            B_v2[n, m] / Arb(2)^m
        end
    end

    #return A_v2, B_v2, C_v2, D_v2

    #####
    # Use that binomial(n, j) <= 2^n to factor out 2^n in B
    #####

    B_v3 = let cs = similar(C_v1)
        for n = 1:N
            for m = 1:K
                cs[n, m] = Arb(2)^n * sum(max(0, n - m):n-1, init = zero(Arb)) do j
                    abs(C_v2[n-j, m-n+j+1])
                end
            end
        end
        cs
    end

    A_v3 = map(1:N) do n
        # Terms for m = 1:3 are all zero
        sum(4:K) do m
            B_v3[n, m] / Arb(2)^m
        end
    end

    G22_2_2_upper_constant_v3 =
        2^3 * sum(2:N) do n
            ((1 + α) / 2)^(n - 1) * sum(0:n-1) do k
                binomial(n, k) * γ^k * A_v3[n-k]
            end
        end

    @show G22_2_2_upper_constant_v3
    #@assert Arblib.overlaps(G22_2_2_upper_constant_v1, G22_2_2_upper_constant_v2)

    #####
    # Factor out 2^n from B all the way
    #####

    B_v4 = let cs = similar(C_v1)
        for n = 1:N
            for m = 1:K
                cs[n, m] = sum(max(0, n - m):n-1, init = zero(Arb)) do j
                    abs(C_v2[n-j, m-n+j+1])
                end
            end
        end
        cs
    end

    A_v4 = map(1:N) do n
        # Terms for m = 1:3 are all zero
        sum(4:K) do m
            B_v4[n, m] / Arb(2)^m
        end
    end

    G22_2_2_upper_constant_v4 =
        2^3 * sum(2:N) do n
            ((1 + α) / 2)^(n - 1) * sum(0:n-1) do k
                binomial(n, k) * γ^k * Arb(2)^(n - k) * A_v4[n-k]
            end
        end

    @show G22_2_2_upper_constant_v4
    @assert Arblib.overlaps(G22_2_2_upper_constant_v3, G22_2_2_upper_constant_v4)

    #return A_v4, B_v4, C_v2, D_v2

    #####
    # Rewrite B in a recursive way
    #####

    B_v5 = let cs = similar(C_v1)
        for n = 1:N
            for m = 1:K
                if n == 1
                    cs[1, m] = abs(C_v2[1, m])
                elseif n <= m
                    cs[n, m] = abs(C_v2[n, m-n+1]) + cs[n-1, m]
                else
                    cs[n, m] = cs[n-1, m]
                end
            end
        end
        cs
    end

    @assert all(Arblib.overlaps.(B_v4, B_v5))

    #####
    # From the recursive definition it is clear that we have B[n, m]
    # <= B[m, m]. We therefore only have to work with B[m, m].
    #####

    B_m_v1 = let cs = Vector{Arb}(undef, size(C_v1, 1))
        for m = 1:min(N, K)
            cs[m] = sum(0:m-1, init = zero(Arb)) do j
                abs(C_v2[m-j, j+1])
            end
        end
        cs
    end

    @assert !any(B_v5 .> B_m_v1')

    #####
    # Give upper bound of A not depending on n.
    #####

    A_bound_v1 = sum(4:min(N, K)) do m
        B_m_v1[m] / Arb(2)^m
    end

    G22_2_2_upper_constant_v5 =
        2^3 * sum(2:N) do n
            ((1 + α) / 2)^(n - 1) * A_bound_v1 * sum(0:n-1) do k
                binomial(n, k) * γ^k * Arb(2)^(n - k)
            end
        end

    @show G22_2_2_upper_constant_v5
    @assert !(G22_2_2_upper_constant_v4 > G22_2_2_upper_constant_v5)

    #####
    # Simplify C[m-j,j+1] in B_m depending on m
    #####

    B_m_v2 = let cs = zeros(Arb, size(C_v1, 1))
        # Only even coefficients are non-zero
        for m = 2:2:min(N, K)
            cs[m] = 2sum(0:m-1, init = zero(Arb)) do j
                abs(D_v2[m-j, j+1])
            end
        end
        cs
    end

    @assert all(Arblib.overlaps.(B_m_v1, B_m_v2))

    #####
    # Remove odd terms in B_b
    #####

    B_m_even_v1 = map(2:2:min(N, K)) do m
        2sum(0:m-1, init = zero(Arb)) do j
            abs(D_v2[m-j, j+1])
        end
    end

    @assert all(Arblib.overlaps.(B_m_v2[2:2:end], B_m_even_v1))

    A_bound_v2 = sum(2:length(B_m_even_v1)) do m
        B_m_even_v1[m] / Arb(4)^(m)
    end

    @assert Arblib.overlaps(A_bound_v1, A_bound_v2)

    #####
    # Numerically B_m_even seems to behave like a^m for a ≈ 2.5 for m
    # >= 2
    #####

    @assert all(2:length(B_m_even_v1)) do m
        Arb(2)^m < B_m_even_v1[m] < Arb(3)^m
    end

    #####
    # Compute bound if we would have that B_m_even is bounded by 3^m
    # for m >= 2
    #####

    A_bound_test_v1 = Arb(9 // 4)

    G22_2_2_upper_constant_test_v6 =
        2^3 *
        A_bound_test_v1 *
        sum(2:N) do n
            ((1 + α) / 2)^(n - 1) * sum(0:n-1) do k
                binomial(n, k) * γ^k * Arb(2)^(n - k)
            end
        end

    @assert A_bound_test_v1 > A_bound_v2

    @show G22_2_2_upper_constant_test_v6
    @assert G22_2_2_upper_constant_test_v6 > G22_2_2_upper_constant_v5

    #####
    # Simplify G22_2_2 in the above case
    #####

    G22_2_2_upper_constant_test_v7 =
        2^3 * A_bound_test_v1 * sum(2:N) do n
            ((1 + α) / 2)^(n - 1) * ((2 + γ)^n - γ^n)
        end

    @assert Arblib.overlaps(G22_2_2_upper_constant_test_v6, G22_2_2_upper_constant_test_v7)

    G22_2_2_upper_constant_test_v8 =
        2^3 * A_bound_test_v1 * (1 + α) / 2 * ((2 + γ)^2 * sum(0:N) do n
            ((1 + α) * (2 + γ) / 2)^n
        end - γ^2 * sum(0:N) do n
            ((1 + α) * γ / 2)^n
        end)

    @assert Arblib.overlaps(G22_2_2_upper_constant_test_v7, G22_2_2_upper_constant_test_v8)

    G22_2_2_upper_constant_test_v9 =
        2^3 * A_bound_test_v1 * (1 + α) / 2 *
        ((2 + γ)^2 / (1 - (1 + α) * (2 + γ) / 2) - γ^2 / (1 - (1 + α) * γ / 2))

    @assert Arblib.overlaps(G22_2_2_upper_constant_test_v8, G22_2_2_upper_constant_test_v9)

    # Check that it works for α overlapping -1
    G22_2_2_upper_constant_test_v10 = let α = union(Arb(-1), α)
        2^3 * A_bound_test_v1 * (1 + α) / 2 *
        ((2 + γ)^2 / (1 - (1 + α) * (2 + γ) / 2) - γ^2 / (1 - (1 + α) * γ / 2))
    end

    @show G22_2_2_upper_constant_test_v10
    @assert Arblib.overlaps(G22_2_2_upper_constant_test_v9, G22_2_2_upper_constant_test_v10)

    #####
    # The above seems to be good enough. So the goal is to prove that
    # B_m_even is bounded by 3^m for m >= 2.
    #####

    #####
    # Numerically we see that all elements of D are positive, we
    # should hence be able to remove the absolute value in B. This
    # should be possible to prove by going back to the definition of D
    # as coefficients in Taylor expansions.
    #####

    @assert all(Arblib.ispositive, D_v2)

    B_m_even_v2 = map(2:2:min(N, K)) do m
        2sum(0:m-1, init = zero(Arb)) do j
            D_v2[m-j, j+1]
        end
    end

    @assert all(Arblib.overlaps.(B_m_even_v1, B_m_even_v2))

    #####
    # Insert definition of B into A
    #####

    A_bound_v3 = 2sum(2:N÷2) do m
        1 / Arb(4)^(m) * sum(0:2m-1, init = zero(Arb)) do j
            D_v2[2m-j, j+1]
        end
    end

    @assert Arblib.overlaps(A_bound_v2, A_bound_v3)

    #####
    # Switch up order of summation. Instead of summing along the skew
    # diagonals we sum along the rows.
    #####

    # Helper method to see which indices are accessed
    D_v3(n, m, accesses) = begin
        accesses[n, m] += 1
        D_v2[n, m]
    end

    # Old version, summing along skew diagonals
    accesses1 = zeros(Int, size(D_v2))
    A_bound_v4 = 2sum(2:N÷2) do m
        1 / Arb(4)^(m) * sum(0:2m-1, init = zero(Arb)) do j
            D_v3(2m - j, j + 1, accesses1)
        end
    end

    @assert Arblib.overlaps(A_bound_v3, A_bound_v4)

    # New version, summing along rows
    accesses2 = zeros(Int, size(D_v2))
    N == K || @warn "method probably needs N == K"
    A_bound_v5 =
        2sum(1:N) do n
            start = n <= 2 ? 2 : 1
            # We need this since the old version doesn't sum below the
            # main skew diagonal.
            stop = K + 1 - n
            if isodd(n)
                sum(D_v3(n, 2m, accesses2) / Arb(4)^m for m = start:(stop+1)÷2) /
                Arb(2)^(n - 1)
            else
                sum(D_v3(n, 2m - 1, accesses2) / Arb(4)^m for m = start:(stop+1)÷2) /
                Arb(2)^(n - 2)
            end
        end

    @assert Arblib.overlaps(A_bound_v4, A_bound_v5)

    #####
    # Slightly modify summation for A_bound
    #####

    A_bound_v6 =
        4sum(1:N) do n
            start = n <= 2 ? 2 : 1
            # We need this since the old version doesn't sum below the
            # main skew diagonal.
            stop = K + 1 - n
            if isodd(n)
                sum(D_v2[n, 2m] / Arb(2)^2m for m = start:(stop+1)÷2) / Arb(2)^n
            else
                sum(D_v2[n, 2m-1] / Arb(2)^(2m - 1) for m = start:(stop+1)÷2) / Arb(2)^n
            end
        end

    @assert Arblib.overlaps(A_bound_v5, A_bound_v6)

    #####
    # Write D in terms of expansion of log(1 - t)
    #####

    D_v4 = let cs = Matrix{Arb}(undef, N, K)
        for n = 1:N
            t = ArbSeries((0, 1), degree = K + n - 1)
            expansion = (-log(1 - t))^n << n
            cs[n, :] = Arblib.coeffs(expansion)
        end
        cs
    end

    @assert all(Arblib.overlaps.(D_v2, D_v4))

    #####
    # Give upper bound of A_bound by summing whole rows
    #####

    A_bound_upper_v1 = 4sum(1:N) do n
        sum(D_v2[n, m] / Arb(2)^m for m = 1:K) / Arb(2)^n
    end

    @assert A_bound_v6 < A_bound_upper_v1

    #####
    # Rewrite in terms of log
    #####

    A_bound_upper_v2 = 4sum(1:N) do n
        (-log(1 - inv(Arb(2))))^n / 2
    end

    @show A_bound_upper_v1 A_bound_upper_v2

    # v2 includes all terms in the m expansion
    @assert A_bound_upper_v1 < A_bound_upper_v2
    @assert isapprox(A_bound_upper_v1, A_bound_upper_v2, atol = 1e-5)

    #####
    # Explicitly compute sum
    #####

    A_bound_upper_v3 = 2 / (1 + log(1 - inv(Arb(2)))) - 2

    @show A_bound_upper_v2 A_bound_upper_v3

    # v3 includes all terms in the n expansion
    @assert A_bound_upper_v2 < A_bound_upper_v3
    @assert isapprox(A_bound_upper_v2, A_bound_upper_v3, atol = 1e-5)

    #####
    # Insert back into G22_2_2
    #####

    G22_2_2_upper_constant_v6 =
        2^3 *
        A_bound_upper_v3 *
        sum(2:N) do n
            ((1 + α) / 2)^(n - 1) * sum(0:n-1) do k
                binomial(n, k) * γ^k * Arb(2)^(n - k)
            end
        end

    @show G22_2_2_upper_constant_v6
    @assert G22_2_2_upper_constant_v5 < G22_2_2_upper_constant_v6

    G22_2_2_upper_constant_v7 =
        2^3 * A_bound_upper_v3 * (1 + α) / 2 * ((2 + γ)^2 * sum(0:N) do n
            ((1 + α) * (2 + γ) / 2)^n
        end - γ^2 * sum(0:N) do n
            ((1 + α) * γ / 2)^n
        end)

    @assert Arblib.overlaps(G22_2_2_upper_constant_v6, G22_2_2_upper_constant_v7)

    G22_2_2_upper_constant_v8 =
        2^3 * A_bound_upper_v3 * (1 + α) / 2 *
        ((2 + γ)^2 / (1 - (1 + α) * (2 + γ) / 2) - γ^2 / (1 - (1 + α) * γ / 2))

    @assert Arblib.overlaps(G22_2_2_upper_constant_v7, G22_2_2_upper_constant_v8)

    # Check that it works for α overlapping -1
    G22_2_2_upper_constant_v9 = let α = union(Arb(-1), α)
        2^3 * A_bound_upper_v3 * (1 + α) / 2 *
        ((2 + γ)^2 / (1 - (1 + α) * (2 + γ) / 2) - γ^2 / (1 - (1 + α) * γ / 2))
    end

    @show G22_2_2_upper_constant_v9
    @assert Arblib.overlaps(G22_2_2_upper_constant_v8, G22_2_2_upper_constant_v9)

    #####
    # Give final closed form formula
    #####

    G22_2_2_upper_constant_final = let α = union(Arb(-1), α)
        2^3 *
        (inv(1 + log(1 - inv(Arb(2)))) - 1) *
        (1 + α) *
        ((2 + γ)^2 / (1 - (1 + α) * (2 + γ) / 2) - γ^2 / (1 - (1 + α) * γ / 2))
    end

    @show G22_2_2_upper_constant_final
    @assert Arblib.overlaps(G22_2_2_upper_constant_v9, G22_2_2_upper_constant_final)

    return G22_2_2_upper_constant_final
end

"""
    _T0_asymptotic_main_2_testing_G22_2_1(x::Arb = Arb(1e-10), α::Arb = Arb(-0.9999), γ::Arb = Arb(0.5), c::Arb = 2Arb(ℯ))

Method for testing the development of [`_T0_asymptotic_main_2`](@ref).
The background for that method is very long and it is easy to make
mistakes. This method tests a lot of the rewrites and simplifications
that is done for [`_T0_asymptotic_main_2`](@ref) to catch any
mistakes.

This looks at evaluating `G22_2_1` from
[`_T0_asymptotic_main_2_testing_remainder`](@ref) in a way that works
for `a = 1` and `b = π / x` with `x` overlapping zero.
"""
function _T0_asymptotic_main_2_testing_G22_2_1(
    x::Arb = Arb(1e-10),
    α::Arb = Arb(-0.9999),
    γ::Arb = Arb(0.5),
    c::Arb = 2Arb(ℯ),
)

    # Enclosure of inv(log(inv(x)))
    invloginvx = if iszero(x)
        zero(x)
    elseif Arblib.contains_zero(x)
        -Arb((inv(log(ubound(Arb, x))), 0))
    else
        -inv(log(x))
    end

    # In the end we want to have a = 1 and b possibly overlapping ∞
    a = Arb(1 + 1e-15)
    b = π / x

    #####
    # Starting version
    #####

    I1_remainder_n1_primitive_v1(t) = (t^2 - 1) * (log(t - 1) + log(t + 1) - 2log(t)) / 2
    I1_remainder_n1_v1 = I1_remainder_n1_primitive_v1(b) - I1_remainder_n1_primitive_v1(a)

    I2_remainder_n1_primitive_v1(t) =
        (
            -3 - Arb(π)^2 / 3 +
            2(log(t - 1) + log(t + 1) + 2log(t)^2) +
            2t^2 * (
                -log(t - 1) - log(t + 1) +
                2log(t) +
                log(t) * (2log(t - 1) + 2log(t + 1) - 4log(t))
            ) +
            2real(polylog(2, Acb(1 - t^2)))
        ) / 8 |> real
    I2_remainder_n1_v1 = I2_remainder_n1_primitive_v1(b) - I2_remainder_n1_primitive_v1(a)

    I3_remainder_n1_v1 = Arb((log(1 + c * x * a), log(1 + c * x * b))) * I1_remainder_n1_v1

    @show I1_remainder_n1_v1 I2_remainder_n1_v1 I3_remainder_n1_v1

    G22_2_1_v1 = -(
        I1_remainder_n1_v1 - invloginvx * I2_remainder_n1_v1 +
        invloginvx * I3_remainder_n1_v1
    )

    @show G22_2_1_v1

    #####
    # Evaluating I1_remainder_n1_primitive_v1 at t = 1 gives zero
    #####

    @assert isapprox(I1_remainder_n1_primitive_v1(Arb(1) + 1e-10), 0, atol = 1e-8)
    @assert isapprox(I1_remainder_n1_primitive_v1(Arb(1) + 1e-20), 0, atol = 1e-16)

    #####
    # Simplify I1_remainder_n1
    #####

    I1_remainder_n1_v2 = I1_remainder_n1_primitive_v1(b)

    @assert isapprox(I1_remainder_n1_v1, I1_remainder_n1_v2, atol = 1e-6)

    #####
    # Rewrite I1_remainder_n1_primitive to allow for evaluation at
    # t=∞. Rewrite it in terms of inv(t). PROVE: that it is decreasing
    # in inv(t) and -0.5 for inv(t) = 0
    #####

    I1_remainder_n1_primitive_inv_t_v1(invt) =
        if iszero(invt)
            -Arb(1 // 2)
        elseif Arblib.contains_zero(invt)
            Arb((
                I1_remainder_n1_primitive_inv_t_v1(zero(invt)),
                I1_remainder_n1_primitive_inv_t_v1(ubound(Arb, invt)),
            ))
        else
            let t = inv(invt)
                (t^2 - 1) * (log(t - 1) + log(t + 1) - 2log(t)) / 2
            end
        end

    @assert Arblib.overlaps(
        I1_remainder_n1_primitive_v1(b),
        I1_remainder_n1_primitive_inv_t_v1(inv(b)),
    )
    @assert isfinite(I1_remainder_n1_primitive_inv_t_v1(inv(b)))
    @assert isapprox(
        I1_remainder_n1_primitive_inv_t_v1(Arb(0)),
        I1_remainder_n1_primitive_inv_t_v1(inv(b)),
        atol = 1e-6,
    )

    # Ensure that it can be evaluated for x overlapping zero
    let x = union(x, Arb(0))
        @show x / π
        @show I1_remainder_n1_primitive_inv_t_v1(x / π)
        @assert isfinite(I1_remainder_n1_primitive_inv_t_v1(x / π))
    end

    #####
    # Use the above to update I1_remainder_n1
    #####

    I1_remainder_n1_v3 = I1_remainder_n1_primitive_inv_t_v1(inv(b))

    @assert isapprox(I1_remainder_n1_v2, I1_remainder_n1_v3, atol = 1e-6)

    #####
    # Write final version that works for x overlapping zero
    #####

    I1_remainder_n1_final = let x = union(x, Arb(0))
        I1_remainder_n1_primitive_inv_t_v1(x / π)
    end

    I3_remainder_n1_final = let x = union(x, Arb(0))
        Arb((log(1 + c * x), log(1 + c * π))) * I1_remainder_n1_final
    end

    @show I1_remainder_n1_final
    @show I3_remainder_n1_final
    @assert isfinite(I1_remainder_n1_final)
    @assert isfinite(I3_remainder_n1_final)
    @assert isapprox(I1_remainder_n1_v1, I1_remainder_n1_final, atol = 1e-6)
    @assert isapprox(
        midpoint(I3_remainder_n1_v1),
        midpoint(I3_remainder_n1_final),
        atol = 1e-6,
    )

    #####
    # Simplify I2_remainder_n1_primitive at t = 1
    #####

    # Insert a directly
    I2_remainder_n1_primitive_one_v1 =
        (
            -3 - Arb(π)^2 / 3 +
            2(log(a - 1) + log(a + 1) + 2log(a)^2) +
            2a^2 * (
                -log(a - 1) - log(a + 1) +
                2log(a) +
                log(a) * (2log(a - 1) + 2log(a + 1) - 4log(a))
            ) +
            2real(polylog(2, Acb(1 - a^2)))
        ) / 8 |> real

    @assert Arblib.overlaps(
        I2_remainder_n1_primitive_v1(a),
        I2_remainder_n1_primitive_one_v1,
    )

    # Insert a = 1 where we get finite terms
    I2_remainder_n1_primitive_one_v2 = let aa = one(a)
        (
            -3 - Arb(π)^2 / 3 +
            2(log(a - 1) + log(aa + 1) + 2log(aa)^2) +
            2aa^2 * (
                -log(a - 1) - log(aa + 1) +
                2log(aa) +
                log(aa) * (2log(a - 1) + 2log(aa + 1) - 4log(aa))
            ) +
            2real(polylog(2, Acb(1 - aa^2)))
        ) / 8 |> real
    end

    @assert isapprox(
        I2_remainder_n1_primitive_one_v1,
        I2_remainder_n1_primitive_one_v2,
        atol = 1e-12,
    )

    # Remove terms that are zero and explicitly insert 1 in other places
    I2_remainder_n1_primitive_one_v3 =
        (-3 - Arb(π)^2 / 3 + 2(log(a - 1) + log(Arb(2))) + 2(-log(a - 1) - log(Arb(2)))) / 8

    @assert Arblib.overlaps(
        I2_remainder_n1_primitive_one_v2,
        I2_remainder_n1_primitive_one_v3,
    )

    # Simplify
    I2_remainder_n1_primitive_one_v4 = -(3 + Arb(π)^2 / 3) / 8

    @assert Arblib.overlaps(
        I2_remainder_n1_primitive_one_v3,
        I2_remainder_n1_primitive_one_v4,
    )

    #####
    # Take out problematic part of I2_remainder_n1_primitive
    #####

    I2_remainder_n1_primitive_part_v1(t) = begin
        r1 = log(t - 1) + log(t + 1) + 2log(t)^2
        r2 = -t^2 * (log(t - 1) + log(t + 1) - 2log(t))
        r3 = 2t^2 * log(t) * (log(t - 1) + log(t + 1) - 2log(t))
        r4 = real(polylog(2, Acb(1 - t^2)))
        r1 + r2 + r3 + r4
    end

    I2_remainder_n1_primitive_v2(t) = begin
        part = I2_remainder_n1_primitive_part_v1(t)

        ((-3 - Arb(π)^2 / 3) / 2 + part) / 4
    end

    @assert isfinite(I2_remainder_n1_primitive_v2(b))
    @assert Arblib.overlaps(
        I2_remainder_n1_primitive_v1(b),
        I2_remainder_n1_primitive_v2(b),
    )

    #####
    # Rewrite in terms of inv(t) to allow for evaluation at infinity
    #####

    I2_remainder_n1_primitive_part_inv_t_v1(invt) =
        let t = inv(invt)
            r1 = log(t - 1) + log(t + 1) + 2log(t)^2
            r2 = -t^2 * (log(t - 1) + log(t + 1) - 2log(t))
            r3 = 2t^2 * log(t) * (log(t - 1) + log(t + 1) - 2log(t))
            r4 = real(polylog(2, Acb(1 - t^2)))
            r1 + r2 + r3 + r4
        end

    I2_remainder_n1_primitive_inv_t_v1(invt) = begin
        part = I2_remainder_n1_primitive_part_inv_t_v1(invt)

        ((-3 - Arb(π)^2 / 3) / 2 + part) / 4
    end

    @assert isfinite(I2_remainder_n1_primitive_inv_t_v1(inv(b)))
    @assert Arblib.overlaps(
        I2_remainder_n1_primitive_v2(b),
        I2_remainder_n1_primitive_inv_t_v1(inv(b)),
    )

    #####
    # PROVE: That I3_remainder_n1_primitive_part is decreasing
    # in t and converges to 1 - π^2 / 6 at t = ∞.
    #####

    I2_remainder_n1_primitive_part_inv_t_v2(invt) =
        if iszero(invt)
            1 - Arb(π)^2 / 6
        elseif Arblib.contains_zero(invt)
            # Increasing in inv_t
            Arb((
                I2_remainder_n1_primitive_part_inv_t_v2(zero(invt)),
                I2_remainder_n1_primitive_part_inv_t_v2(ubound(Arb, invt)),
            ),)
        else
            let t = inv(invt)
                r1 = log(t - 1) + log(t + 1) + 2log(t)^2
                r2 = -t^2 * (log(t - 1) + log(t + 1) - 2log(t))
                r3 = 2t^2 * log(t) * (log(t - 1) + log(t + 1) - 2log(t))
                r4 = real(polylog(2, Acb(1 - t^2)))
                r1 + r2 + r3 + r4
            end
        end

    I2_remainder_n1_primitive_inv_t_v2(invt) = begin
        part = I2_remainder_n1_primitive_part_inv_t_v2(invt)

        ((-3 - Arb(π)^2 / 3) / 2 + part) / 4
    end

    @assert isfinite(I2_remainder_n1_primitive_inv_t_v2(inv(b)))
    @assert Arblib.overlaps(
        I2_remainder_n1_primitive_inv_t_v1(inv(b)),
        I2_remainder_n1_primitive_inv_t_v2(inv(b)),
    )

    # Limit seems correct
    @assert isapprox(
        I2_remainder_n1_primitive_inv_t_v2(Arb(1e-5)),
        I2_remainder_n1_primitive_inv_t_v2(Arb(0)),
        atol = 1e-8,
    )
    # Seems to be increasing in inv(t)
    @assert I2_remainder_n1_primitive_inv_t_v2(Arb(0)) <
            I2_remainder_n1_primitive_inv_t_v2(Arb(1e-5)) <
            I2_remainder_n1_primitive_inv_t_v2(Arb(1e-3))

    #####
    # Insert the above back into I2_remainder_n1
    #####

    I2_remainder_n1_v2 =
        I2_remainder_n1_primitive_inv_t_v2(x / π) - I2_remainder_n1_primitive_one_v4

    @assert isapprox(I2_remainder_n1_v1, I2_remainder_n1_v2, atol = 1e-5)

    #####
    # Simplify
    #####

    I2_remainder_n1_v3 = I2_remainder_n1_primitive_part_inv_t_v2(x / π) / 4

    @assert Arblib.overlaps(I2_remainder_n1_v2, I2_remainder_n1_v3)

    #####
    # Write final version that works for x overlapping zero
    #####

    I2_remainder_n1_final = let x = union(x, zero(x))
        I2_remainder_n1_primitive_part_inv_t_v2(x / π) / 4
    end

    @show I2_remainder_n1_final
    @assert isfinite(I2_remainder_n1_final)
    @assert isapprox(I2_remainder_n1_v1, I2_remainder_n1_final, atol = 1e-5)

    #####
    # Insert everything back into G22_2_1
    #####

    G22_2_1_final = -(
        I1_remainder_n1_final - invloginvx * I2_remainder_n1_final +
        invloginvx * I3_remainder_n1_final
    )

    @show G22_2_1_final
    @assert isapprox(midpoint(G22_2_1_v1), midpoint(G22_2_1_final), atol = 1e-6)
    # They should not be exactly equal but will in practice overlap
    @assert Arblib.overlaps(G22_2_1_v1, G22_2_1_final)
end


"""
    _T0_asymptotic_main_2_testing_G21(x::Arb = Arb(1e-10), α::Arb = Arb(-0.9999), γ::Arb = Arb(0.5), c::Arb = 2Arb(ℯ))

Method for testing the development of [`_T0_asymptotic_main_2`](@ref).
The background for that method is very long and it is easy to make
mistakes. This method tests a lot of the rewrites and simplifications
that is done for [`_T0_asymptotic_main_2`](@ref) to catch any
mistakes.

This looks at evaluating `G21` from
[`_T0_asymptotic_main_2_testing`](@ref) in a way that works for `α`
overlapping `-1` and `x` overlapping zero.
"""
function _T0_asymptotic_main_2_testing_G21(
    x::Arb = Arb(1e-10),
    α::Arb = Arb(-0.9999),
    γ::Arb = Arb(0.5),
    c::Arb = 2Arb(ℯ),
)

    #####
    # Starting point - same as in _T0_asymptotic_main_2_testing
    #####

    G21_v1 =
        let a = one(x),
            b = π / x,
            q0 = (1 + α) * (1 + γ),
            p0 = 1 + α + (1 + α)^2 / 2,
            D = Arb((-log(1 + c * x * b), log(1 + c * x * b)))

            (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) * (
                (D - log(x) - inv(q0)) * (a^(-q0) - b^(-q0)) +
                (log(b) * b^(-q0) - log(a) * a^(-q0))
            )
        end

    @show G21_v1

    #####
    # To make testing easier temporarily take the upper bound of D
    #####

    G21_mid_v1 =
        let a = one(x),
            b = π / x,
            q0 = (1 + α) * (1 + γ),
            p0 = 1 + α + (1 + α)^2 / 2,
            D = abs_ubound(Arb, Arb((-log(1 + c * x * b), log(1 + c * x * b))))

            (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) * (
                (D - log(x) - inv(q0)) * (a^(-q0) - b^(-q0)) +
                (log(b) * b^(-q0) - log(a) * a^(-q0))
            )
        end

    @show G21_mid_v1

    @assert Arblib.overlaps(G21_v1, G21_mid_v1)

    #####
    # Factor out (2 + α) / (1 + γ) and split other factor into part
    # with D and part without.
    #####

    G21_mid_part1_v1 = let q0 = (1 + α) * (1 + γ), p0 = 1 + α + (1 + α)^2 / 2

        (1 - (x / π)^q0) / ((1 - x^p0) * log(inv(x)))
    end

    G21_mid_part2_v1 = let q0 = (1 + α) * (1 + γ), p0 = 1 + α + (1 + α)^2 / 2
        (log(Arb(π)) * (x / π)^q0 - (1 - (x / π)^q0) / q0 - log(x)) / ((1 - x^p0) * log(inv(x)))
    end


    #####
    # Explicitly insert a = 1 and simplify D
    #####

    G21_mid_v2 =
        let b = π / x,
            q0 = (1 + α) * (1 + γ),
            p0 = 1 + α + (1 + α)^2 / 2,
            D = abs_ubound(Arb, Arb((-log(1 + c * π), log(1 + c * π))))

            (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) *
            ((D - log(x) - inv(q0)) * (1 - b^(-q0)) + log(b) * b^(-q0))
        end

    @assert Arblib.overlaps(G21_mid_v1, G21_mid_v2)

    #####
    # Insert b = π / x
    #####

    G21_mid_v3 =
        let q0 = (1 + α) * (1 + γ),
            p0 = 1 + α + (1 + α)^2 / 2,
            D = abs_ubound(Arb, Arb((-log(1 + c * π), log(1 + c * π))))

            (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) *
            ((D - log(x) - inv(q0)) * (1 - (x / π)^q0) + log(π / x) * (x / π)^q0)
        end

    @assert Arblib.overlaps(G21_mid_v2, G21_mid_v3)

    #####
    # Split terms in numerator
    #####

    G21_mid_v4 =
        let q0 = (1 + α) * (1 + γ),
            p0 = 1 + α + (1 + α)^2 / 2,
            D = abs_ubound(Arb, Arb((-log(1 + c * π), log(1 + c * π))))

            (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) * (
                (D - inv(q0)) * (1 - (x / π)^q0) - log(x) +
                log(x) * (x / π)^q0 +
                log(Arb(π)) * (x / π)^q0 - log(x) * (x / π)^q0
            )
        end

    @assert Arblib.overlaps(G21_mid_v3, G21_mid_v4)


    #####
    # Simplify numerator
    #####

    G21_mid_v5 =
        let q0 = (1 + α) * (1 + γ),
            p0 = 1 + α + (1 + α)^2 / 2,
            D = abs_ubound(Arb, Arb((-log(1 + c * π), log(1 + c * π))))

            (2 + α) / (1 + γ) / ((1 - x^p0) * log(inv(x))) *
            ((D - inv(q0)) * (1 - (x / π)^q0) - log(x) + log(Arb(π)) * (x / π)^q0)
        end

    @assert Arblib.overlaps(G21_mid_v4, G21_mid_v5)

    #####
    # Factor out (2 + α) / (1 + γ) and split other factor into part
    # with D and part without.
    #####

    G21_mid_part1_v1 = let q0 = (1 + α) * (1 + γ), p0 = 1 + α + (1 + α)^2 / 2

        (1 - (x / π)^q0) / ((1 - x^p0) * log(inv(x)))
    end

    G21_mid_part2_v1 = let q0 = (1 + α) * (1 + γ), p0 = 1 + α + (1 + α)^2 / 2
        (log(Arb(π)) * (x / π)^q0 - (1 - (x / π)^q0) / q0 - log(x)) / ((1 - x^p0) * log(inv(x)))
    end

    G21_mid_v6 = let D = abs_ubound(Arb, Arb((-log(1 + c * π), log(1 + c * π))))
        (2 + α) / (1 + γ) * (D * G21_mid_part1_v1 + G21_mid_part2_v1)
    end

    @show G21_mid_part1_v1 G21_mid_part2_v1 G21_mid_v6

    @assert Arblib.overlaps(G21_mid_v5, G21_mid_v6)

    #####
    # Write G21_mid_part1 in terms of exp
    #####

    G21_mid_part1_v2 = let q0 = (1 + α) * (1 + γ), p0 = 1 + α + (1 + α)^2 / 2

        (1 - exp(log(x / π) * q0)) / (1 - exp(log(x) * p0)) / log(inv(x))
    end

    @assert Arblib.overlaps(G21_mid_part1_v1, G21_mid_part1_v2)

    #####
    # Write in terms of t = log(x) * p0
    #####

    G21_mid_part1_v3 =
        let p0 = 1 + α + (1 + α)^2 / 2,
            t = log(x) * p0,
            C = (1 + γ) / (1 + (1 + α) / 2) * (1 - log(Arb(π)) / log(x))

            (1 - exp(C * t)) / (1 - exp(t)) / log(inv(x))
        end

    @assert Arblib.overlaps(G21_mid_part1_v2, G21_mid_part1_v3)

    #####
    # Split into two terms by adding and subtracting exp(t) to the
    # numerator
    #####

    G21_mid_part1_v4 =
        let p0 = 1 + α + (1 + α)^2 / 2,
            t = log(x) * p0,
            C = (1 + γ) / (1 + (1 + α) / 2) * (1 - log(Arb(π)) / log(x))

            (1 + exp(t) * (1 - exp((C - 1) * t)) / (1 - exp(t))) / log(inv(x))
        end

    @assert Arblib.overlaps(G21_mid_part1_v3, G21_mid_part1_v4)

    #####
    # Upper bound using that t < 0 so exp(t) < 1
    #####

    G21_mid_part1_upper_v1 =
        let p0 = 1 + α + (1 + α)^2 / 2,
            t = log(x) * p0,
            C = (1 + γ) / (1 + (1 + α) / 2) * (1 - log(Arb(π)) / log(x))

            (1 + (1 - exp((C - 1) * t)) / (1 - exp(t))) / log(inv(x))
        end

    @assert G21_mid_part1_v4 < G21_mid_part1_upper_v1

    #####
    # We want to bound the quotient using that exp((C - 1) * t) >
    # exp(t). For this to hold we must have (C - 1) * t > t. Since t
    # is negative this corresponds to C < 2.
    #####

    G21_mid_part1_upper_v2 =
        let p0 = 1 + α + (1 + α)^2 / 2,
            t = log(x) * p0,
            C = (1 + γ) / (1 + (1 + α) / 2) * (1 - log(Arb(π)) / log(x))

            @assert C < 2

            2 / log(inv(x))
        end

    @assert G21_mid_part1_upper_v1 < G21_mid_part1_upper_v2

    #####
    # However, C < 2 doesn't hold for all values of x. Since 1 + (1 +
    # α) / 2 > 1 an upper bound for C is given by (1 + γ) * (1 -
    # log(Arb(π)) / log(x)). For this to be bounded by 2 we must have
    # x > π^((1 + γ) / (1 - γ)).
    #####

    let x = 1.001 * inv(Arb(π)^((1 + γ) / (1 - γ))),
        C = (1 + γ) / (1 + (1 + α) / 2) * (1 - log(Arb(π)) / log(x))

        @assert !(C < 2)
    end

    let x = 0.999 * inv(Arb(π)^((1 + γ) / (1 - γ))),
        C = (1 + γ) / (1 + (1 + α) / 2) * (1 - log(Arb(π)) / log(x))

        @assert C < 2
    end

    #####
    # Evaluate for x overlapping zero
    #####

    G21_mid_part1_upper_final = let x = union(x, zero(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        2invloginvx
    end

    @show G21_mid_part1_upper_final
    @assert isfinite(G21_mid_part1_upper_final)
    @assert Arblib.overlaps(G21_mid_part1_upper_v1, G21_mid_part1_upper_final)

    #####
    # We want to prove that G21_mid_part2 is bounded by 1. To begin
    # with we check that this seems to be the case and that the bound
    # seems sharp
    #####

    G21_mid_part2_v2(x, α) =
        let q0 = (1 + α) * (1 + γ), p0 = 1 + α + (1 + α)^2 / 2
            (log(Arb(π)) * (x / π)^q0 - (1 - (x / π)^q0) / q0 - log(x)) /
            ((1 - x^p0) * log(inv(x)))
        end

    for x in exp.(range(log(Arb("1e-100000")), log(Arb("1e-1")), 100))
        for α in -1 .+ exp.(range(log(Arb("1e-10")), log(Arb("1e-1")), 100))
            @assert G21_mid_part2_v2(x, α) < 1
        end
    end

    @assert isapprox(G21_mid_part2_v2(Arb("1e-100000"), Arb("1e-10")), 1, atol = 1e-5)

    #####
    # Rewrite expression - multiplying with q0 and using exp
    #####

    G21_mid_part2_v3 = let q0 = (1 + α) * (1 + γ), p0 = 1 + α + (1 + α)^2 / 2
        -(
            q0 * log(Arb(π)) * exp(q0 * log(x) - q0 * log(Arb(π))) - 1 +
            exp(q0 * log(x) - q0 * log(Arb(π))) - q0 * log(x)
        ) / ((1 - exp(p0 * log(x))) * q0 * log(x))
    end

    @assert Arblib.overlaps(G21_mid_part2_v1, G21_mid_part2_v3)

    #####
    # Rewrite in terms of s = q0 * log(π), t = q0 * log(x) and h = p0
    # / q0 = (1 + (1 + α) / 2) / (1 + γ). From α ∈ [-1, 0] and x ∈ [0,
    # 1] we get q0, p0 ∈ [0, 2] and from there s ∈ [0, 2log(π)], t ∈
    # [-∞, 0] and h ∈ [1 / (1 + γ), 2]. Note that these enclosures are way
    # larger than needed, but enough for what we need.
    #####

    G21_mid_part2_v4 =
        let q0 = (1 + α) * (1 + γ),
            s = q0 * log(Arb(π)),
            t = q0 * log(x),
            h = (1 + (1 + α) / 2) / (1 + γ)

            -(s * exp(t - s) - 1 + exp(t - s) - t) / ((1 - exp(h * t)) * t)
        end

    @assert Arblib.overlaps(G21_mid_part2_v3, G21_mid_part2_v4)

    # Check enclosures of s, t and h

    s_enclosure = Arb((0, 2log(Arb(π))))
    h_enclosure = Arb((1 / (1 + γ), 2))

    for x in exp.(range(log(Arb("1e-100000")), log(Arb("1e-1")), 100))
        for α in -1 .+ exp.(range(log(Arb("1e-10")), log(Arb("1e-1")), 100))
            let q0 = (1 + α) * (1 + γ), p0 = 1 + α + (1 + α)^2 / 2
                @assert contains(s_enclosure, q0 * log(Arb(π)))
                @assert Arblib.isnegative(q0 * log(x))
                @assert contains(h_enclosure, p0 / q0)
            end
        end
    end

    #####
    # Simplify slightly
    #####

    G21_mid_part2_v5 =
        let q0 = (1 + α) * (1 + γ),
            s = q0 * log(Arb(π)),
            t = q0 * log(x),
            h = (1 + (1 + α) / 2) / (1 + γ)

            (1 + t - exp(t - s) * (1 + s)) / ((1 - exp(h * t)) * t)
        end

    @assert Arblib.overlaps(G21_mid_part2_v4, G21_mid_part2_v5)

    #####
    # We want to prove that this is bounded by 1. Since denominator is
    # negative we want to check that numerator minus denominator is
    # positive.
    #####

    G21_mid_part2_inequality_v1 =
        let q0 = (1 + α) * (1 + γ),
            s = q0 * log(Arb(π)),
            t = q0 * log(x),
            h = (1 + (1 + α) / 2) / (1 + γ)

            (1 + t - exp(t - s) * (1 + s)) - ((1 - exp(h * t)) * t)
        end

    @show G21_mid_part2_inequality_v1

    @assert G21_mid_part2_inequality_v1 > 0


    #####
    # Simplify expression
    #####

    G21_mid_part2_inequality_v2 =
        let q0 = (1 + α) * (1 + γ),
            s = q0 * log(Arb(π)),
            t = q0 * log(x),
            h = (1 + (1 + α) / 2) / (1 + γ)

            1 + t * exp(h * t) - (1 + s) * exp(t - s)
        end

    @assert Arblib.overlaps(G21_mid_part2_inequality_v1, G21_mid_part2_inequality_v2)

    #####
    # Differentiating w.r.t s we get s * e^(t - s), which is
    # clearly positive.
    #####

    G21_mid_part2_inequality_ds_v1 =
        let q0 = (1 + α) * (1 + γ),
            s = ArbSeries((q0 * log(Arb(π)), 1)),
            t = q0 * log(x),
            h = (1 + (1 + α) / 2) / (1 + γ)

            (1+t*exp(h * t)-(1+s)*exp(t - s))[1]
        end


    G21_mid_part2_inequality_ds_v2 =
        let q0 = (1 + α) * (1 + γ),
            s = q0 * log(Arb(π)),
            t = q0 * log(x),
            h = (1 + (1 + α) / 2) / (1 + γ)

            s * exp(t - s)
        end

    @assert Arblib.overlaps(G21_mid_part2_inequality_ds_v1, G21_mid_part2_inequality_ds_v2)
    @assert Arblib.ispositive(G21_mid_part2_inequality_ds_v2)

    #####
    # It follows that the left hand side for the inequality is lower
    # bounded by the value for s = 0
    #####

    G21_mid_part2_inequality_v3 =
        let q0 = (1 + α) * (1 + γ), t = q0 * log(x), h = (1 + (1 + α) / 2) / (1 + γ)

            1 + t * exp(h * t) - exp(t)
        end

    @show G21_mid_part2_inequality_v3
    @assert G21_mid_part2_inequality_v2 > G21_mid_part2_inequality_v3 > 0

    #####
    # Differentiating w.r.t. t we get exp(t) * ((1 + h * t) * exp((h -
    # 1) * t) - ). We want to check that this is negative.
    #####

    G21_mid_part2_inequality_dt_v1 =
        let q0 = (1 + α) * (1 + γ),
            t = ArbSeries((q0 * log(x), 1)),
            h = (1 + (1 + α) / 2) / (1 + γ)

            (1+t*exp(h * t)-exp(t))[1]
        end


    G21_mid_part2_inequality_dt_v2 =
        let q0 = (1 + α) * (1 + γ), t = q0 * log(x), h = (1 + (1 + α) / 2) / (1 + γ)

            exp(t) * ((1 + h * t) * exp((h - 1) * t) - 1)
        end

    @assert Arblib.overlaps(G21_mid_part2_inequality_dt_v1, G21_mid_part2_inequality_dt_v2)
    @assert Arblib.isnegative(G21_mid_part2_inequality_dt_v2) # Want to show this

    #####
    # If we differentiate the second factor again w.r.t. t we get (2h
    # - 1 + (h - 1) * h * t) * exp((h - 1) * t)
    #####

    G21_mid_part2_inequality_dtt_v1 =
        let q0 = (1 + α) * (1 + γ),
            t = ArbSeries((q0 * log(x), 1)),
            h = (1 + (1 + α) / 2) / (1 + γ)

            ((1+h*t)*exp((h - 1) * t)-1)[1]
        end


    G21_mid_part2_inequality_dtt_v2 =
        let q0 = (1 + α) * (1 + γ), t = q0 * log(x), h = (1 + (1 + α) / 2) / (1 + γ)

            (2h - 1 + h * (h - 1) * t) * exp((h - 1) * t)
        end

    @show G21_mid_part2_inequality_dtt_v1 G21_mid_part2_inequality_dtt_v2
    @assert Arblib.overlaps(
        G21_mid_part2_inequality_dtt_v1,
        G21_mid_part2_inequality_dtt_v2,
    )

    #####
    # We have 2h - 1 = (2 + α - γ) / (1 + γ), which is positive
    #####

    let q0 = (1 + α) * (1 + γ), t = q0 * log(x), h = (1 + (1 + α) / 2) / (1 + γ)

        @assert Arblib.overlaps(2h - 1, (2 + α - γ) / (1 + γ))
        @assert Arblib.ispositive((2 + α - γ) / (1 + γ))
    end

    #####
    # We have h - 1 = ((1 + α) / 2 - γ) / (1 + γ), which is negative
    # if α < 2γ - 1. In particular this always holds for γ = 1 / 2.
    #####


    let q0 = (1 + α) * (1 + γ), t = q0 * log(x), h = (1 + (1 + α) / 2) / (1 + γ)

        @assert Arblib.overlaps(h - 1, ((1 + α) / 2 - γ) / (1 + γ))
        @assert Arblib.isnegative(((1 + α) / 2 - γ) / (1 + γ))
    end

    #####
    # Since t < 0 and h > 0 it follows from the above that 2h - 1 + h
    # * (h - 1) * t is positive. Hence G21_mid_part2_inequality_dt is
    # increasing in t and takes it maximum value at t = 0, where it is
    # zero.
    #####

    G21_mid_part2_inequality_dt_v3 =
        let q0 = (1 + α) * (1 + γ), t = Arb(0), h = (1 + (1 + α) / 2) / (1 + γ)

            exp(t) * ((1 + h * t) * exp((h - 1) * t) - 1)
        end

    @assert iszero(G21_mid_part2_inequality_dt_v3)
    @assert G21_mid_part2_inequality_dt_v2 < G21_mid_part2_inequality_dt_v3

    #####
    # It follows that G21_mid_part2_inequality is decreasing in t and
    # attains attains it minimum at t = 0, which is zero. Hence it is
    # positive for t < 0.
    #####

    G21_mid_part2_inequality_v4 = let q0 = (1 + α) * (1 + γ), t = 0
        h = (1 + (1 + α) / 2) / (1 + γ)

        1 + t * exp(h * t) - exp(t)
    end

    @assert iszero(G21_mid_part2_inequality_v4)
    @assert G21_mid_part2_inequality_v3 > G21_mid_part2_inequality_v4

    #####
    # From the above it follows that G21_mid_part2_upper is upper
    # bounded by 1.
    #####

    let q0 = (1 + α) * (1 + γ),
        s = q0 * log(Arb(π)),
        t = q0 * log(x),
        h = (1 + (1 + α) / 2) / (1 + γ)

        @assert t * exp(h * t) < 0
        @assert -(1 + s) * exp(t - s) < 0
    end

    G21_mid_part2_upper_final = one(x)

    @assert G21_mid_part2_v5 < G21_mid_part2_upper_final

    #####
    # Insert back into G21
    #####

    G21_mid_final = let D = abs_ubound(Arb, Arb((-log(1 + c * π), log(1 + c * π))))
        (2 + α) / (1 + γ) * (D * G21_mid_part1_upper_final + G21_mid_part2_upper_final)
    end

    @assert G21_mid_v6 < G21_mid_final

    G21_final = let D = Arb((-log(1 + c * π), log(1 + c * π)))
        (2 + α) / (1 + γ) * (D * G21_mid_part1_upper_final + G21_mid_part2_upper_final)
    end

    @assert !(G21_v1 > G21_final)

    @show G21_final
end

"""
    _T0_asymptotic_main_2_1(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G21(x) = 1 / log(inv(x)) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
where the integration is taken from `1` to `2`.

For ``t ∈ [1, 2]`` we have
```
log(c + inv(2x)) <= log(c + inv(x * t)) <= log(c + inv(x))
```
allowing us to factor out `log(c + inv(x * t))`, we get
```
G21(x) = log(c + inv(x * Arb((1, 2)))) / log(inv(x)) *
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

**IMPROVE:** This uses slightly different notation compared to the
paper. Partially because it handles general `γ`.
"""
function _T0_asymptotic_main_2_1(α::Arb, γ::Arb, c::Arb)
    # Primitive function of
    # ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) * t
    primitive = let s = -(1 + α)
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

        return logfactor * I
    end
end

"""
    _T0_asymptotic_main_2_2(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G22(x) = 1 / log(inv(x)) *
            ∫ ((t - 1)^(-α - 1) + (t + 1)^(-α - 1) - 2t^(-α - 1)) / (1 + α) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
where the integration is taken from `2` to `π / x`.

This method is meant for non-zero `x`. The integral is computed by
integrating numerically.

To compute better enclosures in the numerical integration we compute a
lower bound by integrating to the lower bound of `π / x` and an upper
bound by integrating to the upper bound of it. We also use that
```
log(c + inv(ubound(x) * t)) <= log(c + inv(x * t)) <= log(c + inv(lbound(x) * t))
```
"""
function _T0_asymptotic_main_2_2(α::Arb, γ::Arb, c::Arb)
    return x::Arb -> begin
        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
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

        return invloginvx * I
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
