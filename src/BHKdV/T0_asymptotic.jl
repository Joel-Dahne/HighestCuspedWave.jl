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
