"""
    _T0_asymptotic_main_1(α::Arb, γ::Arb, c::Arb)

Compute an upper bound of
```
G1(x) = inv((1 - x^p0) * log(inv(x))) *
            ∫ abs((1 - t)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - γ * (1 + α)) * log(c + inv(x * t)) dt
```
where the integration is taken from `0` to `1`, defined in the
asymptotic version of [`T0`](@ref).

This corresponds to `G1` in [`lemma_bhkdv_T0_asymptotic_split`](@ref)
and a bound is given in [`lemma_bhkdv_U0_G1`](@ref). The lemma only
handles the case `γ = 1 / 2` and `c = 2ℯ`, but it is easily
adapted to other values. We give the general proof below.

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
the removable singularity. What remains to handle is
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
                    # Note that t should always be nonnegative in this case
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
        Arblib.ispositive(fx_div_x(s -> (1 - b)^-s + (1 + b)^-s - 2b^-s, -αp1)) ||
        error("integrand must be positive on [b, 1]")

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
                    # Note that t should always be nonnegative in this case
                    let a = Arb(1 // 2), tᵤ = ubound(Arb, t)
                        # Check for monotonicity and return an
                        # indeterminate result otherwise
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
        Arblib.ispositive(fx_div_x(s -> (1 - b)^-s + (1 + b)^-s - 2b^-s, -αp1)) ||
        error("integrand must be positive on [b, 1]")

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
            ∫ ((t - 1)^(-α - 1) + (1 + t)^(-α - 1) - 2t^(-α - 1)) *
                t^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x * t)) dt
```
where the integration is taken from `1` to `π / x`, defined in the
asymptotic version of [`T0`](@ref). This corresponds to `G2` in
[`lemma_bhkdv_T0_asymptotic_split`](@ref).

If `x` does not overlap zero it uses [`_T0_asymptotic_main_2_1`](@ref)
and [`_T0_asymptotic_main_2_2`](@ref) to compute the integral. In that
case we factor out `1 + α` from the integral to get
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
bound. It is split into the lemmas
- [`lemma_bhkdv_U0_G2_split`](@ref)
- [`lemma_bhkdv_U0_G2_M`](@ref)
- [`lemma_bhkdv_U0_G2_R_n1`](@ref)
- [`lemma_bhkdv_U0_G2_R_1`](@ref)
- [`lemma_bhkdv_U0_G2_R_2`](@ref)
We here give an overview of the approach used. Note that the paper
fixes `γ = 1 / 2` and below we also use that assumption. The paper
also fixes `c = 2ℯ` but it is straight forward to adapt the results to
a general `c` and that is the version we give below.

# Split into `G2_M`, `G2_R_1` and `G2_R_2`
By [`lemma_bhkdv_U0_G2_split`](@ref) we have
```
G2(x) = G2_M(x) + G2_R(x)
```
with
```
G2_M(x) = 1 / (1 - x^(1 + α + (1 + α)^2 / 2)) * 1 / (log(1 / x)) * (4 + 2α) / 3 *
    (
        (D(x) - log(x) - 2 / 3 * inv(1 + α)) * (1 - (x / π)^(3 / 2 * (1 + α)))
        + log(π / x) * (x / π)^(3 / 2 * (1 + α))
    )
```
and
```
G2_R(x) = 1 / (1 - x^(1 + α + (1 + α)^2 / 2)) * 1 / log(1 / x) * sum(1:Inf) do n
    (-1)^n * (1 + α)^n / factorial(n) * sum(0:n-1) do k
        binomial(n, k) / 2^k * ∫_1^(π / x) log(t)^k * h(n - k, t) * t * log(c + 1 / (x * t)) dt
    end
end
```
for some `D(x)` satisfying `-log(1 + c * π) <= D(x) <= log(1 + c * π)`
and
```
h(k, t) = log(t - 1)^k + log(t + 1)^k - 2log(t)^k - k * (k - 1 - log(t)) * log(t)^(k - 2) / t^2
```

# Bounding `G2_M(x)`
By [`lemma_bhkdv_U0_G2_M`](@ref) we have the following bound for `G2_M(x)`
```
G2_M(x) <= (4 + 2α) / 3 * (2log(1 + c * π) / log(inv(x)) + 1)
```
Which holds for `x < inv(π^3)`. The factor `log(1 + c * π)` comes from
the upper bound of `D(x)`. In the paper we only give the upper bound
but for the code we instead opt to use the full enclosure of `D(x)`,
this makes it easier to track how good of an approximation the result
is.

# Bounding `G2_R(x)`
The function `G2_R(x)` can be bounded as
```
G2_R(x) <= G2_R_factor(x) * (G2_R_n1(x) + G2_R_1(x) + G2_R_2(x))
```
with
```
G2_R_factor(x) = (1 + α) / (1 - x^(1 + α + (1 + α)^2 / 2)) * (1 + log(1 + c * x) / log(1 / x))
```
```
G2_R_n1 = ∫_1^(π / x) abs(h(1, t)) * t dt
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

## Enclosing `G2_R_factor(x)`
We can enclose `(1 + α) / (1 - x^(1 + α + (1 + α)^2 / 2))` using that
it is increasing in `x`. It is equal to `1 + α` for `x = 0` and for `x
> 0` we can handle the removable singularity.

## Bounding `G2_R_n1(x)`
From [`lemma_bhkdv_U0_G2_R_n1`](@ref) we get that this term is bounded
by `1 / 2`.

## Bounding `G2_R_1(x)`
From [`lemma_bhkdv_U0_G2_R_1`](@ref) we get the following bound
```
G2_R_1(x) <= 2(
    sqrt(ℯ) * (1 + α) / (-α)
    + 4(exp(3(1 + α)) - 3(1 + α) - 1) / 3(1 + α)
    + (2 + α) * exp(3 / 2 * (1 + α))
    - 1
)
```

# Bounding `G2_R_2(x)`
From [`lemma_bhkdv_U0_G2_R_2`](@ref) we get the following bound
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
            # Given by (1 + α) / ((1 - x^(1 + α + (1 + α)^2 / 2))) * (1 + log(1 + c * x) / log(1 / x))
            G2_R_factor = begin
                lower = 1 + α
                upper = let xᵤ = ubound(Arb, x)
                    inv(fx_div_x(αp1 -> 1 - xᵤ^(αp1 + αp1^2 / 2), α + 1, extra_degree = 2))
                end
                Arb((lower, upper)) * (1 + log(1 + c * x) * invloginvx)
            end

            # First term in sum of second factor of remainder term, n = 1
            G2_R_n1 = Arb(1 // 2)

            # Remaining terms in second factor of remainder terms
            # integrated from 1 to 2
            G2_R_1 = begin
                # s = (exp(3(1 + α)) - 3(1 + α) - 1) / 3(1 + α)
                s = fx_div_x(3(1 + α)) do t
                    exp(t) - t - 1
                end

                2(sqrt(Arb(ℯ)) * (1 + α) / (-α) + 9s + (2 + α) * exp(3(1 + α) / 2) - 1)
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


By [`lemma_bhkdv_U0_G21`](@ref) we have for `γ = 1 / 2` and `c = 2ℯ`
```
G2_1(x) = C / log(1 / x) * (1 + 2^-α - 3^-α + α * (-4 + 5 * 2^-α - 2 * 3^-α)) / ((α - 1) * α * (α + 1))
```
for some `C` satisfying
```
2^(-(1 + α) / 2) * log(2ℯ + 1 / 2x) < C < log(2ℯ + 1 / x)
```
For general `γ` and `c` the only change is the value of `C`, and we
get
```
2^(-γ * (1 + α)) * log(c + 1 / 2x) < C < log(c + 1 / x)
```

To enclose `C / log(1 / x)` we note that
```
2^(-γ * (1 + α)) * log(c + 1 / 2x) / log(inv(x)) < C / log(inv(x)) < log(c + 1 / x) / log(inv(x))
```
and
```
log(c + 1 / 2x) / log(inv(x)) = 1 + log(1 / 2 + c * x) / log(1 / x)

log(c + 1 / x) / log(inv(x)) = 1 + log(1 + c * x) / log(1 / x)
```
"""
function _T0_asymptotic_main_2_1(α::Arb, γ::Arb, c::Arb)
    # Enclosure of
    # (1 + 2^-α - 3^-α + α * (-4 + 5 * 2^-α - 2 * 3^-α)) / ((α - 1) * α * (α + 1))
    I = fx_div_x(α + 1) do αp1
        let α = αp1 - 1
            (1 + 2^-α - 3^-α + α * (-4 + 5 * 2^-α - 2 * 3^-α)) / ((α - 1) * α)
        end
    end

    return x::Arb -> begin
        # Enclosure of inv(log(inv(x))) = -inv(log(x))
        invloginvx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            -Arb((inv(log(ubound(Arb, x))), 0))
        else
            -inv(log(x))
        end

        # C / log(inv(x))
        C_div_log_lower = 2^(-γ * (1 + α)) * (1 + log(1 // 2 + c * x) * invloginvx)
        C_div_log_upper = 1 + log(1 + c * x) * invloginvx
        C_div_log = Arb((C_div_log_lower, C_div_log_upper))

        return C_div_log * I
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

The integral is computed by integrating numerically.

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

        # Compute ((t - 1)^-(α + 1) + (t + 1)^-(α + 1) - 2t^-(α + 1)) / (α + 1).
        # It only handles the types of α and t that actually occurs in
        # the integration.
        integrand_part =
            (α, t) -> let
                if α isa ArbSeries && contains(α[0], -1)
                    return fx_div_x(α + 1; extra_degree) do αp1
                        (t - 1)^-αp1 + (t + 1)^-αp1 - 2t^-αp1
                    end
                elseif α isa Arb && contains(α, -1)
                    return fx_div_x(oftype(t, 1 + α)) do αp1
                        (t - 1)^-αp1 + (t + 1)^-αp1 - 2t^-αp1
                    end
                else
                    # Large cancellations so do the computations at a
                    # higher precision
                    return let αp1 = setprecision(α + 1, 2precision(α))
                        ((t - 1)^-αp1 + (t + 1)^-αp1 - 2t^-αp1) / αp1
                    end
                end
            end

        I_lower =
            Arblib.integrate(2, lbound(Arb, π / x), warn_on_no_convergence = false) do t
                if isreal(t)
                    t = real(t)

                    return ArbExtras.enclosure_series(α, degree = 4) do α
                        integrand_part(α, t) * t^(1 - γ * (α + 1))
                    end * log(c + inv(xᵤ * t))
                else
                    return integrand_part(α, Acb(t)) *
                           t^(1 - γ * (α + 1)) *
                           log(c + inv(xᵤ * t))
                end
            end |> real

        I_upper =
            Arblib.integrate(2, ubound(π / x), warn_on_no_convergence = false) do t
                if isreal(t)
                    t = real(t)

                    return ArbExtras.enclosure_series(α, degree = 4) do α
                        integrand_part(α, t) * t^(1 - γ * (α + 1))
                    end * log(c + inv(xₗ * t))
                else
                    return integrand_part(α, Acb(t)) *
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
R(x) = 2sum(1:Inf) do m
    (-1)^m * zeta(-α - 2m) * x^2m / factorial(2m) * sum(0:m-1) do k
        binomial(2m, 2k) * ∫_0^(π / x) t^(2k + 1 - γ * (1 + α)) * log(c + inv(x * t)) dt
    end
end
```
defined in the asymptotic version of [`T0`](@ref), see also
[`lemma_bhkdv_T0_asymptotic_split`](@ref).

To compute the integral we first split the log-factor as
```
log(c + inv(x * t)) = log(1 + c * x * t) - log(x * t)
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
and
```
∫ t^(2k + 1 - γ * (1 + α)) * log(x * t) dt =
    ((2(k + 1) - γ * (1 + α)) * log(π) - 1) * (π / x)^(2(k + 1) - γ * (1 + α)) / (2(k + 1) - γ * (1 + α))^2
    <= log(π) * (π / x)^(2(k + 1)) / (2(k + 1) - γ * (1 + α))
```
Hence
```
∫ t^(2k + 1 - γ * (1 + α)) * log(c + inv(x * t)) dt <=
    log(c + inv(π)) / (2(k + 1) - γ * (1 + α)) * (π / x)^(2(k + 1))
```
Inserting this back into `R(x)` and factoring out `(π / x)^2m` we get
```
R(x) <= 2log(c + inv(π)) * sum(1:Inf) do m
    (-1)^m * zeta(-α - 2m) * π^2m / factorial(2m) * sum(0:m-1) do k
        binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α))
    end
end
```
We split this sum into `m = 1:M-1` and `m = M:Inf`. For the finite
part we sum directly. For the infinite tail we note that
```
sum(0:m-1) do k
    binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α))
end <=
    sum(binomial(2m, 2k) * (1 / π)^(2(m - 1 - k)) / (2k + 1) for k = 0:m-1) <=
    sum(binomial(2m, 2k) * (1 / π)^(2(m - 1 - k)) / (2k + 1) for k = 0:m) =
    π^(2 - 2m) * ((π - 1)^2m + (π + 1)^2m) / 2 <=
    π^(2 - 2m) * (π + 1)^2m <=
    = π^2 * (1 + inv(π))^2m
```
Giving us
```
sum(M:Inf) do m
    (-1)^m * zeta(-α - 2m) * π^2m / factorial(2m) * sum(0:m-1) do k
        binomial(2m, 2k) * (x / π)^(2(m - 1 - k)) / (2(k + 1) - γ * (1 + α))
    end
end <= π^2 * sum((-1)^m * zeta(-α - 2m) * (π * (1 + inv(π)))^2m / factorial(2m) for m = M:Inf)
```
This sum is exactly the same as in the remainder for `clausenc` in
[`clausenc_expansion_remainder`](@ref), substituting `x` for `π * (1 +
inv(π))`, and we can get a bound from that function.
"""
function _T0_asymptotic_remainder(α::Arb, γ::Arb, c::Arb; M = 20)
    # Precompute the factors
    # (-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m))
    # in the finite sum
    factors = [(-1)^m * zeta(-α - 2m) * Arb(π)^2m / factorial(big(2m)) for m = 1:M-1]

    # Compute tail of sum, which is the same as for clausenc_expansion
    # It is valid for x < 1
    tail = let y = π * (1 + inv(Arb(π)))
        Arb(π)^2 * clausenc_expansion_remainder(y, -α, M) * y^2M
    end

    return x::Arb -> begin
        x < 1 || throw(DomainError(x, "must have x < 1"))

        # Sum directly for m = 1:M-1
        main = sum(1:M-1, init = zero(x)) do m
            factors[m] * sum(0:m-1) do k
                binomial(2m, 2k) * abspow(x / π, 2(m - 1 - k)) / (2(k + 1) - γ * (1 + α))
            end
        end

        return 2log(c + inv(Arb(π))) * (main + tail)
    end
end
