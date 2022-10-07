# PROVE: That all occurrences of force = true indeed are zero

"""
    eval_expansion(u0::BHKdVAnsatz, expansion, x)

Evaluate the given expansion.

It requires that `0 <= x < 1`, any negative parts of `x` are ignored.

The terms are stored as `((p, q, i, j, k, l, m), y)`. The parameters `(i, j,
k l, m)` correspond to the term
```
y * x^(i * α + j * p0 - k*u0.v0.α + l*u0.v0.p0 + m)
```
where `α ∈ (-1, -1 + u0.ϵ]` and `p0 = 1 + α + (1 + α)^2 / 2`.

The parameter `p` corresponds to multiplication by the factor
```
a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) * x^-α
```
to the power `p`, which is the part of the expansion for the main term
which is not even powers of `x`. The parameter `q` corresponds to
multiplication by the factor
```
-a0 * (
    gamma(2α) * cospi(α) * x^(-2α) -
    gamma(2α - p0) * cospi((2α - p0) / 2) * x^(-2α + p0) +
    (-zeta(1 - 2α - 2) / 2 + zeta(1 - 2α + p0 - 2) / 2) * x^2
)
```
to the power `q`, which is the part of the expansion for `H` applied
to the main term which is not on the form `x^2m` for `m >= 2`.

# Handling of `p` and `q`
The method currently only support `q = 0` and `p == 0` or `p == 1`.
The other cases are handled specially in [`F0`](@ref) and for that
reason we don't bother implementing them here.

For `x != 0` the factor
```
a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0)
```
has a removable singularity at `α = -1`. To handle this we rewrite it
as
```
a0 * (α + 1) * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) / (α + 1)
```
and use [`finda0αp1`](@ref) to enclose `a0 * (α + 1)`. To evaluate
`gamma(α) * cospi(α / 2)` and `gamma(α - p0) * cospi((α - p0) / 2)` we
rewrite them as
```
gamma(α) * cospi(α / 2) = (cospi(α / 2) / (α + 1)) / (rgamma(α) / (α + 1))
gamma(α - p0) * cospi((α - p0) / 2) = (cospi((α - p0) / 2) / (α + 1)^2) / (rgamma(α - p0) / (α + 1)^2)
```
and use [`fx_div_x`](@ref). We can then evaluate
```
(gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) / (α + 1)
```
using [`fx_div_x`](@ref) as well.

For `x = 0` this doesn't work. For now we use that it converges to
```
(1 - γ - log(x)) / π
```
and that this gives an upper bound.
- **PROVE:** That `(1 - γ - log(x)) / π` gives an upper bound.
  Alternatively we don't use this method for `x = 0` (using
  [`inv_u0_bound`](@ref) instead).

The argument `use_approx_p_and_q` is used for computing approximate
values and is mainly intended for testing.

# Terms not depending on `p` or `q`
The terms for which `p == q == 0` we take out and handle in a separate
sum. We compute an expansion in `α` to get a better enclosure, in most
cases this picks up the monotonicity.
"""
function eval_expansion(
    u0::BHKdVAnsatz{Arb},
    expansion::AbstractDict{NTuple{7,Int},Arb},
    x::Arb;
    use_approx_p_and_q = false,
)
    @assert x < 1

    # We only care about the non-negative part of x
    x = Arblib.nonnegative_part!(zero(x), x)

    # Enclosure of α, α + 1 and a0 * (α + 1)
    α = Arb((-1, -1 + u0.ϵ))
    αp1 = Arblib.nonnegative_part!(zero(u0.ϵ), union(zero(u0.ϵ), u0.ϵ))
    a0αp1 = finda0αp1(α)

    # In-place methods for computing the exponent
    # i * α + j * p0 - k*u0.v0.α + l*u0.v0.p0 + m
    # for both α::Arb and α::ArbSeries

    # Buffers used when computing the exponent for α::Arb
    lower = zero(α)
    upper = zero(α)
    # Precomputed values used for α::Arb: upper bounds of α and p0
    α_upper = -1 + u0.ϵ
    p0_upper = u0.ϵ + u0.ϵ^2 / 2
    _exponent!(exponent::Arb, α::Arb, i, j, k, l, m) = begin
        # Compute the part i * α + j * p0 + m
        # Note that i * α + j * p0 = 3j / 2 + (i + 2j) * α + j * α^2 / 2
        # is increasing in α if i + 2j is non-negative.
        if i + 2j >= 0 && iswide(α)
            # It is increasing in α, evaluated at endpoints
            # We use the full interval of α and not just for the
            # argument given here. In practice this will always be the
            # same, but even if it isn't we still get an enclosure.

            # Lower bound at α = -1 can be done with integers
            # For α = -1 we get
            # i * α + j * (1 + α + (1 + α)^2 // 2) = m - 1
            Arblib.set!(lower, m - i)

            # Upper bound i * α_upper + j * p0_upper + m
            Arblib.set!(upper, m)
            Arblib.addmul!(upper, α_upper, i)
            Arblib.addmul!(upper, p0_upper, j)

            Arblib.union!(exponent, lower, upper)

            # If the lower bound is zero we want to avoid any spurious
            # negative parts
            iszero(lower) && Arblib.nonnegative_part!(exponent, exponent)
        else
            # Note that i * α + j * p0 + m =
            # (i + j) * α + j * (1 + (1 + α)^2 / 2) + m

            # exponent = j * (1 + (1 + α)^2 / 2)
            Arblib.add!(exponent, α, 1)
            Arblib.sqr!(exponent, exponent)
            Arblib.mul_2exp!(exponent, exponent, -1)
            Arblib.add!(exponent, exponent, 1)
            Arblib.mul!(exponent, exponent, j)

            # exponent += (i + j) * α + m
            Arblib.addmul!(exponent, α, i + j)
            Arblib.add!(exponent, exponent, m)
        end

        # Add - k*u0.v0.α + l*u0.v0.p0
        Arblib.submul!(exponent, u0.v0.α, k)
        Arblib.addmul!(exponent, u0.v0.p0, l)

        return exponent
    end

    # IMPROVE: Consider using more buffers and precomputing more
    # values
    _exponent!(exponent::ArbSeries, α::ArbSeries, i, j, k, l, m) = begin
        # Compute the part i * α + j * p0 + m

        # Note that i * α + j * p0 + m = (i + j) * α + j * (1 + (1 + α)^2 / 2) + m

        # exponent = j * (1 + (1 + α)^2 / 2)
        Arblib.add!(exponent, α, 1)
        Arblib.pow_arb_series!(exponent, exponent, Arb(2), length(exponent))
        Arblib.mul_2exp!(exponent, exponent, -1)
        Arblib.add!(exponent, exponent, 1)
        Arblib.mul!(exponent, exponent, Arb(j))

        # exponent += (i * j) * α
        Arblib.add!(exponent, exponent, (i + j) * α)

        exponent0 = Arblib.ref(exponent, 0)
        Arblib.submul!(exponent0, u0.v0.α, k)
        Arblib.addmul!(exponent0, u0.v0.p0, l)
        Arblib.add!(exponent0, exponent0, m)

        return exponent
    end

    # Function for computing sum of terms with p and q zero for a given α
    # It is only called with α::ArbSeries if x is non-zero
    logx = log(x)
    f(α) = begin
        S = zero(α)
        exponent = zero(α)
        term = zero(α)
        buffer1 = zero(Arb)
        buffer2 = zero(α)
        for ((p, q, i, j, k, l, m), y) in expansion
            if !iszero(y) && iszero(p) && iszero(q)
                _exponent!(exponent, α, i, j, k, l, m)
                if α isa Arb
                    term = abspow!(term, x, exponent)
                else
                    Arblib.mul!(term, exponent, logx)
                    Arblib.exp_series!(term, term, length(term))
                end
                Arblib.mul!(term, term, y)
                Arblib.add!(S, S, term)
            end
        end
        S
    end

    # Sum of terms with p and q zero
    res1 = if Arblib.contains_zero(x)
        f(α)
    else
        ArbExtras.enclosure_series(f, α)
    end

    # Compute enclosure of the p-coefficient
    # a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0)
    # Note that this depends on x and is only used when x doesn't contain zero
    p_coefficient = if !Arblib.contains_zero(x)
        # Enclosure of
        # (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) / (α + 1)
        extra_degree = 2
        coefficient_div_α = fx_div_x(αp1, force = true; extra_degree) do r
            if Arblib.contains_zero(r[0])
                # Enclosure of rgamma(α) / (α + 1)
                rgamma1_div_α = fx_div_x(s -> rgamma(s - 1), r; extra_degree)
                # Enclosure of rgamma(α - p0) / (α + 1)^2
                rgamma2_div_α2 = fx_div_x(s -> rgamma(-1 - s^2 / 2), r, 2; extra_degree)
                # Enclosure of cospi(α / 2) / (α + 1)
                cos1_div_α = fx_div_x(s -> cospi((s - 1) / 2), r; extra_degree)
                # Enclosure of cospi((α - p0) / 2) / (α + 1)^2
                cos2_div_α2 =
                    fx_div_x(s -> cospi((-1 - s^2 / 2) / 2), r, 2; extra_degree)

                cos1_div_α / rgamma1_div_α - cos2_div_α2 / rgamma2_div_α2 * x^(r + r^2 / 2)
            else
                gamma(r - 1) * cospi((r - 1) / 2) -
                gamma(-1 - r^2 / 2) * cospi((-1 - r^2 / 2) / 2) * x^(r + r^2 / 2)
            end
        end

        a0αp1 * coefficient_div_α
    else
        indeterminate(α)
    end

    # Sum of terms with either p or q non-zero
    res2 = zero(x)
    exponent = zero(α)
    for ((p, q, i, j, k, l, m), y) in expansion
        if !iszero(y) && !(iszero(p) && iszero(q))
            if isone(p)
                # Add -α to the exponent coming from the p factor
                _exponent!(exponent, α, i - 1, j, k, l, m)

                if iszero(x)
                    if Arblib.ispositive(exponent)
                        term = zero(x)
                    else
                        term = indeterminate(x)
                    end
                elseif Arblib.contains_zero(x)
                    if Arblib.ispositive(exponent)
                        lower = zero(x)
                    else
                        lower = indeterminate(x)
                    end

                    # FIXME: Prove that this gives an upper bound
                    upper = let x = ubound(Arb, x)
                        @warn "Non-rigorous evaluation of expansion" maxlog = 1
                        (1 - Arb(Irrational{:γ}()) - log(x)) / π
                    end

                    term = Arb((lower, upper)) * abspow(x, exponent)
                else
                    term = abspow(x, exponent)
                    Arblib.mul!(term, term, p_coefficient)
                end
            else
                use_approx_p_and_q ||
                    throw(ArgumentError("only p == 0 or p == 1 supported, got p = $p"))

                _exponent!(exponent, α, i, j, k, l, m)

                term = abspow(x, exponent)

                let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, a0 = finda0(α)
                    p_factor =
                        a0 *
                        (
                            gamma(α) * cospi(α / 2) -
                            gamma(α - p0) * cospi((α - p0) / 2) * x^p0
                        ) *
                        x^-α
                    term *= p_factor^p
                end
            end

            # We don't support evaluation of non-zero q
            if !iszero(q)
                use_approx_p_and_q ||
                    throw(ArgumentError("only q == 0 is supported, got q = $q"))
                let α = -1 + u0.ϵ, p0 = 1 + α + (1 + α)^2 / 2, a0 = finda0(α)
                    q_factor =
                        -a0 * (
                            gamma(2α) * cospi(α) * x^(-2α) -
                            gamma(2α - p0) * cospi((2α - p0) / 2) * x^(-2α + p0) +
                            (-zeta(1 - 2α - 2) / 2 + zeta(1 - 2α + p0 - 2) / 2) * x^2
                        )
                    term *= q_factor^q
                end
            end

            #res += y * term
            Arblib.addmul!(res2, y, term)
        end
    end

    return res1 + res2
end

"""
    expansion_ispositive(u0::BHKdVAnsatz, expansion, ϵ)

Attempt to prove that the given expansion is positive on the interval
``(0, ϵ]``. Returns true on success and false on failure. It requires
that `0 < ϵ < 1`

The method can only handle expansions on certain simple forms. More
precisely it requires that all terms have a key of the form
```
(0, 0, 0, 0, 1, l, 0)
```
or
```
(0, 0, 0, 0, 0, 0, m)
```
with `l >= 1` and `m >= 2` Let `L` denote the set of all keys on the
first form and `M` the set of all keys on the second form.

The keys in `L` corresponds to terms of the form
```
a[l] * x^(-u0.v0.α + l*u0.v0.p0)
```
and the keys in `M` to terms of the form
```
b[m] * x^m
```
Let
```
S_L(x) = sum(a[l] * x^(-u0.v0.α + l*u0.v0.p0) for l in L)
S_M(x) = sum(b[m] * x^m for m in M)
```
We want to prove that `S_L(x) + S_M(x)` is positive for all `x` in
``(0, ϵ]``.

To begin with we compute a lower bound of `S_L`. If `a[l]` is negative
then a lower bound for the corresponding term is given by
```
a[l] * x^(-u0.v0.α + u0.v0.p0)
```
The method checks that `a[l]` is negative for all `l >= 2`. This means that
```
S_L(x) >= sum(a[l] * x^(-u0.v0.α + u0.v0.p0) for l in L)
       = sum(a[l] for l in L) * x^(-u0.v0.α + u0.v0.p0)
```
it then checks that `sum(a[l] for l in L) > 0` so that `S_L(x)` is
positive for all `x` in the interval.

The next step is to prove that `abs(S_M(x)) < S_L(x)`. We have
```
abs(S_M(x)) <= sum(abs(b[m]) * x^m for m in M)
```
So it is enough to check
```
sum(abs(b[m]) * x^m for m in M) < sum(a[l] for l in L) * x^(-u0.v0.α + u0.v0.p0)
```
Since `-u0.v0.α + u0.v0.p0 < 2` it is enough to check that this holds
for `x = ϵ`
"""
function expansion_ispositive(
    u0::BHKdVAnsatz{Arb},
    expansion::AbstractDict{NTuple{7,Int},Arb},
    ϵ::Arb,
)
    @assert 0 < ϵ < 1
    @assert 0 < -u0.v0.α + u0.v0.p0 < 2

    # Isolate all keys of the from (0, 0, 0, 0, 1, l, 0)
    expansion_1 = filter(expansion) do ((p, q, i, j, k, l, m), a)
        p == q == i == j == m == 0 && k == 1 && l >= 1
    end

    # Isolate all keys of the from (0, 0, 0, 0, 0, 0, m)
    expansion_2 = filter(expansion) do ((p, q, i, j, k, l, m), b)
        p == q == i == j == k == l == 0 && m >= 2
    end

    # Check that this was all keys
    @assert length(expansion) == length(expansion_1) + length(expansion_2)

    # Extract they key (0, 0, 0, 0, 1, 1, 0) from the expansion
    y₁ = expansion_1[(0, 0, 0, 0, 1, 1, 0)]
    delete!(expansion_1, (0, 0, 0, 0, 1, 1, 0))
    all(Arblib.isnegative, values(expansion_1)) || return false
    Arblib.ispositive(y₁ + sum(values(expansion_1))) || return false

    # This should always hold
    @assert -u0.v0.α + u0.v0.p0 < 2

    sum(abs(b) * ϵ^m for ((p, q, i, j, k, l, m), b) in expansion_2) <
    y₁ + sum(values(expansion_1)) * ϵ^(-u0.v0.α + u0.v0.p0) || return false

    return true
end

"""
    finda0αp1(α)

Compute `finda0(α) * (α + 1)`.

We have that
```
a0 = 2gamma(2α) * cospi(α) / (gamma(α)^2 * cospi(α / 2)^2)
```
Using that
```
gamma(2α) = gamma(2α + 3) / rising(2α, 3) = gamma(2α + 3) / ((2α + 2) * (2α + 1) * 2α)
gamma(α) = gamma(α + 2) / rising(α, 2) = gamma(α + 2) / ((α + 1) * α)
```
We can write this as
```
a0 = 2gamma(2α + 3) / ((2α + 2) * (2α + 1) * 2α) * cospi(α) / (gamma(α + 2)^2 / ((α + 1) * α)^2 * cospi(α / 2)^2)
   = 1 / 2 * gamma(2α + 3) / (2α + 1) * cospi(α) / (gamma(α + 2)^2 / (α + 1) * α * cospi(α / 2)^2)
   = 1 / 2 * cospi(α) * gamma(2α + 3) / gamma(α + 2)^2 * α / (2α + 1) * (α + 1) / cospi(α / 2)^2
```
Multiplying by `α + 1` we have
```
a0 * (α + 1) = 1 / 2 * cospi(α) * gamma(2α + 3) / gamma(α + 2)^2 * α / (2α + 1) * ((α + 1) / cospi(α / 2))^2
```
Using that
```
cospi(α / 2) / (α + 1) = sinc((α + 1) / 2) * π / 2
```
we can rewrite this as
```
a0 * (α + 1) = 2cospi(α) * gamma(2α + 3) / gamma(α + 2)^2 * α / (2α + 1) / (π * sinc((α + 1) / 2))^2
             = 2 / π^2 * α * gamma(2α + 3) * cospi(α) / ((2α + 1) * (gamma(α + 2) * sinc((α + 1) / 2))^2)
```
"""
finda0αp1(α::Arb) =
    (2 / Arb(π)^2) * ArbExtras.enclosure_series(α) do α
        # Use _sinc to allow for ArbSeries overlapping zero
        α * gamma(2α + 3) * cospi(α) / ((2α + 1) * (gamma(α + 2) * _sinc((α + 1) / 2))^2)
    end

finda0αp1(α::ArbSeries) =
    (2 / Arb(π)^2) * α * gamma(2α + 3) * cospi(α) /
    ((2α + 1) * (gamma(α + 2) * _sinc((α + 1) / 2))^2)

"""
    (u0::BHKdVAnsatz)(x, ::Ball)

Evaluate the ansatz `u0` at the point `x`.

The tail term is evaluated directly.

The main term is given by
```
a0 * (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0))
```
We rewrite this as
```
(a0 * (α + 1)) * ((clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) / (α + 1))
```
and enclose the two factors separately.

We enclose `a0 * (α + 1)` using [`finda0αp1`](@ref).

# Enclosing `(clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) / (α + 1)`
If `x::Arb` this is done directly using [`fx_div_x`](@ref). If
`x::ArbSeries` this is not possible since we can't mix expansion in
the argument and the parameter of `clausencmzeta`. Instead we
differentiate manually with respect to `x`. If we let `t = α + 1` then
we can write this as
```
(clausencmzeta(x, 2 - t) - clausencmzeta(x, 2 + t^2 / 2)) / t
```
We then use essentially the same approach as in
`clausencmzeta(x::ArbSeries, s::Arb)`.
"""
function (u0::BHKdVAnsatz{Arb})(x::Union{Arb,ArbSeries}, ::Ball)
    α = Arb((-1, -1 + u0.ϵ))
    αp1 = Arblib.nonnegative_part!(zero(u0.ϵ), union(zero(u0.ϵ), u0.ϵ))

    # Main term

    # Enclosure of a0 * (α + 1)
    a0αp1 = finda0αp1(α)

    # Enclosure of
    # (clausencmzeta(x, 1 - α) - clausencmzeta(x, 1 - α + p0)) / (α + 1)
    extra_degree = 2
    C = let
        if x isa Arb
            fx_div_x(αp1, force = true; extra_degree) do t
                clausencmzeta(x, 2 - t) - clausencmzeta(x, 2 + t^2 / 2)
            end
        elseif x isa ArbSeries
            x₀ = x[0]

            C₀ = zero(x)
            C₀[0] = fx_div_x(αp1, force = true; extra_degree) do t
                clausencmzeta(x₀, 2 - t) - clausencmzeta(x₀, 2 + t^2 / 2)
            end
            for i = 1:Arblib.degree(x)
                if i % 2 == 0
                    C₀[i] = fx_div_x(αp1, force = true; extra_degree) do t
                        (-1)^(i ÷ 2) *
                        (clausenc(x₀, 2 - t - i) - clausenc(x₀, 2 + t^2 / 2 - i)) /
                        factorial(i)
                    end
                else
                    C₀[i] = fx_div_x(αp1, force = true; extra_degree) do t
                        -(-1)^(i ÷ 2) *
                        (clausens(x₀, 2 - t - i) - clausens(x₀, 2 + t^2 / 2 - i)) /
                        factorial(i)
                    end
                end
            end

            x_tmp = copy(x)
            x_tmp[0] = 0

            Arblib.compose(C₀, x_tmp)
        end
    end

    res = a0αp1 * C

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.N0
        s = 1 - u0.v0.α + j * u0.v0.p0
        res += u0.v0.a[j] * clausencmzeta(x, s)
    end

    # Fourier terms
    for n = 1:u0.v0.N1
        res += u0.v0.b[n] * (cos(n * x) - 1)
    end

    return res
end

(u0::BHKdVAnsatz)(x, ::Asymptotic; M::Integer = 3) =
    eval_expansion(u0, u0(x, AsymptoticExpansion(); M), x)

"""
    (u0::BHKdVAnsatz)(x, ::AsymptoticExpansion; M = 3)

Return a dictionary containing the terms in the asymptotic expansion
of `u0` which can then be evaluated with [`eval_expansion`](@ref).

The highest term, `x^2M`, is an error term is which makes sure that
evaluation of the expansion gives an enclosure of the result when
evaluated at `|y| < |x|`.

See [`eval_expansion`](@ref) for more details about how the
coefficients are stored.

For the tail term the expansions are easily computed exactly like for
`BHAnsatz`. For the main term we have to be a bit more careful.

# Main term

## Leading term
The leading term of the expansion of the main term is
```
a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) * x^-α
```
which we don't evaluate at all yet. Instead we store them implicitly
in the expansion.

## Non-leading terms
For the main term the coefficients in front of `x^2m` is given by
```
a0 * (-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / factorial(2m)
```
It has a removable singularity at `α = -1`. To compute an enclosure we
rewrite it as
```
a0 * (α + 1) * (-1)^m * ((zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / (α + 1)) / factorial(2m)
```
and enclose `a0 * (α + 1)` using [`finda0αp1`](@ref) and `(zeta(1 - α
- 2m) - zeta(1 - α + p0 - 2m)) / (α + 1)` using [`fx_div_x`](@ref).

## Remainder term
The remainder term is given by
```
a0 * sum((-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) * x^2m / factorial(2m) for m = M:Inf) / x^2M
```
We want to bound the absolute value of this. We can rewrite it as
```
a0 * (1 + α) * sum((-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / (1 + α) * x^2m / factorial(2m) for m = M:Inf) / x^2M
```
We can enclose `a0 * (1 + α)` using [`finda0αp1`](@ref). What remains
is to bound the absolute value of
```
S = sum((-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / (1 + α) * x^(2m - 2M) / factorial(2m) for m = M:Inf)
```

If we let `t = α + 1` we can write the factor with the removable
singularity as
```
(zeta(2 - 2m - t) - zeta(2 - 2m + t^2 / 2)) / t
= (zeta(2 - 2m - t) - zeta(2 - 2m)) / t - (zeta(2 - 2m + t^2 / 2) - zeta(2 - 2m)) / t
```

We can write the first term as
```
-(zeta(2 - 2m + (-t)) - zeta(2 - 2m)) / (-t)
```
which is on the form `(f(s) - f(0)) / s` and hence given by `f'(ξ)`
for some `ξ` between `0` and `-t`.

For the second term we write it as
```
(zeta(2 - 2m + t^2 / 2) - zeta(2 - 2m)) / t
= t / 2 * (zeta(2 - 2m + t^2 / 2) - zeta(2 - 2m)) / (t^2 / 2)
```
which is given by `t / 2 * f'(ξ)` for some `ξ` between `0` and `t^2 /
2`.

Combining this gives us that
```
(zeta(2 - 2m - t) - zeta(2 - 2m + t^2 / 2)) / t = -dzeta(2 - 2m + ξ1) + t / 2 * dzeta(2 - 2m + ξ2)
```
with `ξ1` between `0` and `-t` and `ξ2` between `0` and `t^2 / 2`. In
terms of interval arithmetic we can write this as (using Arblib
notation)
```
-dzeta(2 - 2m + Arb((-t, 0))) + t / 2 * dzeta(2 - 2m + Arb((0, t^2 / 2)))
```
Thus we can reduce bounding the absolute value of `S` to bounding the
absolute value of
```
S1 = sum((-1)^m * zeta(2 + Arb((-(α + 1), 0)) - 2m) * x^2m / factorial(2m) for m = M:Inf)
S2 = sum((-1)^m * dzeta(2 + Arb((0, (α + 1)^2 / 2)) - 2m) * x^2m / factorial(2m) for m = M:Inf)
```
and combine them as `S1 + (α + 1) / 2 * S2`. These sums are the same
as those appearing in
```
clausenc_expansion_remainder(x, 2 + Arb((-(α + 1), 0)), 1, M)
clausenc_expansion_remainder(x, 2 + Arb((0, (α + 1)^2 / 2)), 1, M)
```
- **PROVE:** This is mostly true, the issue is that we are not
  guaranteed that the value in `Arb((-(α + 1), 0))` taken for the
  argument of the zeta function is the same on each iteration. This
  does not seem to be an issue because we rely on any cancellations
  occurring in the sum. But it would need to be proved.
"""
function (u0::BHKdVAnsatz{Arb})(x, ::AsymptoticExpansion; M::Integer = 3)
    @assert M >= 3

    # Enclosure of α, α + 1 and a0 * (α + 1)
    α = Arb((-1, -1 + u0.ϵ))
    αp1 = Arblib.nonnegative_part!(zero(u0.ϵ), union(zero(u0.ϵ), u0.ϵ))
    a0αp1 = finda0αp1(α)

    res = OrderedDict{NTuple{7,Int},Arb}()

    # Initiate even powers of x
    for m = 1:M
        res[(0, 0, 0, 0, 0, 0, 2m)] = 0
    end

    # Main term

    # Leading term - stored implicitly
    res[(1, 0, 0, 0, 0, 0, 0)] = 1

    # x^2m terms
    for m = 1:M-1
        # Enclosure of
        # (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / (α + 1)
        # = (zeta(2 - s - 2m) - zeta(2 + s^2 / 2 - 2m)) / s
        # with s = α + 1
        zeta_div_α = if m == 1
            # zeta(x::ArbSeries) doesn't handle balls containing
            # zero but centered at a negative number well. For
            # that reason we take a symmetric interval in this
            # case.
            fx_div_x(
                s -> zeta(2 - s - 2m) - zeta(2 + s^2 / 2 - 2m),
                union(-αp1, αp1),
                extra_degree = 2,
                force = true,
            )
        else
            fx_div_x(s -> zeta(2 - s - 2m) - zeta(2 + s^2 / 2 - 2m), αp1, extra_degree = 2)
        end

        coefficient = a0αp1 * (-1)^m * zeta_div_α / factorial(2m)

        res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
    end

    # Remainder term for main term
    remainder =
        a0αp1 * (
            clausenc_expansion_remainder(x, 2 - αp1, 1, M) +
            αp1 / 2 * clausenc_expansion_remainder(x, 2 + αp1^2 / 2, 1, M)
        )
    res[(0, 0, 0, 0, 0, 0, 2M)] += remainder

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.N0
        s = 1 - u0.v0.α + j * u0.v0.p0
        C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)
        res[(0, 0, 0, 0, 1, j, 0)] = C * u0.v0.a[j]
        for m = 1:M-1
            res[(0, 0, 0, 0, 0, 0, 2m)] += p[2m] * u0.v0.a[j]
        end
        res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.a[j]
    end

    # Fourier terms
    if !iszero(u0.v0.N1)
        for m = 1:M-1
            res[(0, 0, 0, 0, 0, 0, 2m)] +=
                (-1)^m * sum(Arb(n)^(2m) * u0.v0.b[n] for n = 1:u0.v0.N1) / factorial(2m)
        end
        Arblib.add_error!(
            res[(0, 0, 0, 0, 0, 0, 2M)],
            sum(Arb(n)^(2M) * abs(u0.v0.b[n]) for n = 1:u0.v0.N1) / factorial(2M),
        )
    end

    return res
end

"""
    H(u0::BHKdVAnsatz, ::Ball)

Returns a function such that `H(u0, Ball())(x)` evaluates
``H^α[u0](x)``.

# Main term
The transform of the main term is given by
```
-a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
```
We rewrite this as
```
-(a0 * (α + 1)) * ((clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0)) / (1 + α))
```
and enclose the two factors separately.

We enclose `a0 * (α + 1)` using [`finda0αp1`](@ref).

## Enclosing `(clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0)) / (α + 1)`
If `x::Arb` this is done directly using [`fx_div_x`](@ref). If
`x::ArbSeries` this is not possible since we can't mix expansion in
the argument and the parameter of `clausencmzeta`. Instead we
differentiate manually with respect to `x`. If we let `t = α + 1` then
we can write this as
```
(clausencmzeta(x, 3 - 2t) - clausencmzeta(x, 3 - t + t^2 / 2)) / t
```
We then use essentially the same approach as in
`clausencmzeta(x::ArbSeries, s::Arb)`.

# Tail
For the tail term we need to make sure that we correctly handle the
fact that the transform depends on the value of `α`. As long as `u0.ϵ`
is sufficiently small we get good enough bounds by just using it as a
ball directly.
"""
function H(u0::BHKdVAnsatz{Arb}, ::Ball)
    α = Arb((-1, -1 + u0.ϵ))
    αp1 = Arblib.nonnegative_part!(zero(u0.ϵ), union(zero(u0.ϵ), u0.ϵ))

    # Enclosure of a0 * (α + 1)
    a0αp1 = finda0αp1(α)

    return x::Union{Arb,ArbSeries} -> begin
        # Main term

        # Enclosure of
        # (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0)) / (α + 1)
        extra_degree = 2
        enclosure_degree = -1
        C = let
            if x isa Arb
                fx_div_x(αp1, force = true; extra_degree, enclosure_degree) do t
                    clausencmzeta(x, 3 - 2t) - clausencmzeta(x, 3 - t + t^2 / 2)
                end
            elseif x isa ArbSeries
                x₀ = x[0]

                C₀ = zero(x)
                C₀[0] =
                    fx_div_x(αp1, force = true; extra_degree, enclosure_degree) do t
                        clausencmzeta(x₀, 3 - 2t) - clausencmzeta(x₀, 3 - t + t^2 / 2)
                    end
                for i = 1:Arblib.degree(x)
                    if i % 2 == 0
                        C₀[i] = fx_div_x(αp1, force = true; extra_degree, enclosure_degree) do t
                            (-1)^(i ÷ 2) * (
                                clausenc(x₀, 3 - 2t - i) -
                                clausenc(x₀, 3 - t + t^2 / 2 - i)
                            ) / factorial(i)
                        end
                    else
                        C₀[i] = fx_div_x(αp1, force = true; extra_degree, enclosure_degree) do t
                            -(-1)^(i ÷ 2) * (
                                clausens(x₀, 3 - 2t - i) -
                                clausens(x₀, 3 - t + t^2 / 2 - i)
                            ) / factorial(i)
                        end
                    end
                end

                x_tmp = copy(x)
                x_tmp[0] = 0

                Arblib.compose(C₀, x_tmp)
            end
        end

        res = -a0αp1 * C

        # Tail term

        let α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α
            # IMPROVE: We might be able to compute better enclosures
            # in α, allowing us to work with larger values of u0.ϵ, by
            # expanding in s. This would require a fair amount of work
            # though, we would probably have to extract the formula
            # from _clausenc_zeta.

            # Clausen terms
            for j = 1:u0.v0.N0
                s = 1 - α - u0.v0.α + j * u0.v0.p0
                res -= u0.v0.a[j] * clausencmzeta(x, s)
            end

            # Fourier terms
            for n = 1:u0.v0.N1
                res -= u0.v0.b[n] * n^α * (cos(n * x) - 1)
            end
        end

        return res
    end
end

function H(
    u0::BHKdVAnsatz,
    ::Asymptotic;
    M::Integer = 3,
    skip_singular_j_until = 0,
    use_approx_p_and_q = false,
)
    f = H(u0, AsymptoticExpansion(); M, skip_singular_j_until)

    return x -> eval_expansion(u0, f(x), x; use_approx_p_and_q)
end

"""
    H(u0::BHKdVAnsatz, ::AsymptoticExpansion; M = 3, skip_singular_j_until::Integer = 0,)

Return a dictionary containing the terms in the asymptotic expansion
of `H(u0)` which can then be evaluated with [`eval_expansion`](@ref).

The highest term, `x^2M`, is an error term is which makes sure that
evaluation of the expansion gives an enclosure of the result when
evaluated at `|y| < |x|`.

See [`eval_expansion`](@ref) for more details about how the
coefficients are stored.

# Main term
For the main term we want the expansion of
```
-a0 * (clausencmzeta(x, 1 - 2α) - clausencmzeta(x, 1 - 2α + p0))
```

## Leading term
The leading term of the expansion of the main term is
```
-a0 * (gamma(2α) * cospi(α) - gamma(2α - p0) * cospi((2α - p0) / 2) * x^p0) * x^(-2α)
```
However, as `α` approaches `-1` this overlaps with the next term in
the expansion, `x^2`, and we need to account for their cancellations.
That term is given by
```
-a0 * (zeta(1 - 2α - 2) / 2 - zeta(1 - 2α + p0 - 2) / 2) * x^2
```
We don't evaluate these terms at all yet. Instead we store the term
```
-a0 * (
    (gamma(2α) * cospi(α) - gamma(2α - p0) * cospi((2α - p0) / 2) * x^p0) * x^(-2α)
    + (zeta(1 - 2α - 2) / 2 - zeta(1 - 2α + p0 - 2) / 2) * x^2
)
```
implicitly in the expansion.

## Non-leading terms
For the main term the coefficients in front of `x^2m` is given by
```
-a0 * (-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / factorial(2m)
```
For `m >= 2` it has a removable singularity at `α = -1`. To compute an
enclosure we rewrite it as
```
-a0 * (α + 1) * (-1)^m * ((zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / (α + 1)) / factorial(2m)
```
and enclose `a0 * (α + 1)` using [`finda0αp1`](@ref) and `(zeta(1 - 2α
- 2m) - zeta(1 - 2α + p0 - 2m)) / (α + 1)` using [`fx_div_x`](@ref).

## Remainder term
The remainder term is given by
```
-a0 * sum((-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) * x^2m / factorial(2m) for m = M:Inf) / x^2M
```
We want to bound the absolute value of this. We can rewrite it as
```
a0 * (1 + α) * sum((-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / (1 + α) * x^2m / factorial(2m) for m = M:Inf) / x^2M
```
We can enclose `a0 * (1 + α)` using [`finda0αp1`](@ref). What remains
is to bound the absolute value of
```
S = sum((-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / (1 + α) * x^(2m - 2M) / factorial(2m) for m = M:Inf)
```

If we let `t = α + 1` we can write the factor with the removable
singularity as
```
(zeta(3 - 2m - 2t) - zeta(3 - 2m - t + t^2 / 2)) / t
= (zeta(3 - 2m - 2t) - zeta(2 - 2m)) / t - (zeta(3 - 2m - t + t^2 / 2) - zeta(2 - 2m)) / t
```

We can write the first term as
```
-2(zeta(3 - 2m + (-2t)) - zeta(3 - 2m)) / (-2t)
```
which is on the form `(f(s) - f(0)) / s` and hence given by `f'(ξ)`
for some `ξ` between `0` and `-2t`.

For the second term we write it as
```
(zeta(3 - 2m - t + t^2 / 2) - zeta(3 - 2m)) / t
= (-1 + t / 2) * (zeta(3 - 2m + t^2 / 2) - zeta(3 - 2m)) / (-t + t^2 / 2)
```
which is given by `(-1 + t / 2) * f'(ξ)` for some `ξ` between `0` and
`-t + t^2 / 2`.

Combining this gives us that
```
(zeta(3 - 2m - 2t) - zeta(3 - 2m -t + t^2 / 2)) / t = -2dzeta(3 - 2m + ξ1) + (-1 + t / 2) * dzeta(3 - 2m + ξ2)
```
with `ξ1` between `0` and `-2t` and `ξ2` between `0` and `-t + t^2 / 2`. In
terms of interval arithmetic we can write this as (using Arblib
notation)
```
-2dzeta(3 - 2m + Arb((-2t, 0))) + (-1 + t / 2) * dzeta(3 - 2m + Arb((-t + t^2 / 2, 0)))
```
Thus we can reduce bounding the absolute value of `S` to bounding the
absolute value of
```
S1 = sum((-1)^m * zeta(3 + Arb((-2(α + 1), 0)) - 2m) * x^2m / factorial(2m) for m = M:Inf)
S2 = sum((-1)^m * dzeta(3 + Arb((0, -(α + 1) + (α + 1)^2 / 2)) - 2m) * x^2m / factorial(2m) for m = M:Inf)
```
and combine them as `2S1 + (1 - (α + 1) / 2) * S2`. These sums are
the same as those appearing in
```
clausenc_expansion_remainder(x, 3 + Arb((-2(α + 1), 0)), 1, M)
clausenc_expansion_remainder(x, 3 + Arb((0, -(α + 1) + (α + 1)^2 / 2)), 1, M)
```
- **PROVE:** This is mostly true, the issue is that we are not
  guaranteed that the value in e.g.g `Arb((-2(α + 1), 0))` taken for
  the argument of the zeta function is the same on each iteration.
  This does not seem to be an issue because we rely on any
  cancellations occurring in the sum. But it would need to be proved.

# Tail
For both the Clausen terms and the Fourier terms we let `α` be a ball.
This gives good enclosures for the Fourier terms. For the Clausen
terms it give good enclosures unless `j` is small. For small values of
`j` the two terms
```
gamma(α + u0.v0.α - j * u0.v0.p0) * cospi((α + u0.v0.α - j * u0.v0.p0) / 2) *
    x^-(α + u0.v0.α - j * u0.v0.p0)
```
and
```
-zeta(-1 - α + u0.v0.α - j * u0.v0.p0) / 2 * x^2
```
have very large cancellations.

If `skip_singular_j_until` is greater than zero then skip the two
above terms for all Clausen functions in the tail from `j = 1` to
`skip_singular_j_until`. This is used in [`F0`](@ref) where these
terms are handled separately.

If `approximate_singular_j_until` is greater than zero then
approximate the two above terms for all Clausen functions in the tail
from `j = 1` to `approximate_singular_j_until`. This is only for
testing. Terms skipped by `skip_singular_j_until` are not
approximated.
"""
function H(
    u0::BHKdVAnsatz{Arb},
    ::AsymptoticExpansion;
    M::Integer = 3,
    skip_singular_j_until::Integer = 0,
    approximate_singular_j_until::Integer = 0,
)
    @assert M >= 3

    skip_singular_j_until > u0.v0.N0 && throw(
        ArgumentError("can't skip more j-terms than there are, j = $j, N0 = $(u0.v0.N0)"),
    )

    # Enclosure of α, α + 1 and a0 * (α + 1)
    α = Arb((-1, -1 + u0.ϵ))
    αp1 = Arblib.nonnegative_part!(zero(u0.ϵ), union(zero(u0.ϵ), u0.ϵ))
    a0αp1 = finda0αp1(α)

    return x -> begin
        res = OrderedDict{NTuple{7,Int},Arb}()

        # Initiate even powers of x
        for m = 1:M
            res[(0, 0, 0, 0, 0, 0, 2m)] = 0
        end

        # Main term

        # Leading term - stored implicitly
        res[(0, 1, 0, 0, 0, 0, 0)] = 1

        # x^2m terms with m >= 2
        for m = 2:M-1
            # Enclosure of
            # (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / (α + 1)
            # = (zeta(3 - 2s - 2m) - zeta(3 - s + s^2 / 2 - 2m)) / s
            # with s = α + 1
            zeta_div_α = fx_div_x(
                s -> zeta(3 - 2s - 2m) - zeta(3 - s + s^2 / 2 - 2m),
                αp1,
                extra_degree = 2,
                force = true,
            )

            coefficient = -a0αp1 * (-1)^m * zeta_div_α / factorial(2m)

            res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
        end

        # Remainder term for main term
        remainder =
            a0αp1 * (
                2clausenc_expansion_remainder(x, 3 - 2αp1, 1, M) +
                (1 - αp1 / 2) * clausenc_expansion_remainder(x, 3 - αp1 + αp1^2 / 2, 1, M)
            )
        res[(0, 0, 0, 0, 0, 0, 2M)] += remainder

        # Tail term

        # Clausen terms
        for j = 1:u0.v0.N0
            s = 1 - α - u0.v0.α + j * u0.v0.p0
            C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

            if j > skip_singular_j_until
                if j <= approximate_singular_j_until
                    C2, _, p2, _ = let s = 1 - (-1 + u0.ϵ) - u0.v0.α + j * u0.v0.p0
                        clausenc_expansion(x, s, M, skip_constant = true)
                    end
                    res[(0, 0, -1, 0, 1, j, 0)] = -C2 * u0.v0.a[j]
                    res[(0, 0, 0, 0, 0, 0, 2)] -= p2[2] * u0.v0.a[j]
                else
                    res[(0, 0, -1, 0, 1, j, 0)] = -C * u0.v0.a[j]
                    res[(0, 0, 0, 0, 0, 0, 2)] -= p[2] * u0.v0.a[j]
                end
            end

            for m = 2:M-1
                res[(0, 0, 0, 0, 0, 0, 2m)] -= p[2m] * u0.v0.a[j]
            end
            res[(0, 0, 0, 0, 0, 0, 2M)] += E * u0.v0.a[j]
        end

        # Fourier terms
        if !iszero(u0.v0.N1)
            for m = 1:M-1
                res[(0, 0, 0, 0, 0, 0, 2m)] -=
                    (-1)^m * sum(n^α * Arb(n)^(2m) * u0.v0.b[n] for n = 1:u0.v0.N1) /
                    factorial(2m)
            end
            Arblib.add_error!(
                res[(0, 0, 0, 0, 0, 0, 2M)],
                sum(n^α * Arb(n)^(2M) * abs(u0.v0.b[n]) for n = 1:u0.v0.N1) / factorial(2M),
            )
        end

        return res
    end
end

function D(u0::BHKdVAnsatz, ::Asymptotic; M::Integer = 3, skip_singular_j_until = 0)
    f = D(u0, AsymptoticExpansion(); M, skip_singular_j_until)

    return x -> eval_expansion(u0, f(x), x)
end

function D(
    u0::BHKdVAnsatz,
    evaltype::AsymptoticExpansion;
    M::Integer = 3,
    skip_singular_j_until = 0,
)
    f = x -> u0(x, evaltype; M)
    g = H(u0, evaltype; M, skip_singular_j_until)

    return x -> begin
        expansion1 = f(x)
        expansion2 = g(x)
        expansion = empty(expansion1)

        # u0^2/2 term
        let expansion1 = collect(expansion1)
            z = zero(x) # Avoid allocating zero multiple times
            for (i, (key1, a1)) in enumerate(expansion1)
                expansion[2 .* key1] = get(expansion, 2 .* key1, z) + a1^2 / 2
                for j = i+1:length(expansion1)
                    (key2, a2) = expansion1[j]
                    key = key1 .+ key2
                    expansion[key] = get(expansion, key, z) + a1 * a2
                end
            end
        end

        # H term
        merge!(+, expansion, expansion2)

        return expansion
    end
end

"""
    F0_bound(u0::BHKdVAnsatz{Arb}, evaltype::Ball = Ball())

Return a function `f` such that the absolute value of `f(x)` is an
upper bound of the absolute value of `F0(u0)(x)`.

More precisely this computes
```
D(u0)(x) / (u0.w(x) * u0.v0(x))
```
Since `u0.v0(x)` gives a lower bound of `u0(x)`, by
[`lemma_bhkdv_monotonicity_alpha`](@ref), this gives a value which has
the same sign as `F0(x)` but is larger in magnitude. This holds as
long as `u0.v0(x)` is positive at least, which is easily checked.

**IMPROVE:** The computation of `u0.v0(x)` and `u0(x)` have many
calculations in common. We could improve performance by using this.
"""
function F0_bound(u0::BHKdVAnsatz{Arb}, evaltype::Ball = Ball())
    g = D(u0, evaltype)

    return x::Union{Arb,ArbSeries} -> begin
        invweight = inv(u0.w(x))

        isfinite(invweight) || return invweight

        invu0v0 = inv(u0.v0(x))

        isfinite(invu0v0) || return invu0v0

        if (x isa Arb && !Arblib.ispositive(invu0v0)) ||
           (x isa ArbSeries && !Arblib.ispositive(Arblib.ref(invu0v0, 0)))
            error("expected u0.v0(x) to be positive, got inv(u0.v0(x)) = $invu0v0")
        end

        return g(x) * invu0v0 * invweight
    end
end


"""
    F0(u0::BHKdVAnsatz{Arb}, ::Asymptotic)

Returns a function such that an **upper bound** of `F0(u0)(x)` is
computed accurately for small values of `x`.

It requires that `0 <= x < 1`, any negative parts of `x` are ignored.

Recall that the expression we are interested in bounding is
```
(u0(x)^2 / 2 + H(u0)(x)) / (u0.w(x) * u0(x))
```
with
```
u0.w(x) = x^(1 - u0.γ * (1 + α)) * log(u0.c + inv(x))
```

# Split into three factors
As a first step we split the expression into three factors which are
all bounded as `x -> 0` that we bound separately.
We write it as
```
(gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x))
* (log(inv(x)) / log(u0.c + inv(x)))
* ((u0(x)^2 / 2 + H(u0)(x)) / (gamma(1 + α) * log(inv(x)) * x^(1 - u0.γ * (1 + α) - α) * (1 - x^p0)))
```
The first factor is bounded using [`inv_u0_bound`](@ref). The second
factor is bounded by noticing that it is `1` at `x = 0` and decreasing
in `x`.

For the third factor we can get an upper bound, in absolute value, by
taking `γ = 0`, giving us
```
F = (u0(x)^2 / 2 + H(u0)(x)) / (gamma(1 + α) * log(inv(x)) * (1 - x^p0) * x^(1 - α))
```
We now explain how to bound `F`

# Bounding `F`
Getting an accurate bound for `F` requires more work.

We use the asymptotic expansion of `u0(x)^2 / 2 + H(u0)(x)`, but some
of the terms in the expansion require extra care. We therefore split
the expansion into three parts
1. The two leading terms, with keys `(2, 0, 0, 0, 0, 0, 0)` and `(0, 1,
    0, 0, 0, 0, 0)`, which we call `P` and `Q` respectively.
2. The leading terms of the Clausen terms in the tail for small values
   of `j`.
3. The remaining terms.
This splits `F` into three terms, which we will call `T1, T2, T3`.

## Handling `T1`: `P` and `Q`
The terms `P` and `Q` are given by
```
P^2 = a0^2 / 2 * (c(α)^2 - 2c(α) * c(α - p0) * x^p0 + c(α - p0)^2 * x^2p0) * x^-2α
```
and
```
Q = -a0 * ((c(2α) - c(2α - p0) * x^p0) * x^-2α - (zeta(-2α - 1) / 2 - zeta(-2α + p0 - 1) / 2) * x^2)
```
where `c(a) = gamma(a) * cospi(a / 2)`. We can handle the removable
singularity at `a = -1` similar to how it is done in
[`inv_u0_bound`](@ref).

By construction `a0` is such that the terms with exponent `x^-2α`
cancel out. This leaves us with
```
P + Q = a0 * (
    (c(2α - p0) - a0 * c(α) * c(α - p0)) * x^(-2α + p0) +
    (zeta(-2α - 1) - zeta(-2α + p0 - 1)) / 2 * x^2 +
    a0 * c(α - p0)^2 / 2 * x^(-2α + 2p0)
)
```
where we can note that `-2α + p0 < 2 < -2α + 2p0`. We are now
interested in bounding
```
(P + Q) / (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0))
```
If we cancel the exponents and reorder the factors slightly we get
```
a0 / gamma(1 + α) * (
    (c(2α - p0) - a0 * c(α) * c(α - p0)) * x^(-α + p0 - 1) / (log(x) * (1 - x^p0)) +
    (zeta(-2α - 1) - zeta(-2α + p0 - 1)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0)) +
    a0 * c(α - p0)^2 / 2 * x^(-α + 2p0 - 1) / (log(x) * (1 - x^p0))
)
```
The factor `a0 / gamma(1 + α)` can be rewritten as
```
a0 * (1 + α) / gamma(2 + α)
```
and enclosed using [`finda0αp1`](@ref).

First we focus on the term
```
T11 = (c(2α - p0) - a0 * c(α) * c(α - p0)) * x^(-α + p0 - 1) / (log(x) * (1 - x^p0))
```
This term is small and fairly stable in `α` (positive and increasing).
To compute an enclosure we split it into three factors
```
T111 = x^(-α + p0 - 1) / log(x)
T112 = (1 + α) / (1 - x^p0)
T113 = (c(2α - p0) - a0 * c(α) * c(α - p0)) / (1 + α)
```
For the first one we can directly get an enclosure using that
```
-α + p0 - 1 = -α + (1 + α + (1 + α)^2 / 2) - 1 = (1 + α)^2 / 2
```
For the second one we note that it is increasing in `x` and takes the
value `1 + α` at `x = 0`. For non-zero `x` we can handle the removable
singularity in `α`.

For the third term, `T113`, we note that `a0` can be written as
```
a0 = 2c(2α) / c(α)^2
```
which allows us to simplify it to
```
T113 = (c(2α - p0) - 2c(2α) * c(α - p0) / c(α)) / (1 + α)
```
Again we handle the removable singularity at `α = -1` as before, this
one requires a bit more work though. To do this we rewrite it as
```
(
    (1 + α) * c(2α - p0)
    - 2((1 + α) * c(2α)) * c(α - p0) / c(α)
) / (1 + α)^2
```

We then consider the two remaining terms together since they mostly
cancel out.
```
T12 = (zeta(-2α - 1) - zeta(-2α + p0 - 1)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0)) +
    a0 * c(α - p0)^2 / 2 * x^(-α + 2p0 - 1) / (log(x) * (1 - x^p0))
    =
    ((zeta(-2α - 1) - zeta(-2α + p0 - 1)) +
    a0 * c(α - p0)^2 * x^(-2α + 2p0 - 2)) / 2 * x^(1 + α) / (log(x) * (1 - x^p0))
```
We treat this in a very similar way as `T11`, by splitting it into
three factors.
```
T121 = x^(1 + α) / 2
T122 = (1 + α) / (1 - x^p0)
T123 = ((zeta(-2α - 1) - zeta(-2α + p0 - 1)) + a0 * c(α - p0)^2 * x^(-2α + 2p0 - 2)) /
    ((1 + α) * log(x))
```
The factor `T121` we can enclose directly. The factor `T122` is the
same as the factor `T112` above. We are hence left enclosing `T123`.
To do that we split it into two terms, letting
```
v(α) = zeta(-2α - 1) - zeta(-2α + p0 - 1)
w(α) = -a0 * c(α - p0)^2 = -2c(2α) * c(α - p0)^2 / c(α)^2
```
and noticing that `-2α + 2p0 - 2 = (1 + α)^2` we can write
```
T123 = (v(α) - w(α) * x^((1 + α)^2)) / ((1 + α) * log(x))
```
Adding and subtracting `w(α)` in the numerator we can split this into
two terms
```
T1231 = (v(α) - w(α)) / ((1 + α) * log(x))
T1232 = w(α) * (1 - x^((1 + α)^2)) / ((1 + α) * log(x))
```
For `T1231` it enough to handle the removable singularity of `(v(α) -
w(α)) / (1 + α)` and then multiply by an enclosure of `inv(log(x))`.

For `T1232` we can rewrite it as
```
T1232 = (1 + α) * w(α) * (1 - x^((1 + α)^2)) / ((1 + α)^2 * log(x))
```
We can compute an enclosure of
```
(1 + α) * w(α) = -2(1 + α) * c(2α) * c(α - p0)^2 / c(α)^2
```
by handling the removable singularity.

We can compute an enclosure of
```
(1 - x^((1 + α)^2) / ((1 + α)^2 * log(x))
```
By letting `t = (1 + α)^2 * log(x)` to write it as
```
(1 - exp(t)) / t
```
We then handle this similarly to how we do it for the tail Clausen
below.

## Handling `T2`: tail Clausen with small `j`
For small values of `j` the two terms
```
gamma(α + u0.v0.α - j * u0.v0.p0) * cospi((α + u0.v0.α - j * u0.v0.p0) / 2) *
    x^-(α + u0.v0.α - j * u0.v0.p0)
```
and
```
-zeta(-1 - α + u0.v0.α - j * u0.v0.p0) / 2 * x^2
```
in the asymptotic expansion of the Clausen functions in the tail have
very large cancellations. It is therefore beneficial to treat them
together to account for the cancellations. Which `j` we treat
separately like this is determined by the argument
`skip_singular_j_until`, it then handles `j = 1:skip_singular_j_until`
separately.

We are thus interested in bounding the first two terms when expanding
the Clausen function in
```
-u0.v0.a[j] * clausenc(x, 1 - α - u0.v0.α + j * u0.v0.p0) /
    (gamma(1 + α) * log(x) * x^(1 - α) * (1 - x^p0))
```
We can get an enclosure of `inv(gamma(1 + α) * (1 - x^p0))` by
noticing that it is increasing in `x` and it is hence enough to
compute at the endpoints of `x`. For `x = 0` it is given by
`inv(gamma(1 + α)) = rgamma(1 + α)`. Otherwise we use the same
approach as in [`inv_u0_bound`](@ref) for enclosing it. We are then
interested in enclosing the rest.

Let `r = -u0.v0.α + j * u0.v0.p0 - 1`, then `r > 0` and for
small values of `j` it is very close to zero. We have `1 - α -
u0.v0.α + j * u0.v0.p0 = 2 - α + r`. The sum of the first two
terms in the asymptotic expansion of the Clausen is then given by
```
gamma(α - 1 - r) * cospi((α - 1 - r) / 2) * x^(1 - α + r) -
    zeta(-α + r) / 2 * x^2
```
Dividing by `log(x) * x^(1 - α)` gives us
```
(gamma(α - 1 - r) * cospi((1 - α + r) / 2) * x^r - zeta(-α + r) / 2 * x^(1 + α)) / log(x)
```
Notice that the order of the terms depends on the value of `α`, in
some cases `x^r` is leading and in some cases `x^(1 + α)`. Adding and
subtracting `zeta(-α + r) / 2 * x^r` we can rewrite this as
```
(gamma(α - 1 - r) * cospi((1 - α + r) / 2) - zeta(-α + r) / 2) * x^r / log(x) +
zeta(-α + r) / 2 * (x^r - x^(1 + α)) / log(x)
```
Let
```
T21 = (gamma(α - 1 - r) * cospi((1 - α + r) / 2) - zeta(-α + r) / 2) * x^r / log(x)
T22 = zeta(-α + r) / 2 * (x^r - x^(1 + α)) / log(x)
```

For `T21` we rewrite it using that
```
gamma(α - 1 - r) = gamma(α + 2 - r) / ((α + 1 - r) * (α - r) * (α - 1 - r))

zeta(-α + r) / 2 = zeta_deflated(-α + r) / 2 - 1 / 2(1 + α - r)
```
giving us
```
T21 = (
    gamma(α + 2 - r) / ((α - r) * (α - 1 - r)) * cospi((1 - α + r) / 2) + 1 // 2
) / (α + 1 - r) - zeta_deflated(-α + r) / 2,
```
This formulation is much more stable and can be accurately enclosed.
For the case when `α + 1 - r` overlaps zero we can handle the
removable singularity.

For `T22` we have to compute an enclosure of
```
zeta(-α + r) / 2 * (x^r - x^(1 + α)) / log(x)
```
Using [`zeta_deflated`](@ref) we have
```
zeta(-α + r) = zeta_deflated(-α + r) - 1 / (1 + α - r)
```
and can split `T22` as
```
zeta_deflated(-α + r) * (x^r - x^(1 + α)) / 2log(x) + (x^r - x^(1 + α)) / (-α + r - 1) / 2log(x)
= T221 + T222
```
The term `T221` can be enclosed directly. For `T222` we rewrite it as
```
inv(2log(x)) * (x^r - x^(1 + α)) / (r - (1 + α))
```
We want to show that this is non-decreasing in `α`. Differentiating
the quotient w.r.t. `α` gives us
```
(x^r - (1 + log(x) * (r - (1 + α))) * x^(1 + α)) / (r - (1 + α))^2
```
The sign depends only on the sign of the numerator, we want to prove
that it is non-negative. Differentiating the numerator w.r.t. `α`
gives us
```
-log(x)^2 * (r - (1 + α)) * x^(1 + α)
```
The sign of this depend only on `r - (1 + α)` and we see that it is
decreasing for `r > 1 + α` and increasing for `r < 1 + α`. The minimum
of
```
x^r - (1 + log(x) * (r - (1 + α))) * x^(1 + α)
```
is hence attained at `r = 1 + α`, where it takes the value
```
x^r - (1 + log(x) * (r - r)) * x^r = 0
```
It follows that it is non-negative for all `α` Hence
```
(x^r - x^(1 + α)) / (r - (1 + α))
```
is non-decreasing in `α` and we can evaluate it at the endpoints.

## Handling `T3`: the remaining terms
Once the terms `P` and `Q` have been taken out from the expansion it
is possible to enclose the remaining ones directly. However the
expansion for the first Clausen function in the tail has very large
cancellations between the first two terms in the expansion and they
are therefore handled separately. The rest of the terms we enclose
directly by using [`eval_expansion`](@ref). If `u0.ϵ` is small this
gives good enough enclosures so that we don't have to handle it in any
more sophisticated way.

Recall that all terms are supposed to be divided by `(gamma(1 + α) *
log(x) * x^(1 - α) * (1 - x^p0))`. The `x^(1 - α)` factor will be
cancelled explicitly. For the remaining part we compute an enclosure.
We are therefore interested in computing an enclosure of
```
inv(log(x)) * inv(gamma(1 + α) * (1 - x^p0))
```
The `inv(log(x))` is easily handled using monotonicity in `x` and the
other factor is the same as in the above section.
"""
function F0(
    u0::BHKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 3,
    ϵ::Arb = Arb(0.5),
    skip_singular_j_until = 100,
)
    @assert ϵ < 1

    # Function for bounding gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x)
    f1 = inv_u0_bound(u0; M, ϵ)

    # Function for enclosing log(inv(x)) / log(u0.c + inv(x))
    f2 = x -> if iszero(x)
        one(x)
    elseif Arblib.contains_zero(x)
        lower = let xᵤ = ubound(Arb, x)
            log(inv(xᵤ)) / log(u0.c + inv(xᵤ))
        end
        upper = one(x)
        Arb((lower, upper))
    else
        log(inv(x)) / log(u0.c + inv(x))
    end

    # Compute the expansion of D(u0), skipping the two leading terms
    # in the expansion of the Clausen term in the tail for j =
    # 1:skip_j_until, which are handled separately
    Du0_expansion = D(u0, AsymptoticExpansion(); M, skip_singular_j_until)(ϵ)
    delete!(Du0_expansion, (2, 0, 0, 0, 0, 0, 0))
    delete!(Du0_expansion, (0, 1, 0, 0, 0, 0, 0))

    # Divide the expansion of D(u0) by x^(1 - α)
    Du0_expansion_div_x_onemα = empty(Du0_expansion)
    for ((p, q, i, j, k, l, m), y) in Du0_expansion
        Du0_expansion_div_x_onemα[(p, q, i + 1, j, k, l, m - 1)] = y
    end

    # Compute enclosures of several values depending only on α, many
    # of them with removable singularities

    # Interval for α
    α = Arb((-1, -1 + u0.ϵ))
    # Interval for α + 1
    αp1 = Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))
    # Use this for computing tighter enclosures
    extra_degree = 2

    # c(a) = gamma(a) * cospi(a / 2) rewritten to handle the removable
    # singularity at a = -1
    c(a) = Arb(π) * gamma(2 + a) * _sinc((1 + a) / 2) / 2a

    # Enclosure of a0 / gamma(1 + α)
    a0_div_gamma_1pα = ArbExtras.enclosure_series(α) do α
        finda0αp1(Arb((-1, -1 + u0.ϵ))) / gamma(2 + α)
    end

    # Enclosure of rgamma(1 + α) / (1 + α) = rgamma(2 + α)
    rgamma_1pα_div_αp1 = rgamma(2 + α)

    # Enclosure of ((1 + α) * c(2α - p0) - 2(1 + α) * c(2α) * c(α - p0) / c(α)) / (1 + α)^2
    T113 = fx_div_x(αp1, 2; extra_degree, force = true) do s
        if Arblib.contains_zero(s[0])
            # Expansion of c(2α - p0) * (1 + α)
            c_2αmp0_mul_α = inv(
                fx_div_x(
                    t -> rgamma(t - 2 - t^2 / 2) / cospi((t - 2 - t^2 / 2) / 2),
                    s;
                    extra_degree,
                ),
            )

            # Enclosure of c(2α) * (1 + α)
            c_2α_mul_α =
                inv(fx_div_x(t -> rgamma(2(t - 1)) / cospi(t - 1), s; extra_degree))

            # Enclosure of c(α - p0)
            c_αmp0 = c(-1 - s^2 / 2)

            # Enclosure of c(α)
            c_α = c(s - 1)

            c_2αmp0_mul_α - 2c_2α_mul_α * c_αmp0 / c_α
        else
            # Evaluate this at higher precision since it is close to
            # the removable singularity
            let α = setprecision(s - 1, 2precision(s)), p0 = (1 + α) + (1 + α)^2 / 2
                setprecision((c(2α - p0) - 2c(2α) * c(α - p0) / c(α)) / (1 + α), precision(s))
            end
        end
    end

    # Enclosure of T1231 * log(x) = (v(α) - w(α)) / (1 + α)
    T1231_mul_logx = fx_div_x(αp1, 2; extra_degree, force = true) do s
        if Arblib.contains_zero(s[0])
            # Enclosure of s * zeta(1 - 2s)
            zeta_1m2s_mul_α = s * zeta_deflated(1 - 2s, one(Arb)) - 1 // 2

            # Enclosure of s * zeta(1 - s + s^2 / 2)
            zeta_1msps2div2_mul_α =
                s * zeta_deflated(1 - s + s^2 / 2, one(Arb)) - 1 / (1 - s / 2)

            v_mul_α = zeta_1m2s_mul_α - zeta_1msps2div2_mul_α

            # Enclosure of c(2α) * (1 + α)
            c_2α_mul_α =
                inv(fx_div_x(t -> rgamma(2(t - 1)) / cospi(t - 1), s; extra_degree))

            # Enclosure of c(α - p0)
            c_αmp0 = c(-1 - s^2 / 2)

            # Enclosure of c(α)
            c_α = c(s - 1)

            w_mul_α = -2c_2α_mul_α * c_αmp0^2 / c_α^2

            v_mul_α - w_mul_α
        else
            # Evaluate this at higher precision since it is close to
            # the removable singularity
            let α = setprecision(s - 1, 2precision(s)), p0 = (1 + α) + (1 + α)^2 / 2
                v = zeta(-2α - 1) - zeta(-2α + p0 - 1)
                w = -2c(2α) * c(α - p0)^2 / c(α)^2
                setprecision(v - w, precision(s))
            end
        end
    end

    # Enclosure of w(α) * (1 + α)
    w_mul_α = ArbExtras.enclosure_series(α) do α
        -finda0αp1(α) * c(-1 + (1 + α)^2 / 2)^2
    end

    # α-factor of T21 for j = 1:skip_singular_j_until
    # Use a very high degree to get a good enclosure.
    T21_α = map(1:skip_singular_j_until) do j
        let r = -u0.v0.α + j * u0.v0.p0 - 1
            ArbExtras.enclosure_series(α, degree = 10) do α
                term1 = let t = -α + r
                    if t isa ArbSeries && is_approx_integer(t[0])
                        # This gives much better enclosures
                        t[0] = union(t[0], Arb(1))
                    end
                    -zeta_deflated(t, one(Arb)) / 2
                end

                term2 = let t = α + 1 - r
                    if t isa ArbSeries && is_approx_integer(t[0])
                        # This gives much better enclosures
                        t[0] = union(t[0], Arb(0))
                    end

                    if (t isa Arb && Arblib.contains_zero(t)) ||
                       (t isa ArbSeries && Arblib.contains_zero(t[0]))
                        fx_div_x(t; extra_degree) do s
                            gamma(s + 1) / ((s - 1) * (s - 2)) * cospi((2 - s) / 2) + 1 // 2
                        end
                    else
                        (gamma(t + 1) / ((t - 1) * (t - 2)) * cospi((t - 2) / 2) + 1 // 2) / t
                    end
                end

                term1 + term2
            end
        end
    end

    # Enclosure of zeta_deflated(-α + r) for j = 1:skip_singular_j_until
    # Use a very high degree to get a good enclosure.
    zeta_deflated_mαpr = map(1:skip_singular_j_until) do j
        let r = -u0.v0.α + j * u0.v0.p0 - 1
            ArbExtras.enclosure_series(α, degree = 10) do α
                if α isa ArbSeries && is_approx_integer(-α[0] + r)
                    # This gives much better enclosures
                    α[0] = union(α[0], r - 1)
                end
                zeta_deflated(-α + r, one(r))
            end
        end
    end

    return x::Arb -> begin
        @assert x <= ϵ

        # Enclosure of inv(log(x))
        invlogx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            Arb((inv(log(ubound(Arb, x))), 0))
        else
            inv(log(x))
        end

        # Enclosure of inv(gamma(1 + α) * (1 - x^p0))
        invgamma1mxp0 = if iszero(x)
            rgamma(1 + α)
        elseif Arblib.contains_zero(x)
            lower = zero(x)
            upper = let xᵤ = ubound(Arb, x)
                # Enclosure of (1 - xᵤ^p0) / (1 + α)
                onemxp0_div_αp1 =
                    fx_div_x(s -> (1 - xᵤ^(s + s^2 / 2)), αp1, extra_degree = 2)
                rgamma_1pα_div_αp1 / onemxp0_div_αp1
            end
            Arb((lower, upper))
        else
            # Enclosure of (1 - x^p0) / (1 + α)
            onemxp0_div_αp1 = fx_div_x(s -> (1 - x^(s + s^2 / 2)), αp1, extra_degree = 2)
            rgamma_1pα_div_αp1 / onemxp0_div_αp1
        end

        # Enclosure for the terms P and Q
        T1 = let
            # Compute an enclosure of T11

            # Enclosure of T111 = x^(-α + p0 - 1) / log(x) using that
            # -α + p0 - 1 = (1 + α)^2 / 2
            T111 = abspow(x, Arblib.nonnegative_part!(zero(x), αp1^2 / 2)) * invlogx

            # Enclosure of T112 = (1 + α) / (1 - x^p0)
            T112 = if iszero(x)
                αp1
            elseif Arblib.contains_zero(x)
                lower = αp1
                upper = let xᵤ = ubound(Arb, x)
                    # Enclosure of inv((1 - xᵤ^p0) / (1 + α))
                    inv(fx_div_x(s -> (1 - xᵤ^(s + s^2 / 2)), αp1; extra_degree))
                end
                Arb((lower, upper))
            else
                # Enclosure of inv((1 - x^p0) / (1 + α))
                inv(fx_div_x(s -> (1 - x^(s + s^2 / 2)), αp1; extra_degree))
            end

            T11 = T111 * T112 * T113

            # Compute an enclosure of T12

            # Enclosure of T121 = x^(1 + α) / 2
            T121 = abspow(x, Arblib.nonnegative_part!(zero(x), αp1)) / 2

            # Enclosure of T122, which is the same as T112
            T122 = T112

            T123 = let
                T1231 = T1231_mul_logx * invlogx

                T1232 = let
                    # Lower and upper bounds of
                    # t = (1 + α)^2 * log(x)
                    tₗ = ubound(Arb, αp1^2) * log(abs_lbound(Arb, x))
                    tᵤ = abs_lbound(Arb, αp1^2) * log(ubound(Arb, x))
                    # Lower and upper bounds of (1 - exp(t)) / t
                    # Using that t <= 0 to handle the singular cases
                    lower = Arblib.isnegative(tᵤ) ? (1 - exp(tᵤ)) / tᵤ : -one(tᵤ)
                    upper = isfinite(tₗ) ? (1 - exp(tₗ)) / tₗ : zero(tₗ)

                    w_mul_α * Arb((lower, upper))
                end

                T1231 + T1232
            end

            T12 = T121 * T122 * T123

            a0_div_gamma_1pα * (T11 + T12)
        end

        # Enclosure of the two leading terms in the expansion of
        # -u0.v0.a[j] * clausenc(x, 1 - α - u0.v0.α + j * u0.v0.p0)
        # for j = 1:skip_singular_j_until
        T2s = map(1:skip_singular_j_until) do j
            let r = -u0.v0.α + j * u0.v0.p0 - 1
                # Enclosure of
                # (gamma(α - 1 - r) * cospi((1 - α + r) / 2) - zeta(-α + r) / 2) * x^r / log(x)
                T21 = T21_α[j] * abspow(x, r) * invlogx

                # Enclosure of zeta_deflated(-α + r) * (x^r - x^(1 + α)) / 2log(x)
                T221 =
                    zeta_deflated_mαpr[j] * (abspow(x, r) - abspow(x, αp1)) * invlogx / 2

                # Enclosure of (x^r - x^(1 + α)) / (r - (1 + α)) / 2log(x)
                # It is non-decreasing in α so we evaluate at the endpoints
                # TODO: Improve this enclosure for wide x, in
                # particular overlapping zero
                T222 =
                    Arb((
                        (abspow(x, r) - 1) / r,
                        (abspow(x, r) - abspow(x, u0.ϵ)) / (r - u0.ϵ),
                    )) * invlogx / 2

                term = T21 + T221 + T222

                -u0.v0.a[j] * term
            end
        end

        T2 = sum(T2s, init = zero(x)) * invgamma1mxp0

        # Enclosure of the remaining terms in the expansion
        T3 = eval_expansion(u0, Du0_expansion_div_x_onemα, x) * invlogx * invgamma1mxp0

        # (u0(x)^2 / 2 + Hu0x) / (log(x) * gamma(1 + α) * x^(1 - α) * (1 - x^p0))
        F = T1 + T2 + T3

        return f1(x) * f2(x) * F
    end
end

"""
    inv_u0_bound(u0::BHKdVAnsatz{Arb})

Return a function `f` such that `f(x)` gives an upper bound of
```
gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x)
```

It assumes that `x` is non-negative, any negative parts of `x` are
ignored.

For non-zero `x` it computes an enclosure of the value. It splits it as
```
gamma(1 + α) * (1 - x^p0) / (u0(x) / x^(-α))
```
and computes an enclosure of
```
gamma(1 + α) * (1 - x^p0) = gamma(2 + α) * (1 - x^p0) / (1 + α)
```
by handling the removable singularity at `α = -1` using
[`fx_div_x`](@ref). For `u0(x) / x^(-α)` it computes an enclosure
using the asymptotic expansion of `u0`.

We now describe how to compute an upper bound when `x` overlaps with
zero.

We start by computing the asymptotic expansion of `u0`. We then split
the expansion into the leading term and a tail. In practice both the
leading term and the tail are positive and the tail is much smaller
than the leading term. Since we are only interested in an upper bound
we can thus just skip the tail. For this to be valid we must ensure
that both the leading term and the tail are positive. For the tail
this is checked using [`expansion_ispositive`](@ref). For the leading
term this is done by checking that the end result is positive (it is
easily seen that the numerator is positive in our case)

What remains is to handle the leading term, it is given by
```
a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0) * x^-α
```
Since we are interested in `gamma(1 + α) * x^-α * (1 - x^p0)` divided
by this we get
```
gamma(1 + α) * (1 - x^p0) /
    (a0 * (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0))
```
We split this further into the two factors
```
F1 = gamma(1 + α) / a0
F2 = (1 - x^p0) /
    (gamma(α) * cospi(α / 2) - gamma(α - p0) * cospi((α - p0) / 2) * x^p0)
```
For `F1` we can write it as
```
F1 = gamma(1 + α) / a0 = gamma(2 + α) / (a0 * (1 + α))
```
and use [`finda0αp1`](@ref).

# `F2`
For `F2` we let `c(a) = gamma(a) * cospi(a / 2)` and then factor it
out, giving us
```
inv(c(α)) * (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
```
We can handle the removable singularity of `c(a)` at `a = -1` by
rewriting it as
```
c(a) = π * gamma(2 + a) * sinc((1 + a) / 2) / 2a
```

For the remaining part we compute an enclosure of the derivative of
`c(α - p0) / c(α)` with respect to `α` and check that this is
negative. This means that `c(α - p0) / c(α)` is decreasing in `α` and
the maximum value is hence attained at `α = -1` where it is `1`. This
means that
```
(1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
```
is upper bounded by `(1 - x^p0) / (1 - x^p0) = 1`. Since `c(α - p0) /
c(α) <= 1` we also get that
```
(1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
```
is non-increasing in `x`. A lower bound can thus be computed by
considering `xᵤ = ubound(x)`. At `xᵤ > 0` we can handle the removable
singularity in `α`
"""
function inv_u0_bound(u0::BHKdVAnsatz{Arb}; M::Integer = 3, ϵ::Arb = Arb(0.5))
    @assert ϵ < 1

    # Interval for α
    α = Arb((-1, -1 + u0.ϵ))
    # Interval for α + 1
    αp1 = Arblib.nonnegative_part!(zero(Arb), Arb((0, u0.ϵ)))
    extra_degree = 2

    # Enclosure of F1
    F1 = ArbExtras.enclosure_series(α) do α
        gamma(2 + α) / finda0αp1(α)
    end

    c(a) = Arb(π) * gamma(2 + a) * _sinc((1 + a) / 2) / 2a

    inv_c_α = ArbExtras.enclosure_series(inv ∘ c, α)

    # Prove that c(α - p0) / c(α) is decreasing in α
    cαmp0_div_cα = let α = ArbSeries((α, 1)), p0 = (1 + α) + (1 + α)^2 / 2
        c(α - p0) / c(α)
    end

    Arblib.isnegative(cαmp0_div_cα[1]) ||
        error("c(α - p0) / c(α) could not be proved to be decreasing")

    # Compute the expansion of u0
    u0_expansion = u0(ϵ, AsymptoticExpansion(); M)

    # Temporarily remove the leading term and check that the rest of
    # the expansion is positive, so that we can remove it from the
    # denominator and still get an upper bound.
    leading_term = u0_expansion[(1, 0, 0, 0, 0, 0, 0)]
    delete!(u0_expansion, (1, 0, 0, 0, 0, 0, 0))
    expansion_ispositive(u0, u0_expansion, ϵ) ||
        error("expansion of u0 could not prove to be positive")
    u0_expansion[(1, 0, 0, 0, 0, 0, 0)] = leading_term

    # Divide the expansion of u0 by x^-α
    u0_expansion_div_x_mα = empty(u0_expansion)
    for ((p, q, i, j, k, l, m), y) in u0_expansion
        u0_expansion_div_x_mα[(p, q, i + 1, j, k, l, m)] = y
    end

    return x::Arb -> begin
        @assert x <= ϵ

        if Arblib.contains_zero(x)
            # Enclose F2
            F2 = if iszero(x)
                # We have (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0) = 1
                inv_c_α
            else
                # Compute lower and upper bounds of
                # (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
                # IMPROVE: Compute tighter enclosures for very small x
                lower = let xᵤ = ubound(Arb, x)
                    # Enclosure of (1 - x^p0) / (α + 1)
                    numerator_div_αp1 =
                        fx_div_x(s -> 1 - xᵤ^(s + s^2 / 2), αp1; extra_degree)
                    # Enclosure of (1 - c(α - p0) / c(α) * x^p0) / (α + 1)
                    denominator_div_αp1 = fx_div_x(
                        s -> 1 - c(-1 - s^2 / 2) / c(s - 1) * xᵤ^(s + s^2 / 2),
                        αp1;
                        extra_degree,
                        force = true,
                    )
                    numerator_div_αp1 / denominator_div_αp1
                end
                upper = one(Arb)

                inv_c_α * Arb((lower, upper))
            end

            F = F1 * F2

            @assert Arblib.ispositive(F) # Should always be the case
        else
            # Enclosure of gamma(1 + α) * (1 - x^p0) = gamma(2 + α) * (1 - x^p0) / (1 + α)
            numerator =
                gamma(2 + α) * fx_div_x(s -> (1 - x^(s + s^2 / 2)), αp1; extra_degree)

            F = numerator / eval_expansion(u0, u0_expansion_div_x_mα, x)
        end

        return F
    end
end
