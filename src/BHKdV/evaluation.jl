"""
    _exponent!(exponent::Arb, α::Arb, i, j, k, l, m, buffer1::Arb, buffer2::Arb, α_upper::Arb, p0_upper::Arb, u0)

Compute the exponent used in [`eval_expansion`](@ref) using inplace
computations. It computes
```
i * α + j * p0 - k * u0.v0.α + l * u0.v0.p0 + m
```

The variables `α_upper` and `p0_upper` should be "global" upper bounds
of `α` and `p0` respectively. Global in the sense that it should hold
for `α = Arb((-1, -1 + u0.ϵ))` and not for the `α` argument given to
the function, which might be a subset of those.

This function is only used internally in [`eval_expansion`](@ref) but
is defined outside of the function for performance reasons. Julia
doesn't optimize properly when it is put inside the function for some
reason.
"""
_exponent!(
    exponent::Arb,
    α::Arb,
    i,
    j,
    k,
    l,
    m,
    α_upper::Arb,
    p0_upper::Arb,
    u0,
    buffer1::Arb,
    buffer2::Arb,
) = begin
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
        Arblib.set!(buffer1, m - i)

        # Upper bound i * α_upper + j * p0_upper + m
        Arblib.set!(buffer2, m)
        Arblib.addmul!(buffer2, α_upper, i)
        Arblib.addmul!(buffer2, p0_upper, j)

        Arblib.union!(exponent, buffer1, buffer2)

        # If the lower bound is zero we want to avoid any spurious
        # negative parts
        iszero(buffer1) && Arblib.nonnegative_part!(exponent, exponent)
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

_exponent!(
    exponent::ArbSeries,
    α::ArbSeries,
    i,
    j,
    k,
    l,
    m,
    _::Arb,
    _::Arb,
    u0,
    buffer1::ArbSeries,
    buffer2::Arb,
) = begin
    # Compute the part i * α + j * p0
    # Note that i * α + j * p0 = (i + j) * α + j * (1 + (1 + α)^2 / 2)

    # exponent = j * (1 + (1 + α)^2 / 2)
    Arblib.add!(buffer1, α, 1)
    Arblib.pow_arb_series!(exponent, buffer1, Arblib.set!(buffer2, 2), length(exponent))
    Arblib.mul_2exp!(exponent, exponent, -1)
    Arblib.add!(exponent, exponent, 1)
    Arblib.mul!(exponent, exponent, Arblib.set!(buffer2, j))

    # exponent += (i * j) * α
    Arblib.mul!(buffer1, α, Arblib.set!(buffer2, i + j))
    Arblib.add!(exponent, exponent, buffer1)

    exponent0 = Arblib.ref(exponent, 0)
    Arblib.submul!(exponent0, u0.v0.α, k)
    Arblib.addmul!(exponent0, u0.v0.p0, l)
    Arblib.add!(exponent0, exponent0, m)
    if iszero(length(exponent.poly)) && !iszero(exponent0)
        # Manually update the degree
        Arblib.cstruct(exponent).length = 1
    end

    return exponent
end

"""
    eval_expansion(u0::BHKdVAnsatz, expansion, x; div_logx = false)

Evaluate the given expansion. If `div_logx` is true then divide the
result by `log(x)`, in case `x` overlaps with zero this allows for
cancellations that improve the enclosure.

It requires that `0 <= x < 1`, any negative parts of `x` are ignored.

The terms are stored as `((p, q, i, j, k, l, m), y)`. The parameters `(i, j,
k l, m)` correspond to the term
```
y * x^(i * α + j * p0 - k*u0.v0.α + l*u0.v0.p0 + m)
```
where `α ∈ (-1, -1 + u0.ϵ]` and `p0 = 1 + α + (1 + α)^2 / 2`.

Let `c(s) = gamma(s) * cospi(s / 2)`. The parameter `p` corresponds to
multiplication by the factor
```
a0 * (c(α) - c(α - p0) * x^p0) * x^-α
```
to the power `p`, which is the part of the expansion for the main term
which is not even powers of `x`. The parameter `q` corresponds to
multiplication by the factor
```
-a0 * (
    c(2α) * x^(-2α) - c(2α - p0) * x^(-2α + p0) +
    (-zeta(1 - 2α - 2) / 2 + zeta(1 - 2α + p0 - 2) / 2) * x^2
)
```
to the power `q`, which is the part of the expansion for `H` applied
to the main term which is not on the form `x^2m` for `m >= 2`.

# Handling of `p` and `q`
The method currently only support `q = 0` and `p == 0` or `p == 1`.
The other cases are handled specially in [`F0`](@ref) and for that
reason we don't bother implementing them here.

For `x > 0` the factor
```
a0 * (c(α) - c(α - p0) * x^p0)
```
has a removable singularity at `α = -1`. To handle this we rewrite it
as
```
a0 * (α + 1) * (c(α) - c(α - p0) * x^p0) / (α + 1)
```
and use [`finda0αp1`](@ref) to compute `a0 * (α + 1)`. For `x > 0` we
can evaluate `(c(α) - c(α - p0) * x^p0) / (α + 1)` using
[`fx_div_x`](@ref).

To better handle `x` overlapping zero we follow in this case
[`lemma_bhkdv_F0_P_factor`](@ref) and add and subtract `c(α) * x^p0`
to the numerator, giving us
```
(c(α) - c(α - p0) * x^p0) / (α + 1) =
    c(α) * (1 - x^p0) / (α + 1) + x^p0 * (c(α) - c(α - p0)) / (1 + α)
```
We can write
```
(1 - x^p0) / (α + 1) = -(1 + (1 + α) / 2) * (x^p0 - 1) / p0
```
to get
```
(c(α) - c(α - p0) * x^p0) / (α + 1) =
    x^p0 * (c(α) - c(α - p0)) / (1 + α) - c(α) * (1 + (1 + α) / 2) * (x^p0 - 1) / p0
```
If we let
```
C1 = a0 * (α + 1) * (c(α) - c(α - p0)) / (1 + α)
C2 = a0 * (α + 1) * c(α) * (1 + (1 + α) / 2)
```
then we have
```
a0 * (α + 1) * (c(α) - c(α - p0) * x^p0) / (α + 1) =
    C1 * x^p0 - C2 * (x^p0 - 1) / p0
```
The factor `(x^p0 - 1) / p0` is not finite near `x = 0` and we
therefore have to take into account the multiplication by `x^s` when
evaluating, we then have
```
x^s * a0 * (α + 1) * (c(α) - c(α - p0) * x^p0) / (α + 1) =
    C1 * x^(s + p0) - C2 * x^s * (x^p0 - 1) / p0
```

We can use [`x_pow_s_x_pow_t_m1_div_t`](@ref) to evaluate `x^s * (x^p0
- 1) / p0`, however this gives very poor results when `s` is very
small. If `div_logx` is true this can be handled in a much better way.
From the documentation of [`x_pow_s_x_pow_t_m1_div_t`](@ref) we have
that `x^s * (x^p0 - 1) / p0` is increasing in `p0`. A lower bound is
hence given by letting `p0 = 0`, for which we get `x^s * log(x)`.
Since `x^p0 - 1 <= 0` for `0 <= x < 1` and `p0 > 0` a trivial upper
bound is given by `0`. We thus have the enclosure `[x^s * log(x), 0]`
and when `div_logx` is true we can explicitly divide by `log(x)` to
get the enclosure `[0, x^s]`.

# Terms not depending on `p` or `q`
The terms for which `p == q == 0` we take out and handle in a separate
sum. We compute an expansion in `α` to get a better enclosure, in most
cases this picks up the monotonicity in `α`.
"""
function eval_expansion(
    u0::BHKdVAnsatz{Arb},
    expansion::AbstractDict{NTuple{7,Int},Arb},
    x::Arb;
    div_logx = false,
)
    @assert x < 1

    # We only care about the non-negative part of x
    x = Arblib.nonnegative_part!(zero(x), x)

    # Enclosure of α, α + 1, p0 and a0 * (α + 1)
    α = Arb((-1, -1 + u0.ϵ))
    αp1 = Arblib.nonnegative_part!(zero(u0.ϵ), union(zero(u0.ϵ), u0.ϵ))
    p0 = Arblib.nonnegative_part!(zero(α), αp1 + αp1^2 / 2)

    # Precomputed values used in _exponents!, upper bounds of α and p0
    α_upper = -1 + u0.ϵ
    p0_upper = u0.ϵ + u0.ϵ^2 / 2

    # Function for computing sum of terms with p = q = 0 for a given α
    # It is only called with α::ArbSeries if x is non-zero
    logx = log(x)
    f(α) =
        let exponent = zero(α), term = zero(α), buffer1 = zero(α), buffer2 = zero(Arb)
            S = zero(α)
            for ((p, q, i, j, k, l, m), y) in expansion
                if !iszero(y) && iszero(p) && iszero(q)
                    _exponent!(
                        exponent,
                        α,
                        i,
                        j,
                        k,
                        l,
                        m,
                        α_upper,
                        p0_upper,
                        u0,
                        buffer1,
                        buffer2,
                    )
                    if α isa Arb
                        abspow!(term, x, exponent)
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
    res1 = Arblib.contains_zero(x) ? f(α) : ArbExtras.enclosure_series(f, α)

    # c(s) = gamma(s) * cospi(s / 2) = π * gamma(s + 2) * sinc(s + 2) / 2s
    c(s) = begin
        sinc_sp2 = _sinc((1 + s) / 2)

        if s isa ArbSeries && round(Float64(s[0])) == -1
            # _sinc performs poorly for arguments close to zero, it is
            # often better to slightly widen the argument to include
            # zero
            t = (1 + s) / 2
            t[0] = union(t[0], Arb(0))
            sinc_sp2_wide = _sinc(t)

            for i = 0:Arblib.degree(s)
                sinc_sp2[i] = intersect(sinc_sp2[i], sinc_sp2_wide[i])
            end
        end

        return π * gamma(s + 2) * sinc_sp2 / 2s
    end

    # Compute enclosure of the p-coefficient
    # a0 * (c(α) - c(α - p0) * x^p0)
    # Note that this depends on x and is only used when x doesn't contain zero
    p_coefficient = if !Arblib.contains_zero(x)
        ArbExtras.enclosure_series(α) do α
            # Enclosure of (c(α) - c(α - p0) * x^p0) / (α + 1)
            coefficient_div_α =
                if (α isa Arb && contains(α, -1)) || (α isa ArbSeries && contains(α[0], -1))
                    fx_div_x(α + 1, force = true, extra_degree = 2) do r::ArbSeries
                        let α = r - 1, p0 = r + r^2 / 2
                            c(α) - c(α - p0) * x^p0
                        end
                    end
                else
                    coefficient_div_α = let p0 = 1 + α + (1 + α)^2 / 2
                        (c(α) - c(α - p0) * x^p0) / (1 + α)
                    end
                end

            finda0αp1(α) * coefficient_div_α
        end
    else
        indeterminate(α)
    end

    # Compute enclosures of C1 and C2, these are only used if x
    # overlaps zero.
    if Arblib.contains_zero(x)
        # Enclosure of a0 * (α + 1) * (c(α) - c(α - p0)) / (1 + α)
        C1 = ArbExtras.enclosure_series(α) do α
            # Enclosure of (c(α) - c(α - p0)) / (α + 1)
            coefficient_div_α =
                if (α isa Arb && contains(α, -1)) || (α isa ArbSeries && contains(α[0], -1))
                    fx_div_x(α + 1, force = true, extra_degree = 2) do r::ArbSeries
                        let α = r - 1, p0 = r + r^2 / 2
                            c(α) - c(α - p0)
                        end
                    end
                else
                    coefficient_div_α = let p0 = 1 + α + (1 + α)^2 / 2
                        (c(α) - c(α - p0)) / (1 + α)
                    end
                end

            finda0αp1(α) * coefficient_div_α
        end

        # Enclosure of a0 * (α + 1) * c(α) * (1 + (1 + α) / 2)
        C2 = ArbExtras.enclosure_series(α) do α
            finda0αp1(α) * c(α) * (1 + (1 + α) / 2)
        end
    else
        C1, C2 = indeterminate(α), indeterminate(α)
    end

    # Sum of terms with either p or q non-zero
    res2 = zero(x)
    let term = zero(α), exponent = zero(α), buffer1 = zero(α), buffer2 = zero(α)
        for ((p, q, i, j, k, l, m), y) in expansion
            if !iszero(y) && !(iszero(p) && iszero(q))
                # We don't implement any other cases since they are
                # handled specially by F0.
                iszero(p) ||
                    isone(p) ||
                    throw(ArgumentError("only p == 0 or p == 1 supported, got p = $p"))
                iszero(q) || throw(ArgumentError("only q == 0 supported, got q = $q"))

                # Can assume that p = 1 and q = 0

                # Add -α to the exponent coming from the p factor
                _exponent!(
                    exponent,
                    α,
                    i - 1,
                    j,
                    k,
                    l,
                    m,
                    α_upper,
                    p0_upper,
                    u0,
                    buffer1,
                    buffer2,
                )

                if Arblib.contains_zero(x)
                    # term = C1 * abspow(x, exponent + p0)
                    Arblib.add!(term, exponent, p0)
                    abspow!(term, x, term)
                    Arblib.mul!(term, term, C1)

                    # term -= C2 * x_pow_s_x_pow_t_m1_div_t(x, exponent, p0)
                    # Possibly dividing by log(x) if div_logx is true.
                    D = if div_logx
                        # Perform the division by log(x) directly
                        Arb((zero(x), abspow(x, exponent)))
                    else
                        Arb((logabspow(x, 1, exponent), zero(x)))
                    end
                    Arblib.submul!(term, C2, D)
                else
                    abspow!(term, x, exponent)
                    Arblib.mul!(term, term, p_coefficient)
                end

                #res += y * term
                Arblib.addmul!(res2, y, term)
            end
        end
    end

    if div_logx
        # Enclosure of inv(log(x))
        invlogx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            Arb((inv(log(ubound(Arb, x))), 0))
        else
            inv(log(x))
        end

        if Arblib.contains_zero(x)
            # In this case the division by log(x) is done directly
            # when computing res2
            return res1 * invlogx + res2
        else
            return (res1 + res2) * invlogx
        end
    else
        return res1 + res2
    end
end

"""
    expansion_ispositive(u0::BHKdVAnsatz, expansion, ϵ)

Attempt to prove that the given expansion is positive on the interval
``(0, ϵ]``. Returns true on success and false on failure. It requires
that `0 < ϵ < 1`. This is used in [`inv_u0_normalised`](@ref) for
computing an enclosure.

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
then a lower bound for the corresponding term is given by taking `l =
1` in the exponent, i.e.
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
    expansion_1 = filter(expansion) do ((p, q, i, j, k, l, m), _)
        p == q == i == j == m == 0 && k == 1 && l >= 1
    end

    # Isolate all keys of the from (0, 0, 0, 0, 0, 0, m)
    expansion_2 = filter(expansion) do ((p, q, i, j, k, l, m), _)
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

Compute `u0(x)`. It is based on [`equation_bhkdv_u0`](@ref) but some
care has to be taken when evaluating it.

The two sums are evaluated directly.

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

            ArbExtras.compose_zero(C₀, x)
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

Compute an expansion of `u0` around zero with remainder terms valid
on the interval ``[-abs(x), abs(x)]``.

The expansion is returned as a dictionary, which can then be evaluated
with [`eval_expansion`](@ref). See [`eval_expansion`](@ref) for more
details about how the coefficients are stored.

The expansion is based on [`lemma_bhkdv_asymptotic_expansion`](@ref).
In the lemma expansions from the sum of Clausen terms and the sum of
Fourier terms are combined. For the computation we compute the
expansions for these two sums separately and then combine them.

The infinite sum in the Lemma is truncated and the tail is combined
into one remainder term with the coefficient `x^2M`.

For the terms coming from `u0.v0` the expansions are easily computed
exactly like for `BHAnsatz`. For the main term we have to be a bit
more careful.

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
We can rewrite it as
```
a0 * (1 + α) * sum((-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / (1 + α) * x^2m / factorial(2m) for m = M:Inf) / x^2M
```
We can enclose `a0 * (1 + α)` using [`finda0αp1`](@ref). What remains
is to bound the absolute value of
```
S = sum((-1)^m * (zeta(1 - α - 2m) - zeta(1 - α + p0 - 2m)) / (1 + α) * x^(2m - 2M) / factorial(2m) for m = M:Inf)
```
From [`lemma_clausen_remainder_bhkdv`](@ref) we get
```
abs(S) <= abs(S1) + (1 + α) / 2 * abs(S2)
```
with
```
S1 = sum((-1)^m * dzeta(s1[m] - 2m) * x^2m / factorial(2m) for m = M:Inf)
S2 = sum((-1)^m * dzeta(s2[m] - 2m) * x^2m / factorial(2m) for m = M:Inf)
```
with `1 - α <= s1[m] <= 2` and `2 <= s2[m] <= 2 + (1 + α)^2 / 2`.

From [`lemma_clausen_derivative_remainder_interval`](@ref) we get that
`S1` and `S2` can be bounded in the same way as the sums bounded by
[`clausenc_expansion_remainder`](@ref). The only difference being that
we have to give interval arguments for `s` that enclose `s1[m]` and
`s2[m]` respectively. This is straight forward since we have uniform
bounds for `s1[m]` and `s2[m]`.
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
        zeta_div_αp1 = if m == 1
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

        coefficient = a0αp1 * (-1)^m * zeta_div_αp1 / factorial(2m)

        res[(0, 0, 0, 0, 0, 0, 2m)] += coefficient
    end

    # Remainder term for main term
    # s1 and s2 are enclosures of s1[m] and s2[m] for all m >= M.
    remainder = let s1 = 2 - αp1, s2 = 2 + αp1^2 / 2
        a0αp1 * (
            clausenc_expansion_remainder(x, s1, 1, M) +
            αp1 / 2 * clausenc_expansion_remainder(x, s2, 1, M)
        )
    end
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

Return a function for computing `H(u0)(x)` as given in
[`equation_bhkdv_Hu0`](@ref).

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

                ArbExtras.compose_zero(C₀, x)
            end
        end

        res = -a0αp1 * C

        # Tail term

        let α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α
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

# H(u0, Asymptotic()) doesn't exist because eval_expansion doesn't
# implement evaluation of the main term of H directly.

"""
    H(u0::BHKdVAnsatz, ::AsymptoticExpansion; M = 3, skip_singular_j_until::Integer = 0,)

Return a function such that `H(u0)(x)` computes an expansion of
`H(u0)` around zero with remainder terms valid on the interval
``[-abs(x), abs(x)]``.

The expansion is returned as a dictionary, which can then be evaluated
with [`eval_expansion`](@ref). See [`eval_expansion`](@ref) for more
details about how the coefficients are stored.

The expansion is based on [`lemma_bhkdv_asymptotic_expansion`](@ref).
In the lemma expansions from the sum of Clausen terms and the sum of
Fourier terms are combined. For the computation we compute the
expansions for these two sums separately and then combine them.

The infinite sum in the Lemma is truncated and the tail is combined
into one remainder term with the coefficient `x^2M`.

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
We can rewrite it as
```
a0 * (1 + α) * sum((-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / (1 + α) * x^2m / factorial(2m) for m = M:Inf) / x^2M
```
We can enclose `a0 * (1 + α)` using [`finda0αp1`](@ref). What remains
is to bound the absolute value of
```
S = sum((-1)^m * (zeta(1 - 2α - 2m) - zeta(1 - 2α + p0 - 2m)) / (1 + α) * x^(2m - 2M) / factorial(2m) for m = M:Inf)
```
From [`lemma_clausen_remainder_bhkdv`](@ref) we get
```
abs(S) <= 2abs(S1) + (1 - (1 + α) / 2) * abs(S2)
```
with
```
S1 = sum((-1)^m * dzeta(s3[m] - 2m) * x^2m / factorial(2m) for m = M:Inf)
S2 = sum((-1)^m * dzeta(s4[m] - 2m) * x^2m / factorial(2m) for m = M:Inf)
```
with `1 - 2α <= s1[m] <= 3` and `2 - α + (1 + α)^2 / 2 <= s2[m] <= 3`.

From [`lemma_clausen_derivative_remainder_interval`](@ref) we get that
`S1` and `S2` can be bounded in the same way as the sums bounded by
[`clausenc_expansion_remainder`](@ref). The only difference being that
we have to give interval arguments for `s` that enclose `s3[m]` and
`s4[m]` respectively. This is straight forward since we have uniform
bounds for `s3[m]` and `s4[m]`.

# Terms from `u0.v0`
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
"""
function H(
    u0::BHKdVAnsatz{Arb},
    ::AsymptoticExpansion;
    M::Integer = 3,
    skip_singular_j_until::Integer = 0,
)
    @assert M >= 3

    skip_singular_j_until > u0.v0.N0 && throw(
        ArgumentError(
            "can't skip more j-terms than there are, skip_singular_j_until = $skip_singular_j_until, N0 = $(u0.v0.N0)",
        ),
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
        # s3 and s4 are enclosures of s3[m] and s4[m] for all m >= M.
        remainder = let s3 = 3 - 2αp1, s4 = 3 - αp1 + αp1^2
            a0αp1 * (
                2clausenc_expansion_remainder(x, s3, 1, M) +
                (1 - αp1 / 2) * clausenc_expansion_remainder(x, s4, 1, M)
            )
        end
        res[(0, 0, 0, 0, 0, 0, 2M)] += remainder

        # Tail term

        # Clausen terms
        for j = 1:u0.v0.N0
            s = 1 - α - u0.v0.α + j * u0.v0.p0
            C, _, p, E = clausenc_expansion(x, s, M, skip_constant = true)

            if j > skip_singular_j_until
                res[(0, 0, -1, 0, 1, j, 0)] = -C * u0.v0.a[j]
                res[(0, 0, 0, 0, 0, 0, 2)] -= p[2] * u0.v0.a[j]
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

# defect(u0, Asymptotic()) doesn't exist because eval_expansion
# doesn't implement evaluation of the main term of H directly.

"""
    defect(u0::FractionalKdVAnsatz{Arb}, ::AsymptoticExpansion; M = 5, skip_singular_j_until = 0)

Return a function such that `defect(u0, AsymptoticExpansion())(x)`
compute an expansion of
```
H(u0)(x) + u0(x)^2 / 2
```
around zero with a remainder term valid on the interval `[-abs(x),
abs(x)]`.

The expansion is returned as a dictionary, which can then be evaluated
with [`eval_expansion`](@ref).

The expansion is computed by first computing the expansions of `u0(x)`
and `H(u0)(x)`. The expansion for `u0(x)^2 / 2` is then computed from
the expansion of `u0(x)` by computing the products between all terms.
Finally the expansion of `H(u0)(x) + u0(x)^2 / 2` is given by simply
adding the two expansions together.
"""
function defect(
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

        # u0^2 / 2 term
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

Return a function that computes a value that is larger than
`F0(u0)(x)` in magnitude and has the same sign. The absolute value of
this then gives an upper bound of `abs(F0(u0)(x))`.

It uses [`lemma_bhkdv_main_term_limit`](@ref) which implies that
`u0.v0(x)` gives a lower bound of `u0(x)`. Hence
```
defect(u0)(x) / (u0.w(x) * u0.v0(x))
```
this gives a value which has the same sign as `F0(x)` but is larger in
magnitude. This holds as long as `u0.v0(x)` is positive at least,
which is easily checked.

**IMPROVE:** The computation of `u0.v0(x)` and `u0(x)` have many
calculations in common. We could improve performance by using this.
"""
function F0_bound(u0::BHKdVAnsatz{Arb}, evaltype::Ball = Ball())
    # Assert that the lemma holds
    @assert lemma_bhkdv_main_term_limit(u0)

    g = defect(u0, evaltype)

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

Returns a function that computes a value which absolute value is an
**upper bound** of `abs(F0(u0)(x))`. It is computed in a way that is
works well for small values of `x`.

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
all bounded as `x -> 0` that we bound separately. We write it as
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
See [`equation_bhkdv_F0_factor`](@ref) where the factor is introduced
in the paper. We now explain how to bound `F`.

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

For `T111` one we can directly get an enclosure using that
```
-α + p0 - 1 = -α + (1 + α + (1 + α)^2 / 2) - 1 = (1 + α)^2 / 2
```

For `T112` we use that it is non-decreasing in both `x` and `α`. For
`x` this is immediate. For `α` we let `t = 1 + α`, giving us
```
T112 = t / (1 - x^(t + t^2 / 2))
```
Differentiation w.r.t. `t` gives
```
(1 - x^(t + t^2 / 2) + t * (1 + t) * log(x) * x^(t + t^2 / 2)) / (1 - x^(t + t^2 / 2))^2
```
The sign depends only on the numerator
```
1 + (t * (1 + t) * log(x) - 1) * x^(t + t^2 / 2)
```
Differentiating the numerator we get
```
t * log(x) * (1 + (t + 1)^2 * log(x)) * x^(t + t^2 / 2)
```
This is non-negative as long as `1 + (t + 1)^2 * log(x)` is negative,
which holds whenever `x < exp(-1)`. So, for `x < exp(-1)` the
numerator
```
1 + (t * (1 + t) * log(x) - 1) * x^(t + t^2 / 2)
```
is non-decreasing in `t`, a lower bound is hence given at `t = 0`
where we get `0`. It follows that the derivative of `T112` w.r.t. `α`
is non-negative and hence it is non-decreasing in `α`. A lower bound
for `T112` is thus given by `(1 + α) / (1 - x^p0)` evaluated at `α =
-1`, where it can be seen to be equal to `-inv(log(x))`. We get an
upper bound by evaluating at an upper bound for `x` and `α`. See also
[`x_pow_s_x_pow_t_m1_div_t`](@ref).

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
We can get an enclosure of `inv(gamma(1 + α) * (1 - x^p0))` by writing
it as `inv(gamma(2 + α)) * (1 + α) / (1 - x^p0)` and using the method
for enclosing `T112`. We are then interested in enclosing the rest.

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
The term `T221` can be enclosed directly.

For `T222` we rewrite it as
```
inv(2log(x)) * (x^r - x^(1 + α)) / (r - (1 + α))
```
We want to show that this is decreasing in `α`. Differentiating w.r.t.
`α` gives us
```
inv(2log(x)) * (x^r - (1 + log(x) * (r - (1 + α))) * x^(1 + α)) / (r - (1 + α))^2
```
The sign depends only on the sign of the numerator
```
x^r - (1 + log(x) * (r - (1 + α))) * x^(1 + α)
```
we want to prove that it is non-negative. Differentiating the
numerator w.r.t. `α` gives us
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
It follows that it is non-negative for all `α`. Since `log(x)` is
negative we get that the full derivative is non-positive and hence
```
T222 = inv(2log(x)) * (x^r - x^(1 + α)) / (r - (1 + α))
```
is non-increasing in `α` and we can evaluate it at the endpoints.
Furthermore the upper bound is increasing in `x`, the upper bound is
given by
```
inv(2log(x)) * (x^r - 1) / r
= inv(2) * (x^r - 1) / (r * log(x))
= inv(2) * (exp(r * log(x)) - 1) / (r * log(x))
```
which is increasing in `x`.

## Handling `T3`: the remaining terms
For the remaining terms we factor out `inv(gamma(1 + α) * (1 - x^p0))`
and enclose it separately, in the same way as for `T2` above.

For the terms in the expansion we note that they all have an exponent
greater than `1 - α` and we can therefore explicitly cancel the
division by `x^(1 - α)`. The remaining part of the denominator is
`log(1 / x)`, which is handled internally by [`eval_expansion`](@ref)
by giving it the argument `div_logx = true`. This makes use of
[`lemma_bhkdv_F0_P_factor`](@ref) for improving the computed
enclosures when `x` overlaps zero.
"""
function F0(
    u0::BHKdVAnsatz{Arb},
    ::Asymptotic;
    M::Integer = 3,
    ϵ::Arb = Arb(0.5),
    skip_singular_j_until = u0.v0.N0,
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

    # Compute the expansion of defect(u0), skipping the two leading
    # terms in the expansion of the Clausen term in the tail for j =
    # 1:skip_j_until, which are handled separately
    defect_u0_expansion = defect(u0, AsymptoticExpansion(); M, skip_singular_j_until)(ϵ)
    delete!(defect_u0_expansion, (2, 0, 0, 0, 0, 0, 0))
    delete!(defect_u0_expansion, (0, 1, 0, 0, 0, 0, 0))

    # Divide the expansion of defect(u0) by x^(1 - α)
    defect_u0_expansion_div_x_onemα = empty(defect_u0_expansion)
    for ((p, q, i, j, k, l, m), y) in defect_u0_expansion
        defect_u0_expansion_div_x_onemα[(p, q, i + 1, j, k, l, m - 1)] = y
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

        xᵤ = ubound(Arb, x)

        # Enclosure of inv(log(x))
        invlogx = if iszero(x)
            zero(x)
        elseif Arblib.contains_zero(x)
            Arb((inv(log(ubound(Arb, x))), 0))
        else
            inv(log(x))
        end

        # Compute an enclosure of (1 + α) / (1 - x^p0)
        αp1_div_onemxp0 = if iszero(x)
            αp1
        else
            if x < exp(Arb(-1)) # Use monotonicity on α and x
                lower = if Arblib.contains_zero(x)
                    αp1
                else
                    -invlogx
                end
                upper = let αp1ᵤ = ubound(Arb, αp1)
                    αp1ᵤ / (1 - xᵤ^(αp1ᵤ + αp1ᵤ^2 / 2))
                end
                Arb((lower, upper))
            else
                inv(fx_div_x(s -> (1 - x^(s + s^2 / 2)), αp1, extra_degree = 2))
            end
        end

        # Enclosure of inv(gamma(1 + α) * (1 - x^p0))
        # = gamma(2 + α) * (1 + α) / (1 - x^p0)
        invgamma1mxp0 = gamma(2 + α) * αp1_div_onemxp0

        # Enclosure for the terms P and Q
        T1 = let
            # Compute an enclosure of T11

            # Enclosure of T111 = x^(-α + p0 - 1) / log(x) using that
            # -α + p0 - 1 = (1 + α)^2 / 2
            T111 = abspow(x, Arblib.nonnegative_part!(zero(x), αp1^2 / 2)) * invlogx

            # Enclosure of T112 = (1 + α) / (1 - x^p0)
            T112 = αp1_div_onemxp0

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
        x_pow_u0ϵ = abspow(x, u0.ϵ)
        x_pow_αp1 = abspow(x, αp1)
        invlogx_div_2 = Arblib.mul_2exp!(zero(x), invlogx, -1)
        T2s = map(1:skip_singular_j_until) do j
            let r = -u0.v0.α + j * u0.v0.p0 - 1
                x_pow_r = abspow(x, r)

                # Enclosure of
                # (gamma(α - 1 - r) * cospi((1 - α + r) / 2) - zeta(-α + r) / 2) * x^r / log(x)
                T21 = T21_α[j] * x_pow_r * invlogx

                # Enclosure of zeta_deflated(-α + r) * (x^r - x^(1 + α)) / 2log(x)
                T221 = zeta_deflated_mαpr[j] * (x_pow_r - x_pow_αp1) * invlogx_div_2

                # Enclosure of (x^r - x^(1 + α)) / (r - (1 + α)) / 2log(x)
                # It is non-increasing in α so we evaluate at the endpoints
                T222_lower = if Arblib.contains_zero(x)
                    zero(x)
                else
                    # IMPROVE: This can be improved for wide x
                    (x_pow_r - x_pow_u0ϵ) / (r - u0.ϵ) * invlogx_div_2
                end

                T222_upper = if iszero(x)
                    zero(x)
                else
                    (abspow(xᵤ, r) - 1) / (2r * log(xᵤ))
                end
                T222 = Arb((T222_lower, T222_upper))
                # The enclosure for the lower bound is in general
                # worse and sometimes it gets bigger than the upper
                # bound. It is therefore to take the minimum with the
                # enclosure and the upper bound.
                T222 = Arblib.min!(T222, T222, T222_upper)

                # term = -u0.v0.a[j] * (T21 + T221 + T222)
                term = T21 + T221
                Arblib.add!(term, term, T222)
                Arblib.mul!(term, term, u0.v0.a[j])
                Arblib.neg!(term, term)
                term
            end
        end

        T2 = sum(T2s, init = zero(x)) * invgamma1mxp0

        # Enclosure of the remaining terms in the expansion
        T3 =
            eval_expansion(u0, defect_u0_expansion_div_x_onemα, x, div_logx = true) *
            invgamma1mxp0

        # (u0(x)^2 / 2 + Hu0x) / (log(x) * gamma(1 + α) * x^(1 - α) * (1 - x^p0))
        F = T1 + T2 + T3

        return f1(x) * f2(x) * F
    end
end

"""
    inv_u0_bound(u0::BHKdVAnsatz{Arb})

Return a function that computes an upper bound of
```
gamma(1 + α) * x^(-α) * (1 - x^p0) / u0(x)
```

The function bounds the expression appearing in
[`lemma_bhkdv_inv_u0_normalised`](@ref).

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
using the same method as for enclosing `T112` in [`F0`](@ref).

We now describe how to compute an upper bound when `x` overlaps with
zero. This is similar to what is done in the proof of
[`lemma_bhkdv_inv_u0_normalised`](@ref).

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
considering an upper bound of `x`, due to instabilities for very small
`x` we either take `xᵤ = ubound(x)` or `xᵤ = 1e-1000`, whichever is
larger. At `xᵤ > 0` we can handle the removable singularity in `α`.
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

        # If x is very small it is better to use the version for x
        # containing zero
        if x < Arb("1e-1000")
            x = union(x, Arb(0))
        end

        if Arblib.contains_zero(x)
            # Enclose F2
            F2 = if iszero(x)
                # We have (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0) = 1
                inv_c_α
            else
                # Compute lower and upper bounds of
                # (1 - x^p0) / (1 - c(α - p0) / c(α) * x^p0)
                lower = let xᵤ = ubound(Arb, x)
                    # We don't want x to be too small due to instabilities
                    xᵤ = max(xᵤ, Arb("1e-1000"))
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

            # If F is not positive we are not guaranteed that we have
            # an enclosure
            Arblib.ispositive(F) || return indeterminate(x)
        else
            # Compute an enclosure of inv((1 + α) / (1 - x^p0))
            inv_αp1_div_onemxp0 = if x < exp(Arb(-1)) # Use monotonicity in α and x
                # The lower and upper bounds refer to
                # (1 + α) / (1 - x^p0)
                # Since we want the inverse we invert them and switch
                # there order
                lower = -inv(log(x))
                upper = let xᵤ = ubound(Arb, x), αp1ᵤ = ubound(Arb, αp1)
                    αp1ᵤ / (1 - xᵤ^(αp1ᵤ + αp1ᵤ^2 / 2))
                end
                Arb((inv(upper), inv(lower)))
            else
                inv(fx_div_x(s -> (1 - x^(s + s^2 / 2)), αp1, extra_degree = 2))
            end

            # Enclosure of gamma(1 + α) * (1 - x^p0) = gamma(2 + α) * (1 - x^p0) / (1 + α)
            numerator = gamma(2 + α) * inv_αp1_div_onemxp0

            F = numerator / eval_expansion(u0, u0_expansion_div_x_mα, x)
        end

        return F
    end
end
