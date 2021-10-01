"""
    (u0::BHKdVAnsatz)(x, ::Ball)

Evaluate the ansatz `u0` at the point `x` using direct ball arithmetic
(not an asymptotic approach).

The tail term is evaluated directly.

To evaluate the main term, given by
```
a0 * (Ci(x, 1 - α) - Ci(x, 1 - α + p0) - (zeta(1 - α) - zeta(1 - α + p0)))
```
we make use of the fact that this converges to
```
u0.v0.a0 * (Ci(x, 2, 1) - zeta(2, d = 1))
```
, which is the main term for `BHAnsatz`, as `α -> -1`. We therefore
evaluate this function and add precomputed \$L^\\infty\$ bounds for
```
a0 * (Ci(x, 1 - α) - Ci(x, 1 - α + p0) - (zeta(1 - α) - zeta(1 - α + p0))) - u0.v0.a0 * (Ci(x, 2, 1) - zeta(2, d = 1))
```
valid for the entire range `α ∈ (-1, -1 + u0.ϵ]`.

**TODO:** Compute rigorous \$L^\\infty\$ bounds for the above
  expression. We currently use a heuristic value.
"""
function (u0::BHKdVAnsatz{Arb})(x, ::Ball)
    # Main term

    # Approximation
    res = 2 / Arb(π)^2 * (Ci(x, 2, 1) - zeta(Arb(2), d = 1))

    # Add error bounds
    @warn "L^∞ bounds not rigorously computed - using heuristic values" maxlog = 1
    if u0.ϵ <= 0.0003
        Arblib.add_error!(res, Mag(2e-4))
    else
        throw(ErrorException("no L^∞ bounds for ϵ = $ϵ"))
    end

    # Tail term

    # Clausen terms
    for j = 1:u0.v0.v0.N0
        s = 1 - u0.v0.v0.α + j * u0.v0.v0.p0
        res += u0.v0.v0.a[j] * (Ci(x, s) - zeta(s))
    end

    # Fourier terms
    for n = 1:u0.v0.N
        res += u0.v0.b[n] * (cos(n * x) - 1)
    end

    return res
end

"""
    H(u0::BHKdVAnsatz, ::Ball)

Returns a function such that `H(u0, Ball())(x)` gives an enclosure of
`H^-α[u0](x)` for all values of `α ∈ (-1, -1 + u0.ϵ]`. It uses direct
ball arithmetic (not an asymptotic approach).

The transform of the main term is given by
```
-a0 * (Ci(x, 1 - 2α) - Ci(x, 1 - 2α + p0) - (zeta(1 - 2α) - zeta(1 - 2α + p0)))
```
we make use of the fact that this converges to
```
-u0.v0.a0 * (Ci(x, 3, 1) - zeta(3, d = 1))
```
, which is the main term for `BHAnsatz`, as `α -> -1`. We therefore
evaluate this function and add precomputed \$L^\\infty\$ bounds for
```
-a0 * (Ci(x, 1 - 2α) - Ci(x, 1 - 2α + p0) - (zeta(1 - 2α) - zeta(1 - 2α + p0))) + u0.v0.a0 * (Ci(x, 3, 1) - zeta(3, d = 1))
```
valid for the entire range `α ∈ (-1, -1 + u0.ϵ]`.

**TODO:** Compute rigorous \$L^\\infty\$ bounds for the above
  expression. We currently use a heuristic value.

For the tail term we need to make sure that we correctly handle the
fact that the transform depends on the value of `α`.

For the Fourier terms we do this directly, the transformation takes
`cos(n * x)` to `-n^α * cos(n * x)` and in this case we just let `α`
be a ball containing `(-1, -1 + u0.ϵ]`.

For the Clausen functions we have to be a bit more careful. The
transformation takes `Ci(x, s)` to `-Ci(x, s - α)` but putting `α` as
a ball doesn't give good enclosures. It gives a large overestimations,
in particular when `s - α` is close to an integer, and fails when it
overlaps with an integer. Instead we use the approximation given by
taking the Hilbert transform and bounding the error. The Hilbert
transform of `Ci(x, s)` is given by `-Ci(x, s + 1)`. For the error we
currently use a heuristic, for `s = 2` and `ϵ = 0.0003` the error
maximum error in `x` is around `3e-4` and for larger values of `s` we
get a smaller error.

**TODO:** How to properly bound the error when using `-Ci(x, s + 1)`
instead of `-Ci(x, s - α)`?

**TODO:** We will probably have to improve on the enclosure to get a
sufficiently good defect in the end.

"""
function H(u0::BHKdVAnsatz{Arb}, ::Ball)


    @warn "L^∞ bounds for main term not rigorously computed - using heuristic values"
    @warn "L^∞ bounds for Clausen tail term not rigorously computed - using heuristic values"

    return x -> begin

        # Main term

        # Approximation
        res = -u0.v0.a0 * (Ci(x, 3, 1) - zeta(Arb(3), d = 1))

        # Add error bounds
        if u0.ϵ <= 0.0003
            Arblib.add_error!(res, Mag(1e-4))
        else
            throw(ErrorException("no L^∞ bounds for ϵ = $ϵ"))
        end

        # Tail term

        # Clausen terms
        for j = 1:u0.v0.v0.N0
            s = 2 - u0.v0.v0.α + j * u0.v0.v0.p0
            term = Ci(x, s) - zeta(s)
            if s >= 2 && u0.ϵ <= 0.0003
                Arblib.add_error!(term, Mag(3e-4))
            else
                # Error bounds only valid for s >= 2
                throw(ErrorException("no L^∞ bounds for s = $s and ϵ = $ϵ"))
            end
            res -= u0.v0.v0.a[j] * term
        end

        # Fourier terms
        let α = Arb((-1, -1 + u0.ϵ)) # Ball containing the range of α
            for n = 1:u0.v0.N
                res -= u0.v0.b[n] * n^α * (cos(n * x) - 1)
            end
        end

        return res
    end
end
