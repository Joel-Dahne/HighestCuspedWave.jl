"""
    n0_bound(u0::KdVZeroAnsatz; rtol, verbose)

Enclose the value of `n₀`. This is the supremum of
```
N(x) = u0.w(x) / 2u0(x)
```
for `0 < x < π`.

It uses an asymptotic expansion on the whole interval. The expansion
is given by computing the asymptotic expansion of `u0`. The terms in
the expansion are then evaluated in `α`, giving us an expansion in `x`
where the coefficients are just balls. We factor out `x^-α` from this
expansion and compute
```
x^(1 + α) / 2(x^α * u0(x))
```
"""
function n0_bound(u0::KdVZeroAnsatz; rtol = Arb(1e-3), verbose = false)
    # Compute expansion of u0 and evaluate the terms in α
    expansion = let expansion = u0(Arb(3.2), AsymptoticExpansion())
        res = empty(expansion, Arb)
        for (key, value) in expansion
            res[key] = value(u0.α)
        end
        res
    end

    # Function for computing x^(1 + α) / 2(x^α * u0(x))
    g(x) = abspow(x, 1 + u0.α) / 2eval_expansion(u0, expansion, x, offset_i = -1)

    # This will typically not reach the tolerance because of the
    # overestimations in g. It simply terminates when it reaches the
    # maximum number of evaluations. Since g is very fast to evaluate
    # this is not an issue.
    m = ArbExtras.maximum_enclosure(
        g,
        Arf(0),
        ubound(Arb(π)),
        abs_value = true,
        degree = 0;
        rtol,
        verbose,
    )

    return m
end
