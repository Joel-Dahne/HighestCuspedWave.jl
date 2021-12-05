"""
    alpha0(u0::KdVZeroAnsatz; verbose)

Compute an upper bound of
```
u0.w(x) / 2u0(x)
```
for `0 < x < π`, valid for the full range of `α`.

Using that `u0.w(x) = x` we get
```
x / 2u0(x)
```
Next using the method [`u0_div_xmα`](@ref) we can write `u0(x)` as
```
u0(x) = x^-α * f(x)
```
giving us
```
x / (2x^-α * f(x)) = x^(1 + α) / 2f(x)
```
"""
function alpha0(u0::KdVZeroAnsatz{Arb}; verbose = false)
    f = u0_div_xmα(u0, Asymptotic(), ϵ = Arb(3.5), M = 5)

    # Function for computing x^(1 + α) / 2f(x)
    g(x) = abspow(x, 1 + u0.α) / 2f(x)

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
        verbose,
    )

    return m
end
