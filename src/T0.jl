export T0

"""
    T0(u0::AbstractAnsatz, evaltype = Ball())
Returns a function such that T0(u0)(x) computes the integral T_0
from the paper. This is the function whose supremum on [0, π] gives
`D₀`.

In general it returns an **enclosure** but in some cases it only
returns an **upper bound**.

The strategy for evaluation depends on evaltype and the precise
ansatz.
"""
T0(u0::AbstractAnsatz; kwargs...) = T0(u0, Ball(); kwargs...)

function T0(u0::AbstractAnsatz, evaltype::Asymptotic; nonasymptotic_u0 = false)
    f = T01(u0, evaltype; nonasymptotic_u0)
    g = T02(u0, evaltype; nonasymptotic_u0)
    return x -> f(x) + g(x)
end

"""
    T01(u0::AbstractAnsatz, evaltype = Ball(); δ1, δ2)
Returns a function such that T01(u0; δ1, δ2)(x) computes the integral
T_{0,1} from the paper.

In general it returns an **enclosure** but in some cases it only
returns an **upper bound**.

The strategy for evaluation depends on evaltype and the precise
ansatz.
"""
T01(u0::AbstractAnsatz; kwargs...) = T01(u0, Ball(); kwargs...)

"""
    T02(u0::AbstractAnsatz, evaltype = Ball(); δ2, ϵ)
Returns a function such that T02(u0; δ2, ϵ)(x) computes the integral
T_{0,2} from the paper.

In general it returns an **enclosure** but in some cases it only
returns an **upper bound**.

The strategy for evaluation depends on evaltype and the precise
ansatz.
"""
T02(u0::AbstractAnsatz; kwargs...) = T02(u0, Ball(); kwargs...)
