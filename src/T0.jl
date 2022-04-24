export T0

"""
    T0(u0::AbstractAnsatz, evaltype = Ball())

Return a function such that `T0(u0)(x)` computes the function ``𝒯(x)``
from the paper. This is the function whose supremum on ``[0, π]``
gives `D₀`.

In general it only returns an **upper bound**, but in some cases it
returns an **enclosure**.
"""
T0(u0::AbstractAnsatz; kwargs...) = T0(u0, Ball(); kwargs...)

"""
    T01(u0::AbstractAnsatz, evaltype = Ball())

Return a function such that `T01(u0)(x)` computes the part of
`T0(u0)(x)` that corresponds to the integration from `0` to `x`.

In general it only returns an **upper bound**, but in some cases it
returns an **enclosure**.
"""
T01(u0::AbstractAnsatz; kwargs...) = T01(u0, Ball(); kwargs...)

"""
    T02(u0::AbstractAnsatz, evaltype = Ball())

Return a function such that `T02(u0)(x)` computes the part of
`T0(u0)(x)` that corresponds to the integration from `x` to `π`.

In general it only returns an **upper bound**, but in some cases it
returns an **enclosure**.
"""
T02(u0::AbstractAnsatz; kwargs...) = T02(u0, Ball(); kwargs...)
