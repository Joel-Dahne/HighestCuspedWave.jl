export n0_bound, delta0_bound, D0_bound, D0_bounded_by

"""
    n0_bound(u0::AbstractAnsatz{Arb})

Compute upper bound of `n₀` from the paper.
"""
function n0_bound end

"""
    delta0_bound(u0::AbstractAnsatz{Arb})

Compute upper bound of `δ₀` from the paper.
"""
function delta0_bound end

"""
    D0_bound(u0::AbstractAnsatz{Arb})

Compute upper bound of `D₀` from the paper.
"""
function D0_bound end

"""
    D0_bounded_by(u0::AbstractAnsatz{Arb}, C::Arb)

Return true if `D0(u0)` is bounded by `C`. Return false if it is not
bounded by `C` or if proving the bound fails.
"""
function D0_bounded_by end
