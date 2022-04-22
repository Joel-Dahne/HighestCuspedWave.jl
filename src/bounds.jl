export n0_bound, delta0, CB, CB_bounded_by

"""
    n0_bound(u0::AbstractAnsatz{Arb})

Compute upper bound of `n₀` from the paper.
"""
function n0_bound end

"""
    delta0(u0::AbstractAnsatz{Arb})

Compute upper bound of `δ₀` from the paper.
"""
function delta0 end

"""
    CB(u0::AbstractAnsatz{Arb})

Compute upper bound of ``C_B`` from the paper.
"""
function CB end

"""
    CB_bounded_by(u0::AbstractAnsatz{Arb}, C::Arb)

Return true if `CB(u0)` is bounded by `C`. Return false if it is not
bounded by `C` or if proving the bound fails.
"""
function CB_bounded_by end
