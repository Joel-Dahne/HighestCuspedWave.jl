export alpha0, delta0, CB, CB_bounded_by

"""
    alpha0(u0::AbstractAnsatz{Arb})

Compute upper bound of `α₀` from the paper.
"""
function alpha0 end

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
