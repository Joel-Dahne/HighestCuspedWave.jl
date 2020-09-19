function prove(u0::FractionalKdVAnsatz{arb};
               verbose = false,
               )

    α₀ = alpha0(u0)
    δ₀ = delta0(u0)

    if verbose
        @show α₀ δ₀
    end

    C_B = CB_estimate(u0)

    β = 1/(1 - C_B)
    C = 1/(4α₀*β^2)
    D = 1 - 2sqrt(α₀*δ₀)

    if verbose
        @show C_B C
        println("Must have C_B ≤ 1 - 2√(α₀δ₀) = $D")
    end

    if !(δ₀ < C)
        return false
    end

    res = CB_bounded_by(u0, D, show_trace = true)

    return res
end
