# Contains error terms used in the asymptotic expansions

"""
    E(u0::FractionalKdVAnsatz{Arb}, M::Integer)

Compute an enclosure of ``E_{u_0}(x)`` from the paper, skipping the
`x^2M` factor.
"""
function E(u0::FractionalKdVAnsatz{Arb}, M::Integer)
    π = convert(Arb, pi)

    return x::Arb -> begin
        # Compute error bounds for the Clausen functions
        E_bound1 = zero(u0.α)
        for j = 0:u0.N0
            E_bound1 += (2π)^(j * u0.p0) * abs(zeta(2M + u0.α - j * u0.p0) * u0.a[j])
        end
        E_bound1 *= 2(2π)^(2 - u0.α - 2M) / (4π^2 - x^2)

        # Compute error bounds for the Fourier terms
        E_bound2 = zero(u0.α)
        for n = 1:u0.N1
            E_bound2 += Arb(n)^(2M) * abs(u0.b[n])
        end
        E_bound2 /= factorial(2M)

        return Arblib.add_error!(zero(u0.α), E_bound1 + E_bound2)
    end
end

"""
    EH(u0::FractionalKdVAnsatz{Arb}, M::Integer)

Compute an enclosure of ``E_{H^{-\alpha}u_0}(x)`` from the paper,
skipping the `x^2M` factor.
"""
function EH(u0::FractionalKdVAnsatz{Arb}, M::Integer)
    π = convert(Arb, pi)

    return x -> begin
        # Compute error bounds for the Clausians
        E_bound1 = zero(u0.α)
        for j = 0:u0.N0
            E_bound1 += (2π)^(j * u0.p0) * abs(zeta(2M + 2u0.α - j * u0.p0) * u0.a[j])
        end
        E_bound1 *= 2 * (2π)^(2 - 2u0.α - 2M) / (4π^2 - x^2)

        # Compute error bounds for the Fourier terms
        E_bound2 = zero(u0.α)
        for n = 1:u0.N1
            E_bound2 += n^(2M + u0.α) * abs(u0.b[n])
        end
        E_bound2 /= factorial(2M)

        return Arblib.add_error!(zero(u0.α), E_bound1 + E_bound2)
    end
end
