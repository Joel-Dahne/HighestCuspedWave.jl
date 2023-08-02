### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(
                Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                "AbstractPlutoDingetjes",
            )].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 52d3f93a-701a-11ec-1841-c365d4d26244
begin
    using Pkg
    Pkg.activate("../", io = devnull)
    using Arblib, Folds, HighestCuspedWave, LaTeXStrings, Plots, PlutoUI

    pgfplotsx()

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ 8af96418-5c46-44bc-8cd6-37c052d1cc0c
md"""
# Figures
This notebook is responsible for generating the figures in the paper not directly related to the computer assisted proof.
"""

# ╔═╡ b12eb977-09e6-4437-948c-c94ee2202a83
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save CheckBox(default = false))
"""

# ╔═╡ 98bed951-46c9-43cf-bcb0-cc88276ae22e
md"""
## Behaviour of $a_{\alpha,0}$ and $p_\alpha$
"""

# ╔═╡ dc7a7f9b-d050-45c8-a1eb-193198e416ee
let
    # Compute data
    αs = collect(range(-1, 0, length = 100)[2:end-1])
    a0s = Float64.(HighestCuspedWave.finda0.(Arb.(αs)))
    # The limit at α = 0 is 0
    push!(αs, 0)
    push!(a0s, 0)

    # Sanity checks
    @assert all(isfinite, a0s)
    @assert all(<=(0), a0s)

    # Make plot
    pl = plot(
        αs,
        a0s,
        xlabel = L"\alpha",
        ylabel = L"a_{\alpha,0}",
        ylims = (-5, 0),
        guidefontsize = 18,
        tickfontsize = 18,
        legend = :none,
    )

    save && savefig(pl, "../figures/a_0.pdf")

    pl
end

# ╔═╡ ccd152d4-074d-462e-bd69-50b98ed9256c
let
    # Compute data
    αs = collect(range(-1, 0, length = 100)[2:end-1])
    p0s = Float64.(HighestCuspedWave.findp0.(Arb.(αs)))
    # The limit at α = -1 is 0
    pushfirst!(αs, -1)
    pushfirst!(p0s, 0)
    # Compute the limit at α = 0
    push!(αs, 0)
    push!(p0s, HighestCuspedWave.expansion_p0(KdVZeroAnsatz, Arb(0), Arb(0)).p[0])

    # Sanity checks
    @assert all(isfinite, p0s)
    @assert all(>=(0), p0s)

    # Make plot
    pl = plot(
        αs,
        p0s,
        xlabel = L"\alpha",
        ylabel = L"p_\alpha",
        guidefontsize = 18,
        tickfontsize = 18,
        legend = :none,
    )

    save && savefig(pl, "../figures/p_alpha.pdf")

    pl
end

# ╔═╡ ac7f6b0a-7b7e-4720-80d6-11e76c7d5eb4
md"""
## Behaviour of $\mathcal{T}_{\alpha}(0)$
"""

# ╔═╡ f382326b-c95c-4bd4-aba6-e2790c0f4d16
let
    # Compute data
    αs = collect(range(-1, 0, length = 100)[2:end-1])
    T0s = Folds.map(αs, basesize = 1) do α
        u0 = FractionalKdVAnsatz(Arb(α), p = Arb(1), use_bhkdv = false)
        Float64(T0(u0, Asymptotic(), M = 5)(Arb(0)))
    end

    # Sanity checks
    @assert all(isfinite, T0s)
    @assert T0s[end] < 1 < T0s[1]
    @assert T0s[length(T0s)÷2] < 1 # Midpoint less than 1

    # Make plot
    pl = plot(
        αs,
        T0s,
        xlabel = L"\alpha",
        ylabel = L"\mathcal{T}_{\alpha}(0)",
        guidefontsize = 18,
        tickfontsize = 18,
        legend = :none,
    )

    save && savefig(pl, "../figures/T_alpha_0.pdf")

    pl
end

# ╔═╡ Cell order:
# ╟─8af96418-5c46-44bc-8cd6-37c052d1cc0c
# ╠═52d3f93a-701a-11ec-1841-c365d4d26244
# ╟─b12eb977-09e6-4437-948c-c94ee2202a83
# ╟─98bed951-46c9-43cf-bcb0-cc88276ae22e
# ╟─dc7a7f9b-d050-45c8-a1eb-193198e416ee
# ╟─ccd152d4-074d-462e-bd69-50b98ed9256c
# ╟─ac7f6b0a-7b7e-4720-80d6-11e76c7d5eb4
# ╟─f382326b-c95c-4bd4-aba6-e2790c0f4d16
