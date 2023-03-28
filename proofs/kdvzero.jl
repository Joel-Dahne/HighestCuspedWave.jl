### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ 28ec11ff-acb6-4be8-9cdd-45c42cfb839d
begin
    using Pkg
    Pkg.activate("../", io = devnull)
    using Arblib, ArbExtras, Folds, HighestCuspedWave, LaTeXStrings, Plots, PlutoUI

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ 4e3eaa14-583c-11ec-0fa5-e12df02e492f
md"""
# Fractional KdV Equation for $\alpha$ close to zero

This notebook contains the computer assisted part of the proof of existence of a $2\pi$-periodic highest cusped wave for the fractional KdV equations for $\alpha$ in the interval $(\epsilon, 0)$. More precisely it contains the proof of Lemma 12.9, 12.10 and 12.11.
"""

# ╔═╡ d1c5d171-c504-4c49-a6cb-7a707f3c7cc3
md"""
This is the $\epsilon$ we use for the computations.
"""

# ╔═╡ f27e0113-5e77-44af-958a-4b9e4a13a40b
ϵ = Arb(-1.2e-3)

# ╔═╡ 27785eb7-54c7-47f7-8029-d37887bcdf1e
md"""
## Compute the approximation $u_\alpha$
"""

# ╔═╡ 34009ff0-dadd-4464-9187-890723d94a3e
u0 = KdVZeroAnsatz(Arb((ϵ, 0)))

# ╔═╡ c5e2a657-dae0-4243-951c-a90df8bf1aaf
md"""
We work on the interval $[0, \pi]$ where $u_\alpha$ looks like this.
"""

# ╔═╡ 4fee88af-73dd-4b47-86b9-501c0eb0f0e3
let xs = range(Arb(0), π, length = 100)
    ys = Folds.map(x -> u0(x)(u0.α), xs)
    plot(
        xs,
        ys,
        ribbon = radius.(Arb, ys),
        xlabel = L"x",
        ylabel = L"u_\alpha(x)",
        label = "",
        m = :circle,
        ms = 1,
    )
end

# ╔═╡ 5a06011d-4fd2-47fd-adc0-e8828c0696eb
md"""
## Compute constants
We are now ready to compute upper bounds of $n_\alpha$, $\Delta_\delta$ and $\Delta_D$. The precise upper bounds given in the paper are produced in the next part.
"""

# ╔═╡ b85e1c39-99b2-4bc7-af04-854abf3de52b
md"""
### Bound $n_\alpha$
This corresponds to Lemma 12.9. We compute an enclosure of $n_\alpha$ and plot it together with $N_\alpha(x)$ for $x \in [0, \pi]$.
"""

# ╔═╡ 13373308-42e3-44ac-bf99-2a5145e6a5df
@time n0_time = @elapsed n0 = n0_bound(u0, verbose = true)

# ╔═╡ fa901388-6f0a-4de5-8323-d9b4832062f6
n0_xs, n0_ys = let xs = range(Arb(0), π, length = 100)
    N = x -> u0.w(x) / 2u0(x)(u0.α)
    ys = Folds.map(N, xs)
    ys[1] = 0 # It is zero at x = 0
    xs, ys
end

# ╔═╡ af90b425-59cb-4229-873a-bb7448713dbe
let pl = plot(legend = :bottomright)
    plot!(pl, n0_xs, n0_ys, ribbon = radius.(Arb, n0_ys), label = "", m = :circle, ms = 2)
    hline!(pl, [n0], ribbon = [radius(Arb, n0)], color = :green, label = "n0 bound")
    pl
end

# ╔═╡ d977785c-751a-4e2a-8634-9ecd1039e3f8
md"""
### Bound $\Delta_\delta$
This corresonds to Lemma 12.10. For a fixed $x$ we can compute a Taylor model around $\alpha = 0$ of $F_\alpha(x)$ corresonding to the defect for the given $x$.
"""

# ╔═╡ 054103ee-ac90-4310-afa3-55291712829d
F0(u0)(Arb(2.5)) # Note that the Taylor model is printed with the variable x and not α

# ╔═╡ 38e338b0-f7e7-45e6-b84d-45c6946736b3
md"""
Note that the polynomial for the Taylor model is zero. We bound the remainder term and also plot it as a function of $x$.
"""

# ╔═╡ 581ee14a-acac-4a3e-b854-f46a6f26dcc3
@time Δδ_time = @elapsed Δδ = delta0_bound(u0, verbose = true).p[2]

# ╔═╡ 51603714-bb7e-4691-b80d-7b18cf159b94
Δδ_xs, Δδ_ys = let xs = range(Arb(0), π, length = 200)
    f = let F0_nonasym = F0(u0), F0_asym = F0(u0, Asymptotic(), ϵ = Arb(3.2), M = 10)
        # We only use the asymptotic evaluation
        x -> begin
            expansion = F0_asym(x)
            @assert iszero(expansion.p[0]) && iszero(expansion.p[1])
            expansion.p[2]
        end
    end
    ys = Folds.map(f, xs)
    xs, ys
end

# ╔═╡ 03cf0016-f828-4049-90ae-5ce560ab644b
let pl = plot(legend = :bottomright)
    plot!(pl, Δδ_xs, Δδ_ys, ribbon = radius.(Arb, Δδ_ys), label = "", m = :circle, ms = 1)
    hline!([Δδ], ribbon = [radius(Arb, Δδ)], color = :green, label = "Δδ bound")
    hline!([-Δδ], ribbon = [radius(Arb, Δδ)], color = :green, label = "")
    pl
end

# ╔═╡ ce0c9835-7448-4344-86d8-83e89181208b
md"""
### Bound $Δ_D$
This corresonds to Lemma 12.11. Similar to for $Δ_F$ we can compute a Taylor model around $\alpha = 0$ of $T_\alpha$ for a fixed $x$.
"""

# ╔═╡ 353a6fb2-054c-4398-9c6f-2eafda12e4fa
T0(u0)(Arb(2.5)) # Note that the expansion is printed with the variable x and not α

# ╔═╡ ce8d3cb4-ee97-4f9a-8466-e1709b478bd5
md"""
Note that the first coefficient in the expansion is $1$. We can bound the remainder term and also plot it as a function of $x$.
"""

# ╔═╡ 39ce2082-2b36-4ce5-bcda-809c9c5214ec
@time ΔD_time = @elapsed ΔD = D0_bound(u0, verbose = true).p[1]

# ╔═╡ 9a568a65-91d0-4758-9cbb-15cc088f9aef
ΔD_xs, ΔD_ys = let xs = range(Arb(0), π, length = 200)[2:end]
    f = let T0_nonasym = T0(u0), T0_asym = T0(u0, Asymptotic())
        x -> begin
            if x < 0.1
                expansion = T0_asym(x)
            else
                expansion = T0_nonasym(x)
            end
            @assert isone(expansion.p[0])
            expansion.p[1]
        end
    end
    ys = Folds.map(f, xs)
    xs, ys
end

# ╔═╡ feb38bc2-3963-4143-8e1c-9077a967fba4
let pl = plot(legend = :bottomright)
    plot!(pl, ΔD_xs, ΔD_ys, ribbon = radius.(Arb, ΔD_ys), label = "", m = :circle, ms = 1)
    hline!([ΔD], ribbon = [radius(Arb, ΔD)], color = :green, label = "ΔD bound")
    pl
end

# ╔═╡ 65291862-b961-4410-a020-23f7abe358a2
md"""
## Check inequality
In the end the inequality we want to be satisfied is

$$\delta_\alpha < \frac{(1 - D_\alpha)^2}{4n_\alpha}$$

This is equivalent to

$$\Delta_\delta < \frac{\Delta_D^2}{4n_\alpha}$$

which we can check
"""

# ╔═╡ 4abcc0c4-0c90-4e6f-9fdb-9416eb166a8a
Δδ < ΔD^2 / 4n0

# ╔═╡ 93ae49fd-f041-4475-9dc0-823b0aa01178
md"""
## Prepare for publishing
We compute rounded values of the upper bounds that are given in the paper and check that they satisfy the required inequality. We also produce figures for the paper.
"""

# ╔═╡ dc035349-06ff-42c9-8956-932144c703c0
proved, n0_rounded, Δδ_rounded, ΔD_rounded = HighestCuspedWave.round_for_publishing_kdvzero(
    n0,
    Δδ,
    ΔD,
    sigdigits_n₀ = 3,
    sigdigits_Δδ = 2,
    sigdigits_ΔD = 3,
)

# ╔═╡ c9a21c38-e27c-47b7-b95b-8ab737408f6f
if proved
    @info "Inequality holds! 🥳" n0_rounded Δδ_rounded ΔD_rounded
else
    @error "Inequality doesn't hold!" n0_rounded Δδ_rounded ΔD_rounded
end

# ╔═╡ ee0c1ff9-baed-416f-b0db-7740ebdd288e
md"""
When giving the enclosures in the paper we want to make sure that the printed version is less than the given upper bound. Here we check that this is the case. Note that this is not important for the result to hold, just for the printed values to look nice.
"""

# ╔═╡ 80c61213-b801-4aef-9c0e-7a9f74f880f3
n0_printed = string(n0)

# ╔═╡ a0376989-d0e2-4be5-b9f8-c694aafe83d4
Arb(n0_printed) < n0_rounded

# ╔═╡ cdd7a87e-70d6-4c93-b65c-b6e35afbd721
Δδ_printed = string(Δδ)

# ╔═╡ a03d3f5f-e9e6-45bf-9b7a-c77a38e7d05d
Arb(Δδ_printed) < Δδ_rounded

# ╔═╡ 83078a21-fd16-4b4d-9cba-e92c66f3dbeb
md"""
For $\Delta_D$ we need to print a few more digits to get the enclosures correct.
"""

# ╔═╡ ca217129-82ba-4062-a387-67f4f78b87c1
ΔD_printed = Arblib.string_nice(ΔD, 2, UInt(1))

# ╔═╡ 4f1d46e1-2301-43a3-a2ad-897aaf16b2b1
Arb(ΔD_printed) > ΔD_rounded

# ╔═╡ 14a37967-e894-4571-8aa6-8733c544dee2
md"""
The plots in the paper were produced with `pgfplotsx`. However the plots look nicer in the notebook with `gr` so this is not the default.
"""

# ╔═╡ c25e8495-28fe-4454-a9c4-10b5ed43b917
#pgfplotsx()

# ╔═╡ 3994a705-3b25-45be-b67a-92f835ace2b8
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save CheckBox(default = false))
"""

# ╔═╡ 534cf54b-9044-43ca-8d47-d8917e2856cb
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"N_\alpha(x)",
        guidefontsize = 18,
        tickfontsize = 18,
    )
    plot!(
        pl,
        Float64.(n0_xs),
        Float64.(n0_ys),
        ribbon = Float64.(radius.(Arb, n0_ys)),
        m = :circle,
        ms = 1,
    )
    hline!(pl, Float64[n0_rounded], color = :green, linestyle = :dash)
    save && savefig(pl, "../figures/KdVZero-N.pdf")
    pl
end

# ╔═╡ 750101e6-7433-406a-8776-5914f81d0634
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"\Delta_{M_T(x)}",
        guidefontsize = 18,
        tickfontsize = 18,
    )
    plot!(
        pl,
        Float64.(ΔD_xs),
        Float64.(ΔD_ys),
        ribbon = Float64.(radius.(Arb, ΔD_ys)),
        m = :circle,
        ms = 1,
    )
    hline!(pl, Float64[ΔD_rounded], color = :green, linestyle = :dash)
    save && savefig(pl, "../figures/KdVZero-Delta-D.pdf")
    pl
end

# ╔═╡ 939b4965-b0bb-45ef-966b-aa8ecd3e3f3a
Δδ_rounded_goal = ΔD_rounded^2 / 4n0_rounded

# ╔═╡ 52d4822c-9c20-484f-9eda-63dc4eaeda4a
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"\Delta_{M_F(x)}",
        guidefontsize = 24,
        tickfontsize = 24,
    )
    plot!(
        pl,
        Float64.(Δδ_xs),
        Float64.(Δδ_ys),
        ribbon = Float64.(radius.(Arb, Δδ_ys)),
        m = :circle,
        ms = 1,
    )
    hline!(Float64[-Δδ_rounded, Δδ_rounded], color = :green, linestyle = :dash)
    hline!(Float64[-Δδ_rounded_goal, Δδ_rounded_goal], color = :red, linestyle = :dot)
    save && savefig(pl, "../figures/KdVZero-Delta-delta.pdf")
    pl
end

# ╔═╡ Cell order:
# ╠═28ec11ff-acb6-4be8-9cdd-45c42cfb839d
# ╟─4e3eaa14-583c-11ec-0fa5-e12df02e492f
# ╟─d1c5d171-c504-4c49-a6cb-7a707f3c7cc3
# ╠═f27e0113-5e77-44af-958a-4b9e4a13a40b
# ╟─27785eb7-54c7-47f7-8029-d37887bcdf1e
# ╠═34009ff0-dadd-4464-9187-890723d94a3e
# ╟─c5e2a657-dae0-4243-951c-a90df8bf1aaf
# ╟─4fee88af-73dd-4b47-86b9-501c0eb0f0e3
# ╟─5a06011d-4fd2-47fd-adc0-e8828c0696eb
# ╟─b85e1c39-99b2-4bc7-af04-854abf3de52b
# ╠═13373308-42e3-44ac-bf99-2a5145e6a5df
# ╟─fa901388-6f0a-4de5-8323-d9b4832062f6
# ╟─af90b425-59cb-4229-873a-bb7448713dbe
# ╟─d977785c-751a-4e2a-8634-9ecd1039e3f8
# ╠═054103ee-ac90-4310-afa3-55291712829d
# ╟─38e338b0-f7e7-45e6-b84d-45c6946736b3
# ╠═581ee14a-acac-4a3e-b854-f46a6f26dcc3
# ╟─51603714-bb7e-4691-b80d-7b18cf159b94
# ╟─03cf0016-f828-4049-90ae-5ce560ab644b
# ╟─ce0c9835-7448-4344-86d8-83e89181208b
# ╠═353a6fb2-054c-4398-9c6f-2eafda12e4fa
# ╟─ce8d3cb4-ee97-4f9a-8466-e1709b478bd5
# ╠═39ce2082-2b36-4ce5-bcda-809c9c5214ec
# ╟─9a568a65-91d0-4758-9cbb-15cc088f9aef
# ╟─feb38bc2-3963-4143-8e1c-9077a967fba4
# ╟─65291862-b961-4410-a020-23f7abe358a2
# ╠═4abcc0c4-0c90-4e6f-9fdb-9416eb166a8a
# ╟─93ae49fd-f041-4475-9dc0-823b0aa01178
# ╠═dc035349-06ff-42c9-8956-932144c703c0
# ╟─c9a21c38-e27c-47b7-b95b-8ab737408f6f
# ╟─ee0c1ff9-baed-416f-b0db-7740ebdd288e
# ╠═80c61213-b801-4aef-9c0e-7a9f74f880f3
# ╠═a0376989-d0e2-4be5-b9f8-c694aafe83d4
# ╠═cdd7a87e-70d6-4c93-b65c-b6e35afbd721
# ╠═a03d3f5f-e9e6-45bf-9b7a-c77a38e7d05d
# ╟─83078a21-fd16-4b4d-9cba-e92c66f3dbeb
# ╠═ca217129-82ba-4062-a387-67f4f78b87c1
# ╠═4f1d46e1-2301-43a3-a2ad-897aaf16b2b1
# ╟─14a37967-e894-4571-8aa6-8733c544dee2
# ╠═c25e8495-28fe-4454-a9c4-10b5ed43b917
# ╟─3994a705-3b25-45be-b67a-92f835ace2b8
# ╟─534cf54b-9044-43ca-8d47-d8917e2856cb
# ╟─750101e6-7433-406a-8776-5914f81d0634
# ╠═939b4965-b0bb-45ef-966b-aa8ecd3e3f3a
# ╟─52d4822c-9c20-484f-9eda-63dc4eaeda4a
