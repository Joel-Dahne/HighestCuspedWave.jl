### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 28ec11ff-acb6-4be8-9cdd-45c42cfb839d
begin
    using Pkg, Revise
    Pkg.activate("../")
    using Arblib, ArbExtras, Folds, HighestCuspedWave, LaTeXStrings, Plots, PlutoUI

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ 4e3eaa14-583c-11ec-0fa5-e12df02e492f
md"""
# Fractional KdV Equation for $\alpha$ close to zero

This notebook contains the computer assisted part of the proof of existence of a $2\pi$-periodic highest cusped wave for the fractional KdV equations for $\alpha$ in the interval $(\epsilon, 0)$. The fractional KdV equation is given by

$f_t + f f_x = |D|^\alpha f_x$

For traveling waves the ansatz $f(x, t) = \varphi(x - ct)$ reduces the equation to

$-c\varphi' + \varphi \varphi' = |D|^\alpha \varphi'$

where $c$ is the wave speed.
"""

# ╔═╡ c7f82e62-33c7-4214-813d-1e63d7fb1abc
md"""
For details about the proof see the paper, we here give a very short overview of the reduction to a fixed point problem.

With the ansatz $\varphi(x) = c - u(x)$ and integrating we can further reduce it to the equation

$\frac{1}{2}u^2 = -\mathcal{H}^\alpha[u]$

where $\mathcal{H}^\alpha[u](x)$ is the operator

$\mathcal{H}^\alpha[u](x) = |D|^\alpha u(x) - |D|^\alpha u(0).$

If we now make the ansatz

$u(x) = u_\alpha(x) + w_\alpha(x)v(x)$

with $u_\alpha$ an approximate solution and $w_\alpha$ a fixed weight function we can, with some work, reduce the problem of proving the existence of a solution to the equation to proving existence of a fixed point of the operator

$G_\alpha[v] = (I - T_\alpha)^{-1}(-F_\alpha - N_\alpha v^2)$

Here

$N_\alpha(x) = \frac{w_\alpha(x)}{2u_\alpha(x)}$

$F_\alpha(x) = \frac{1}{w_\alpha(x)u_\alpha(x)}\left(\mathcal{H}^\alpha[u_\alpha](x) + \frac{1}{2}u_\alpha(x)^2\right)$

$T_\alpha[v](x) = -\frac{1}{w_\alpha(x) u_\alpha(x)}\mathcal{H}^\alpha[wv](x)$

Let $n_\alpha = \|N_\alpha\|_{L^\infty}$, $\delta_\alpha = \|F_\alpha\|_{L^\infty}$ and $D_\alpha = \|T_\alpha\|$. Proving the existence of a fixed point for $G_\alpha$ reduces to checking the inequality

$\delta_\alpha < \frac{(1 - D_\alpha)^2}{4n_\alpha}.$
"""

# ╔═╡ 6055dad1-da62-4048-a6fb-e8da26b5ef9a
md"""
To be able to handle all $\alpha \in (\epsilon, 0)$ we have to understand the behaviour of $n_\alpha$, $\delta_\alpha$ and $D_\alpha$ as $\alpha$ converges to $0$.

The value $n_0$ converges to a finite non-zero value as $\alpha \to 0$ and it is enough to compute an enclosure of $n_\alpha$ valid on the full interval $(\epsilon, 0)$.

As $\alpha \to 0$ we have that $\delta_\alpha$ converges to zero, which is good in terms of satisfying the inequality. However $D_\alpha \to 1$ from below, which is bad in terms of satisfying the inequality. It is therefore not possible to compute any uniform enclosure of $\delta_\alpha$ or $D_\alpha$ such that the inequality holds on the whole interval.

Since a uniform enclosure is not enough we will instead compute an enclosure depending on $\alpha$. More preciely we will find compute intervals $\Delta_F$ and $\Delta_B$ such that

$\delta_\alpha \in \Delta_\delta \alpha^2 \text{ and } D_\alpha \in 1 - \Delta_D \alpha$

for every $\alpha \in (\epsilon, 0)$. The inequality then reduces to checking

$\Delta_\delta \alpha^2 \leq \frac{\Delta_D \alpha^2}{4n_\alpha} \iff \Delta_\delta \leq \frac{\Delta_D^2}{4n_\alpha}.$

In what follows we compute enclosures of $n_\alpha$, $\Delta_\delta$ and $\Delta_D$.
"""

# ╔═╡ d1c5d171-c504-4c49-a6cb-7a707f3c7cc3
md"""
This is the $\epsilon$ we use for the computations.
"""

# ╔═╡ f27e0113-5e77-44af-958a-4b9e4a13a40b
ϵ = Arb(-1.2e-3)

# ╔═╡ e2f6fded-8452-4e17-a969-a040a362a26f
md"""
**NOTE:** in the code we don't use the subscript $\alpha$ but in most cases use `0`, for example $u_\alpha$ becomes `u0`.
"""

# ╔═╡ 27785eb7-54c7-47f7-8029-d37887bcdf1e
md"""
## Compute approximation
"""

# ╔═╡ 34009ff0-dadd-4464-9187-890723d94a3e
u0 = KdVZeroAnsatz(Arb((ϵ, 0)))

# ╔═╡ c5e2a657-dae0-4243-951c-a90df8bf1aaf
md"""
We will work on the interval $[0, \pi]$ where $u_\alpha$ looks like this.
"""

# ╔═╡ 4fee88af-73dd-4b47-86b9-501c0eb0f0e3
let xs = range(Arb(0), π, length = 100)
    ys = Folds.map(x -> u0(x)(u0.α), xs)
    plot(xs, ys, ribbon = radius.(Arb, ys), label = "", m = :circle, ms = 1)
end

# ╔═╡ 5a06011d-4fd2-47fd-adc0-e8828c0696eb
md"""
## Computing constants
"""

# ╔═╡ b85e1c39-99b2-4bc7-af04-854abf3de52b
md"""
### Bound $n_\alpha$
"""

# ╔═╡ d465c80c-9e34-4198-93b2-854243aea6d4
md"""
We start by computing an enclosure of $n_\alpha$ and plot it together with $N_\alpha(x)$ for $x \in [0, \pi]$.
"""

# ╔═╡ 13373308-42e3-44ac-bf99-2a5145e6a5df
n0_time = @elapsed n0 = n0_bound(u0, verbose = true)

# ╔═╡ fa901388-6f0a-4de5-8323-d9b4832062f6
n0_xs, n0_ys = let xs = range(Arb(0), π, length = 100)[2:end]
    N = x -> u0.w(x) / 2u0(x)(u0.α)
    ys = Folds.map(N, xs)
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
## Bound $\Delta_\delta$
"""

# ╔═╡ 073a34c8-a2e9-434c-b047-d1f1fffd590b
md"""
Next we compute an enclosure of $Δ_\delta$. For a fixed $x$ we can compute a Taylor model around $\alpha = 0$ of $F_\alpha$ corresonding to the defect for the given $x$.
"""

# ╔═╡ 054103ee-ac90-4310-afa3-55291712829d
F0(u0)(Arb(2.5)) # Note that the Taylor model is printed with the variable x and not α

# ╔═╡ 38e338b0-f7e7-45e6-b84d-45c6946736b3
md"""
Note that the polynomial for the Taylor model is zero, we can plot the remainder term as a function of $x$.
"""

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

# ╔═╡ 616c4ebf-b2ca-4c5d-a300-30f4c1918223
Δδ_xs_asym, Δδ_ys_asym =
    let xs = exp.(range(log(Arb("1e-10")), log(Arb("1e-1")), length = 200))
        f = let F0 = F0(u0, Asymptotic(), ϵ = one(Arb))
            x -> begin
                expansion = F0(x)
                @assert iszero(expansion.p[0]) && iszero(expansion.p[1])
                expansion.p[2]
            end
        end
        ys = Folds.map(f, xs)
        xs, ys
    end

# ╔═╡ 581ee14a-acac-4a3e-b854-f46a6f26dcc3
Δδ_time = @elapsed Δδ =
	delta0_bound(u0, maxevals = 10000, verbose = true, atol = Arb(0.0025)).p[2]

# ╔═╡ 03cf0016-f828-4049-90ae-5ce560ab644b
let pl = plot(legend = :bottomright)
    plot!(pl, Δδ_xs, Δδ_ys, ribbon = radius.(Arb, Δδ_ys), label = "", m = :circle, ms = 1)
    hline!([Δδ], ribbon = [radius(Arb, Δδ)], color = :green, label = "Δδ bound")
    hline!([-Δδ], ribbon = [radius(Arb, Δδ)], color = :green, label = "")
    pl
end

# ╔═╡ 17860367-ec0a-4dd4-8650-7bd99c5320c6
let pl = plot(legend = :bottomright)
    plot!(
        pl,
        Δδ_xs_asym,
        Δδ_ys_asym,
        ribbon = radius.(Arb, Δδ_ys_asym),
        label = "",
        m = :circle,
        ms = 1,
        xaxis = :log10,
    )
    hline!([Δδ], ribbon = [radius(Arb, Δδ)], color = :green, label = "Δδ bound")
    hline!([-Δδ], ribbon = [radius(Arb, Δδ)], color = :green, label = "")
end

# ╔═╡ ce0c9835-7448-4344-86d8-83e89181208b
md"""
## Bound $Δ_D$
"""

# ╔═╡ e12220bf-4434-424d-afcc-bf204fd81411
md"""
Finally we compute an enclosure of $Δ_D$. Similar to for $Δ_F$ we can compute a Taylor model around $\alpha = 0$ of $T_\alpha$ for a fixed $x$.
"""

# ╔═╡ 353a6fb2-054c-4398-9c6f-2eafda12e4fa
T0(u0)(Arb(2.5)) # Note that the expansion is printed with the variable x and not α

# ╔═╡ ce8d3cb4-ee97-4f9a-8466-e1709b478bd5
md"""
Note that the first coefficient in the expansion is $1$. We can plot the remainder term as a function of $x$.
"""

# ╔═╡ 9a568a65-91d0-4758-9cbb-15cc088f9aef
ΔD_xs, ΔD_ys = let xs = range(Arb(0), π, length = 200)[2:end]
    f = let T0_nonasym = T0(u0), T0_asym = T0(u0, Asymptotic())
        x -> begin
            if x < 1
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

# ╔═╡ f01b115e-0f4a-403a-807c-92d876c7e7c6
ΔD_xs_asym, ΔD_ys_asym = let xs = exp.(range(log(Arb("1e-10")), log(Arb(0.9)), length = 200))
    f = let T0 = T0(u0, Asymptotic(), ϵ = one(Arb))
        x -> begin
            expansion = T0(x)
            @assert isone(expansion.p[0])
            expansion.p[1]
        end
    end
    ys = Folds.map(f, xs)
    xs, ys
end

# ╔═╡ 39ce2082-2b36-4ce5-bcda-809c9c5214ec
ΔD_time = @elapsed ΔD = D0_bound(u0, verbose = true).p[1]

# ╔═╡ feb38bc2-3963-4143-8e1c-9077a967fba4
let pl = plot(legend = :bottomright)
    plot!(pl, ΔD_xs, ΔD_ys, ribbon = radius.(Arb, ΔD_ys), label = "", m = :circle, ms = 1)
    hline!([ΔD], ribbon = [radius(Arb, ΔD)], color = :green, label = "ΔD bound")
    pl
end

# ╔═╡ 89599cbf-2571-4475-a360-bc4e165d71d8
let pl = plot(legend = :topleft)
    plot!(
        pl,
        ΔD_xs_asym,
        ΔD_ys_asym,
        ribbon = radius.(Arb, ΔD_ys_asym),
        label = "",
        m = :circle,
        ms = 1,
        xaxis = :log10,
    )
    hline!([ΔD], ribbon = [radius(Arb, ΔD)], color = :green, label = "ΔD bound")
end

# ╔═╡ 0ba5ede7-9c13-46c4-b288-47ef301ad092
lbound(ΔD)

# ╔═╡ 11e831ac-4e18-4822-bc53-ee99394ab64f
md"""
The inequality we want to check is $\Delta_\delta \leq \frac{\Delta_D^2}{4n_\alpha}$
"""

# ╔═╡ e21e7ba1-7bd9-47c6-aa75-cebc4e413c4d
Δδ

# ╔═╡ b142ad21-e34e-42b5-b73f-8a32e71feb26
ΔD^2 / 4n0

# ╔═╡ 04ea9588-f613-4715-9f82-d23b786cc0ce
Δδ < ΔD^2 / 4n0

# ╔═╡ 93ae49fd-f041-4475-9dc0-823b0aa01178
md"""
# Prepare for publishing
We compute rounded values of the upper bounds that are given in the paper, as well as produce figures for the paper.
"""

# ╔═╡ 8289dd0e-c400-444d-b6df-aa525e41b3a8
md"""
**TODO:** Prepare values of bounds to be printed in the paper.
"""

# ╔═╡ 47d4faa9-ecdd-4a04-b485-7e75333632ab
n0_rounded = round(Arblib.get_d(ubound(n0), RoundUp), RoundUp, sigdigits = 5)

# ╔═╡ 9d432aa8-ecb4-44d7-9180-595c9167cf09
Δδ_rounded = round(Arblib.get_d(ubound(Δδ), RoundUp), RoundUp, sigdigits = 5)

# ╔═╡ 3ebb3003-a281-43ba-9122-e75e20e84041
ΔD_rounded = round(Arblib.get_d(lbound(ΔD), RoundDown), RoundDown, sigdigits = 5)

# ╔═╡ 14a37967-e894-4571-8aa6-8733c544dee2
md"""
The plots in the paper were produced with `pgfplotsx`, this is not a dependency so is not enabled by default. The plots look a bit weird with the default backend.
"""

# ╔═╡ c25e8495-28fe-4454-a9c4-10b5ed43b917
#pgfplotsx()

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
    savefig(pl, "../figures/publication/KdVZero-N.pdf")
    pl
end

# ╔═╡ 750101e6-7433-406a-8776-5914f81d0634
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"\Delta_D(x)",
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
    savefig(pl, "../figures/publication/KdVZero-Delta-D.pdf")
    pl
end

# ╔═╡ 939b4965-b0bb-45ef-966b-aa8ecd3e3f3a
Δδ_rounded_goal = ΔD_rounded^2 / 4n0_rounded

# ╔═╡ 52d4822c-9c20-484f-9eda-63dc4eaeda4a
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"\Delta_\delta(x)",
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
    savefig(pl, "../figures/publication/KdVZero-Delta-delta.pdf")
    pl
end

# ╔═╡ Cell order:
# ╠═28ec11ff-acb6-4be8-9cdd-45c42cfb839d
# ╟─4e3eaa14-583c-11ec-0fa5-e12df02e492f
# ╟─c7f82e62-33c7-4214-813d-1e63d7fb1abc
# ╟─6055dad1-da62-4048-a6fb-e8da26b5ef9a
# ╟─d1c5d171-c504-4c49-a6cb-7a707f3c7cc3
# ╠═f27e0113-5e77-44af-958a-4b9e4a13a40b
# ╟─e2f6fded-8452-4e17-a969-a040a362a26f
# ╟─27785eb7-54c7-47f7-8029-d37887bcdf1e
# ╠═34009ff0-dadd-4464-9187-890723d94a3e
# ╟─c5e2a657-dae0-4243-951c-a90df8bf1aaf
# ╠═4fee88af-73dd-4b47-86b9-501c0eb0f0e3
# ╟─5a06011d-4fd2-47fd-adc0-e8828c0696eb
# ╟─b85e1c39-99b2-4bc7-af04-854abf3de52b
# ╟─d465c80c-9e34-4198-93b2-854243aea6d4
# ╠═13373308-42e3-44ac-bf99-2a5145e6a5df
# ╟─fa901388-6f0a-4de5-8323-d9b4832062f6
# ╟─af90b425-59cb-4229-873a-bb7448713dbe
# ╟─d977785c-751a-4e2a-8634-9ecd1039e3f8
# ╟─073a34c8-a2e9-434c-b047-d1f1fffd590b
# ╠═054103ee-ac90-4310-afa3-55291712829d
# ╟─38e338b0-f7e7-45e6-b84d-45c6946736b3
# ╟─51603714-bb7e-4691-b80d-7b18cf159b94
# ╟─616c4ebf-b2ca-4c5d-a300-30f4c1918223
# ╠═581ee14a-acac-4a3e-b854-f46a6f26dcc3
# ╟─03cf0016-f828-4049-90ae-5ce560ab644b
# ╟─17860367-ec0a-4dd4-8650-7bd99c5320c6
# ╟─ce0c9835-7448-4344-86d8-83e89181208b
# ╟─e12220bf-4434-424d-afcc-bf204fd81411
# ╠═353a6fb2-054c-4398-9c6f-2eafda12e4fa
# ╟─ce8d3cb4-ee97-4f9a-8466-e1709b478bd5
# ╟─9a568a65-91d0-4758-9cbb-15cc088f9aef
# ╟─f01b115e-0f4a-403a-807c-92d876c7e7c6
# ╠═39ce2082-2b36-4ce5-bcda-809c9c5214ec
# ╟─feb38bc2-3963-4143-8e1c-9077a967fba4
# ╟─89599cbf-2571-4475-a360-bc4e165d71d8
# ╠═0ba5ede7-9c13-46c4-b288-47ef301ad092
# ╟─11e831ac-4e18-4822-bc53-ee99394ab64f
# ╠═e21e7ba1-7bd9-47c6-aa75-cebc4e413c4d
# ╠═b142ad21-e34e-42b5-b73f-8a32e71feb26
# ╠═04ea9588-f613-4715-9f82-d23b786cc0ce
# ╟─93ae49fd-f041-4475-9dc0-823b0aa01178
# ╟─8289dd0e-c400-444d-b6df-aa525e41b3a8
# ╠═47d4faa9-ecdd-4a04-b485-7e75333632ab
# ╠═9d432aa8-ecb4-44d7-9180-595c9167cf09
# ╠═3ebb3003-a281-43ba-9122-e75e20e84041
# ╟─14a37967-e894-4571-8aa6-8733c544dee2
# ╠═c25e8495-28fe-4454-a9c4-10b5ed43b917
# ╟─534cf54b-9044-43ca-8d47-d8917e2856cb
# ╟─750101e6-7433-406a-8776-5914f81d0634
# ╠═939b4965-b0bb-45ef-966b-aa8ecd3e3f3a
# ╟─52d4822c-9c20-484f-9eda-63dc4eaeda4a
