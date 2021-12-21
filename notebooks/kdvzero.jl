### A Pluto.jl notebook ###
# v0.17.3

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
    using Pkg, Revise
    Pkg.activate("../")
    using Arblib, ArbExtras, HighestCuspedWave, Plots, PlutoUI

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ 4e3eaa14-583c-11ec-0fa5-e12df02e492f
md"""
# Fractional KdV Equations for $\alpha$ close to zero

This notebook contains illustrates the computer assisted part of the proof of existence of a highest cusped traveling wave for the fractional KdV equations for $\alpha$ in the interval $(\epsilon, 0)$. The equation we are interested in is given by

$u_t + \left(\frac{u^2}{2}\right)_x = H^{-\alpha}[u]$

where $H^{-\alpha}[u]$ is the fractional Hilbert transform of $u$. For the details of the proof see the paper, this notebook focuses on the computer assisted parts.
"""

# ╔═╡ fa84cc84-d9af-4f76-b01b-a088572626bf
md"""
What we actually prove is that there is $2\pi$ periodic, even solution which at $x = 0$ behaves like 

$u(x) = |x|^{-\alpha} + \mathcal{O}(|x|).$

The general idea is to construct an approximate solution, $u_0$, of the equation and then prove that there is a true solution

$u(x) = u_0(x) + |x|u_e(x)$

where $u_e \in L^\infty(\mathbb{T})$ is an error term which is bounded on the interval. To do this we use a fixed point argument and for it to go through we need to control three different values that depend on the choice of $u_0$.
- $\alpha_0 = \sup_{x \in \mathbb{T}} \left|\frac{w(x)}{2u_0(x)}\right|;$
- The defect, $\delta_0$, how far our approximate solution is from satisfying the differential equation;
- A number $C_B$ which depends on the norm of an operator that occurs in the fixed point equation.
They need to satisfy the **inequality**

$\delta \leq \frac{1}{4\alpha_0 \beta^2}$

where $\beta = \frac{1}{1 - C_B}$ and we require that $C_B < 1$.
"""

# ╔═╡ 6055dad1-da62-4048-a6fb-e8da26b5ef9a
md"""
To be able to handle the full interval $(\epsilon, 0)$ we have to understand the behaviour of $\alpha_0$, $\delta_0$ and $\beta$ as $\alpha$ converges to $0$.

The value $\alpha_0$ converges to a finite non-zero value as $\alpha \to 0$ and it is enough to compute an enclosure of $\alpha_0$ valid on the full interval $(\epsilon, 0)$.

As $\alpha \to 0$ we have that $\delta_0$ converges to zero, which is good in terms of satisfying the inequality. However $C_B \to 1$ from below and hence $\beta \to \infty$ which is bad in terms of satisfying the inequality. It is therefore not possible to compute any uniform enclosure of $\delta_0$ or $C_B$ such that the inequality holds on the whole interval.

Since a uniform enclosure is not enough we will instead compute an enclosure which is paratemtric in $\alpha$. More preciely we will find intervals $A$ and $B$ such that 

$\delta_0 \in A \alpha^2 \text{ and } C_B \in 1 - B \alpha$

for every $\alpha \in (\epsilon, 0)$. The inequality then reduces to checking

$A\alpha^2 \leq \frac{1}{4\alpha_0 \frac{1}{B^2 \alpha^2}} = \frac{B^2}{4\alpha_0}\alpha^2 \iff A \leq \frac{B^2}{4\alpha_0}.$

In what follows we compute enclosures of $\alpha_0$, $A$ and $B$.
"""

# ╔═╡ f27e0113-5e77-44af-958a-4b9e4a13a40b
ϵ = Arb(-1e-8)

# ╔═╡ 34009ff0-dadd-4464-9187-890723d94a3e
u0 = KdVZeroAnsatz(Arb((ϵ, 0)))

# ╔═╡ 5a06011d-4fd2-47fd-adc0-e8828c0696eb
md"""
## Computing constants
"""

# ╔═╡ f101b712-69a4-4f38-add8-f3a43a43a5d6
md"""
The code can either compute rigorous error bounds for the required constants or use estimates. The estimates are given by simply evaluating the corresponding functions on a number of points and taking the maximum. For the defect `δ0` we also make sure to use points asymptically close to `0` since that's where the largest defect is found. Check the constants to use rigorous error bounds for
-  `δ₀` $(@bind use_rigorous_bounds_δ0 CheckBox(default = false))
-  `C_B` $(@bind use_rigorous_bounds_C_B CheckBox(default = false))
Notice that the rigorous error bounds take **significantly** longer time to compute with.
"""

# ╔═╡ c6207967-17ea-41ff-8049-5e940d62374e
md"""
The code uses
- **$(ifelse(use_rigorous_bounds_δ0, "rigorous bounds", "estimates"))** for `δ₀`
- **$(ifelse(use_rigorous_bounds_C_B, "rigorous bounds", "estimates"))** for `C_B`
"""

# ╔═╡ d465c80c-9e34-4198-93b2-854243aea6d4
md"""
We start by computing an enclosure of $\alpha_0$. Recall that it is given by

$\alpha_0 = \sup_{x \in \mathbb{T}} \left|\frac{w(x)}{2u_0(x)}\right|.$

We can plot $\frac{w(x)}{2u_0(x)}$ for $x \in (0, pi)$ to see how it behaves.
"""

# ╔═╡ fa901388-6f0a-4de5-8323-d9b4832062f6
α0_xs, α0_ys = let xs = range(Arb(0), π, length = 200)[2:end]
    ys = similar(xs)
    f(x) = u0.w(x) / 2u0(x)(u0.α)
    Threads.@threads for i in eachindex(xs)
        ys[i] = f(xs[i])
    end
    xs, ys
end

# ╔═╡ 1a9363ba-a8bb-4ee9-b4ed-83c31e196b49
md"""
Compute an enclosure of $\alpha_0$
"""

# ╔═╡ 13373308-42e3-44ac-bf99-2a5145e6a5df
α0 = alpha0(u0, verbose = true)

# ╔═╡ 28f1459b-ac86-4316-9523-cb50b0b3a692
let pl = plot(legend = :bottomright)
    plot!(pl, α0_xs, α0_ys, ribbon = radius.(Arb, α0_ys), label = "", m = :dot, ms = 1)
    hline!(pl, [α0], ribbon = [radius(Arb, α0)], color = :green, label = "α₀")
end

# ╔═╡ 073a34c8-a2e9-434c-b047-d1f1fffd590b
md"""
Next we compute an enclosure of $A$. For a fixed $x$ we can compute an expansion around $\alpha = 0$ of $F_0(u_0)$ corresonding to the defect for the given $x$.
"""

# ╔═╡ 054103ee-ac90-4310-afa3-55291712829d
F0(u0)(Arb(2.5)) # Note that the expansion is printed with the variable x and not α

# ╔═╡ 38e338b0-f7e7-45e6-b84d-45c6946736b3
md"""
We can now plot the coefficient in front of $\alpha^2$ as a function of $x$.
"""

# ╔═╡ 51603714-bb7e-4691-b80d-7b18cf159b94
A_xs, A_ys = let xs = range(Arb("1e-1"), π, length = 200)
    ys = similar(xs)
    f = let F0 = F0(u0)
        x -> begin
            expansion = F0(x)
            @assert iszero(expansion[0]) && iszero(expansion[1])
            expansion[2]
        end
    end
    Threads.@threads for i in eachindex(xs)
        ys[i] = f(xs[i])
    end
    xs, ys
end

# ╔═╡ 616c4ebf-b2ca-4c5d-a300-30f4c1918223
A_xs_asym, A_ys_asym =
    let xs = exp.(range(log(Arb("1e-10")), log(Arb("1e-1")), length = 200))

        ys = similar(xs)
        f = let F0 = F0(u0, Asymptotic(), ϵ = one(Arb))
            x -> begin
                expansion = F0(x)
                @assert iszero(expansion[0]) && iszero(expansion[1])
                expansion[2]
            end
        end
        Threads.@threads for i in eachindex(xs)
            ys[i] = f(xs[i])
        end
        xs, ys
    end

# ╔═╡ 0044f472-bfd7-462a-aee8-256913beaad3
A = if use_rigorous_bounds_δ0
    delta0(u0, verbose = true, maxevals = 10000)[2]
else
    max(maximum(abs.(A_ys)), maximum(abs.(A_ys_asym)))
end

# ╔═╡ 03cf0016-f828-4049-90ae-5ce560ab644b
let pl = plot(legend = :bottomright)
    plot!(pl, A_xs, A_ys, ribbon = radius.(Arb, A_ys), label = "", m = :dot, ms = 1)
    hline!([A], ribbon = [radius(Arb, A)], color = :green, label = "A bound")
    hline!([-A], ribbon = [radius(Arb, A)], color = :green, label = "")
end

# ╔═╡ 17860367-ec0a-4dd4-8650-7bd99c5320c6
let pl = plot(legend = :bottomright)
    plot!(
        pl,
        A_xs_asym,
        A_ys_asym,
        ribbon = radius.(Arb, A_ys_asym),
        label = "",
        m = :dot,
        ms = 1,
        xaxis = :log10,
    )
    hline!([A], ribbon = [radius(Arb, A)], color = :green, label = "A bound")
    hline!([-A], ribbon = [radius(Arb, A)], color = :green, label = "")
end

# ╔═╡ e12220bf-4434-424d-afcc-bf204fd81411
md"""
Finally we compute an enclosure of $B$. Similar to for $A$ we can compute an expansion around $\alpha = 0$ of $T_0(u_0)$ for a fixed $x$.
"""

# ╔═╡ 353a6fb2-054c-4398-9c6f-2eafda12e4fa
T0(u0)(Arb(2.5)) # Note that the expansion is printed with the variable x and not α

# ╔═╡ ce8d3cb4-ee97-4f9a-8466-e1709b478bd5
md"""
Note that the first coefficient in the expansion is $1$. We can plot the linear coefficient as a function of $x$.
"""

# ╔═╡ 9a568a65-91d0-4758-9cbb-15cc088f9aef
B_xs, B_ys = let xs = range(Arb(0), π, length = 200)[2:end]
    ys = similar(xs)
    f = let T0 = T0(u0)
        x -> begin
            expansion = T0(x)
            @assert isone(expansion[0])
            expansion[1]
        end
    end
    Threads.@threads for i in eachindex(xs)
        ys[i] = f(xs[i])
    end
    xs, ys
end

# ╔═╡ f01b115e-0f4a-403a-807c-92d876c7e7c6
B_xs_asym, B_ys_asym =
    let xs = exp.(range(log(Arb("1e-10")), log(Arb("1e-1")), length = 200))
        ys = similar(xs)
        f = let T0 = T0(u0, Asymptotic(), ϵ = one(Arb))
            x -> begin
                expansion = T0(x)
                @assert isone(expansion[0])
                expansion[1]
            end
        end
        Threads.@threads for i in eachindex(xs)
            ys[i] = f(xs[i])
        end
        xs, ys
    end

# ╔═╡ 17755368-cba0-446a-9d64-edcd32b11be2
B = if use_rigorous_bounds_C_B
    CB(u0, verbose = true)[1]
else
    min(minimum(B_ys), minimum(B_ys_asym))
end

# ╔═╡ feb38bc2-3963-4143-8e1c-9077a967fba4
let pl = plot(legend = :bottomright)
    plot!(pl, B_xs, B_ys, ribbon = radius.(Arb, B_ys), label = "", m = :dot, ms = 1)
    hline!([B], ribbon = [radius(Arb, B)], color = :green, label = "A bound")
end

# ╔═╡ 89599cbf-2571-4475-a360-bc4e165d71d8
let pl = plot(legend = :topleft)
    plot!(
        pl,
        B_xs_asym,
        B_ys_asym,
        ribbon = radius.(Arb, B_ys_asym),
        label = "",
        m = :dot,
        ms = 1,
        xaxis = :log10,
    )
    hline!([B], ribbon = [radius(Arb, B)], color = :green, label = "A bound")
end

# ╔═╡ 11e831ac-4e18-4822-bc53-ee99394ab64f
md"""
The inequality we want to check is $A \leq \frac{B^2}{4\alpha_0}$
"""

# ╔═╡ e21e7ba1-7bd9-47c6-aa75-cebc4e413c4d
A

# ╔═╡ b142ad21-e34e-42b5-b73f-8a32e71feb26
B^2 / 4α0

# ╔═╡ 04ea9588-f613-4715-9f82-d23b786cc0ce
A < B^2 / 4α0

# ╔═╡ Cell order:
# ╟─28ec11ff-acb6-4be8-9cdd-45c42cfb839d
# ╟─4e3eaa14-583c-11ec-0fa5-e12df02e492f
# ╟─fa84cc84-d9af-4f76-b01b-a088572626bf
# ╟─6055dad1-da62-4048-a6fb-e8da26b5ef9a
# ╠═f27e0113-5e77-44af-958a-4b9e4a13a40b
# ╠═34009ff0-dadd-4464-9187-890723d94a3e
# ╟─5a06011d-4fd2-47fd-adc0-e8828c0696eb
# ╟─f101b712-69a4-4f38-add8-f3a43a43a5d6
# ╟─c6207967-17ea-41ff-8049-5e940d62374e
# ╟─d465c80c-9e34-4198-93b2-854243aea6d4
# ╠═fa901388-6f0a-4de5-8323-d9b4832062f6
# ╟─1a9363ba-a8bb-4ee9-b4ed-83c31e196b49
# ╠═13373308-42e3-44ac-bf99-2a5145e6a5df
# ╠═28f1459b-ac86-4316-9523-cb50b0b3a692
# ╟─073a34c8-a2e9-434c-b047-d1f1fffd590b
# ╠═054103ee-ac90-4310-afa3-55291712829d
# ╟─38e338b0-f7e7-45e6-b84d-45c6946736b3
# ╟─51603714-bb7e-4691-b80d-7b18cf159b94
# ╟─616c4ebf-b2ca-4c5d-a300-30f4c1918223
# ╠═0044f472-bfd7-462a-aee8-256913beaad3
# ╟─03cf0016-f828-4049-90ae-5ce560ab644b
# ╟─17860367-ec0a-4dd4-8650-7bd99c5320c6
# ╟─e12220bf-4434-424d-afcc-bf204fd81411
# ╠═353a6fb2-054c-4398-9c6f-2eafda12e4fa
# ╟─ce8d3cb4-ee97-4f9a-8466-e1709b478bd5
# ╟─9a568a65-91d0-4758-9cbb-15cc088f9aef
# ╟─f01b115e-0f4a-403a-807c-92d876c7e7c6
# ╠═17755368-cba0-446a-9d64-edcd32b11be2
# ╟─feb38bc2-3963-4143-8e1c-9077a967fba4
# ╟─89599cbf-2571-4475-a360-bc4e165d71d8
# ╟─11e831ac-4e18-4822-bc53-ee99394ab64f
# ╠═e21e7ba1-7bd9-47c6-aa75-cebc4e413c4d
# ╠═b142ad21-e34e-42b5-b73f-8a32e71feb26
# ╠═04ea9588-f613-4715-9f82-d23b786cc0ce
