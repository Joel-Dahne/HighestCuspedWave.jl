### A Pluto.jl notebook ###
# v0.19.11

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

# ╔═╡ a0ab3d57-b420-43c2-b69b-c403dde1f3ad
begin
    using Pkg, Revise
    Pkg.activate("../")
    using Arblib, ArbExtras, Folds, HighestCuspedWave, LaTeXStrings, Plots, PlutoUI

    setprecision(Arb, 128)

    nothing
end

# ╔═╡ 3426f2ac-f96f-11eb-22b0-2b3f9ccb38b9
md"""
# Fractional KdV Equation for $\alpha$ close to one

This notebook contains the computer assisted part of the proof of existence of a $2\pi$-periodic highest cusped wave for the fractional KdV equations for $\alpha$ in the interval $(-1, -1 + \epsilon)$. The fractional KdV equation is given by

$f_t + f f_x = |D|^\alpha f_x$

For traveling waves the ansatz $f(x, t) = \varphi(x - ct)$ reduces the equation to

$-c\varphi' + \varphi \varphi' = |D|^\alpha \varphi'$

where $c$ is the wave speed.
"""

# ╔═╡ d9575be4-0776-46d2-ad22-9fc8bc614a49
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

# ╔═╡ 5119010f-c9b1-41f5-9dac-cbb6f2b3a24d
md"""
This is the ϵ we use for the computations.
"""

# ╔═╡ 8747f80d-82f6-4e2f-bfba-d87280649950
ϵ = Arb(1e-4)

# ╔═╡ 52ef575d-cb95-45f9-8eb3-d0d97b625e77
md"""
**NOTE:** in the code we don't use the subscript $\alpha$ but in most cases use `0`, for example $u_\alpha$ becomes `u0`.
"""

# ╔═╡ 4767de5e-bbfb-4e3a-9f83-8036de0916cc
md"""
## Compute approximation
"""

# ╔═╡ 69dde124-8988-4f70-837d-01e940d199e4
u0 = BHKdVAnsatz(ϵ, BHAnsatz{Arb}())

# ╔═╡ 0a7c70da-f6d6-4484-baf5-0ae51ef3e349
md"""
We will work on the interval $[0, \pi]$ where $u_\alpha$ looks like this.
"""

# ╔═╡ 1b9e2283-03f9-4f5a-9143-85984586d77c
let xs = range(Arb(0), π, length = 100)
    ys = Folds.map(x -> u0(x), xs)
    ys[1] = 0
    plot(xs, ys, ribbon = radius.(Arb, ys), label = "", m = :circle, ms = 1)
end

# ╔═╡ 43ff127c-f7fa-4ff0-9827-36fc9507fb0b
md"## Computing constants
"

# ╔═╡ f0baf2ec-3f73-4d55-9ce4-754d94d7f3ce
md"""
The code can either compute rigorous error bounds for $\delta_\alpha$ and $D_\alpha$ or use estimates. The estimates are given by simply evaluating the corresponding functions on a number of points and taking the maximum. For the defect $\delta_\alpha$ we also make sure to use points asymptically close to $0$ since that's where the largest defect is found. Check the constants to use rigorous error bounds for
-  $δ_\alpha$ $(@bind use_rigorous_bounds_δ0 CheckBox(default = false))
-  $D_\alpha$ $(@bind use_rigorous_bounds_D0 CheckBox(default = false))
Notice that the rigorous error bounds take **significantly** longer time to compute with. On the order of 30 minutes using 12 threads.
"""

# ╔═╡ 3e6b7582-bb9f-46be-84de-f568dec6780e
md"""
The code uses
- **$(ifelse(use_rigorous_bounds_δ0, "rigorous bounds", "estimates"))** for $δ_\alpha$
- **$(ifelse(use_rigorous_bounds_D0, "rigorous bounds", "estimates"))** for $D_\alpha$
"""

# ╔═╡ d6c88a1f-94e2-4e01-bc6a-6f201f93c379
md"""
### Bound $n_\alpha$
"""

# ╔═╡ b24c29d2-1933-49c1-94e4-7a765acf355e
n0_time, n0 = begin
    n0_time = @elapsed n0 = n0_bound(u0)
    n0_time, n0
end

# ╔═╡ f1dce520-a035-43e6-9e08-4696a14c5a54
n0_xs, n0_ys = let xs = range(Arb(0), π, length = 100)[2:end]
    N = x -> u0.w(x) / 2u0(x)
    ys = Folds.map(N, xs)
    xs, ys
end

# ╔═╡ 22573416-3918-4bb9-9ec2-06b87d923c1d
n0_asym_xs, n0_asym_ys =
    let xs = exp.(range(log(Arb("1e-100")), log(Arb("1e-1")), length = 100))
        N = x -> u0.w(x) / 2u0(x, Asymptotic())
        ys = Folds.map(N, xs)
        xs, ys
    end

# ╔═╡ ab9d59df-3488-4f8a-a321-d23aab7e01d4
let pl = plot(legend = :bottomright)
    plot!(pl, n0_xs, n0_ys, ribbon = radius.(Arb, n0_ys), label = "", m = :circle, ms = 1)
    hline!(pl, [n0], ribbon = [radius(Arb, n0)], color = :green, label = "n0 bound")
    pl
end

# ╔═╡ 681e590a-a297-4dcc-8528-2f24e05df30f
let pl = plot(legend = :bottomleft, xaxis = :log10)
    plot!(
        pl,
        n0_asym_xs,
        n0_asym_ys,
        ribbon = radius.(Arb, n0_asym_ys),
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!(pl, [n0], ribbon = [radius(Arb, n0)], color = :green, label = "n bound")
    pl
end

# ╔═╡ af8899b0-eac1-442d-90ef-9d399aeb170c
md"""
### Bound $\delta_\alpha$
"""

# ╔═╡ cf99b665-f81a-4e17-8b9f-9357904cf676
md"""
The code uses **$(ifelse(use_rigorous_bounds_δ0, "rigorous bounds", "estimates"))** for $\delta_\alpha$.
"""

# ╔═╡ 3d72c6ae-5b44-491f-bfee-c8ea23224ea9
md"""
For the defect we do three different plots. One non-asymptotic plot on $[0.1, \pi]$, one asymptotic plot on $[10^{-5}, 0.1]$ and when even more asymptotic plot on $[10^{-100}, 10^{-5}]$.
"""

# ╔═╡ b8c5ba34-748e-4c4b-be9c-135240287351
δ0_xs, δ0_ys = let xs = range(Arb(0.1), Arb(π), length = 100)
    ys = Folds.map(F0(u0), xs)
    xs, ys
end

# ╔═╡ f88046ed-d4c4-4ce4-a424-96c87dc87711
δ0_asym_xs, δ0_asym_ys =
    let xs = exp.(range(log(Arb("1e-5")), log(Arb("1e-1")), length = 200))
        ys = Folds.map(F0(u0, Asymptotic(), ϵ = 2xs[end]), xs)
        xs, ys
    end

# ╔═╡ 1bb84607-4f0f-4e7b-a24f-258b4e581c2c
δ0_very_asym_xs, δ0_very_asym_ys =
    let xs = exp.(range(log(Arb("1e-100")), log(Arb("1e-5")), length = 200))
        ys = Folds.map(F0(u0, Asymptotic(), ϵ = 2xs[end]), xs)
        xs, ys
    end

# ╔═╡ 150a963b-03e2-404e-98e4-0fa2cd516dd3
δ0_bound, δ0_time = if use_rigorous_bounds_δ0
    δ0_time = @elapsed δ0_bound = delta0_bound(u0, verbose = true)
    δ0_bound, δ0_time
else
    missing, missing
end

# ╔═╡ fa5c8c0b-a5b4-4438-a219-1604c065897e
δ0 = if use_rigorous_bounds_δ0
    δ0_bound
else
    max(maximum(abs.(δ0_ys)), maximum(abs.(δ0_asym_ys)), maximum(abs.(δ0_very_asym_ys)))
end

# ╔═╡ 6dfdb4b5-fbd1-4b0d-96d9-8d4985c7b0dd
md"""
### Bound $D_\alpha$
"""

# ╔═╡ cd53f6b7-6a4c-478e-9cd6-c140b7e56d92
md"""
The code uses **$(ifelse(use_rigorous_bounds_D0, "rigorous bounds", "estimates"))** for $D_\alpha$
"""

# ╔═╡ baae040c-58f6-4657-b92f-6ca96740cbc8
D0_xs, D0_ys = let xs = range(Arb(0.1), π, length = 100)
    ys = Folds.map(T0(u0, Ball()), xs)
    xs, ys
end

# ╔═╡ 26345012-4f0c-454c-b3f4-ee9dbffb53d3
D₀_bound, D0_time = if use_rigorous_bounds_D0
    D0_time = @elapsed D₀_bound = D0_bound(u0, verbose = true)
    D₀_bound, D0_time
else
    missing, missing
end

# ╔═╡ b0577d0f-77ba-4035-9d3b-ae4d6e5c624f
D0 = if use_rigorous_bounds_D0
    D₀_bound
else
    maximum(D0_ys)
end

# ╔═╡ 89d54f92-8b3a-4913-875d-31068856fb62
let pl = plot(legend = :bottomright, xlims = (0, NaN))
    plot!(
        pl,
        D0_xs,
        D0_ys,
        ribbon = Arblib.radius.(Arb, D0_ys),
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!(pl, [D0], ribbon = [Arblib.radius(Arb, D0)], color = :green, label = "D0 bound")
    pl
end

# ╔═╡ 34ca657f-e539-49c4-8e44-8f897c816e5f
md"""
The inequality we want to check is $\delta_\alpha < \frac{(1 - D_\alpha)^2}{4n_\alpha}$
"""

# ╔═╡ 27854d39-0483-4db1-9653-d3f0761f1205
δ0

# ╔═╡ bacd05e6-7d35-4f8d-97e6-80cfd6c231d6
δ0_goal = (1 - D0)^2 / 4n0

# ╔═╡ ac03e920-25ad-4127-ad87-00e907701da3
let pl = plot()
    plot!(
        pl,
        δ0_xs,
        δ0_ys,
        ribbon = Arblib.radius.(Arb, δ0_ys),
        m = :circle,
        ms = 1,
        label = "",
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "δ0 bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    hline!([δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "δ0 goal")
    hline!([-δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "")
    pl
end

# ╔═╡ 9e18768c-40d0-45a8-b1fe-01d092b50f52
let pl = plot()
    plot!(
        pl,
        δ0_asym_xs,
        δ0_asym_ys,
        ribbon = Arblib.radius.(Arb, δ0_asym_ys),
        m = :circle,
        ms = 1,
        label = "",
        xaxis = :log10,
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "δ0 bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    hline!([δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "δ0 goal")
    hline!([-δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "")
    pl
end

# ╔═╡ 15e21dde-fb44-47d8-83d8-9f5ffffab74d
let pl = plot(legend = :bottomright)
    plot!(
        pl,
        δ0_very_asym_xs,
        δ0_very_asym_ys,
        ribbon = Arblib.radius.(Arb, δ0_very_asym_ys),
        m = :circle,
        ms = 1,
        label = "",
        xaxis = :log10,
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "δ0 bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    hline!([δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "δ0 goal")
    hline!([-δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "")
end

# ╔═╡ c3ce65c8-be0a-46fc-ab17-725aa62c60d6
δ0 < δ0_goal

# ╔═╡ fe5b35d4-e831-4d18-a73a-b6faac1d51a2
md"""
# Prepare for publishing
We compute rounded values of the upper bounds that are given in the paper, as well as produce figures for the paper.
"""

# ╔═╡ 9b641a34-21e0-4690-acdc-fec29f4603d9
md"""
**TODO:** This should possibly be done in a `round_for_publishing` method.
"""

# ╔═╡ 7ec5d2ec-de6d-4c2b-abd4-e7a4e60f92ff
n0_rounded = round(Arblib.get_d(ubound(n0), RoundUp), RoundUp, sigdigits = 3)

# ╔═╡ 8ed0be92-8c35-4495-880f-9c653013333a
δ0_rounded = round(Arblib.get_d(ubound(δ0), RoundUp), RoundUp, sigdigits = 2)

# ╔═╡ 45a99eb3-a82c-456a-ae4d-b1068141b3de
D0_rounded = round(Arblib.get_d(ubound(D0), RoundUp), RoundUp, sigdigits = 3)

# ╔═╡ f630741a-ce2c-4152-9aac-374129ce5022
md"""
Check that the inequality holds after rounding.
"""

# ╔═╡ d3f87464-24d2-4b0d-b176-5d5722a18ce2
inequality_holds_rounded = Arb(δ0_rounded) < (1 - Arb(D0_rounded))^2 / 4Arb(n0_rounded)

# ╔═╡ 7a5071ea-47bd-4a02-a9c0-8d23bfdfc9b8
md"""
When giving the enclosures in the paper we want to make sure that the printed version is less than the given upper bound. Here we check that this is the case. Note that this is not important for the result to hold, just for the printed values to look nice.
"""

# ╔═╡ 147e1d3b-da36-4ac5-84d5-59d37595cfd5
n0_printed = string(n0)

# ╔═╡ 14659193-b02c-47f3-9742-e0946e1c8b4f
Arb(n0_printed) < n0_rounded

# ╔═╡ 2bc90f20-9ef2-40a2-b4b4-32b7c049bd74
δ0_printed = string(δ0)

# ╔═╡ 6b534294-6c31-4420-9817-41907802c7bc
Arb(δ0_printed) < δ0_rounded

# ╔═╡ a590ffc8-3a0c-49e5-9692-fcb4ad7b67cb
D0_printed = string(D0)

# ╔═╡ 3c2078dd-6ab0-442d-a2b2-eaa12e21bf7f
Arb(D0_printed) < D0_rounded

# ╔═╡ d47bef6f-f292-463a-bb88-34e3e062d59d
md"""
The plots in the paper were produced with `pgfplotsx`, this is not a dependency so is not enabled by default. The plots look a bit weird with the default backend.
"""

# ╔═╡ bedfbb2f-c4bd-4d0c-9d32-0c893902286d
#pgfplotsx()

# ╔═╡ 1a469179-78d9-4a23-a67c-2e66e14f2791
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
    savefig(pl, "../figures/publication/BHKdV-N.pdf")
    pl
end

# ╔═╡ 85272438-efd7-43aa-9e7a-017623c93f50
let pl = plot(
        xaxis = :log10,
        legend = :none,
        xlabel = L"x",
        ylabel = L"N_\alpha(x)",
        guidefontsize = 18,
        tickfontsize = 18,
    )
    plot!(
        pl,
        Float64.(n0_asym_xs),
        Float64.(n0_asym_ys),
        ribbon = Float64.(radius.(Arb, n0_asym_ys)),
        m = :circle,
        ms = 1,
    )
    hline!(pl, Float64[n0_rounded], color = :green, linestyle = :dash)
    savefig(pl, "../figures/publication/BHKdV-N-asymptotic.pdf")
    pl
end

# ╔═╡ d3d93fd5-13b6-42b6-9275-4dec2e893910
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"\mathcal{T}_\alpha(x)",
        guidefontsize = 18,
        tickfontsize = 18,
    )
    plot!(
        pl,
        Float64.(D0_xs),
        Float64.(D0_ys),
        ribbon = Float64.(radius.(Arb, D0_ys)),
        m = :circle,
        ms = 1,
    )
    hline!(pl, Float64[D0_rounded], color = :green, linestyle = :dash)
    savefig(pl, "../figures/publication/BHKdV-T.pdf")
    pl
end

# ╔═╡ cb8d73a3-c025-4d94-b8a7-b952c6f7b264
δ0_rounded_goal = (1 - Arb(D0_rounded))^2 / 4Arb(n0_rounded)

# ╔═╡ a1e3b94f-77ce-496b-96ff-705ed410801b
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"F_\alpha(x)",
        guidefontsize = 18,
        tickfontsize = 18,
        xlims = (0, NaN),
    )
    plot!(
        pl,
        Float64.(δ0_xs),
        Float64.(δ0_ys),
        ribbon = Float64.(radius.(Arb, δ0_ys)),
        m = :circle,
        ms = 1,
    )
    hline!(Float64[-δ0_rounded, δ0_rounded], color = :green, linestyle = :dash)
    hline!(Float64[-δ0_rounded_goal, δ0_rounded_goal], color = :red, linestyle = :dot)
    savefig(pl, "../figures/publication/BHKdV-F.pdf")
    pl
end

# ╔═╡ 9117aa37-3fe1-47ee-b168-836e47ed74c1
let pl = plot(
        xaxis = :log10,
        legend = :none,
        xlabel = L"x",
        ylabel = L"F_\alpha(x)",
        guidefontsize = 18,
        tickfontsize = 18,
    )
    plot!(
        pl,
        Float64.(δ0_asym_xs),
        Float64.(δ0_asym_ys),
        ribbon = Float64.(radius.(Arb, δ0_asym_ys)),
        m = :circle,
        ms = 1,
    )
    hline!(Float64[-δ0_rounded, δ0_rounded], color = :green, linestyle = :dash)
    hline!(Float64[-δ0_rounded_goal, δ0_rounded_goal], color = :red, linestyle = :dot)
    savefig(pl, "../figures/publication/BHKdV-F-asymptotic.pdf")
    pl
end

# ╔═╡ Cell order:
# ╠═a0ab3d57-b420-43c2-b69b-c403dde1f3ad
# ╟─3426f2ac-f96f-11eb-22b0-2b3f9ccb38b9
# ╟─d9575be4-0776-46d2-ad22-9fc8bc614a49
# ╟─5119010f-c9b1-41f5-9dac-cbb6f2b3a24d
# ╠═8747f80d-82f6-4e2f-bfba-d87280649950
# ╟─52ef575d-cb95-45f9-8eb3-d0d97b625e77
# ╟─4767de5e-bbfb-4e3a-9f83-8036de0916cc
# ╠═69dde124-8988-4f70-837d-01e940d199e4
# ╟─0a7c70da-f6d6-4484-baf5-0ae51ef3e349
# ╟─1b9e2283-03f9-4f5a-9143-85984586d77c
# ╟─43ff127c-f7fa-4ff0-9827-36fc9507fb0b
# ╟─f0baf2ec-3f73-4d55-9ce4-754d94d7f3ce
# ╟─3e6b7582-bb9f-46be-84de-f568dec6780e
# ╟─d6c88a1f-94e2-4e01-bc6a-6f201f93c379
# ╠═b24c29d2-1933-49c1-94e4-7a765acf355e
# ╟─f1dce520-a035-43e6-9e08-4696a14c5a54
# ╟─22573416-3918-4bb9-9ec2-06b87d923c1d
# ╟─ab9d59df-3488-4f8a-a321-d23aab7e01d4
# ╟─681e590a-a297-4dcc-8528-2f24e05df30f
# ╟─af8899b0-eac1-442d-90ef-9d399aeb170c
# ╟─cf99b665-f81a-4e17-8b9f-9357904cf676
# ╟─3d72c6ae-5b44-491f-bfee-c8ea23224ea9
# ╟─b8c5ba34-748e-4c4b-be9c-135240287351
# ╟─f88046ed-d4c4-4ce4-a424-96c87dc87711
# ╟─1bb84607-4f0f-4e7b-a24f-258b4e581c2c
# ╠═150a963b-03e2-404e-98e4-0fa2cd516dd3
# ╠═fa5c8c0b-a5b4-4438-a219-1604c065897e
# ╟─ac03e920-25ad-4127-ad87-00e907701da3
# ╟─9e18768c-40d0-45a8-b1fe-01d092b50f52
# ╟─15e21dde-fb44-47d8-83d8-9f5ffffab74d
# ╟─6dfdb4b5-fbd1-4b0d-96d9-8d4985c7b0dd
# ╟─cd53f6b7-6a4c-478e-9cd6-c140b7e56d92
# ╟─baae040c-58f6-4657-b92f-6ca96740cbc8
# ╠═26345012-4f0c-454c-b3f4-ee9dbffb53d3
# ╠═b0577d0f-77ba-4035-9d3b-ae4d6e5c624f
# ╟─89d54f92-8b3a-4913-875d-31068856fb62
# ╟─34ca657f-e539-49c4-8e44-8f897c816e5f
# ╠═27854d39-0483-4db1-9653-d3f0761f1205
# ╠═bacd05e6-7d35-4f8d-97e6-80cfd6c231d6
# ╠═c3ce65c8-be0a-46fc-ab17-725aa62c60d6
# ╟─fe5b35d4-e831-4d18-a73a-b6faac1d51a2
# ╟─9b641a34-21e0-4690-acdc-fec29f4603d9
# ╠═7ec5d2ec-de6d-4c2b-abd4-e7a4e60f92ff
# ╠═8ed0be92-8c35-4495-880f-9c653013333a
# ╠═45a99eb3-a82c-456a-ae4d-b1068141b3de
# ╟─f630741a-ce2c-4152-9aac-374129ce5022
# ╠═d3f87464-24d2-4b0d-b176-5d5722a18ce2
# ╟─7a5071ea-47bd-4a02-a9c0-8d23bfdfc9b8
# ╠═147e1d3b-da36-4ac5-84d5-59d37595cfd5
# ╠═14659193-b02c-47f3-9742-e0946e1c8b4f
# ╠═2bc90f20-9ef2-40a2-b4b4-32b7c049bd74
# ╠═6b534294-6c31-4420-9817-41907802c7bc
# ╠═a590ffc8-3a0c-49e5-9692-fcb4ad7b67cb
# ╠═3c2078dd-6ab0-442d-a2b2-eaa12e21bf7f
# ╟─d47bef6f-f292-463a-bb88-34e3e062d59d
# ╠═bedfbb2f-c4bd-4d0c-9d32-0c893902286d
# ╟─1a469179-78d9-4a23-a67c-2e66e14f2791
# ╟─85272438-efd7-43aa-9e7a-017623c93f50
# ╟─d3d93fd5-13b6-42b6-9275-4dec2e893910
# ╠═cb8d73a3-c025-4d94-b8a7-b952c6f7b264
# ╟─a1e3b94f-77ce-496b-96ff-705ed410801b
# ╟─9117aa37-3fe1-47ee-b168-836e47ed74c1
