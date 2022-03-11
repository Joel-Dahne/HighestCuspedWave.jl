### A Pluto.jl notebook ###
# v0.18.1

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
    using Arblib, ArbExtras, Folds, HighestCuspedWave, Plots, PlutoUI

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ 3426f2ac-f96f-11eb-22b0-2b3f9ccb38b9
md"""
# Burgers-Hilberts Equation
This notebook contains the computer assisted part of the proof for the existence of a $2\pi$-periodic highest cusped traveling wave for the Burgers-Hilbert equation. The Buregers-Hilbert equation is given by

$f_t + f f_x = H[f]$

where $H[f]$ is the Hilbert transform of $f$. For traveling waves the ansatz $f(x, t) = \varphi(x - ct)$ reduces it to the equation

$-c\varphi' + \varphi\varphi' = H[\varphi],$

where $c$ is the wave speed.
"""

# ╔═╡ a75ba6b9-fe0b-4d3a-a88b-cc8ef9492404
md"""
For details about the proof see the paper, we here give a very short overview of the reduction to a fixed point problem. 

With the ansatz $\varphi(x) = c - u(x)$ and integrating the we can further reduce it to the equation

$\frac{1}{2}u^2 = -\mathcal{H}[u]$

where

$\mathcal{H}[u](x) = \int H[u](x) - \int H[u](0).$

If we now make the ansatz

$u(x) = u_0(x) + w(x)v(x)$

with $u_0$ an approximate solution and $w$ a fixed weight function we can, with some work, reduce the problem of proving the existence of a solution to the equation to proving existence of a fixed point of the operator

$G[v] = (I - T)^{-1}(-F - Nv^2)$

Here

$N(x) = \frac{w(x)}{2u_0(x)}$

$F(x) = \frac{1}{w(x)u_0(x)}\left(\mathcal{H}[u_0](x) + \frac{1}{2}u_0(x)^2\right)$

$T[v](x) = -\frac{1}{w(x) u_0(x)}\mathcal{H}[wv](x)$

Let $n_0 = \|N\|_{L^\infty}$, $\delta_0 = \|F\|_{L^\infty}$ and $D_0 = \|T\|$. Proving the existence of a fixed point for $G$ reduces to checking the inequality

$\delta_0 \leq \frac{(1 - D_0)^2}{4n_0}.$

This notebook is concerned with bounding these three values and checking this inequality.
"""

# ╔═╡ 66888021-535e-4f26-86c0-0db989a84be1
md"""
The first step is to compute the approximate solution $u_0$. The ansatz consists of three parts
- a leading Clausen term,
- a tail of Clausen terms to get a small defect close to $x = 0$;
- a tail of Fourier terms to get a small defect globally;
We use $N_0 = 1929$ Clausen terms in the tail and they are determined by considering the fractional KdV equation with $α = -0.9997$. We take $N_1 = 16$ Fourier terms. For more details about the construction see the paper.
"""

# ╔═╡ a063a9a2-c2c2-4c99-9df1-9fce888baad2
u0 = BHAnsatz{Arb}(Arb(-0.9997), 1929, 16)

# ╔═╡ 0a7c70da-f6d6-4484-baf5-0ae51ef3e349
md"""
We will work on the interval $[0, \pi]$, where $u_0$ looks like this.
"""

# ╔═╡ 1b9e2283-03f9-4f5a-9143-85984586d77c
let xs = range(Arb(0), π, length = 100)
    ys = Folds.map(u0, xs)
    plot(xs, ys, ribbon = Arblib.radius.(Arb, ys), label = "", m = :circle, ms = 1)
end

# ╔═╡ 1e209f4f-4fc3-4f03-8b8f-d9aa977d25ff
md"""
We can also plot what would be the traveling wave solution $\varphi(x)$ corrsponding to the approximation $u_0$.
"""

# ╔═╡ 73ae2ee7-d722-4ad8-8fc7-a57781180d35
let xs = Arb.(range(-4, 4, length = 101))
    ys = Folds.map(xs) do x
        abs(x) < 0.01 ? u0(x, Asymptotic()) : u0(x)
    end
    ys = -ys
    plot(xs, ys, linewidth = 2, axis = ([], false), legend = :none)
end

# ╔═╡ 43ff127c-f7fa-4ff0-9827-36fc9507fb0b
md"""
## Bounding constants
To prove the result we need to bound three different value $n_0$, $\delta_0$ and $D_0$, given by the supremum on the interval $[0, \pi]$ for the functions $N(x)$, $F(x)$ and $T(x)$ respectively. Note that the notation $T(x)$ is a slight abuse of notation, it denotes the function which supremum gives the norm of the operator $T$, see the paper for more details about it.
"""

# ╔═╡ 821589e9-ae94-445c-9bf4-34b8a587a6f5
md"""
For historical reasons the code uses a **different notation** than the paper. The names for the functions in the library are therefore not the same as those used in the paper. We have that `alpha0` corresponds to $n_0$, `delta0` to $\delta_0$ and `CB` to $D_0$.
"""

# ╔═╡ f0baf2ec-3f73-4d55-9ce4-754d94d7f3ce
md"""
The code can either compute rigorous error bounds for the required constants or use estimates. The estimates are given by simply evaluating the corresponding functions on a number of points and taking the maximum. For the defect $\delta_0$ we also make sure to use points asymptically close to $0$ since that's where the largest defect is found. Check the constants to use rigorous error bounds for
-  $n_0$ $(@bind use_rigorous_bounds_n0 CheckBox(default = false))
-  $\delta_0$ $(@bind use_rigorous_bounds_δ0 CheckBox(default = false))
-  $D_0$ $(@bind use_rigorous_bounds_D0 CheckBox(default = false))
Notice that the rigorous error bounds take longer time to compute with. On an AMD Ryzen 9 5900X with 12 cores the computation of $n_0$ takes around 7 seconds, the computation of $\delta_0$ around 15 minutes and the computation of $D_0$ around 2 minutes. The computation of $\delta_0$ also uses a **large amount of memory**, around 20GB.
"""

# ╔═╡ 4c4cbf2a-3aec-4257-9fac-d8a0418d12d7
md"""
### Bound $n_0$
"""

# ╔═╡ 28978335-5797-4ddc-bfb4-9b04a0f9c4ac
md"""
The code uses **$(ifelse(use_rigorous_bounds_n0, "rigorous bounds", "estimates"))** for $n_0$.
"""

# ╔═╡ f1dce520-a035-43e6-9e08-4696a14c5a54
n0_xs, n0_ys = let xs = range(Arb(0), π, length = 100)[2:end]
    N = x -> u0.w(x) / 2u0(x)
    ys = Folds.map(N, xs)
    xs, ys
end

# ╔═╡ f1ca794b-e1fc-481b-9b75-1e42a0b48a58
n0_bound, n0_time = if use_rigorous_bounds_n0
    n0_time = @elapsed n0_bound = alpha0(u0, verbose = true)
    n0_bound, n0_time
else
    missing, missing
end

# ╔═╡ 61151255-15d4-45ec-a3ad-573c46d34d93
n0 = if use_rigorous_bounds_n0
    n0_bound
else
    maximum(abs.(n0_ys))
end

# ╔═╡ ab9d59df-3488-4f8a-a321-d23aab7e01d4
let pl = plot(legend = :bottomright)
    plot!(pl, n0_xs, n0_ys, ribbon = radius.(Arb, n0_ys), label = "", m = :circle, ms = 2)
    hline!(pl, [n0], ribbon = [radius(Arb, n0)], color = :green, label = "n₀")
    pl
end

# ╔═╡ 36b6ce98-6b64-4eee-a7d7-613772e4b79f
md"""
### Bound $\delta_0$
"""

# ╔═╡ 237d130e-d0c0-400d-b8bf-42e88a5a889d
md"""
The code uses **$(ifelse(use_rigorous_bounds_δ0, "rigorous bounds", "estimates"))** for $\delta_0$.
"""

# ╔═╡ 67aa36b0-b77c-4531-a248-f7d474ffd47d
md"""
For the defect we do three different plots. one non-asymptotic plot on $[0.1, \pi]$, one asymptotic plot on $[10^{-100}, 0.1]$ and when even more asymptotic plot on $[10^{-20000}, 10^{-100}]$.
"""

# ╔═╡ 67786c2c-a101-44db-a1b6-f191d2703bb0
md"""
We want to have the inequality $\delta_0 \leq \frac{(1 - D_0)^2}{4n_0}$. We therefore include the value of the right hand side in this inequality in the plots, this is the goal that we want the defect to be smaller than.
"""

# ╔═╡ b8c5ba34-748e-4c4b-be9c-135240287351
δ0_xs, δ0_ys = let xs = range(Arb(1e-3), Arb(π), length = 100)
    ys = Folds.map(F0(u0), xs)
    xs, ys
end

# ╔═╡ 1bb84607-4f0f-4e7b-a24f-258b4e581c2c
δ0_asym_xs, δ0_asym_ys =
    let xs = exp.(range(log(Arb("1e-100")), log(Arb("1e-1")), length = 200))
        ys = Folds.map(F0(u0, Asymptotic(), ϵ = 2xs[end]), xs)
        xs, ys
    end

# ╔═╡ fb6c12ad-3391-4623-a201-412335742930
δ0_very_asym_xs, δ0_very_asym_ys =
    let xs = exp.(range(log(Arb("1e-20000")), log(Arb("1e-100")), length = 200))
        ys = Folds.map(F0(u0, Asymptotic(), ϵ = 2xs[end]), xs)
        xs, ys
    end

# ╔═╡ 664df5e6-152a-4b9a-adf8-352bdef59055
δ0_bound, δ0_time, δ0_subintervals = if use_rigorous_bounds_δ0
    δ0_time = @elapsed δ0_bound, δ0_subintervals... =
        delta0(u0, return_subresults = true, verbose = true)
    δ0_bound, δ0_time, δ0_subintervals
else
    missing, missing, (Arb(NaN), Arb(NaN), Arb(NaN), Arb(NaN), Arb(NaN))
end

# ╔═╡ 150a963b-03e2-404e-98e4-0fa2cd516dd3
δ0 = if use_rigorous_bounds_δ0
    δ0_bound
else
    maximum(abs.([δ0_ys; δ0_asym_ys; δ0_very_asym_ys]))
end

# ╔═╡ 61e8da72-69c0-4af9-bab0-454442dad253
md"""
### Bound $D_0$
"""

# ╔═╡ f4ede9a7-c57e-44ce-adae-4345205fa4e4
md"""
The code uses **$(ifelse(use_rigorous_bounds_D0, "rigorous bounds", "estimates"))** for $D_0$
"""

# ╔═╡ de4546e1-4a9f-4d37-b59c-ee4509d09868
D0_xs, D0_ys = let xs = range(Arb(1e-1), π, length = 100)
    ys = Folds.map(HighestCuspedWave.T0(u0, Ball()), xs)
    xs, ys
end

# ╔═╡ de5bed6f-7079-40a9-a5eb-f315abc20ddf
D0_bound, D0_time = if use_rigorous_bounds_D0
    D0_time = @elapsed D0_bound = CB(u0, verbose = true)
    D0_bound, D0_time
else
    missing, missing
end

# ╔═╡ b0577d0f-77ba-4035-9d3b-ae4d6e5c624f
D0 = if use_rigorous_bounds_D0
    D0_bound
else
    maximum(abs.(D0_ys))
end

# ╔═╡ 358b2691-3b1b-4733-a173-acf987a221ea
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
        label = "defect",
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "δ₀ bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    hline!([δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "δ₀ goal")
    hline!([-δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "")
    pl
end

# ╔═╡ 15e21dde-fb44-47d8-83d8-9f5ffffab74d
let pl = plot()
    plot!(
        pl,
        δ0_asym_xs,
        δ0_asym_ys,
        ribbon = Arblib.radius.(Arb, δ0_asym_ys),
        m = :circle,
        ms = 1,
        label = "defect",
        xaxis = :log10,
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "δ₀ bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    hline!([δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "δ₀ goal")
    hline!([-δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "")
    pl
end

# ╔═╡ f5ee2d3b-4a4c-411b-9e46-35164f6e3e83
let pl = plot()
    plot!(
        pl,
        log.(δ0_very_asym_xs) ./ log(Arb(10)),
        δ0_very_asym_ys,
        ribbon = Arblib.radius.(Arb, δ0_very_asym_ys),
        m = :circle,
        ms = 1,
        label = "defect",
        xlabel = "log10(x)",
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "δ₀ bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    hline!([δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "δ₀ goal")
    hline!([-δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "")
    pl
end

# ╔═╡ 89d54f92-8b3a-4913-875d-31068856fb62
let pl = plot(legend = :bottomright)
    plot!(
        pl,
        D0_xs,
        D0_ys,
        ribbon = Arblib.radius.(Arb, D0_ys),
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!(pl, [D0], ribbon = [radius(Arb, D0)], color = :green, label = "D₀")
    pl
end

# ╔═╡ 3f30187c-0bdb-4cde-9c50-a9ccaaf8bb4b
md"""
# Prepare for publishing
We compute rounded values of the upper bounds that are given in the paper, as well as figures for the paper.
"""

# ╔═╡ 7aab4811-7e37-4828-8e5a-18423f9330e1
proved, n0_rounded, δ0_rounded, D0_rounded =
    HighestCuspedWave.round_for_publishing(n0, δ0, D0, sigdigits = 5)

# ╔═╡ 4d5054ff-6001-4d34-913b-a8029017d217
if proved
    if use_rigorous_bounds_n0 && use_rigorous_bounds_D0 && use_rigorous_bounds_δ0
        "The result holds! 🎉🎉🎉"
    else
        "The inequality holds but the bounds are not rigorous"
    end
else
    "Could not prove the result 😥"
end

# ╔═╡ f25bcedf-cf4f-45f6-929a-abba66e7730c
md"""
When giving the enclosures in the paper we want to make sure that the printed version is less than the given upper bound. Here we check that this is the case. Note that this is not important for the result to hold, just for the printed values to look nice.
"""

# ╔═╡ ff9462b6-3ff1-4fd9-af14-eebce04f349c
Arb(string(n0)) < n0_rounded

# ╔═╡ 7fc04a7b-1082-4af7-ac87-426e73276c18
all([Arb(string(δ)) < δ0_rounded for δ in (δ0, δ0_subintervals...)])

# ╔═╡ 121874d0-43ef-4463-8046-f91985cfc6fe
Arb(string(D0)) < D0_rounded

# ╔═╡ f0110b12-a782-4840-abbf-ed623175c276
let pl = plot(legend = :none, xlabel = "\$x\$", ylabel = "\$N(x)\$")
    plot!(pl, n0_xs, n0_ys, ribbon = radius.(Arb, n0_ys), m = :circle, ms = 2)
    hline!(pl, [n0_rounded], color = :green)
    savefig(pl, "../figures/publication/BH-N.pdf")
    pl
end

# ╔═╡ 16b4c640-6b4c-46a3-a505-e0eb7bbb2632
let pl = plot(legend = :none, xlabel = "\$x\$", ylabel = "\$T(x)\$")
    plot!(pl, D0_xs, D0_ys, ribbon = radius.(Arb, D0_ys), m = :circle, ms = 2)
    hline!(pl, [D0_rounded], color = :green, label = "C_B")
    savefig(pl, "../figures/publication/BH-T.pdf")
    pl
end

# ╔═╡ fe196d34-c822-4da0-8ed1-955d2fb6ffa6
δ0_goal_rounded = (1 - D0_rounded)^2 / 4n0_rounded

# ╔═╡ 4654b37d-fcb5-43ab-8b67-5f946637943d
let pl = plot(legend = :none, xlabel = "\$x\$", ylabel = "\$F(x)\$")
    plot!(pl, δ0_xs, δ0_ys, ribbon = radius.(Arb, δ0_ys), m = :circle, ms = 2)
    hline!([-δ0_rounded, δ0_rounded], color = :green)
    hline!([-δ0_goal_rounded, δ0_goal_rounded], color = :red)
    savefig(pl, "../figures/publication/BH-F.pdf")
    pl
end

# ╔═╡ 1c2d23d0-9e30-4830-a923-3ce15e8deb04
let pl = plot(legend = :none, xlabel = "\$x\$", ylabel = "\$F(x)\$", xaxis = :log10)
    plot!(
        pl,
        δ0_asym_xs,
        δ0_asym_ys,
        ribbon = radius.(Arb, δ0_asym_ys),
        m = :circle,
        ms = 2,
    )
    hline!([-δ0_rounded, δ0_rounded], color = :green)
    hline!([-δ0_goal_rounded, δ0_goal_rounded], color = :red)
    savefig(pl, "../figures/publication/BH-F-asymptotic-1.pdf")
    pl
end

# ╔═╡ 7e1e6379-4cde-4a6f-9390-41618a989490
let pl = plot(legend = :none, xlabel = "\$\\log_{10}(x)\$", ylabel = "\$F(x)\$")
    plot!(
        pl,
        log.(δ0_very_asym_xs) ./ log(Arb(10)),
        δ0_very_asym_ys,
        ribbon = Arblib.radius.(Arb, δ0_very_asym_ys),
        m = :circle,
        ms = 2,
    )
    hline!([-δ0_rounded, δ0_rounded], color = :green)
    hline!([-δ0_goal_rounded, δ0_goal_rounded], color = :red)
    savefig(pl, "../figures/publication/BH-F-asymptotic-2.pdf")
    pl
end

# ╔═╡ Cell order:
# ╟─a0ab3d57-b420-43c2-b69b-c403dde1f3ad
# ╟─3426f2ac-f96f-11eb-22b0-2b3f9ccb38b9
# ╟─a75ba6b9-fe0b-4d3a-a88b-cc8ef9492404
# ╟─66888021-535e-4f26-86c0-0db989a84be1
# ╠═a063a9a2-c2c2-4c99-9df1-9fce888baad2
# ╟─0a7c70da-f6d6-4484-baf5-0ae51ef3e349
# ╟─1b9e2283-03f9-4f5a-9143-85984586d77c
# ╟─1e209f4f-4fc3-4f03-8b8f-d9aa977d25ff
# ╟─73ae2ee7-d722-4ad8-8fc7-a57781180d35
# ╟─43ff127c-f7fa-4ff0-9827-36fc9507fb0b
# ╟─821589e9-ae94-445c-9bf4-34b8a587a6f5
# ╟─f0baf2ec-3f73-4d55-9ce4-754d94d7f3ce
# ╟─4c4cbf2a-3aec-4257-9fac-d8a0418d12d7
# ╟─28978335-5797-4ddc-bfb4-9b04a0f9c4ac
# ╟─f1dce520-a035-43e6-9e08-4696a14c5a54
# ╠═f1ca794b-e1fc-481b-9b75-1e42a0b48a58
# ╠═61151255-15d4-45ec-a3ad-573c46d34d93
# ╟─ab9d59df-3488-4f8a-a321-d23aab7e01d4
# ╟─36b6ce98-6b64-4eee-a7d7-613772e4b79f
# ╟─237d130e-d0c0-400d-b8bf-42e88a5a889d
# ╟─67aa36b0-b77c-4531-a248-f7d474ffd47d
# ╟─67786c2c-a101-44db-a1b6-f191d2703bb0
# ╠═358b2691-3b1b-4733-a173-acf987a221ea
# ╟─b8c5ba34-748e-4c4b-be9c-135240287351
# ╟─1bb84607-4f0f-4e7b-a24f-258b4e581c2c
# ╟─fb6c12ad-3391-4623-a201-412335742930
# ╠═664df5e6-152a-4b9a-adf8-352bdef59055
# ╠═150a963b-03e2-404e-98e4-0fa2cd516dd3
# ╟─ac03e920-25ad-4127-ad87-00e907701da3
# ╟─15e21dde-fb44-47d8-83d8-9f5ffffab74d
# ╟─f5ee2d3b-4a4c-411b-9e46-35164f6e3e83
# ╟─61e8da72-69c0-4af9-bab0-454442dad253
# ╟─f4ede9a7-c57e-44ce-adae-4345205fa4e4
# ╟─de4546e1-4a9f-4d37-b59c-ee4509d09868
# ╠═de5bed6f-7079-40a9-a5eb-f315abc20ddf
# ╠═b0577d0f-77ba-4035-9d3b-ae4d6e5c624f
# ╟─89d54f92-8b3a-4913-875d-31068856fb62
# ╟─3f30187c-0bdb-4cde-9c50-a9ccaaf8bb4b
# ╠═7aab4811-7e37-4828-8e5a-18423f9330e1
# ╟─4d5054ff-6001-4d34-913b-a8029017d217
# ╟─f25bcedf-cf4f-45f6-929a-abba66e7730c
# ╠═ff9462b6-3ff1-4fd9-af14-eebce04f349c
# ╠═7fc04a7b-1082-4af7-ac87-426e73276c18
# ╠═121874d0-43ef-4463-8046-f91985cfc6fe
# ╟─f0110b12-a782-4840-abbf-ed623175c276
# ╟─16b4c640-6b4c-46a3-a505-e0eb7bbb2632
# ╠═fe196d34-c822-4da0-8ed1-955d2fb6ffa6
# ╟─4654b37d-fcb5-43ab-8b67-5f946637943d
# ╟─1c2d23d0-9e30-4830-a923-3ce15e8deb04
# ╟─7e1e6379-4cde-4a6f-9390-41618a989490
