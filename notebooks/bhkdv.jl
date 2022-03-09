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
    using Arblib, ArbExtras, HighestCuspedWave, Plots, PlutoUI

    setprecision(Arb, 128)

    nothing
end

# ╔═╡ 3426f2ac-f96f-11eb-22b0-2b3f9ccb38b9
md"""
# From Burgers Hilberts to KdV
This notebook contains the computer assisted part of the proof for $\alpha \in (-1, -1 + \epsilon)$

$u_t + \left(\frac{u}{2}\right)_x = H^\alpha[u]$

The approximation of the wave looks like this
"""

# ╔═╡ 6c2a16c7-2bcf-4a2b-9466-38309a36b937
md"""
What we actually prove is that there is $2\pi$ periodic, even solution which at $x = 0$ behaves like

$u(x) = |x|^{-\alpha} + \mathcal{O}(|x|\log|x|).$

The general idea is to construct an approximate solution, $u_0$, of the equation and then prove that there is a true solution

$u(x) = u_0(x) + |x|\log(c + 1 / |x|)u_e(x)$

where $u_e \in L^\infty(\mathbb{T})$ is an error term which is bounded on the interval. To do this we use a fixed point argument and for it to go through we need to control three different values that depend on the choice of $u_0$.
- $\alpha_0 = \sup_{x \in \mathbb{T}} \left|\frac{w(x)}{2u_0(x)}\right|;$
- The defect, $\delta_0$, how far our approximate solution is from satisfying the differential equation;
- A number $C_B$ which depends on the norm of an operator that occurs in the fixed point equation.
They need to satisfy the **inequality**

$\delta \leq \frac{1}{4\alpha_0 \beta^2}$

where $\beta = \frac{1}{1 - C_B}$ and we require that $C_B < 1$.

The construction of an approximate solution with the right asymptotic behaviour is one of the most complicated parts of the work. It relies heavily on the so called Clausen functions

$C_s(x) = \sum_{n = 1}^\infty \frac{\cos(nx)}{n^s},\ S_s(x) = \sum_{n = 1}^\infty \frac{\sin(nx)}{n^s}$

and there derivatives with respect to the order $s$.
"""

# ╔═╡ 66888021-535e-4f26-86c0-0db989a84be1
md"The first step is to compute the approximate solution $u_0$. The ansatz consists of three parts
- the main term which is a sum of two Clausen functions
- the tail of the solution used for $\alpha = -1$.
"

# ╔═╡ ccc30af6-7eea-4fcf-932a-d5cb8bccda99
α = Arb(-0.9997)

# ╔═╡ 8747f80d-82f6-4e2f-bfba-d87280649950
ϵ = 1e-5 * (1 + α)

# ╔═╡ 69dde124-8988-4f70-837d-01e940d199e4
u0 = BHKdVAnsatz(midpoint(Arb, ϵ), BHAnsatz{Arb}(α, 1929, 16), γ = Arb(0.5))

# ╔═╡ 73ae2ee7-d722-4ad8-8fc7-a57781180d35
let xs = Arb.(range(-4, 4, length = 101))
    ys = similar(xs)
    Threads.@threads for i in eachindex(xs)
        ys[i] = u0(xs[i])
    end
    ys = u0(zero(Arb), Asymptotic()) .- ys
    plot(
        xs,
        ys,
        ribbon = Arblib.radius.(Arb, ys),
        label = "",
        linewidth = 2,
        axis = ([], false),
    )
end

# ╔═╡ 0a7c70da-f6d6-4484-baf5-0ae51ef3e349
md"We can plot the solution on $[0, \pi]$."

# ╔═╡ 1b9e2283-03f9-4f5a-9143-85984586d77c
let xs = range(Arb(0), π, length = 100)
    ys = similar(xs)
    Threads.@threads for i in eachindex(xs)
        ys[i] = u0(xs[i])
    end
    plot(xs, ys, ribbon = Arblib.radius.(Arb, ys), label = "", m = :circle, ms = 1)
end

# ╔═╡ 43ff127c-f7fa-4ff0-9827-36fc9507fb0b
md"## Computing constants
To prove the result we need to bound three different constants, `α0`, `C_B` and `δ0`. They are all three given by the supremum for `x ∈ [0, π]` using the following three functions
- For `α0` the function is `u0.w(x) / 2u0(x)`
- For `C_B` the function is `T0(u0)`
- For `δ0` the function is `F0(u0)`
For the result to hold we need the constants to satisfy
```
δ0 <= 1 / (4α0 * β^2)
```
where `β = 1 / (1 - C_B)`.

The values of `α0` and `C_B` are more or less the same for any sufficiently good approximation whereas `δ0`, the defect, depends heavily on how good of an approximation the currest solution is. For that reason we start by computing `α₀` and `C_B` and from that we can then get an upper bound for what `δ0` should be. We the compute `δ0` and see if it satisfies this bound.
"

# ╔═╡ f0baf2ec-3f73-4d55-9ce4-754d94d7f3ce
md"""
The code can either compute rigorous error bounds for the required constants or use estimates. The estimates are given by simply evaluating the corresponding functions on a number of points and taking the maximum. For the defect `δ0` we also make sure to use points asymptically close to `0` since that's where the largest defect is found. Check the constants to use rigorous error bounds for
-  `α₀` $(@bind use_rigorous_bounds_α0 CheckBox(default = false))
-  `C_B` $(@bind use_rigorous_bounds_C_B CheckBox(default = false))
-  `δ₀` $(@bind use_rigorous_bounds_δ0 CheckBox(default = false))
Notice that the rigorous error bounds take **significantly** longer time to compute with.
"""

# ╔═╡ 3e6b7582-bb9f-46be-84de-f568dec6780e
md"""
The code uses
- **$(ifelse(use_rigorous_bounds_α0, "rigorous bounds", "estimates"))** for `α₀`
- **$(ifelse(use_rigorous_bounds_C_B, "rigorous bounds", "estimates"))** for `C_B`
- **$(ifelse(use_rigorous_bounds_δ0, "rigorous bounds", "estimates"))** for `δ₀`
"""

# ╔═╡ f1dce520-a035-43e6-9e08-4696a14c5a54
α0_xs, α0_ys = let xs = range(Arb(0), π, length = 100)[2:end]
    ys = similar(xs)
    f(x) = u0.w(x) / 2u0(x)
    Threads.@threads for i in eachindex(xs)
        ys[i] = f(xs[i])
    end
    xs, ys
end

# ╔═╡ 22573416-3918-4bb9-9ec2-06b87d923c1d
α0_asym_xs, α0_asym_ys =
    let xs = exp.(range(log(Arb("1e-100")), log(Arb("1e-1")), length = 100))
        ys = similar(xs)
        f(x) = u0.w(x) / 2u0(x, Asymptotic())
        Threads.@threads for i in eachindex(xs)
            ys[i] = f(xs[i])
        end
        xs, ys
    end

# ╔═╡ 61151255-15d4-45ec-a3ad-573c46d34d93
α0 = if use_rigorous_bounds_α0
    @time alpha0(u0, verbose = true)
else
    # The value at x = 0 is π / 2
    max(maximum(α0_ys), Arb(π) / 2)
end

# ╔═╡ ab9d59df-3488-4f8a-a321-d23aab7e01d4
let pl = plot(legend = :bottomright)
    plot!(pl, α0_xs, α0_ys, ribbon = radius.(Arb, α0_ys), label = "", m = :circle, ms = 1)
    hline!(pl, [α0], ribbon = [radius(Arb, α0)], color = :green, label = "α₀")
    savefig(pl, "../figures/publication/BHKdV-N.pdf")
    pl
end

# ╔═╡ 681e590a-a297-4dcc-8528-2f24e05df30f
let pl = plot(legend = :bottomleft, xaxis = :log10)
    plot!(
        pl,
        α0_asym_xs,
        α0_asym_ys,
        ribbon = radius.(Arb, α0_asym_ys),
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!(pl, [α0], ribbon = [radius(Arb, α0)], color = :green, label = "α₀")
    savefig(pl, "../figures/publication/BHKdV-N-asymptotic.pdf")
    pl
end

# ╔═╡ 87b3712f-1cf7-4d02-b1d1-d68d7f4464d6
md"Next we want to determine the norm, `C_B`, given by the maximum of `T0(u0)` on the interval `[0, π]`. Again we can plot its value on a few points on the interval and take the maximum, giving us a non-rigorous estimate."

# ╔═╡ de4546e1-4a9f-4d37-b59c-ee4509d09868
C_B_xs, C_B_ys = let xs = range(Arb(0), π, length = 100)[2:end]
    ys = similar(xs)
    f = T0(u0)
    Threads.@threads for i in eachindex(xs)
        ys[i] = f(xs[i])
    end
    xs, ys
end

# ╔═╡ b0577d0f-77ba-4035-9d3b-ae4d6e5c624f
C_B = if use_rigorous_bounds_C_B
    @time CB(u0, verbose = true)
else
    maximum(C_B_ys)
end

# ╔═╡ 89d54f92-8b3a-4913-875d-31068856fb62
let pl = plot(legend = :bottomright)
    plot!(
        pl,
        C_B_xs,
        C_B_ys,
        ribbon = Arblib.radius.(Arb, C_B_ys),
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!(pl, [C_B], ribbon = [Arblib.radius(Arb, C_B)], color = :green, label = "C_B")
    savefig(pl, "../figures/publication/BHKdV-T.pdf")
    pl
end

# ╔═╡ 5bc3a7c1-9acd-49c4-afd4-e0a94b3b02d7
md"Now we need the defect, `δ₀`, to be smaller than `1 / (4α₀ * β^2)` where `β = 1 / (1 - C_B)`"

# ╔═╡ 25927755-307c-40b7-adc9-4577ab8c2f61
β = 1 / (1 - C_B)

# ╔═╡ 358b2691-3b1b-4733-a173-acf987a221ea
δ0_goal = 1 / (4α0 * β^2)

# ╔═╡ 67aa36b0-b77c-4531-a248-f7d474ffd47d
md"""
We can now plot the defect on the interval `[0, π]`. We do three different plots. one non-asymptotic plot on `[0.1, π]`, one asymptotic plot on `[1e-5, 0.1]` and when even more asymptotic plot on `[1e-100, 1e-5]`. For the last interval we only compute an upper bound of the absolute value, instead of an enclosure.
"""

# ╔═╡ b8c5ba34-748e-4c4b-be9c-135240287351
δ0_xs, δ0_ys = let xs = range(Arb(0), π, length = 200)[2:end]
    ys = similar(xs)
    f = F0(u0)
    Threads.@threads for i in eachindex(xs)
        ys[i] = f(xs[i])
    end
    xs, ys
end

# ╔═╡ f88046ed-d4c4-4ce4-a424-96c87dc87711
δ0_little_asym_xs, δ0_little_asym_ys =
    let xs = exp.(range(log(Arb("1e-5")), log(Arb("1e-1")), length = 100))
        ys = similar(xs)
        f = F0(u0)
        Threads.@threads for i in eachindex(xs)
            ys[i] = f(xs[i])
        end
        xs, ys
    end

# ╔═╡ 1bb84607-4f0f-4e7b-a24f-258b4e581c2c
δ0_asym_xs, δ0_asym_ys =
    let xs = exp.(range(log(Arb("1e-100")), log(Arb("1e-5")), length = 100))
        ys = similar(xs)
        f = F0(u0, Asymptotic())
        Threads.@threads for i in eachindex(xs)
            ys[i] = f(xs[i])
        end
        xs, ys
    end

# ╔═╡ 150a963b-03e2-404e-98e4-0fa2cd516dd3
δ0 = if use_rigorous_bounds_δ0
    @time delta0(u0, verbose = true)
else
    max(maximum(abs.(δ0_ys)), maximum(abs.(δ0_little_asym_ys)), maximum(abs.(δ0_asym_ys)))
end

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
    savefig(pl, "../figures/publication/BHKdV-F.pdf")
    pl
end

# ╔═╡ 9e18768c-40d0-45a8-b1fe-01d092b50f52
let pl = plot()
    plot!(
        pl,
        δ0_little_asym_xs,
        δ0_little_asym_ys,
        ribbon = Arblib.radius.(Arb, δ0_little_asym_ys),
        m = :circle,
        ms = 1,
        label = "defect",
        xaxis = :log10,
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "δ₀ bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    hline!([δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "δ₀ goal")
    hline!([-δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "")
    savefig(pl, "../figures/publication/BHKdV-F-asymptotic.pdf")
    pl
end

# ╔═╡ 15e21dde-fb44-47d8-83d8-9f5ffffab74d
let pl = plot(legend = :bottomright)
    plot!(
        pl,
        δ0_asym_xs,
        δ0_asym_ys,
        ribbon = Arblib.radius.(Arb, δ0_asym_ys),
        m = :circle,
        ms = 1,
        label = "defect upper bound",
        xaxis = :log10,
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "δ₀ bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    hline!([δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "δ₀ goal")
    hline!([-δ0_goal], ribbon = [radius(Arb, δ0_goal)], color = :red, label = "")
end

# ╔═╡ Cell order:
# ╟─a0ab3d57-b420-43c2-b69b-c403dde1f3ad
# ╟─3426f2ac-f96f-11eb-22b0-2b3f9ccb38b9
# ╟─73ae2ee7-d722-4ad8-8fc7-a57781180d35
# ╟─6c2a16c7-2bcf-4a2b-9466-38309a36b937
# ╟─66888021-535e-4f26-86c0-0db989a84be1
# ╠═ccc30af6-7eea-4fcf-932a-d5cb8bccda99
# ╠═8747f80d-82f6-4e2f-bfba-d87280649950
# ╠═69dde124-8988-4f70-837d-01e940d199e4
# ╟─0a7c70da-f6d6-4484-baf5-0ae51ef3e349
# ╟─1b9e2283-03f9-4f5a-9143-85984586d77c
# ╟─43ff127c-f7fa-4ff0-9827-36fc9507fb0b
# ╟─f0baf2ec-3f73-4d55-9ce4-754d94d7f3ce
# ╟─3e6b7582-bb9f-46be-84de-f568dec6780e
# ╟─f1dce520-a035-43e6-9e08-4696a14c5a54
# ╠═22573416-3918-4bb9-9ec2-06b87d923c1d
# ╠═61151255-15d4-45ec-a3ad-573c46d34d93
# ╟─ab9d59df-3488-4f8a-a321-d23aab7e01d4
# ╟─681e590a-a297-4dcc-8528-2f24e05df30f
# ╟─87b3712f-1cf7-4d02-b1d1-d68d7f4464d6
# ╠═de4546e1-4a9f-4d37-b59c-ee4509d09868
# ╠═b0577d0f-77ba-4035-9d3b-ae4d6e5c624f
# ╟─89d54f92-8b3a-4913-875d-31068856fb62
# ╟─5bc3a7c1-9acd-49c4-afd4-e0a94b3b02d7
# ╠═25927755-307c-40b7-adc9-4577ab8c2f61
# ╠═358b2691-3b1b-4733-a173-acf987a221ea
# ╟─67aa36b0-b77c-4531-a248-f7d474ffd47d
# ╟─b8c5ba34-748e-4c4b-be9c-135240287351
# ╠═f88046ed-d4c4-4ce4-a424-96c87dc87711
# ╟─1bb84607-4f0f-4e7b-a24f-258b4e581c2c
# ╠═150a963b-03e2-404e-98e4-0fa2cd516dd3
# ╟─ac03e920-25ad-4127-ad87-00e907701da3
# ╟─9e18768c-40d0-45a8-b1fe-01d092b50f52
# ╟─15e21dde-fb44-47d8-83d8-9f5ffffab74d
