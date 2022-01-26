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

# ╔═╡ 906f72a8-f912-11eb-15d0-ef321ee21aa5
begin
    using Pkg, Revise
    Pkg.activate("../")
    using Arblib, ArbExtras, HighestCuspedWave, Plots, PlutoUI

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ bf052495-c3a5-4c75-9fb6-24e113c9aedf
md"""
# Fractional KdV Equations for fixed $\alpha$
This notebook contains illustrates the computer assisted part of the proof of existence of a highest cusped traveling wave for the fractional KdV equations for a fixed $\alpha$. Let $\alpha \in (-1, 0)$, the equation we are interested in is given by

$u_t + \left(\frac{u^2}{2}\right)_x = H^{-\alpha}[u]$

where $H^{-\alpha}[u]$ is the fractional Hilbert transform of $u$. For the details of the proof see the paper, this notebook focuses on the computer assisted parts.
"""

# ╔═╡ b0ed36f7-7252-4c50-b840-ec09760691fa
md"""
What we actually prove is that there is $2\pi$ periodic, even solution which at $x = 0$ behaves like

$u(x) = |x|^{-\alpha} + \mathcal{O}(|x|^p).$

The general idea is to construct an approximate solution, $u_0$, of the equation and then prove that there is a true solution

$u(x) = u_0(x) + |x|^{-\alpha}u_e(x)$

where $u_e \in L^\infty(\mathbb{T})$ is an error term which is bounded on the interval. To do this we use a fixed point argument and for it to go through we need to control three different values that depend on the choice of $u_0$.
- $\alpha_0 = \sup_{x \in \mathbb{T}} \left|\frac{w(x)}{2u_0(x)}\right|;$
- The defect, $\delta_0$, how far our approximate solution is from satisfying the differential equation;
- A number $C_B$ which depends on the norm of an operator that occurs in the fixed point equation.
They need to satisfy the **inequality**

$\delta \leq \frac{1}{4\alpha_0 \beta^2}$

where $\beta = \frac{1}{1 - C_B}$ and we require that $C_B < 1$.
"""

# ╔═╡ 4854d9df-785c-484d-96c5-ca5942a01a69
md"""
The first step is to set the value of $\alpha$ that we use.
"""

# ╔═╡ b376a940-bb33-4e6d-bd26-5406ecee51b8
α = Arblib.add_error!(Arb(-0.3), Arb(1e-6))

# ╔═╡ 3e54d065-f8d7-4660-b2c2-4b53e5ade4f5
md"""
We then compute the ansatz we will use and plot it.
"""

# ╔═╡ e5b763ed-ab57-4324-848a-cbfb84000ecb
u0 = FractionalKdVAnsatz(α)

# ╔═╡ f59706aa-f167-43f7-b154-2f8fb5f4eb77
let xs = range(0, π, length = 200)
    ys = u0.(xs)
    ys[1] = 0
    plot(xs, ys, ribbon = Arblib.radius.(Arb, ys), label = "")
end

# ╔═╡ cef61488-2c3c-43a5-bcc3-c33fee24419c
md"""
## Computing constants
To prove the result we need to bound three different constants, $\alpha_0$, $\delta_0$ and $C_B$. They are all three given by the supremum for $x \in [0, \pi]$ using the following three functions
- For $\alpha_0$ the function is $\frac{w(x)}{2u_0(x)}$
- For $\delta_0$ the function is `F0(u0)`
- For $C_B$ the function is `T0(u0)`
For the result to hold we need the constants to satisfy

$$\delta_0 \leq \frac{1}{4\alpha_0 \beta^2}$$

where $\beta = \frac{1}{1 - C_B}$.

Computing $\alpha_0$ and $\delta_0$ is very quick whereas computing $C_B$ is more costly. For that reason we first compute enclosures of $\alpha_0$ and $\delta_0$ and then we check so that $C_B$ satisfies

$$C_B \leq 1 - 2\sqrt{\alpha_0\delta_0}$$

ensuring that the required inequality holds.
"""

# ╔═╡ 1e1392ee-0439-4a72-b112-5a2795a06fbc
md"""
We compute an enclosure of $\alpha_0$ and also plot the value of the corresponding function depending on $x$.
"""

# ╔═╡ 36d31c09-9882-4add-9a51-2c1d43c6d64a
α₀ = alpha0(u0, verbose = true)

# ╔═╡ 35cc4d3c-2ec8-4e08-b721-cca697141276
let pl = plot(legend = :bottomright)
    xs = range(Arb(0), π, length = 100)
    ys = map(x -> u0.w(x) / 2u0(x), xs)
    plot!(pl, xs, ys, ribbon = Arblib.radius.(Arb, ys), label = "")
    hline!(pl, [α₀], ribbon = [Arblib.radius(Arb, α₀)], color = :green, label = "α₀")
end

# ╔═╡ 28915e75-308c-4a93-aedd-e1ed1955603f
md"""
Next we compute an enclosure of $\delta_0$ and plot the defect as a function of $x$. Since the maximum of the defect is attained close to zero we make two plots, one on the interval $[0, \pi]$ and one on the interval $[0, 0.1]$ using a log-scaling.
"""

# ╔═╡ 141836ea-b6cd-4696-a86d-7dc7d1964131
δ₀ = HighestCuspedWave.delta0(u0, verbose = true)

# ╔═╡ 44d3ee74-e7e3-45a9-9d09-e0d3565dd4e5
let pl = plot(legend = :bottomright)
    _, xs, ys = delta0_estimate(u0, n = 500, return_values = true, ϵ = 1e-2)
    plot!(pl, xs, ys, ribbon = Arblib.radius.(Arb, ys), label = "")
    hline!(pl, [δ₀], ribbon = Arblib.radius.(Arb, [δ₀]), color = :green, label = "δ₀")
    hline!(pl, [-δ₀], ribbon = Arblib.radius.(Arb, [-δ₀]), color = :green, label = "")
end

# ╔═╡ 131ee414-0b04-4fc2-b4c4-c15d78ddce0b
let pl = plot(legend = :bottomleft, xaxis = :log10)
    xs = exp.(range(log(Arb("1e-10")), log(Arb("1e-1")), length = 200))
    ys = map(F0(u0, Asymptotic()), xs)
    plot!(pl, xs, ys, ribbon = Arblib.radius.(Arb, ys), label = "")
    hline!(pl, [δ₀], ribbon = Arblib.radius.(Arb, [δ₀]), color = :green, label = "δ₀")
    hline!(pl, [-δ₀], ribbon = Arblib.radius.(Arb, [-δ₀]), color = :green, label = "")
end

# ╔═╡ 9b1a28fa-23ed-4012-a86f-f4bb0cc7d279
md"""
Now we can compute the required bound for $C_B$.
"""

# ╔═╡ 9f7fbf89-e36b-412d-9bfc-9e51537747ec
C_B_goal = 1 - 2Arblib.sqrtpos!(zero(α₀), α₀ * δ₀)

# ╔═╡ c3427f8e-9c4c-424a-adc1-0c378fb2613a
norm_zero = T0(u0, Asymptotic())(Arb(0))

# ╔═╡ ca251b0a-ee3c-4754-bc84-dab89a46c8da
norm_zero < C_B_goal

# ╔═╡ a60681bc-905c-4d9e-89b1-1e28b8ae270f
md"""
We then prove that $C_B$ is bounded by this value and also plot $C_B$ together with this bound.

Since the bound takes some time to compute you first need to check this box for it to run

Prove bound for $C_B$: $(@bind bound_C_B CheckBox(default = false))
"""

# ╔═╡ b2a15239-cedf-4143-bacc-e7a0fb53943e
bound_holds = if bound_C_B
    @time CB_bounded_by(u0, lbound(C_B_goal), verbose = true)
else
    missing
end

# ╔═╡ 8ec4333e-7eee-4dc3-a5eb-1154866bbaf4
if ismissing(bound_holds)
    md"Bound not checked"
elseif bound_holds
    md"The bound holds! So the traveling wave exists! 🎉🎉🎉"
else
    md"The bound could not be proved to hold 😦"
end

# ╔═╡ da412f46-bb81-44a1-a12e-98a542db8320
C_B_xs, C_B_ys = let xs = range(Arb(0), π, length = 100)[2:end]
    ys = similar(xs)
    f = T0(u0, Ball())
    Threads.@threads for i in eachindex(xs)
        ys[i] = f(xs[i])
    end

    xs, ys
end

# ╔═╡ 58f1233b-30d1-4f8c-abf4-078df6ccbee4
let pl = plot(legend = :bottomright, ylims = (NaN, 1))
    plot!(
        pl,
        C_B_xs,
        C_B_ys,
        ribbon = radius.(Arb, C_B_ys),
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!(
        pl,
        [C_B_goal],
        ribbon = [radius(Arb, C_B_goal)],
        color = :green,
        label = "C_B goal",
    )
end

# ╔═╡ Cell order:
# ╠═906f72a8-f912-11eb-15d0-ef321ee21aa5
# ╟─bf052495-c3a5-4c75-9fb6-24e113c9aedf
# ╟─b0ed36f7-7252-4c50-b840-ec09760691fa
# ╟─4854d9df-785c-484d-96c5-ca5942a01a69
# ╠═b376a940-bb33-4e6d-bd26-5406ecee51b8
# ╟─3e54d065-f8d7-4660-b2c2-4b53e5ade4f5
# ╟─e5b763ed-ab57-4324-848a-cbfb84000ecb
# ╟─f59706aa-f167-43f7-b154-2f8fb5f4eb77
# ╟─cef61488-2c3c-43a5-bcc3-c33fee24419c
# ╟─1e1392ee-0439-4a72-b112-5a2795a06fbc
# ╠═36d31c09-9882-4add-9a51-2c1d43c6d64a
# ╠═35cc4d3c-2ec8-4e08-b721-cca697141276
# ╟─28915e75-308c-4a93-aedd-e1ed1955603f
# ╠═141836ea-b6cd-4696-a86d-7dc7d1964131
# ╠═44d3ee74-e7e3-45a9-9d09-e0d3565dd4e5
# ╟─131ee414-0b04-4fc2-b4c4-c15d78ddce0b
# ╟─9b1a28fa-23ed-4012-a86f-f4bb0cc7d279
# ╠═9f7fbf89-e36b-412d-9bfc-9e51537747ec
# ╠═c3427f8e-9c4c-424a-adc1-0c378fb2613a
# ╠═ca251b0a-ee3c-4754-bc84-dab89a46c8da
# ╟─a60681bc-905c-4d9e-89b1-1e28b8ae270f
# ╠═b2a15239-cedf-4143-bacc-e7a0fb53943e
# ╟─8ec4333e-7eee-4dc3-a5eb-1154866bbaf4
# ╠═da412f46-bb81-44a1-a12e-98a542db8320
# ╠═58f1233b-30d1-4f8c-abf4-078df6ccbee4
