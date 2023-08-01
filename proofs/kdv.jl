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

# ╔═╡ 906f72a8-f912-11eb-15d0-ef321ee21aa5
begin
    using Pkg
    Pkg.activate("../", io = devnull)
    using Arblib, ArbExtras, Folds, HighestCuspedWave, LaTeXStrings, Plots, PlutoUI

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ bf052495-c3a5-4c75-9fb6-24e113c9aedf
md"""
# Fractional KdV Equations for fixed $\alpha$

This notebook illustrates the computer assisted part of the proof of existence of a $2\pi$-periodic highest cusped traveling wave for the fractional KdV equations for a fixed $\alpha \in (-1, 0)$. It is related to Lemma 12.5, 12.6 and 12.7 but only treats one fixed $\alpha$ and not the full interval. It is not a part of the proof in the paper but only serves as an illustration.
"""

# ╔═╡ 4854d9df-785c-484d-96c5-ca5942a01a69
md"""
This is the $\alpha$ we use for the computations.
"""

# ╔═╡ b376a940-bb33-4e6d-bd26-5406ecee51b8
α = setball(Arb, -0.6, 1e-10)

# ╔═╡ 3e54d065-f8d7-4660-b2c2-4b53e5ade4f5
md"""
## Compute the approximation $u_\alpha$
"""

# ╔═╡ e5b763ed-ab57-4324-848a-cbfb84000ecb
u0 = FractionalKdVAnsatz(α)

# ╔═╡ 3c366570-463c-4481-b4fd-5cab7da1a007
md"""
We work on the interval $[0, \pi]$ where $u_\alpha$ looks like this.
"""

# ╔═╡ f59706aa-f167-43f7-b154-2f8fb5f4eb77
let xs = range(0, π, length = 200)
    ys = u0.(xs)
    ys[1] = 0
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

# ╔═╡ cef61488-2c3c-43a5-bcc3-c33fee24419c
md"""
## Computing constants
We are now ready to compute upper bounds of $n_\alpha$, $\delta_\alpha$ and $D_\alpha$.
"""

# ╔═╡ ab97adab-dd1d-4826-82a7-227bec87b2ed
md"""
### Bound $n_\alpha$
This corresponds to Lemma 12.5. We compute an enclosure of $n_\alpha$ and plot it together with $N_\alpha(x)$ for $x \in [0, \pi]$.
"""

# ╔═╡ 36d31c09-9882-4add-9a51-2c1d43c6d64a
@time n0_time = @elapsed n0 = n0_bound(u0, verbose = true)

# ╔═╡ f8afd857-4153-4c52-aeeb-1fff14fad24f
n0_xs, n0_ys = let xs = range(Arb(0), π, length = 100)
    N = x -> u0.w(x) / 2u0(x)
    ys = Folds.map(N, xs)
    ys[1] = 0 # It is zero at x = 0
    xs, ys
end

# ╔═╡ 35cc4d3c-2ec8-4e08-b721-cca697141276
let pl = plot(legend = :bottomright)
    plot!(
        pl,
        n0_xs,
        n0_ys,
        ribbon = radius.(Arb, n0_ys),
        xlabel = L"x",
        ylabel = L"N_\alpha(x)",
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!(pl, [n0], ribbon = [radius(Arb, n0)], color = :green, label = L"n_0" * " bound")
    pl
end

# ╔═╡ 179dee23-66a9-4cd6-a08a-c4400bd15440
md"""
### Bound $\delta_\alpha$
This corresonds to Lemma 12.6. We compute an enclosure of $\delta_\alpha$ and plot it together with $F_\alpha(x)$ for $x \in [0, \pi]$.
"""

# ╔═╡ 141836ea-b6cd-4696-a86d-7dc7d1964131
@time δ0_time = @elapsed δ0 = delta0_bound(u0, verbose = true)

# ╔═╡ e2afde10-d8d6-4bcf-9478-e9d6c6d07150
δ0_xs, δ0_ys = let xs = range(Arb(0.1), Arb(π), length = 200)
    ys = Folds.map(F0(u0), xs)
    xs, ys
end

# ╔═╡ c6133c35-b4fe-49f3-b226-2d49d092484d
δ0_asym_xs, δ0_asym_ys =
    let xs = exp.(range(log(Arb("1e-5")), log(Arb("1e-1")), length = 200))
        ys = Folds.map(F0(u0, Asymptotic()), xs)
        xs, ys
    end

# ╔═╡ efd0a5b9-c585-4fda-94e1-00247ea9b97a
md"""
In the end we want to have the inequality $\delta_\alpha \leq \frac{(1 - D_\alpha)^2}{4n_\alpha}$. In the next section we compute estimates of $D_\alpha$ and we can from that get an estimate for the right hand side this is the goal that we want the defect to be smaller than. In general the bound is much smaller than the goal, making it hard to distinguish the different scales. By default we therefore don't include the goal.
- Include goal $(@bind include_δ0_goal CheckBox(default = false))
"""

# ╔═╡ 48586f5e-77df-4a10-8610-dd24a986eca4
md"""
### Bound $D_\alpha$
This corresonds to Lemma 12.7. In this case we don't compute an enclosure of $D_\alpha$, but only prove that it is bounded by

$$1 - 2\sqrt{n_\alpha \delta_\alpha}$$

We then plot this bound together with $\mathcal{T}_\alpha(x)$.
"""

# ╔═╡ 9f7fbf89-e36b-412d-9bfc-9e51537747ec
D0_bound = 1 - 2sqrt(n0 * δ0)

# ╔═╡ 2d5cb115-da9d-40d3-af1c-6cd2333a808d
@time D0_time =
    @elapsed D0_bound_holds = D0_bounded_by(u0, lbound(D0_bound), verbose = true)

# ╔═╡ 8ec4333e-7eee-4dc3-a5eb-1154866bbaf4
if D0_bound_holds
    @info "The bound holds! So the traveling wave exists! 🎉🎉🎉"
else
    @error "The bound could not be proved to hold 😦"
end

# ╔═╡ 87a6990e-7be5-489e-812a-ca5034061784
D0_xs, D0_ys = let xs = collect(range(Arb(0), π, length = 100)[2:end])
    ys = Folds.map(xs) do x
        if x < 0.1
            T0(u0, Asymptotic())(x)
        else
            T0(u0, Ball())(x)
        end
    end
    pushfirst!(xs, 0)
    pushfirst!(ys, T0(u0, Asymptotic())(Arb(0)))
    xs, ys
end

# ╔═╡ 58f1233b-30d1-4f8c-abf4-078df6ccbee4
let pl = plot(legend = :bottomright)
    plot!(
        pl,
        D0_xs,
        D0_ys,
        ribbon = radius.(Arb, D0_ys),
        xlabel = L"x",
        ylabel = L"\mathcal{T}_\alpha(x)",
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!(
        pl,
        [D0_bound],
        ribbon = [radius(Arb, D0_bound)],
        color = :green,
        label = L"D_\alpha" * " bound",
    )
    pl
end

# ╔═╡ 16c74532-55ab-4c2f-af5a-2bab8b5f9668
md"""
From the values used in the plot we can compute an estimate of $D_\alpha$. This estimate allows us to compute an estimated goal for $\delta_\alpha$. This estimated goal is used in the figures for $F_\alpha(x)$ above.
"""

# ╔═╡ 8af1970a-c99c-44d0-b09d-3031cd520ca5
D0_estimate = maximum(D0_ys)

# ╔═╡ 123a5e6f-cbc4-46f8-958a-1b4ae9989186
δ0_goal_estimate = (1 - D0_estimate)^2 / 4n0

# ╔═╡ 60f43fcb-c8fb-4bbc-a6ac-50658d8947a3
let pl = plot()
    plot!(
        pl,
        δ0_xs,
        δ0_ys,
        ribbon = Arblib.radius.(Arb, δ0_ys),
        xlabel = L"x",
        ylabel = L"F_\alpha(x)",
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = L"\delta_0" * " bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    if include_δ0_goal
        hline!(
            [δ0_goal_estimate],
            ribbon = [radius(Arb, δ0_goal_estimate)],
            color = :red,
            label = L"\delta_0" * " goal",
        )
        hline!(
            [-δ0_goal_estimate],
            ribbon = [radius(Arb, δ0_goal_estimate)],
            color = :red,
            label = "",
        )
    end
    pl
end

# ╔═╡ 84eb673b-5d56-4d53-971d-60d385b45894
let pl = plot(xaxis = :log10, legend = :bottomleft)
    plot!(
        pl,
        δ0_asym_xs,
        δ0_asym_ys,
        ribbon = radius.(Arb, δ0_asym_ys),
        xlabel = L"x",
        ylabel = L"F_\alpha(x)",
        label = "",
        m = :circle,
        ms = 1,
    )
    hline!([δ0], ribbon = [radius(Arb, δ0)], color = :green, label = L"\delta_0" * " bound")
    hline!([-δ0], ribbon = [radius(Arb, δ0)], color = :green, label = "")
    if include_δ0_goal
        hline!(
            [δ0_goal_estimate],
            ribbon = [radius(Arb, δ0_goal_estimate)],
            color = :red,
            label = L"\delta_0" * " goal",
        )
        hline!(
            [-δ0_goal_estimate],
            ribbon = [radius(Arb, δ0_goal_estimate)],
            color = :red,
            label = "",
        )
    end
    pl
end

# ╔═╡ Cell order:
# ╠═906f72a8-f912-11eb-15d0-ef321ee21aa5
# ╟─bf052495-c3a5-4c75-9fb6-24e113c9aedf
# ╟─4854d9df-785c-484d-96c5-ca5942a01a69
# ╠═b376a940-bb33-4e6d-bd26-5406ecee51b8
# ╟─3e54d065-f8d7-4660-b2c2-4b53e5ade4f5
# ╟─e5b763ed-ab57-4324-848a-cbfb84000ecb
# ╟─3c366570-463c-4481-b4fd-5cab7da1a007
# ╟─f59706aa-f167-43f7-b154-2f8fb5f4eb77
# ╟─cef61488-2c3c-43a5-bcc3-c33fee24419c
# ╟─ab97adab-dd1d-4826-82a7-227bec87b2ed
# ╠═36d31c09-9882-4add-9a51-2c1d43c6d64a
# ╟─f8afd857-4153-4c52-aeeb-1fff14fad24f
# ╟─35cc4d3c-2ec8-4e08-b721-cca697141276
# ╟─179dee23-66a9-4cd6-a08a-c4400bd15440
# ╠═141836ea-b6cd-4696-a86d-7dc7d1964131
# ╟─e2afde10-d8d6-4bcf-9478-e9d6c6d07150
# ╟─c6133c35-b4fe-49f3-b226-2d49d092484d
# ╟─efd0a5b9-c585-4fda-94e1-00247ea9b97a
# ╟─60f43fcb-c8fb-4bbc-a6ac-50658d8947a3
# ╟─84eb673b-5d56-4d53-971d-60d385b45894
# ╟─48586f5e-77df-4a10-8610-dd24a986eca4
# ╠═9f7fbf89-e36b-412d-9bfc-9e51537747ec
# ╠═2d5cb115-da9d-40d3-af1c-6cd2333a808d
# ╟─8ec4333e-7eee-4dc3-a5eb-1154866bbaf4
# ╟─87a6990e-7be5-489e-812a-ca5034061784
# ╟─58f1233b-30d1-4f8c-abf4-078df6ccbee4
# ╟─16c74532-55ab-4c2f-af5a-2bab8b5f9668
# ╠═8af1970a-c99c-44d0-b09d-3031cd520ca5
# ╠═123a5e6f-cbc4-46f8-958a-1b4ae9989186
