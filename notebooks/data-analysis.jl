### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ c5afc656-9737-11ed-2d9d-95fffa9fe4ef
begin
    using Pkg, Revise
    Pkg.activate("../")
    using Arblib,
        DataFrames, Folds, HighestCuspedWave, LaTeXStrings, Statistics, Plots, PlutoUI

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ 6de44628-7021-4b4b-b8d1-6514ab6c0dd6
md"""
## Load the data
We start by reading in the data from the computations on Dardel.
"""

# ╔═╡ ad318833-410f-436c-ac72-2acb02d95c88
directories = [
    "../PDC/data/proof/proof-2023-01-09T14:54:40",
    "../PDC/data/proof/proof-2023-01-10T09:26:04",
    "../PDC/data/proof/proof-2023-01-11T02:43:33",
    "../PDC/data/proof/proof-2023-01-10T09:27:26",
    "../PDC/data/proof/proof-2023-01-16T10:59:09",
    "../PDC/data/proof/proof-2023-01-16T11:28:54",
]

# ╔═╡ 05d571c5-82b5-49ab-8514-a611ab31aa19
data_dardel = HighestCuspedWave.read_proof_data_dir(directories)

# ╔═╡ 86807f3f-d55b-4d18-9b78-09e6a618f066
HighestCuspedWave.check_proof_data(data_dardel) # Run some sanity checks on the data

# ╔═╡ 3c022eff-6c19-4947-b113-adaec80f2d9c
md"""
For most of the subintervals the proof should have succeeded. However there is typically a few for which it didn't succeed.
"""

# ╔═╡ 0914811f-0b45-4347-97c6-652aa31de73f
filter(:proved => !, data_dardel)

# ╔═╡ 06c6637c-ead0-4f96-bf6d-90268fbf7ce4
md"""
The subintervals that failed are bisected and the calculation is run on each part. This can be done using the following function. However, to save on computation time this has already been precompued and we just load the data.
"""

# ╔═╡ f906e0ab-0b2d-4b71-b2b2-b6e3e4faa78e
function rerun_failed(data::DataFrame)
    failed = filter(:proved => !, data)
    αs_bisected = collect(
        Iterators.flatten([
            (Arb((lbound(α), midpoint(α))), Arb((midpoint(α), ubound(α)))) for α in failed.α
        ]),
    )

    HighestCuspedWave.prove(αs_bisected)
end

# ╔═╡ 4390863a-9692-4983-8e3f-a9638340c6ea
data_rerun = HighestCuspedWave.read_proof_data("../data/proof-rerun-2023-01-18.csv")

# ╔═╡ 5d937e1e-c4a5-4c9a-b602-ad8a767ed1e2
md"""
Finally we take the data to be the one from Dardel that succeeded plus the bisected ones.
"""

# ╔═╡ 525e5539-82f1-4d37-809f-338bf32e8a5a
data = sort!(
    vcat(filter(:proved => identity, data_dardel), data_rerun),
    :α,
    by = α -> midpoint(α),
)

# ╔═╡ 376f8999-6033-46f6-9b6b-e1eeb139a1ec
@assert all(data.proved)

# ╔═╡ 1dc6509c-2217-43b6-a467-ccbf9fade853
@assert all(data.proved)

# ╔═╡ ff4b9389-8fa3-4bc9-b6fb-23504a58a576
md"""
We double check that it actually covers the interval claimed in the paper.
"""

# ╔═╡ 097efc47-79f9-4412-a059-9f7bd6be967b
begin
    # The first and last values are as claimed
    @assert contains(data.α[1], Arb("-0.9999"))
    @assert contains(data.α[end], Arb("-0.0012"))

    for i = 1:length(data.α)-1
        @assert Arblib.overlaps(data.α[i], data.α[i+1])
    end
end

# ╔═╡ 4a36cbeb-12af-47c2-b9af-b1300b06d8a6
md"""
## Analyse the data
"""

# ╔═╡ 9c0b90d8-e6d8-4374-bf89-de304b344205
md"""
### Plots
"""

# ╔═╡ 237ed8a1-c221-41de-bd33-3bdc53aa80fa
save = false

# ╔═╡ 0dd1a8fe-34ac-416d-b2e4-3d1eaf801318
save_dir = "../figures/publication/"

# ╔═╡ a42d3a45-551b-4fff-a934-8be95e60672d
let pl = plot(xlabel = L"\alpha", ylabel = L"p", legend = :none)
    scatter!(pl, data.α, data.p, msw = 0, ms = 1)
    save && savefig(pl, joinpath(save_dir, "interval-p.png"))
    pl
end

# ╔═╡ f68bebf9-bc64-4f2a-ad23-6795054cecbd
let pl = plot(xlabel = L"\alpha", ylabel = L"n_\alpha", ylims = (0, NaN))
    scatter!(pl, data.α, data.n₀_rounded, msw = 0, ms = 1, label = "Bound")
    save && savefig(pl, joinpath(save_dir, "interval-n.png"))
    pl
end

# ╔═╡ e7a9ab6b-a529-4327-bf40-165aaf52854b
let pl = plot(xlabel = L"\alpha", ylabel = L"\delta_\alpha")
    scatter!(pl, data.α, data.δ₀_rounded, msw = 0, ms = 1, label = "Bound")
    goal_estimate = @. (1 - data.D₀_estimate)^2 / 4data.n₀_rounded
    scatter!(pl, data.α, goal_estimate, msw = 0, ms = 1, label = "Estimated goal")
    save && savefig(pl, joinpath(save_dir, "interval-delta.png"))
    pl
end

# ╔═╡ e089b20a-03b2-442b-9806-ad53f277c8aa
let pl = plot(xlabel = L"\alpha", ylabel = L"D_\alpha")
    scatter!(pl, data.α, data.D₀_rounded, msw = 0, ms = 1, label = "Bound")
    scatter!(pl, data.α, data.D₀_estimate, msw = 0, ms = 1, label = "Estimate")
    save && savefig(pl, joinpath(save_dir, "interval-D.png"))
    pl
end

# ╔═╡ be553dc7-1ddc-4239-8dea-0f857ffa5e28
let pl = plot(xlabel = L"\alpha")
    scatter!(pl, data.α, data.u0_N0, msw = 0, ms = 1, label = L"N_{\alpha,0}")
    scatter!(pl, data.α, data.u0_N1, msw = 0, ms = 1, label = L"N_{\alpha,1}")
    save && savefig(pl, joinpath(save_dir, "interval-N0-N1.png"))
    pl
end

# ╔═╡ 8e6f4347-e7ad-448e-8109-4eb7a29ee434
md"""
#### Asymptotic plots near $\alpha = -1$
"""

# ╔═╡ 00feca50-1d56-48ba-8104-fd477f5a5464
data_near_m1 = filter(:α => <(-0.9), data)

# ╔═╡ efe296dd-25ad-4ffd-aad7-2f1abf3112a5
let pl = plot(
        xlabel = L"1 + \alpha",
        ylabel = L"n_\alpha",
        legend = :none,
        ylims = (0, NaN),
        xaxis = :log10,
    )
    scatter!(pl, 1 .+ data_near_m1.α, data_near_m1.n₀_rounded, msw = 0, ms = 1)
    pl
end

# ╔═╡ d07b4c54-c6ec-4b59-b577-9c8b218b3e7b
let pl = plot(
        xlabel = L"1 + \alpha",
        ylabel = L"\delta_\alpha",
        legend = :none,
        xaxis = :log10,
        ylims = (0, NaN),
    )
    scatter!(pl, 1 .+ data_near_m1.α, data_near_m1.δ₀_rounded, msw = 0, ms = 1)
    goal_estimate = @. (1 - data_near_m1.D₀_estimate)^2 / 4data_near_m1.n₀_rounded
    scatter!(pl, 1 .+ data_near_m1.α, goal_estimate, msw = 0, ms = 1)
    pl
end

# ╔═╡ 72b00cfe-7f2a-4043-adbf-9a393aea782a
let pl = plot(
        xlabel = L"1 + \alpha",
        ylabel = L"D_\alpha",
        legend = :none,
        xaxis = :log10,
        ylims = (NaN, 1),
    )
    scatter!(pl, 1 .+ data_near_m1.α, data_near_m1.D₀_rounded, msw = 0, ms = 1)
    scatter!(pl, 1 .+ data_near_m1.α, data_near_m1.D₀_estimate, msw = 0, ms = 1)
    pl
end

# ╔═╡ e39d9c7e-8d49-457e-be7b-d2f9d6d070be
md"""
### Time consumption
"""

# ╔═╡ e89cc470-90da-4576-b2ef-39e2c2f36249
let pl = plot(xlabel = L"\alpha", ylabel = "seconds", legend = :none)
    scatter!(
        pl,
        1 .+ data.α,
        [data.n₀_time data.δ₀_time data.D₀_time],
        msw = 0,
        ms = 1,
        xaxis = :log10,
    )
    pl
end

# ╔═╡ 829e9bf7-0455-4325-a1fc-1d0e1f41ef33
extrema(data.n₀_time), mean(data.n₀_time)

# ╔═╡ efe8fff0-be11-4259-aace-fe75c4940df9
extrema(data.δ₀_time), mean(data.δ₀_time)

# ╔═╡ cd119999-843c-4a06-b7ee-2efe080339c1
extrema(data.D₀_time), (mean(data.D₀_time))

# ╔═╡ ec55a995-a331-4ec4-aa2f-5ecd3839a877
md"""
Core hours used for different parts of the computations. Based on the assumption that all computations were run with 4 threads.
"""

# ╔═╡ 123fa81c-f257-4bfb-9434-e8f8c2f4c3aa
total_core_hours_dardel =
    sum(
        data_dardel.u0_time +
        data_dardel.n₀_time +
        data_dardel.δ₀_time +
        data_dardel.D₀_estimate_time +
        map(x -> isfinite(x) ? x : zero(x), data_dardel.D₀_time),
    ) / 3600 * 4

# ╔═╡ ca68456d-8aab-4fb0-acb8-6363456ada1c
u0_core_hours_dardel = sum(data_dardel.u0_time) / 3600 * 4

# ╔═╡ aa8a7689-89e4-4dc3-af81-d8a60dfe630b
n₀_core_hours_dardel = sum(data_dardel.n₀_time) / 3600 * 4

# ╔═╡ 5faca96a-6094-4911-99bb-08efc6563b05
δ₀_core_hours_dardel = sum(data_dardel.δ₀_time) / 3600 * 4

# ╔═╡ 33121b6b-c9a6-4902-bdd1-970b42c89968
D₀_estimate_core_hours_dardel = sum(data_dardel.D₀_estimate_time) / 3600 * 4

# ╔═╡ f187c0dd-19e9-435e-acb5-05cc23956412
D₀_core_hours_dardel = sum(filter(isfinite, data_dardel.D₀_time)) / 3600 * 4

# ╔═╡ Cell order:
# ╠═c5afc656-9737-11ed-2d9d-95fffa9fe4ef
# ╟─6de44628-7021-4b4b-b8d1-6514ab6c0dd6
# ╠═ad318833-410f-436c-ac72-2acb02d95c88
# ╠═05d571c5-82b5-49ab-8514-a611ab31aa19
# ╠═86807f3f-d55b-4d18-9b78-09e6a618f066
# ╟─3c022eff-6c19-4947-b113-adaec80f2d9c
# ╠═0914811f-0b45-4347-97c6-652aa31de73f
# ╟─06c6637c-ead0-4f96-bf6d-90268fbf7ce4
# ╠═f906e0ab-0b2d-4b71-b2b2-b6e3e4faa78e
# ╠═4390863a-9692-4983-8e3f-a9638340c6ea
# ╠═376f8999-6033-46f6-9b6b-e1eeb139a1ec
# ╟─5d937e1e-c4a5-4c9a-b602-ad8a767ed1e2
# ╠═525e5539-82f1-4d37-809f-338bf32e8a5a
# ╠═1dc6509c-2217-43b6-a467-ccbf9fade853
# ╟─ff4b9389-8fa3-4bc9-b6fb-23504a58a576
# ╠═097efc47-79f9-4412-a059-9f7bd6be967b
# ╟─4a36cbeb-12af-47c2-b9af-b1300b06d8a6
# ╟─9c0b90d8-e6d8-4374-bf89-de304b344205
# ╠═237ed8a1-c221-41de-bd33-3bdc53aa80fa
# ╠═0dd1a8fe-34ac-416d-b2e4-3d1eaf801318
# ╠═a42d3a45-551b-4fff-a934-8be95e60672d
# ╠═f68bebf9-bc64-4f2a-ad23-6795054cecbd
# ╠═e7a9ab6b-a529-4327-bf40-165aaf52854b
# ╠═e089b20a-03b2-442b-9806-ad53f277c8aa
# ╠═be553dc7-1ddc-4239-8dea-0f857ffa5e28
# ╟─8e6f4347-e7ad-448e-8109-4eb7a29ee434
# ╠═00feca50-1d56-48ba-8104-fd477f5a5464
# ╠═efe296dd-25ad-4ffd-aad7-2f1abf3112a5
# ╠═d07b4c54-c6ec-4b59-b577-9c8b218b3e7b
# ╠═72b00cfe-7f2a-4043-adbf-9a393aea782a
# ╟─e39d9c7e-8d49-457e-be7b-d2f9d6d070be
# ╠═e89cc470-90da-4576-b2ef-39e2c2f36249
# ╠═829e9bf7-0455-4325-a1fc-1d0e1f41ef33
# ╠═efe8fff0-be11-4259-aace-fe75c4940df9
# ╠═cd119999-843c-4a06-b7ee-2efe080339c1
# ╟─ec55a995-a331-4ec4-aa2f-5ecd3839a877
# ╠═123fa81c-f257-4bfb-9434-e8f8c2f4c3aa
# ╠═ca68456d-8aab-4fb0-acb8-6363456ada1c
# ╠═aa8a7689-89e4-4dc3-af81-d8a60dfe630b
# ╠═5faca96a-6094-4911-99bb-08efc6563b05
# ╠═33121b6b-c9a6-4902-bdd1-970b42c89968
# ╠═f187c0dd-19e9-435e-acb5-05cc23956412
