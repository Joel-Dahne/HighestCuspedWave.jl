### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 864c4f00-8bfc-11ec-093c-efa21e3a2e46
begin
    using Pkg, Revise
    Pkg.activate("../")
    using Arblib, DataFrames, Folds, HighestCuspedWave, Plots, PlutoUI
    setprecision(Arb, 128)
    nothing
end

# ╔═╡ 2615beec-85f6-4d1c-9f75-4a777ac6b025
md"""
Set up the values of $α$ to consider.

The width of the ball for $α$ that we can use depends on the value of $α$. We here give a list of subintervals of $(-1, 0)$ and the number of parts each subinterval needs to be split into.
"""

# ╔═╡ 13b37502-1104-4ed4-a6f3-715afabd509e
md"""
We set all parameters used in one go so that we can easily change several at once without having to recompute several times.
"""

# ╔═╡ a3596ccc-e650-4fba-8993-f9dc5215101a
md"""
Intervals to look closer at:
"""

# ╔═╡ e27b528e-7326-4265-badf-e7ca97b77f0c
begin
    interval_parameters = [
        (Arb((-0.95, -0.9)), 8000), # Check
        (Arb((-0.9, -0.85)), 8000), # Check
        (Arb((-0.85, -0.8)), 4000), # Check
        (Arb((-0.8, -0.7)), 2000), # Good for estimation
        (Arb((-0.7, -0.6)), 1000), # Good for estimation
        (Arb((-0.6, -0.5)), 1500), # Good for estimation
        (Arb((-0.5, -0.41)), 1500), # Good for estimation
        (Arb((-0.41, -0.33)), 2000), # Good for estimation
        (Arb((-0.33, -0.16)), 25000),
        (Arb((-0.16, -0.1)), 1000),
        (Arb((-0.1, -0.05)), 1000),
        (Arb((-0.05, -0.025)), 1000),
        (Arb((-0.025, -0.0125)), 2000),
        (Arb((-0.0125, -0.00625)), 4000),
        (Arb((-0.00625, -0.003125)), 8000),
        (Arb((-0.003125, -0.0015625)), 16000),
        (Arb((-0.0015625, -0.0012)), 16000),
    ]

    # Index for interval to use
    interval_index = 4

    # Determine how many of its subintervals to compute for
    use_fraction = false
    fraction = 0.01
    number = 100

    failsafe = 2000 # To avoid getting too many intervals by mistake

    thin = false # If true then set radius of αs to zero

    # If this is true it only estimates the norm
    only_estimate_D0 = true
end

# ╔═╡ d3bf0151-cdbb-4f93-a989-4e569b930146
αs = let
    # Pick interval to use
    interval, n = interval_parameters[interval_index]

    αs = HighestCuspedWave.mince(interval, n)

    if use_fraction
        m = round(Int, n * fraction)
    else
        m = number
    end
    m = min(m, n, failsafe)

    indices = round.(Int, range(1, n, m))

    # The left endpoint of the interval sometimes picks up the parameters for the interval to the left. To avoid this discontinuity we prioritise using the subinterval with index 2 instead of with index 1 if it is not already used.
    if !(2 ∈ indices)
        indices[1] = 2
    end

    αs = αs[indices]

    if thin
        αs = midpoint.(Arb, αs)
    end

    αs
end

# ╔═╡ a8d1346e-dee6-451c-8e26-3bcc1c6432f0
md"""
Loop over all values of $α$ and prove the theorem for each one, collecting data about the proof.
"""

# ╔═╡ 3f251c0a-b873-48c8-8ec8-188438960ebe
md"""
**TODO:**
- Progress logging
"""

# ╔═╡ ca05c284-dde3-4275-9a9a-301d75d9a3c2
data = Folds.map(αs, basesize = 1) do α
    proof_data = HighestCuspedWave.prove(
        α,
        M = 10,
        threaded = false,
        verbose = false,
        extra_verbose = false;
        only_estimate_D0,
    )

    (α = α, proof_data...)
end

# ╔═╡ 6d0c7fe1-003c-448c-9e74-098a9298cc14
md"""
# Process data from proofs
"""

# ╔═╡ f25e728d-20c0-4abd-9ed0-81f9c54b80e7
df = DataFrame(data)

# ╔═╡ b2bf455c-9bf5-40fe-8764-d346e9bd8289
df_failed = filter(:proved_estimate => !identity, df)

# ╔═╡ 7c08744f-cdef-4232-a261-87cf98dafa0a
Arblib.dump_string.(df_failed.α)

# ╔═╡ 78391786-747a-4bc7-ba97-eac74a1f5c60
HighestCuspedWave.round_for_publishing.(df.n₀, df.δ₀, df.D₀)

# ╔═╡ e97195d0-ed8d-46bd-a1fb-da9c6d92cd94
let bound = (1 .- df.D₀_estimate) .^ 2 ./ 4df.n₀
    plot(
        df.α,
        [df.δ₀ bound],
        ribbon = radius.(Arb, [df.δ₀ bound]),
        m = :circle,
        labels = ["δ₀" "Required bound"],
        xaxis = "α",
        legend = :bottomright,
    )
end

# ╔═╡ 7a2ce48c-2e95-40ba-980b-5afc0923b2cf
# ╠═╡ disabled = true
#=╠═╡
let bound = (1 .- df.D₀_estimate) .^ 2 ./ 4df.n₀
    plot(
        df.α,
        [df.δ₀ bound],
        ribbon = radius.(Arb, [df.δ₀ bound]),
        m = :circle,
        labels = ["δ₀" "Required bound"],
        xaxis = "α",
        legend = :bottomright,
    )
end
  ╠═╡ =#

# ╔═╡ 3c2c1f6a-bde5-4e22-81f3-c0c1e5232d6a
let bound = 1 .- 2Arblib.sqrtpos.(df.n₀ .* df.δ₀)
    plot(
        df.α,
        [df.D₀_estimate bound],
        ribbon = radius.(Arb, [df.D₀_estimate bound]),
        m = :circle,
        labels = ["D₀_estimate" "Required bound"],
        xaxis = "α",
        legend = :bottomright,
    )
end

# ╔═╡ 65e473ad-e479-4c7e-879d-5f455f21f1ec
# ╠═╡ disabled = true
#=╠═╡
let bound = 1 .- 2Arblib.sqrtpos.(df.n₀ .* df.δ₀)
    plot(
        df.α,
        [df.D₀_estimate bound],
        ribbon = radius.(Arb, [df.D₀_estimate bound]),
        m = :circle,
        labels = ["D₀_estimate" "Required bound"],
        xaxis = "α",
        legend = :bottomright,
    )
end
  ╠═╡ =#

# ╔═╡ 2b14fd7b-d55f-4480-b9dc-47f3d6dcbebf
plot(
    df.α,
    [df.proved df.proved_estimate],
    m = :circle,
    ylims = (-0.1, 1.1),
    xaxis = "α",
    labels = ["Proved" "Estimate holds"],
)

# ╔═╡ 674db1b6-092c-42fb-a42d-d1f84f608382
plot(df.α, df.p, m = :circle, ylims = (-0.1, 1.1), xaxis = "α", yaxis = "p", legend = :none)

# ╔═╡ 0ea0346e-e6b0-4516-b07f-a9f21828ff2c
total_time = if only_estimate_D0
    df.n₀_time + df.δ₀_time + df.D₀_estimate_time + df.u0_time
else
    df.n₀_time + df.δ₀_time + df.D₀_estimate_time + df.D₀_time + df.u0_time
end

# ╔═╡ 99e8152f-b87f-4060-ad44-083bba641fcc
plot(
    df.α,
    [df.n₀_time df.δ₀_time df.D₀_estimate_time df.D₀_time df.u0_time total_time],
    m = :circle,
    ms = 2,
    labels = ["α₀" "δ₀" "D₀ est" "D₀" "u0" "total"],
    title = "Time to compute bound",
    xaxis = "α",
    yaxis = "Seconds",
    legend = :outerright,
)

# ╔═╡ 25340bd8-75b8-4dc4-b59b-ecf733ad65a9
mean_time = sum(total_time) / length(total_time)

# ╔═╡ 9eca2429-1d03-46fa-9610-9193845868f6
md"""
Estimated time to compute for the full interval, without parallelisation.
"""

# ╔═╡ 8796c4f5-269c-4852-b68b-5c4e1388405c
time_estimate_hours = mean_time * interval_parameters[interval_index][2] / 3600

# ╔═╡ 0aae8959-60a5-41cc-ba1a-499dd5ac37e0
time_estimate_days = time_estimate_hours / 24

# ╔═╡ Cell order:
# ╠═864c4f00-8bfc-11ec-093c-efa21e3a2e46
# ╟─2615beec-85f6-4d1c-9f75-4a777ac6b025
# ╟─13b37502-1104-4ed4-a6f3-715afabd509e
# ╟─a3596ccc-e650-4fba-8993-f9dc5215101a
# ╠═e27b528e-7326-4265-badf-e7ca97b77f0c
# ╠═d3bf0151-cdbb-4f93-a989-4e569b930146
# ╟─a8d1346e-dee6-451c-8e26-3bcc1c6432f0
# ╟─3f251c0a-b873-48c8-8ec8-188438960ebe
# ╠═ca05c284-dde3-4275-9a9a-301d75d9a3c2
# ╟─6d0c7fe1-003c-448c-9e74-098a9298cc14
# ╠═f25e728d-20c0-4abd-9ed0-81f9c54b80e7
# ╠═b2bf455c-9bf5-40fe-8764-d346e9bd8289
# ╠═7c08744f-cdef-4232-a261-87cf98dafa0a
# ╠═78391786-747a-4bc7-ba97-eac74a1f5c60
# ╠═e97195d0-ed8d-46bd-a1fb-da9c6d92cd94
# ╠═7a2ce48c-2e95-40ba-980b-5afc0923b2cf
# ╟─3c2c1f6a-bde5-4e22-81f3-c0c1e5232d6a
# ╠═65e473ad-e479-4c7e-879d-5f455f21f1ec
# ╠═2b14fd7b-d55f-4480-b9dc-47f3d6dcbebf
# ╟─674db1b6-092c-42fb-a42d-d1f84f608382
# ╠═0ea0346e-e6b0-4516-b07f-a9f21828ff2c
# ╠═99e8152f-b87f-4060-ad44-083bba641fcc
# ╠═25340bd8-75b8-4dc4-b59b-ecf733ad65a9
# ╟─9eca2429-1d03-46fa-9610-9193845868f6
# ╠═8796c4f5-269c-4852-b68b-5c4e1388405c
# ╠═0aae8959-60a5-41cc-ba1a-499dd5ac37e0
