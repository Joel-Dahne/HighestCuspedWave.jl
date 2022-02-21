### A Pluto.jl notebook ###
# v0.17.5

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

# ╔═╡ e27b528e-7326-4265-badf-e7ca97b77f0c
interval_parameters = [
    (Arb((-1 // 6, -0.1)), 1000),
    (Arb((-0.1, -0.05)), 1000),
    (Arb((-0.05, -0.025)), 1000),
    (Arb((-0.025, -0.0125)), 2000),
    (Arb((-0.0125, -0.00625)), 4000),
    (Arb((-0.00625, -0.003125)), 8000),
    (Arb((-0.003125, -0.0015625)), 16000),
    (Arb((-0.0015625, -0.0012)), 16000),
]

# ╔═╡ d3bf0151-cdbb-4f93-a989-4e569b930146
αs = let
	# Pick interval to use
    interval, n = interval_parameters[1]

	# Determine how many of its subintervals to compute for
    use_fraction = false
    fraction = 0.01
    number = 24

    failsafe = 100 # To avoid getting too many intervals by mistake

    αs = HighestCuspedWave.mince(interval, n)
    if use_fraction
        m = round(Int, n * fraction)
    else
        m = number
    end

    m = min(m, n, failsafe)
	
	indices = round.(Int, range(1, n, m))
    αs[indices]
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
        only_estimate_CB = false,
        threaded = false,
        verbose = false,
        extra_verbose = false,
    )

    (α = α, proof_data..., prec = precision(α))
end

# ╔═╡ 6d0c7fe1-003c-448c-9e74-098a9298cc14
md"""
# Process data from proofs
"""

# ╔═╡ f25e728d-20c0-4abd-9ed0-81f9c54b80e7
df = DataFrame(data)

# ╔═╡ 78391786-747a-4bc7-ba97-eac74a1f5c60
HighestCuspedWave.format_for_publishing.(df.α₀, df.δ₀, df.C_B)

# ╔═╡ 3c2c1f6a-bde5-4e22-81f3-c0c1e5232d6a
plot(
    df.α,
    [df.C_B_estimate (1 .- 2sqrt.(df.α₀ .* df.δ₀))],
    ribbon = radius.(Arb, [df.C_B_estimate (1 .- 2sqrt.(df.α₀ .* df.δ₀))]),
    m = :circle,
    labels = ["C_B_estimate" "Required bound"],
    xaxis = "α",
    legend = :bottomright,
)

# ╔═╡ 674db1b6-092c-42fb-a42d-d1f84f608382
plot(
	df.α, 
	df.p, 
	m = :circle,
	ylims = (-0.1, 1.1),
	xaxis = "α",
	yaxis = "p",
	legend = :none,
)

# ╔═╡ 2b14fd7b-d55f-4480-b9dc-47f3d6dcbebf
plot(
	df.α, 
	df.proved,
	m = :circle,
	ylims = (-0.1, 1.1),
	xaxis = "α",
	title = "Proved",
	legend = :none,
)

# ╔═╡ 99e8152f-b87f-4060-ad44-083bba641fcc
plot(
    df.α,
    [df.α₀_time df.δ₀_time df.C_B_time],
    m = :circle,
    labels = ["α₀" "δ₀" "C_B"],
    title = "Time to compute bound",
    xaxis = "α",
    yaxis = "Seconds",
)

# ╔═╡ Cell order:
# ╠═864c4f00-8bfc-11ec-093c-efa21e3a2e46
# ╟─2615beec-85f6-4d1c-9f75-4a777ac6b025
# ╠═e27b528e-7326-4265-badf-e7ca97b77f0c
# ╠═d3bf0151-cdbb-4f93-a989-4e569b930146
# ╟─a8d1346e-dee6-451c-8e26-3bcc1c6432f0
# ╠═3f251c0a-b873-48c8-8ec8-188438960ebe
# ╠═ca05c284-dde3-4275-9a9a-301d75d9a3c2
# ╟─6d0c7fe1-003c-448c-9e74-098a9298cc14
# ╠═f25e728d-20c0-4abd-9ed0-81f9c54b80e7
# ╠═78391786-747a-4bc7-ba97-eac74a1f5c60
# ╟─3c2c1f6a-bde5-4e22-81f3-c0c1e5232d6a
# ╟─674db1b6-092c-42fb-a42d-d1f84f608382
# ╟─2b14fd7b-d55f-4480-b9dc-47f3d6dcbebf
# ╟─99e8152f-b87f-4060-ad44-083bba641fcc
