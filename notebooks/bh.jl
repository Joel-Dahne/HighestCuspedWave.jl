### A Pluto.jl notebook ###
# v0.19.14

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

    setprecision(Arb, 100)

    nothing
end

# ╔═╡ 3426f2ac-f96f-11eb-22b0-2b3f9ccb38b9
md"""
# Burgers-Hilberts Equation
This notebook contains the computer assisted part of the proof for the existence of a $2\pi$-periodic highest cusped traveling wave for the Burgers-Hilbert equation. The Burgers-Hilbert equation is given by

$f_t + f f_x = H[f]$

where $H[f]$ is the Hilbert transform of $f$. For traveling waves the ansatz $f(x, t) = \varphi(x - ct)$ reduces it to the equation

$-c\varphi' + \varphi\varphi' = H[\varphi],$

where $c$ is the wave speed.
"""

# ╔═╡ a75ba6b9-fe0b-4d3a-a88b-cc8ef9492404
md"""
For details about the proof see the paper, we here give a very short overview of the reduction to a fixed point problem.

With the ansatz $\varphi(x) = c - u(x)$ and integrating we can further reduce it to the equation

$\frac{1}{2}u^2 = -\mathcal{H}[u]$

where $\mathcal{H}[u](x)$ is the integral of the Hilbert transform with the constant of integration taken such that $\mathcal{H}[u](0) = 0$. If we now make the ansatz

$u(x) = u_0(x) + w(x)v(x)$

with $u_0$ an approximate solution and $w$ a fixed weight function we can, with some work, reduce the problem of proving the existence of a solution to the equation to proving existence of a fixed point of the operator

$G[v] = (I - T)^{-1}(-F - Nv^2)$

Here

$N(x) = \frac{w(x)}{2u_0(x)}$

$F(x) = \frac{1}{w(x)u_0(x)}\left(\mathcal{H}[u_0](x) + \frac{1}{2}u_0(x)^2\right)$

$T[v](x) = -\frac{1}{w(x) u_0(x)}\mathcal{H}[wv](x)$

Let $n_0 = \|N\|_{L^\infty}$, $\delta_0 = \|F\|_{L^\infty}$ and $D_0 = \|T\|$. Proving the existence of a fixed point for $G$ reduces to checking the inequality

$\delta_0 < \frac{(1 - D_0)^2}{4n_0}.$

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
To prove the result we need to bound three different value $n_0$, $\delta_0$ and $D_0$, given by the supremum on the interval $[0, \pi]$ for the functions $N(x)$, $F(x)$ and $\mathcal{T}(x)$ respectively. Here $\mathcal{T}(x)$ denotes the function which supremum gives the norm of the operator $T$, see the paper for more details about it.
"""

# ╔═╡ f0baf2ec-3f73-4d55-9ce4-754d94d7f3ce
md"""
Computing rigorous bounds for $\delta_0$ and $D_0$ is fairly time consuming. For that reason the code can either compute rigorous error bounds or use estimates. The estimates are given by simply evaluating the corresponding functions on a number of points and taking the maximum. For the defect $\delta_0$ we also make sure to use points asymptically close to $0$ since that's where the largest defect is found. Check the constants to use rigorous error bounds for
-  $\delta_0$ $(@bind use_rigorous_bounds_δ0 CheckBox(default = !isdefined(Main, :PlutoRunner)))
-  $D_0$ $(@bind use_rigorous_bounds_D0 CheckBox(default = !isdefined(Main, :PlutoRunner)))
Notice that the rigorous error bounds take longer time to compute with. On an AMD Ryzen 9 5900X the computation of $\delta_0$ around 15 minutes and the computation of $D_0$ around 2 minutes when using 12 threads. Note that for $n_0$ we always compute rigorous bounds since it doesn't take much time.
"""

# ╔═╡ 4c4cbf2a-3aec-4257-9fac-d8a0418d12d7
md"""
### Bound $n_0$
"""

# ╔═╡ 9c42cbca-805d-4a76-9d7d-68d40129f10a
md"""
We start by computing an enclosure of $n_0$ and plot it together with $N_0(x)$ for $x \in [0, \pi]$.
"""

# ╔═╡ f1ca794b-e1fc-481b-9b75-1e42a0b48a58
@time n0_time = @elapsed n0 = n0_bound(u0, verbose = true)

# ╔═╡ f1dce520-a035-43e6-9e08-4696a14c5a54
n0_xs, n0_ys = let xs = range(Arb(0), π, length = 100)[2:end]
    N = x -> u0.w(x) / 2u0(x)
    ys = Folds.map(N, xs)
    xs, ys
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

# ╔═╡ 05740cb8-a7fd-4add-a4c7-95b617508f28
md"""
Next we compute bounds of $\delta_0$. The code is set to use **$(ifelse(use_rigorous_bounds_δ0, "rigorous bounds", "estimates"))** for $\delta_0$.
"""

# ╔═╡ 664df5e6-152a-4b9a-adf8-352bdef59055
δ0_bound, δ0_time, δ0_subintervals = if use_rigorous_bounds_δ0
    @time δ0_time = @elapsed δ0_bound, δ0_subintervals... =
        delta0_bound(u0, return_subresults = true, verbose = true)
    δ0_bound, δ0_time, δ0_subintervals
else
    missing, missing, (Arb(NaN), Arb(NaN), Arb(NaN), Arb(NaN), Arb(NaN))
end

# ╔═╡ 67aa36b0-b77c-4531-a248-f7d474ffd47d
md"""
We do three different plots. One non-asymptotic plot on $[0.1, \pi]$, one asymptotic plot on $[10^{-100}, 0.1]$ and when even more asymptotic plot on $[10^{-20000}, 10^{-100}]$.
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

# ╔═╡ 67786c2c-a101-44db-a1b6-f191d2703bb0
md"""
We want to have the inequality $\delta_0 \leq \frac{(1 - D_0)^2}{4n_0}$. We therefore include the value of the right hand side in this inequality in the plots, this is the goal that we want the defect to be smaller than.
"""

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

# ╔═╡ 77d49bd5-2a4e-468d-a19a-b2e16a132f78
md"""
Next we compute bounds of $D_0$. The code is set to use **$(ifelse(use_rigorous_bounds_D0, "rigorous bounds", "estimates"))** for $D_0$.
"""

# ╔═╡ de5bed6f-7079-40a9-a5eb-f315abc20ddf
D0_enclosure, D0_time = if use_rigorous_bounds_D0
    @time D0_time = @elapsed D0_enclosure = D0_bound(u0, verbose = true)
    D0_enclosure, D0_time
else
    missing, missing
end

# ╔═╡ de4546e1-4a9f-4d37-b59c-ee4509d09868
D0_xs, D0_ys = let xs = range(Arb(1e-3), π, length = 100)
    ys = Folds.map(T0(u0, Ball()), xs)
    xs, ys
end

# ╔═╡ b0577d0f-77ba-4035-9d3b-ae4d6e5c624f
D0 = if use_rigorous_bounds_D0
    D0_enclosure
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
    if use_rigorous_bounds_D0 && use_rigorous_bounds_δ0
        @info "The result holds! 🎉🎉🎉"
    else
        @info "The inequality holds but the bounds are not rigorous"
    end
else
    @info "Could not prove the result 😥"
end

# ╔═╡ f25bcedf-cf4f-45f6-929a-abba66e7730c
md"""
When giving the enclosures in the paper we want to make sure that the printed version is less than the given upper bound. Here we check that this is the case. Note that this is not important for the result to hold, just for the printed values to look nice.
"""

# ╔═╡ 8d00fdce-e433-464e-adcd-88c21552439d
n0_printed = string(n0)

# ╔═╡ ff9462b6-3ff1-4fd9-af14-eebce04f349c
Arb(n0_printed) < n0_rounded

# ╔═╡ a80ae53d-c3ec-4f51-9866-1b2f7a23c85f
δ0_printed = string(δ0)

# ╔═╡ c940e613-622f-43e0-939a-bf30c033b810
δ0_subintervals_printed = string.(δ0_subintervals)

# ╔═╡ 7fc04a7b-1082-4af7-ac87-426e73276c18
all([
    Arb(δ_printed) < δ0_rounded for (δ_printed, δ) in
    zip((δ0_printed, δ0_subintervals_printed...), (δ0, δ0_subintervals...))
])

# ╔═╡ 74392233-bb8a-44d4-8f7f-ff8b5090142e
D0_printed = ArbExtras.format_interval(getinterval(D0)...)

# ╔═╡ 2494dd2b-088d-402f-8c28-29b3614a5fc6
md"""
The plots in the paper were produced with `pgfplotsx`, this is not a dependency so is not enabled by default. The plots look a bit weird with the default backend.
"""

# ╔═╡ 34858269-7714-4210-a678-2b91357c3eb7
#pgfplotsx()

# ╔═╡ f0110b12-a782-4840-abbf-ed623175c276
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"N(x)",
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
    savefig(pl, "../figures/publication/BH-N.pdf")
    pl
end

# ╔═╡ 16b4c640-6b4c-46a3-a505-e0eb7bbb2632
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"\mathcal{T}(x)",
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
    savefig(pl, "../figures/publication/BH-T.pdf")
    pl
end

# ╔═╡ fe196d34-c822-4da0-8ed1-955d2fb6ffa6
δ0_goal_rounded = (1 - D0_rounded)^2 / 4n0_rounded

# ╔═╡ 4654b37d-fcb5-43ab-8b67-5f946637943d
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"F(x)",
        guidefontsize = 24,
        tickfontsize = 24,
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
    hline!(Float64[-δ0_goal_rounded, δ0_goal_rounded], color = :red, linestyle = :dot)
    savefig(pl, "../figures/publication/BH-F.pdf")
    pl
end

# ╔═╡ 1c2d23d0-9e30-4830-a923-3ce15e8deb04
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"F(x)",
        xaxis = :log10,
        guidefontsize = 24,
        tickfontsize = 24,
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
    hline!(Float64[-δ0_goal_rounded, δ0_goal_rounded], color = :red, linestyle = :dot)
    savefig(pl, "../figures/publication/BH-F-asymptotic-1.pdf")
    pl
end

# ╔═╡ 7e1e6379-4cde-4a6f-9390-41618a989490
let pl = plot(
        legend = :none,
        xlabel = L"x",
        ylabel = L"F(x)",
        guidefontsize = 24,
        tickfontsize = 24,
        xticks = ([-2e4, -1e4, 0], [L"10^{-20000}", L"10^{-10000}", L"0"]),
    )
    plot!(
        pl,
        Float64.(log.(δ0_very_asym_xs) ./ log(Arb(10))),
        Float64.(δ0_very_asym_ys),
        ribbon = Float64.(Arblib.radius.(Arb, δ0_very_asym_ys)),
        m = :circle,
        ms = 1,
    )
    hline!(Float64[-δ0_rounded, δ0_rounded], color = :green, linestyle = :dash)
    hline!(Float64[-δ0_goal_rounded, δ0_goal_rounded], color = :red, linestyle = :dot)
    savefig(pl, "../figures/publication/BH-F-asymptotic-2.pdf")
    pl
end

# ╔═╡ 8a9cd9da-1719-4b2d-98bf-4e13de8a5ca8
md"""
### Computing the first Fourier coefficient
We can compute the first Fourier coefficient of the solution $u(x) = u_0(x) + w(x)v(x)$ by computing them for $u_0(x)$ and using the bounds for $v(x)$.
"""

# ╔═╡ 64a709b3-0e29-4c03-8bb0-c8bd173ea5a2
md"""
All terms in the approximation $u_0$ are given in terms of their Fourier series, it is therefore easy to compute the first Fourier coefficient. We have

$C_2^{(1)}(x) = \sum_{n = 1}^\infty \frac{\cos(nx)}{n^2}\log(n)$

which for $n = 1$ gives us the coefficient $0$. Next we have

$C_s(x) = \sum_{n = 1}^\infty \frac{\cos(nx)}{n^s}$

which means the first Fourier coefficient is $1$, hence the Fourier coefficient of the Clausen tail is the sum of the coefficients. Finally for the Fourier tail the first Fourier coefficient is simply the first coefficient in the tail.
"""

# ╔═╡ 589218e4-b34f-4f92-b5f4-ae875f300e0e
u0_fourier1 = sum(u0.a) + u0.b[1]

# ╔═╡ d568ad46-2eb7-4b3f-a226-87e8927348ea
md"""
To bound we first fourier coefficient of $w(x)v(x)$ we want to bound

$\frac{1}{\pi}\int_{-\pi}^\pi w(x)v(x)\cos(x)\ dx = \frac{2}{\pi}\int_{0}^\pi w(x)v(x)\cos(x)\ dx.$

Taking the absolute value we get that an upper bound is given by

$\frac{2\|v\|_{L^\infty}}{\pi}\int_{0}^\pi x\sqrt{\log(1 + 1/x)}|\cos(x)|\ dx$
"""

# ╔═╡ cb47f946-fc9d-46c7-9808-e1d4c0b87e5d
vnorm = (1 - D0 - sqrt((1 - D0)^2 - 4δ0 * n0)) / 4n0

# ╔═╡ 790d36ed-3629-4812-bd07-c8e0763205ae
wv_fourier1 = let
    integrand(x; analytic) =
        if Arblib.contains_zero(x)
            analytic && return Arblib.indeterminate!(zero(x))
            x = real(x)
            xᵤ = ubound(Arb, x)

            # Enclosure of x * sqrt(log(1 + inv(x)))
            xsqrtlogx = Arb((0, xᵤ * sqrt(log(1 + inv(xᵤ)))))

            return xsqrtlogx * abs(cos(x))
        else
            x *
            Arblib.sqrt_analytic!(zero(x), log(1 + inv(x)), analytic) *
            Arblib.real_abs!(zero(x), cos(x), analytic)
        end

    bound = real(2vnorm / π * Arblib.integrate(integrand, 0, π, check_analytic = true))

    Arblib.add_error!(Arb(0), bound)
end

# ╔═╡ 54cd2d1c-940b-494a-9546-8c2d0834a2c3
u_fourier1 = u0_fourier1 + wv_fourier1

# ╔═╡ cb6add10-d51a-474e-b708-0c65e896f99a
ArbExtras.format_interval(getinterval(u_fourier1)...)

# ╔═╡ 857b5dfe-89ba-498b-a6b8-130df6cd5008
md"""
# Computing the mean

We can get the wavespeed for the mean-zero solution to the Burgers-Hilbert equation by computing the mean of the solution $u(x) = u_0(x) + w(x)v(x)$. We can compute the mean of $u_0(x)$ directly and bound the mean for $w(x)v(x)$ using the bounds for $v(x)$ given above.
"""

# ╔═╡ da547c40-d829-4c2e-bbdb-b9452b18c8b9
md"""
The mean of $u_0$ can be computed by noticing that the mean of $\tilde{C}_s^{(1)}(x)$ is $-\zeta'(s)$, the mean of $\tilde{C}_s(x)$ is $-\zeta(s)$ and the mean of $\cos(nx) - 1$ is $-1$.

For $w(x)v(x)$ the mean is upper bounded by the mean of $w(x)$ times the bound for $v(x)$ and lower bounded by the negation of this.
"""

# ╔═╡ 6f7c14aa-fbc8-4241-be80-eff13fe1fe88
mean_u0 = let
    # Main term
    res = -u0.a0 * HighestCuspedWave.dzeta(Arb(2))

    # Clausen terms
    for j = 1:u0.N0
        s = 1 - u0.α + j * u0.p0
        res -= u0.a[j] * HighestCuspedWave.zeta(s)
    end

    # Fourier terms
    for n = 1:u0.N1
        res -= u0.b[n]
    end

    res
end

# ╔═╡ ea9e4964-383f-4334-a60f-b8acdc224bca
mean_w = let
    integrand(x; analytic) =
        if Arblib.contains_zero(x)
            analytic && return Arblib.indeterminate!(zero(x))

            # Use monotonicity
            let xᵤ = ubound(Arb, real(x))
                Acb((0, xᵤ * sqrt(log(1 + inv(xᵤ)))))
            end
        else
            x * Arblib.sqrt_analytic!(zero(x), log(1 + inv(x)), analytic)
        end

    real(Arblib.integrate(integrand, 0, π, check_analytic = true)) / π
end

# ╔═╡ 811cdad5-ff5d-4e45-a78f-ef68e4acea1b
mean_wv = mean_w * Arb((-vnorm, vnorm))

# ╔═╡ cbdfcb3a-722f-4ba7-b23f-0230cf7b0631
mean_u = mean_u0 + mean_wv

# ╔═╡ f0ffa9e5-5a53-4d65-9535-22aae7a50186
ArbExtras.format_interval(getinterval(mean_u)...)

# ╔═╡ 87b4508a-2b23-4c3e-b50c-2141f58013b0
md"""
### Plot of approximation with error bounds
"""

# ╔═╡ aa12b03c-9236-4f04-acb5-6eac04234bc4
let xs = range(Arb(0), π, length = 200)
    ys = Folds.map(u0, xs)
    error = Folds.map(x -> u0.w(x) * vnorm, xs)
    error[1] = 0
    xs = Float64[-reverse(xs[2:end]); xs]
    ys = Float64[reverse(ys[2:end]); ys]
    error = [reverse(error[2:end]); error]
    pl = plot(
        xs,
        ys,
        ribbon = Float64.(error),
        xlabel = L"x",
        ylabel = "",
        legend = :none,
        xticks = ([-π, 0, π], [L"-\pi", L"0", L"\pi"]),
        guidefontsize = 16,
        tickfontsize = 16,
    )
    savefig(pl, "../figures/publication/BH-u.pdf")
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
# ╟─f0baf2ec-3f73-4d55-9ce4-754d94d7f3ce
# ╟─4c4cbf2a-3aec-4257-9fac-d8a0418d12d7
# ╟─9c42cbca-805d-4a76-9d7d-68d40129f10a
# ╠═f1ca794b-e1fc-481b-9b75-1e42a0b48a58
# ╟─f1dce520-a035-43e6-9e08-4696a14c5a54
# ╟─ab9d59df-3488-4f8a-a321-d23aab7e01d4
# ╟─36b6ce98-6b64-4eee-a7d7-613772e4b79f
# ╟─05740cb8-a7fd-4add-a4c7-95b617508f28
# ╠═664df5e6-152a-4b9a-adf8-352bdef59055
# ╟─67aa36b0-b77c-4531-a248-f7d474ffd47d
# ╟─b8c5ba34-748e-4c4b-be9c-135240287351
# ╟─1bb84607-4f0f-4e7b-a24f-258b4e581c2c
# ╟─fb6c12ad-3391-4623-a201-412335742930
# ╟─67786c2c-a101-44db-a1b6-f191d2703bb0
# ╠═358b2691-3b1b-4733-a173-acf987a221ea
# ╠═150a963b-03e2-404e-98e4-0fa2cd516dd3
# ╟─ac03e920-25ad-4127-ad87-00e907701da3
# ╟─15e21dde-fb44-47d8-83d8-9f5ffffab74d
# ╟─f5ee2d3b-4a4c-411b-9e46-35164f6e3e83
# ╟─61e8da72-69c0-4af9-bab0-454442dad253
# ╟─77d49bd5-2a4e-468d-a19a-b2e16a132f78
# ╠═de5bed6f-7079-40a9-a5eb-f315abc20ddf
# ╟─de4546e1-4a9f-4d37-b59c-ee4509d09868
# ╠═b0577d0f-77ba-4035-9d3b-ae4d6e5c624f
# ╟─89d54f92-8b3a-4913-875d-31068856fb62
# ╟─3f30187c-0bdb-4cde-9c50-a9ccaaf8bb4b
# ╠═7aab4811-7e37-4828-8e5a-18423f9330e1
# ╟─4d5054ff-6001-4d34-913b-a8029017d217
# ╟─f25bcedf-cf4f-45f6-929a-abba66e7730c
# ╠═8d00fdce-e433-464e-adcd-88c21552439d
# ╠═ff9462b6-3ff1-4fd9-af14-eebce04f349c
# ╠═a80ae53d-c3ec-4f51-9866-1b2f7a23c85f
# ╠═c940e613-622f-43e0-939a-bf30c033b810
# ╠═7fc04a7b-1082-4af7-ac87-426e73276c18
# ╠═74392233-bb8a-44d4-8f7f-ff8b5090142e
# ╟─2494dd2b-088d-402f-8c28-29b3614a5fc6
# ╠═34858269-7714-4210-a678-2b91357c3eb7
# ╟─f0110b12-a782-4840-abbf-ed623175c276
# ╟─16b4c640-6b4c-46a3-a505-e0eb7bbb2632
# ╠═fe196d34-c822-4da0-8ed1-955d2fb6ffa6
# ╟─4654b37d-fcb5-43ab-8b67-5f946637943d
# ╟─1c2d23d0-9e30-4830-a923-3ce15e8deb04
# ╟─7e1e6379-4cde-4a6f-9390-41618a989490
# ╟─8a9cd9da-1719-4b2d-98bf-4e13de8a5ca8
# ╟─64a709b3-0e29-4c03-8bb0-c8bd173ea5a2
# ╠═589218e4-b34f-4f92-b5f4-ae875f300e0e
# ╟─d568ad46-2eb7-4b3f-a226-87e8927348ea
# ╠═cb47f946-fc9d-46c7-9808-e1d4c0b87e5d
# ╠═790d36ed-3629-4812-bd07-c8e0763205ae
# ╠═54cd2d1c-940b-494a-9546-8c2d0834a2c3
# ╠═cb6add10-d51a-474e-b708-0c65e896f99a
# ╟─857b5dfe-89ba-498b-a6b8-130df6cd5008
# ╟─da547c40-d829-4c2e-bbdb-b9452b18c8b9
# ╠═6f7c14aa-fbc8-4241-be80-eff13fe1fe88
# ╠═ea9e4964-383f-4334-a60f-b8acdc224bca
# ╠═811cdad5-ff5d-4e45-a78f-ef68e4acea1b
# ╠═cbdfcb3a-722f-4ba7-b23f-0230cf7b0631
# ╠═f0ffa9e5-5a53-4d65-9535-22aae7a50186
# ╟─87b4508a-2b23-4c3e-b50c-2141f58013b0
# ╟─aa12b03c-9236-4f04-acb5-6eac04234bc4
