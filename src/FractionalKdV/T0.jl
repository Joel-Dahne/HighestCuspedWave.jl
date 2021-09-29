# TODO: There is a lot more tuning to be done!
function T0(
    u0::FractionalKdVAnsatz{arb},
    evaltype::Ball;
    δ0::arb = ifelse(isone(u0.p), parent(u0.α)(1e-4), parent(u0.α)(1e-3)),
    δ1::arb = ifelse(isone(u0.p), parent(u0.α)(1e-4), parent(u0.α)(1e-3)),
    δ2::arb = parent(u0.α)(1e-2),
    ϵ::arb = 1 + u0.α,
    rtol = -1.0,
    atol = -1.0,
    show_trace = false,
)
    f = T01(u0, evaltype; δ0, δ1, rtol, atol, show_trace)
    g = T02(u0, evaltype; δ2, ϵ, rtol, atol, show_trace)

    return x -> begin
        ## Integral on [0, x] - Change to t = y/x
        part1 = f(x)

        # Short circuit on a non-finite result
        isfinite(part1) || return part1

        ## Integral on [x, π]
        part2 = g(x)
        return part1 + part2
    end
end

function T0(
    u0::FractionalKdVAnsatz{Arb},
    evaltype::Ball;
    δ0::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    δ1::Arf = ifelse(isone(u0.p), Arf(1e-4), Arf(1e-3)),
    δ2::Arf = Arf(1e-2),
    ϵ::Arb = 1 + u0.α,
    skip_div_u0 = false,
)
    f = T01(u0, evaltype, skip_div_u0 = true; δ0, δ1)
    g = T02(u0, evaltype, skip_div_u0 = true; δ2, ϵ)

    return x -> begin
        ## Integral on [0, x] - Change to t = y/x
        part1 = f(x)

        # Short circuit on a non-finite result
        isfinite(part1) || return part1

        ## Integral on [x, π]
        part2 = g(x)

        isfinite(part2) || return part2

        if skip_div_u0
            return part1 + part2
        else
            return (part1 + part2) / u0(x)
        end
    end
end
