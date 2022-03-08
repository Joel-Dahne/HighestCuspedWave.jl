function T0(
    u0::BHAnsatz,
    evaltype::Ball;
    δ1::Arb = Arb(1e-5),
    δ2::Arb = Arb(1e-5),
    skip_div_u0 = false,
)
    f = T01(u0, evaltype, skip_div_u0 = true; δ1)
    g = T02(u0, evaltype, skip_div_u0 = true; δ2)

    return x -> begin
        ## Integral on [0, x]
        part1 = f(x)

        # Short circuit on non-finite value
        isfinite(part1) || return part1

        ## Integral on [x, π]
        part2 = g(x)

        isfinite(part2) || return part2 # Avoid computing u0(x)

        if skip_div_u0
            return part1 + part2
        else
            return (part1 + part2) / u0(x)
        end
    end
end

function T0(
    u0::BHAnsatz,
    evaltype::Asymptotic;
    non_asymptotic_u0 = false,
    ϵ::Arb = Arb(2e-1),
)
    f = T01(u0, evaltype; non_asymptotic_u0, ϵ)
    g = T02(u0, evaltype; non_asymptotic_u0, ϵ)

    return x -> f(x) + g(x)
end
