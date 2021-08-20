function T0(
    u0::BHAnsatz,
    evaltype::Ball;
    δ0::Arb = Arb(1e-10),
    δ1::Arb = Arb(1e-10),
    δ2::Arb = Arb(1e-10),
    skip_div_u0 = false,
)
    f = T01(u0, evaltype, skip_div_u0 = true; δ0, δ1)
    g = T02(u0, evaltype, skip_div_u0 = true; δ2)

    return x -> begin
        ## Integral on [0, x]
        part1 = f(x)

        # Short circuit on non-finite value
        isfinite(part1) || return part1

        ## Integral on [x, π]
        part2 = g(x)

        if skip_div_u0
            return part1 + part2
        else
            return (part1 + part2) / u0(x)
        end
    end
end
