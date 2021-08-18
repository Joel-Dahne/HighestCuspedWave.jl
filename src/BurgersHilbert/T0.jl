function T0(u0::BHAnsatz, evaltype::Ball; δ0 = 1e-10, δ1 = 1e-10, δ2 = 1e-10)
    f = T01(u0, evaltype; δ0, δ1)
    g = T02(u0, evaltype; δ2)

    return x -> begin
        ## Integral on [0, x]
        part1 = f(x)

        # Short circuit on non-finite value
        isfinite(part1) || return part1

        ## Integral on [x, π]
        part2 = g(x)

        return part1 + part2
    end
end
