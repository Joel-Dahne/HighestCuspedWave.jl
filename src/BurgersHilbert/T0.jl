function T0(
    u0::BHAnsatz,
    evaltype::Ball;
    δ0 = 1e-10,
    δ1 = 1e-10,
    δ2 = 1e-10,
    rtol = -1, # Not used
    atol = -1, # Not used
    show_trace = false, # Not used
)
    f = T01(u0, evaltype; δ0, δ1)
    g = T02(u0, evaltype; δ2)

    return x -> begin
        ## Integral on [0, x]
        part1 = f(x)

        if isnan(part1)
            # Short circuit on NaN
            return part1
        end

        ## Integral on [x, π]
        part2 = g(x)
        return 2(part1 + part2)
    end
end
