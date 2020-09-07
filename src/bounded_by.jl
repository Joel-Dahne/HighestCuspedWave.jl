export bounded_by

"""
    bounded_by(f, a, b, C)
Return true if `|f|` is bounded by `C` on the interval `[a, b]`, i.e.
`|f(x)| <= C ∀ x ∈ [a, b].

If `use_taylor` is true then compute the maximum on each subinterval
by computing the maximum of the Taylor expansion.
"""
function bounded_by(f,
                    a::arb,
                    b::arb,
                    C::arb;
                    show_trace = false,
                    show_evaluations = false,
                    use_taylor = false,
                    )
    if a > b
        # Empty interval, always true
        return true
    end
    if a == b
        # Thin interval, check the only existing point
        return f(a) < C
    end
    # Only works for finite values of a and b
    @assert isfinite(a) && isfinite(b)

    intervals = [(a, b)]

    iterations = 0
    max_value = parent(a)(NaN)

    if show_trace
        @printf "%6s %11s %s\n" "Iter" "Intervals" "Bound"
    end

    while !isempty(intervals)
        iterations += 1

        if show_trace
            @printf "%6d %11d %s\n" iterations length(intervals) string(getinterval(max_value))
        end

        res = similar(intervals, arb)
        Threads.@threads for i in eachindex(intervals)
            if use_taylor
                res[i] = ArbTools.maximumtaylor(f, intervals[i], 8, absmax = true)
            else
                res[i] = abs(f(setinterval(intervals[i]...)))
            end

            if show_evaluations
                @show res[i]
            end
        end

        next_intervals = Vector{eltype(intervals)}()
        max_value = parent(a)(-Inf)

        for (i, (c, d)) in enumerate(intervals)
            y = res[i]

            max_value = max(max_value, y)
            if y > C
                # If f([c, d]) is greater than C then C is not a bound
                # and we return false
                @show (c, d) ArbTools.getinterval(y)
                return false
            elseif !(y <= C)
                # If we cannot determine if f([c, d]) is less than or
                # equal to C then bisect it
                midpoint = 0.5*(c + d)
                push!(next_intervals, (c, midpoint))
                push!(next_intervals, (midpoint, d))
            end
        end

        intervals = next_intervals
    end

    return true
end
