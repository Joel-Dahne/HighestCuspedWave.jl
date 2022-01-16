export polylog

"""
    polylog(s, z)

Compute the polylogarithm ``Li_s(z)``.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.
"""
function polylog(s::Union{Acb,Integer}, z::Acb)
    if iswide(s) # If this is true then s is always an Acb
        # Degree of Taylor expansion, could possibly be tuned
        degree = 2

        s_mid = Acb(Arblib.midref(Arblib.realref(s)), Arblib.midref(Arblib.imagref(s)))

        # Compute the rest term of the Taylor expansion
        w = Arblib.polylog_series!(
            AcbSeries(degree = degree + 1, prec = precision(z)),
            AcbSeries([s, 1]),
            z,
            degree + 2,
        )

        restterm = (s - s_mid)^(degree + 1) * w[degree+1]

        # Compute the Taylor polynomial at the midpoint of x
        w_mid = Arblib.polylog_series!(
            AcbSeries(prec = precision(z); degree),
            AcbSeries([s_mid, 1]),
            z,
            degree + 1,
        )

        # Evaluate the Taylor polynomial on s - s_mid and add the rest
        # term
        res = w_mid(s - s_mid) + restterm

        # If the resulting enclosure is not contained in the enclosure
        # coming from w[0] then take their intersection. Notice that
        # they will always intersect, so taking the intersection is
        # always okay.
        if !Arblib.contains(w[0], res)
            Arblib.intersection!(
                Arblib.realref(res),
                Arblib.realref(res),
                Arblib.realref(w[0]),
            )
            Arblib.intersection!(
                Arblib.imagref(res),
                Arblib.imagref(res),
                Arblib.imagref(w[0]),
            )
        end

        return res
    end

    return Arblib.polylog!(zero(z), s, z)
end

polylog(s::AcbSeries, z::Acb) = Arblib.polylog_series!(zero(s), s, z, length(s))
