# Contains coefficients used in the asymptotic expansions

"""
    a0(u0::FractionalKdVAnsatz, j::Integer)

Compute ``a_j^0`` from Lemma 3.2.
"""
function a0(u0::FractionalKdVAnsatz, j::Integer)
    s = 1 - u0.α + j * u0.p0
    f(s) = gamma(1 - s) * sinpi(s / 2)

    if iswide(s)
        return u0.a[j] * ArbExtras.enclosure_series(f, s, degree = 2)
    end

    return u0.a[j] * f(s)
end

"""
    K0(u0::FractionalKdVAnsatz, m::Integer)

Compute ``K_m^0`` from Lemma 3.2.
"""
function K0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for j = 0:u0.N0
        s = 1 - u0.α + j * u0.p0 - 2m
        if iswide(s)
            res += u0.a[j] * ArbExtras.enclosure_series(zeta, s, degree = 2)
        else
            res += u0.a[j] * zeta(s)
        end
    end

    # We don't expect 2m to be greater than 20 so we can do this with
    # Int64
    res /= (-1)^m * factorial(2m)

    return res
end

"""
    termK0(u0::FractionalKdVAnsatz, m::Integer)

Compute `K0(u0, m) + (-1)^m * sum(n^(2m) * b_n for n in 1:u0.N1) /
factorial(2m)` occurring in Lemma 3.2.

"""
function termK0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for n = 1:u0.N1
        res += oftype(res, n)^(2m) * u0.b[n]
    end

    # We don't expect 2m to be greater than 20 so we can do this with
    # Int64
    res /= (-1)^m * factorial(2m)

    return res + K0(u0, m)
end

"""
    A0(u0::FractionalKdVAnsatz, j::Integer)

Compute ``A_j^0`` from Lemma 3.2.
"""
function A0(u0::FractionalKdVAnsatz, j::Integer)
    s = 1 - 2u0.α + j * u0.p0

    f(s) = gamma(1 - s) * sinpi(s / 2)

    if iswide(s)
        return u0.a[j] * ArbExtras.enclosure_series(f, s, degree = 2)
    end

    return u0.a[j] * f(s)
end

"""
    L0(u0::FractionalKdVAnsatz, m::Integer)

Compute ``L_m^0`` from Lemma 3.2.
"""
function L0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for j = 0:u0.N0
        s = 1 - 2 * u0.α + j * u0.p0 - 2m
        if iswide(s)
            res += u0.a[j] * ArbExtras.enclosure_series(zeta, s, degree = 2)
        else
            res += u0.a[j] * zeta(s)
        end
    end

    # We don't expect 2m to be greater than 20 so we can do this with
    # Int64
    res /= (-1)^m * factorial(2m)

    return res
end

"""
    termL0(u0::FractionalKdVAnsatz, m::Integer)

Compute `L0(u0, m) + (-1)^m * sum(n^(2m + α) * b_n for n in 1:N1) /
factorial(2m)` occurring in Lemma 3.2.
"""
function termL0(u0::FractionalKdVAnsatz, m::Integer)
    res = zero(u0.α)

    for n = 1:u0.N1
        res += oftype(res, n)^(2m + u0.α) * u0.b[n]
    end

    # We don't expect 2m to be greater than 20 so we can do this with
    # Int64
    res /= (-1)^m * factorial(2m)

    return res + L0(u0, m)
end
