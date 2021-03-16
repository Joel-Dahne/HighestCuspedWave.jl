"""
        eval_expansion(u0::BHAnsatz, expansion, x)
    Evaluate the given expansion. The term `((i, m), y)` is evaluated
    to `y*log(abs(x))^i*abs(x)^m` and then they are all summed.

    In general x needs to be given both when computing the expansion and
    when evaluating it.
"""
function eval_expansion(u0::BHAnsatz{T},
                        expansion::AbstractDict{NTuple{2, Int}, T},
                        x,
                        ) where {T}
    res = zero(u0.a0)

    for ((i, m), y) in expansion
        if i > 0
            # TODO: Handle the case when 0 ∈ x, i > 0, m > 0. Then this
            # should be finite but will not work as given now.
            res += y*log(abs(x))^i*abspow(x, m)
        else
            res += y*abspow(x, m)
        end
    end

    return res
end

function (u0::BHAnsatz{T})(x, ::Ball) where {T}
    conv = T == arb ? parent(u0.a0) : a -> convert(T, a)
    x = conv(x)
    res = u0.a0*(Ci(x, 2, 1) - zeta(conv(2), d = 1))
    res += u0.a1*(Ci(x, 2) - zeta(conv(2)))

    for n in 1:u0.N
        res += u0.b[n]*(cos(n*x) - 1)
    end

    return res
end

function (u0::BHAnsatz)(x, ::Asymptotic; M::Integer = 3)
    return eval_expansion(u0, u0(x, AsymptoticExpansion(); M), x)
end

"""
    (u0::AbstractAnsatz)(x, ::AsymptoticExpansion; M = 3)
Return a dictionary containing the terms in the asymptotic expansion
of `u0`. The element `((i, m), a)` in the dictionary corresponds to
the term `a*log(abs(x))^i*abs(x)^m`. The highest term is an error term
which makes sure that evaluation of the expansion gives an enclosure
of the result.
"""
function (u0::BHAnsatz{arb})(x, ::AsymptoticExpansion; M::Integer = 3)
    # TODO: Check this
    #@assert M >= 3

    RR = parent(u0.a0)
    π = RR(pi)
    γ = RR(stieltjes(arb, 0))

    res = OrderedDict{NTuple{2, Int}, arb}()

    res[(1, 1)] = -u0.a0*π/2
    res[(0, 1)] = -u0.a0*(γ - 1)*π/2 - u0.a1*π/2

    for m in 1:M-1
        term = u0.a0*zeta(RR(2 - 2m), d = 1) + u0.a1*zeta(RR(2 - 2m))
        for n in 1:u0.N
            term += fmpz(n)^(2m)*u0.b[n]
        end
        res[(0, 2m)] = (-1)^m*term/factorial(fmpz(2m))
    end

    #TODO: Add error term

    return res
end

function (u0::BHAnsatz{Arb})(x, ::AsymptoticExpansion; M::Integer = 3)
    # TODO: Check this
    #@assert M >= 3

    π = Arb(Irrational{:π}())
    γ = Arb(Irrational{:γ}())

    res = OrderedDict{NTuple{2, Int}, Arb}()

    res[(1, 1)] = -u0.a0*π/2
    res[(0, 1)] = -u0.a0*(γ - 1)*π/2 - u0.a1*π/2

    for m in 1:M-1
        term = u0.a0*zeta(Arb(2 - 2m), d = 1) + u0.a1*zeta(Arb(2 - 2m))
        for n in 1:u0.N
            term += Arb(n)^(2m)*u0.b[n]
        end
        res[(0, 2m)] = (-1)^m*term/factorial(Arb(2m))
    end

    #TODO: Add error term

    return res
end

function H(u0::BHAnsatz{T}, ::Ball) where {T}
    conv = T == arb ? parent(u0.a0) : a -> convert(T, a)

    return x -> begin
        x = conv(x)
        res = -u0.a0*(Ci(x, 3, 1) - zeta(conv(3), d = 1))
        res += -u0.a1*(Ci(x, 3) - zeta(conv(3)))

        for n in 1:u0.N
            res -= u0.b[n]/n*(cos(n*x) - 1)
        end

        return res
    end
end

function H(u0::BHAnsatz, ::Asymptotic; M::Integer = 3)
    f = H(u0, AsymptoticExpansion(); M)

    return x -> eval_expansion(u0, f(x), x)
end

"""
    H(u0::AbstractAnsatz, ::AsymptoticExpansion; M = 3)
Returns a function such that `H(u0, AsymptoticExpansion)(x)` returns a
dictionary containing the terms in the asymptotic expansion of
`H(u0)`. The element `((i, m), a)` in the dictionary corresponds to
the term `a*log(abs(x))^i*abs(x)^m`. The highest term is an error term
which makes sure that evaluation of the expansion gives an enclosure
of the result.
"""
function H(u0::BHAnsatz{arb}, ::AsymptoticExpansion; M::Integer = 3)
    # TODO: Check this
    #@assert M >= 3

    RR = parent(u0.a0)
    π = RR(pi)
    γ = RR(stieltjes(arb, 0))
    γ₁ = RR(stieltjes(arb, 1))

    return x -> begin
        res = OrderedDict{NTuple{2, Int}, arb}()

        res[(2, 2)] = -u0.a0/4
        res[(1, 2)] = u0.a0*(3//4 - γ/2) - u0.a1/2

        for m in 1:M - 1
            if m == 1
                term = -u0.a0*(3γ - γ^2 - 2γ₁ - 7//2 + (π^2)/12)/2 - u0.a1*3//2
            else
                term = -u0.a0*zeta(RR(3 - 2m), d = 1) - u0.a1*zeta(RR(3 - 2m))
            end
            for n in 1:u0.N
                term -= fmpz(n)^(2m - 1)*u0.b[n]
            end
            res[(0, 2m)] = (-1)^m*term/factorial(fmpz(2m))
        end

        # TODO: Add error term

        return res
    end
end

function H(u0::BHAnsatz{Arb}, ::AsymptoticExpansion; M::Integer = 3)
    # TODO: Check this
    #@assert M >= 3

    π = Arb(Irrational{:π}())
    γ = Arb(Irrational{:γ}())
    γ₁ = stieltjes(Arb, 1)

    return x -> begin
        res = OrderedDict{NTuple{2, Int}, Arb}()

        res[(2, 2)] = -u0.a0/4
        res[(1, 2)] = u0.a0*(3//4 - γ/2) - u0.a1/2

        for m in 1:M - 1
            if m == 1
                term = -u0.a0*(3γ - γ^2 - 2γ₁ - 7//2 + (π^2)/12)/2 - u0.a1*3//2
            else
                term = -u0.a0*zeta(Arb(3 - 2m), d = 1) - u0.a1*zeta(Arb(3 - 2m))
            end
            for n in 1:u0.N
                term -= Arb(n)^(2m - 1)*u0.b[n]
            end
            res[(0, 2m)] = (-1)^m*term/factorial(Arb(2m))
        end

        # TODO: Add error term

        return res
    end
end

function D(u0::BHAnsatz, ::Asymptotic; M::Integer = 3)
    f = D(u0, AsymptoticExpansion(); M)

    return x -> eval_expansion(u0, f(x), x)
end

function D(u0::BHAnsatz, evaltype::AsymptoticExpansion; M::Integer = 3)
    f = x -> u0(x, evaltype, M = M)
    g = H(u0, evaltype, M = M)

    return x -> begin
        expansion1 = f(x)
        expansion2 = g(x)
        expansion = empty(expansion1)

        # u0^2/2 term
        for ((i1, m1), a1) in expansion1
            for ((i2, m2), a2) in expansion1
                key = (i1 + i2, m1 + m2)
                expansion[key] = get(expansion, key, zero(x)) + a1*a2/2
            end
        end

        # H term
        merge!(+, expansion, expansion2)

        # TODO: Will have to handle terms which are supposed to be
        # identically equal to zero.

        return expansion
    end
end

"""
    D(u0::BHAnsatz, xs::AbstractVector)
Returns a function such that `D(u0, xs)(b)` computes `D(u0)(x)` on the
points `x ∈ xs` with `u0.b` set to the given values. Does this in an
efficient way by precomputing as much as possible.

NOTE: This is **not** rigorous!
"""
function D(u0::BHAnsatz, xs::AbstractVector)
    u0_xs_a0_precomputed = zeros(length(xs))
    u0_xs_b_precomputed = zeros(length(xs), u0.N)
    Hu0_xs_a0_precomputed = zeros(length(xs))
    Hu0_xs_b_precomputed = zeros(length(xs), u0.N)

    for i in eachindex(xs)
        x = xs[i]

        u0_xs_a0_precomputed[i] = u0.a0*(Ci(x, 2, 1) - zeta(2, d = 1)) + u0.a1*(Ci(x, 2) - zeta(2))
        Hu0_xs_a0_precomputed[i] = -u0.a0*(Ci(x, 3, 1) - zeta(3, d = 1)) - u0.a1*(Ci(x, 3) - zeta(3))

        for n in 1:u0.N
            u0_xs_b_precomputed[i, n] = cos(n*x) - 1
            Hu0_xs_b_precomputed[i, n] = -(cos(n*x) - 1)/n
        end
    end

    return b -> begin
        return (
            (u0_xs_a0_precomputed .+ u0_xs_b_precomputed*b).^2 ./ 2
            .+ (Hu0_xs_a0_precomputed .+ Hu0_xs_b_precomputed*b)
        )
    end
end
