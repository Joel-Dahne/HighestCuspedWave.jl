export T0

"""
    T0(u0::FractionalKdVAnsatz, evaltype)
Returns a function such that T0(u0)(x) is the function whose supremum
on [0, π] gives C_B. The strategy for evaluation depends on type of
evaltype.

TODO: There is a lot more tuning to be done!
"""
T0(u0::FractionalKdVAnsatz; kwargs...) = T0(u0, Ball(); kwargs...)

function T0(u0::FractionalKdVAnsatz{arb},
            evaltype::Ball;
            rtol = -1.0,
            atol = -1.0,
            show_trace = false,
            δ0::arb = ifelse(isone(u0.p), parent(u0.α)(1e-4), parent(u0.α)(1e-3)),
            δ1::arb = ifelse(isone(u0.p), parent(u0.α)(1e-4), parent(u0.α)(1e-3)),
            δ2::arb = parent(u0.α)(1e-2),
            ϵ::arb = 1 + u0.α,
            )
    # Set up parameters for T01
    # δ0 and δ1 are given as arguments

    # Set up parameters for T02
    # δ2 and ϵ are given as arguments

    f = T01(u0, evaltype; δ0, δ1, rtol, atol, show_trace)
    g = T02(u0, evaltype; δ2, ϵ, rtol, atol, show_trace)

    return x -> begin
        ## Integral on [0, x] - Change to t = y/x
        part1 = f(x)

        if isnan(part1)
            # Short circuit on NaN
            return part1
        end

        ## Integral on [x, π]
        part2 = g(x)
        return part1 + part2
    end
end

function T0(u0::FractionalKdVAnsatz{T}, evaltype::Asymptotic) where {T}
    @warn "Asymptotic evaluation of T0 yet not implemented"
    return x -> begin
        return zero(u0.α)

        ## Integral on [0, x]
        part1 = T01(u0, evaltype)(x)

        if isnan(part1)
            # Short circuit on NaN
            return part1
        end

        ## Integral on [x, π]
        part2 = T02(u0, evaltype)(x)

        return part1 + part2
    end
end

"""
    T01(u0::FractionalKdVAnsatz; δ1, δ2)
Returns a function such that T01(u0, δ1 = δ1, δ2 = δ2)(x) computes the integral
T_{0,1} from the paper.
"""
T01(u0; kwargs...) = T01(u0, Ball(); kwargs...)

function T01(u0::FractionalKdVAnsatz{arb},
             evaltype::Ball;
             δ0::arb = parent(u0.α)(1e-2),
             δ1::arb = parent(u0.α)(1e-2),
             rtol = -1.0,
             atol = -1.0,
             show_trace = false,
             )
    f = T011(u0, evaltype; δ0)
    g = T012(u0, evaltype; δ0, δ1, rtol, atol, show_trace)
    h = T013(u0, evaltype; δ1)
    return x -> begin
        return f(x) + g(x) + h(x)
    end
end

"""
    T01(u0, Asymptotic())
Returns a function such that `T01(u0, Asymptotic())(x)` computes the
integral T_{0,1} from the paper using an evaluation strategy that
works asymptotically as `x` goes to 0.

It splits the Clausians into the main singular part and the analytic
expansion. The integral of the singular part is computed by finding
where the integrand is positive respectively negative and then
integrating explicitly. The expansion is integrated term wise and the
resulting sum is bounded.

FIXME: There is something wrong somewhere. The result doesn't match up
completely with the non-asymptotic version for small values of `x`. I
don't know which one is the correct one...
```
u0 = FractionalKdVAnsatz(RR(-0.6), pp = RR(1))
x = RR(1e-6)
a = HighestCuspedWave.T01(u0, Ball())(x)
b = HighestCuspedWave.T01(u0, Asymptotic())(x)
overlaps(a, b)
```
"""
function T01(u0::FractionalKdVAnsatz{arb},
             ::Asymptotic,
             )
    @warn "T01(u0, Asymptotic()) is not yet complete"
    Γ = Nemo.gamma
    α = u0.α
    p = u0.p
    RR = parent(α)
    CC = ComplexField(prec(RR))
    π = RR(pi)

    return x -> begin
        if p == 1
            # TODO: Bound the tail - the terms go to zero extremely
            # fast so it should be negligible
            c_ϵ = zero(α)
            for m in 1:10
                m = fmpz(m)
                c_ϵ += (-one(α))^m/(factorial(fmpz(2m)))*zeta(-α - 2m)*
                    m*(RR(4)^m - 1)/((m + 1)*(2m + 1))*abspow(x, RR(2m - 2))
            end
            c_ϵ *= 2
        else
            # TODO: Compute c_ϵ from the sum
            c_ϵ = one(α)
        end
        # Compute c_α = ∫0^1 |(1 - t)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)|tᵖ dt
        # Find the unique zero of the integrand on [0, 1]
        # PROVE: That there is at most one zero on [0, 1]
        f = t -> (1 - t)^(-1 - α) + (1 + t)^(-1 - α) - 2t^(-1 - α)
        roots, flags = isolateroots(
            f,
            parent(α)(0.1),
            parent(α)(0.9),
            refine = true,
            evaltype = :taylor,
        )
        @assert only(flags)
        s = setunion(only(roots)...)

        # PROVE: That the imaginary parts cancel out
        c_α = let p = CC(p), α = CC(α), s = CC(s), m1 = CC(-1)
            m1^(1 - p)*(beta_inc(1 + p, -α, m1) - 2beta_inc(1 + p, -α, -s)) -
                2beta_inc(1 + p, -α, s)
        end
        @assert contains_zero(imag(c_α))
        c_α = real(c_α) + (2 - 4s^(p - α))/(α - p) + Γ(-α)*Γ(1 + p)/Γ(1 - α + p)

        # Version without asymptotically expanding u0(x)
        #res = abs(Γ(1 + α)*sinpi(α/2))*c_α*abspow(x, p - α) + ball(zero(c_ϵ), c_ϵ)*abspow(x, 4)
        #return res/(π*abspow(x, p)*u0(x))

        res = abs(Γ(1 + α)*sinpi(α/2))*c_α + ball(zero(c_ϵ), c_ϵ)*abspow(x, 3 + α)
        # Ball containing 1 + hat(u0)(x)
        L = ball(parent(α)(1), c(u0, ArbTools.abs_ubound(x))*abspow(x, u0.p0))

        return L/(π*a0(u0, 0))*res
    end
end

"""
    T011(u0::FractionalKdVAnstaz{arb}; δ0)
Returns a function such that T011(u0, δ0 = δ0)(x) computes the
integral T_{0,1,1} from the paper.

The strategy for evaluation is to compute an expansion for the
integrand which is then integrated termwise on the interval. The first
two terms inside the absolute value are analytic and their Taylor
expansions of degree `N - 1` are computed at `t = 0` and the error
term is enclose. The last term inside the absolute value is not
analytic and instead we use the expansion for Clausians from the paper
and also here enclose the error term.

The value inside the absolute value has constant sign so we can remove
it. Switching integration and summation gives us terms of the form
`∫_0^δ0 t^s*t^p dt` where `s` depends on the term and `p = u0.p`. This
is easily calculated to be `δ0^(s + p + 1*/(s + p + 1)`. The error
handled as constant values which are just multiplied by the length of
the interval.
"""
T011(u0::FractionalKdVAnsatz{arb}; kwargs...) = T011(u0, Ball(); kwargs...)

function T011(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ0::arb = parent(u0.α)(1e-2),
              N::Integer = 3,
              )
    Γ = Nemo.gamma
    α = u0.α

    PP = ArbPolyRing(parent(α), :x)

    M = div(N, 2) + 1

    return x -> begin
        # Analytic terms
        (P, P_E) = taylor_with_error(zero(α), setunion(zero(α), δ0), N) do t
            Ci(x*(1 - t), -α) + Ci(x*(1 + t), -α)
        end
        P_restterm = ball(zero(α), P_E*δ0^N)

        # Singular term
        (C, e, P2, P2_E) = Ci_expansion(x*δ0, -α, M)
        C *= x^e
        for m = 1:M-1
            P2[2m] *= x^(2m)
        end
        P2_restterm = P2_E*(x*δ0)^(2M)

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        res -= 2C*δ0^(e + u0.p + 1)/(e + u0.p + 1)

        # Integrate the analytic terms
        full_series = P - 2P2
        for i = 0:N-1
            res += full_series[i]*δ0^(i + u0.p + 1)/(i + u0.p + 1)
        end

        # Add the error term
        res += δ0*(P_restterm - P2_restterm)

        # Prove: that the expression inside the absolute value of the
        # integrand is negative
        return -res*x/(parent(α)(π)*u0(x))
    end
end

"""
    T012(u0::FractionalKdVAnsatz{arb}; δ0, δ1)
Returns a function such that T012(u0, δ0 = δ0, δ1 = δ1)(x) computes
the integral T_{0,1,2} from the paper.
"""
T012(u0::FractionalKdVAnsatz{arb}; kwargs...) = T012(u0, Ball(); kwargs...)

function T012(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ0::arb = parent(u0.α)(1e-2),
              δ1::arb = parent(u0.α)(1e-2),
              rtol = -1.0,
              atol = -1.0,
              show_trace = false,
              )
    α = u0.α
    CC = ComplexField(prec(parent(α)))
    mα = CC(-α)
    a = CC(δ0)
    b = CC(1 - δ1)

    return x -> begin
        # PROVE: That there are no branch cuts that interact with the
        # integral
        F(t, analytic = false) = begin
            if isreal(t)
                res = CC(Ci(x*(1 - real(t)), -α) + Ci(x*(1 + real(t)), -α) - 2Ci(x*real(t), -α))
            else
                res = Ci(x*(1 - t), mα) + Ci(x*(1 + t), mα) - 2Ci(x*t, mα)
            end

            ArbTools.real_abs(res, analytic = analytic)*t^u0.p
        end
        res = real(ArbTools.integrate(CC, F, a, b,
                                      rel_tol = rtol,
                                      abs_tol = atol,
                                      eval_limit = 3000,
                                      verbose = Int(show_trace),
                                      ))

        return res*x/(parent(u0.α)(π)*u0(x))
    end
end

"""
    T013(u0::FractionalKdVAnstaz{arb}; δ1)
Returns a function such that T013(u0, δ0 = δ0)(x) computes the integral
T_{0,1,3} from the paper.

The strategy for evaluation is the same as for T011 except that the
first term is singular and the last two are analytic and their Taylor
expansion is computed at `t = 1`.

The integral that needs to be computed in this case is `∫_(1 - δ1)^1
(1 - t)^s*t^p dt` which is given by `Γ(1 + s)*Γ(1 + p)/Γ(2 + s + p) -
B(1 + p, 1 + s; 1 - δ1)` where B(a, b; z) is the incomplete
Beta-function.

If `x` is equal or very close to π (determined by `ϵ`) then the Taylor
expansion gives a very poor approximation for `Ci(x*(t + 1), -α)`. In
this case we make use of the fact that it's 2π periodic and even, so
that `Ci(x*(t + 1), -α) = Ci(x*(t + 1) - 2π, -α) = Ci(2π - x*(t + 1),
-α)`, to be able to use the asymptotic expansion instead. That gives
us the integral ```∫_(1 - δ1)^1 (2π - x*(t + 1))^s*t^p dt` which is
given by `(2π - x)^(1 + p + s)*x^(-1 - p)*(B(1 + p, 1 + s, x/(2π - x))
- B(1 + p, 1 + s, (x - δ1*x)/(2π - x)))```. The value of `x/(2π - x)`
will always be less than or equal to 1 for `x` less than or equal to
π, however due to overestimation the enclosing ball might contain
values greater than one, we therefore have to use `beta_inc_zeroone`
to be able to get finite results in that case.

TODO: We could precompute some of the values, in particular the
beta_inc functions can be precomputed.
"""
T013(u0::FractionalKdVAnsatz{arb}; kwargs...) = T013(u0, Ball(); kwargs...)

function T013(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ1::arb = parent(u0.α)(1e-2),
              ϵ::arb = parent(u0.α)(1e-2),
              N::Integer = 3,
              )
    Γ = Nemo.gamma
    α = u0.α
    π = parent(α)(pi)

    PP = ArbPolyRing(parent(α), :x)

    M = div(N, 2) + 1

    return x -> begin
        # Determine if the asymptotic expansion or the Taylor
        # expansion should be used for the second term
        use_asymptotic = π - x < ϵ

        # Analytic terms
        (P, E) = taylor_with_error(one(α), setunion(1 - δ1, one(α)), N) do t
            if !use_asymptotic
                return Ci(x*(1 + t), -α) - 2Ci(x*t, -α)
            else
                return -2Ci(x*t, -α)
            end
        end
        P_restterm = ball(zero(α), E*δ1^N)

        # Singular term
        (C, e, P2, P2_E) = Ci_expansion(x*δ1, -α, M)
        C *= x^e
        for m = 1:M-1
            P2[2m] *= x^(2m)
        end
        P2_restterm = P2_E*(x*δ1)^(2M)

        # Compute the integral
        res = zero(α)
        # Integrate the singular term
        # Using ∫_(1 - δ1)^1 |t - 1|^s*t^(p) dt = ∫_(1 - δ1)^1 (1 - t)^s*t^(p) dt
        res += C*(Γ(1 + e)*Γ(1 + u0.p)/Γ(2 + e + u0.p) - beta_inc(1 + u0.p, 1 + e, 1 - δ1))

        # Integrate the analytic part
        # Using ∫_(1-δ1)^1 (t-1)^i*t^(p) dt = (-1)^i ∫_(1-δ1)^1 (t-1)^i*t^(p) dt
        full_series = P + P2
        for i = 0:N-1
            res += full_series[i] * (-1)^i * (
                Γ(parent(α)(1 + i))*Γ(1 + u0.p)/Γ(2 + i + u0.p)
                - beta_inc(1 + u0.p, parent(α)(1 + i), 1 - δ1)
            )
        end

        # Add the error term
        res += δ1*(P_restterm + P2_restterm)

        if use_asymptotic
            # Handle asymptotic expansion of Ci(x*(t + 1), -α)
            (C, e, P3, P3_E) = Ci_expansion(2π - x*(2 - δ1), -α, M)
            P3_restterm = P2_E*(2π - x*(2 - δ1))^(2M)

            # Add the singular part to the integral
            res += C*(2π - x)^(1 + u0.p + e)*x^(-1 - u0.p)*(
                beta_inc_zeroone(1 + u0.p, 1 + e, x/(2π - x))
                - beta_inc_zeroone(1 + u0.p, 1 + e, (x - δ1*x)/(2π - x))
            )

            for i = 0:2:N-1
                # Only even terms
                res += P3[i]*(2π - x)^(1 + u0.p + i)*x^(-1 - u0.p)*(
                beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), x/(2π - x))
                - beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), (x - δ1*x)/(2π - x))
                )
            end

            # Add error term
            res += δ1*P3_restterm
        end

        # Prove: that the expression inside the absolute value of the
        # integrand is positive
        return res*x/(parent(u0.α)(π)*u0(x))
    end
end

"""
    T02(u0::FractionalKdVAnsatz; δ2)
Returns a function such that T02(u0, δ2 = δ2, ϵ = ϵ)(x) computes the
integral T_{0,2} from the paper.

If `u0.p == 1` it uses a closed form expression for the integral.
Otherwise it computes the integral directly for most of the interval
and uses and asymptotic expansion for `y < x + δ2`.

If `x` is close to π (`π - x < ϵ`) then use only the asymptotic
expansion for the full integral.

The choice of both `δ2` and `ϵ` can be tuned a lot. For `δ2` it
depends on both `x` and `α`, whereas for `ϵ` it only depends on `α`.
For `δ2` it should be larger when both `x` and `α` are larger. In
theory we could compute with several values and take the best result,
but that would likely be to costly. For `ϵ` it's mainly a question of
cost, we don't want to compute the expensive integral when we believe
the asymptotic expansion will be the best anyway, larger values of `α`
should give a lower value of `ϵ`. The choice of `ϵ` is partially based
on /figures/optimal-epsilon-choice.png, but can likely be tuned
further.

TODO: Look closer at computing with the asymptotic expansion and using
the best result. Consider rewriting `T021` and `T022` to be more like
the other methods here.
"""
T02(u0; kwargs...) = T02(u0, Ball(); kwargs...)

function T02(u0::FractionalKdVAnsatz{arb},
             ::Ball;
             δ2::arb = parent(u0.α)(1e-2),
             ϵ::arb = 1 + u0.α,
             rtol = -1.0,
             atol = -1.0,
             show_trace = false,
             )
    π = parent(u0.α)(pi)

    if u0.p == 1
        # Use the closed form expression
        α = u0.α
        p = u0.p
        return x -> begin
            # TODO: Handle the case when x contains π. Then Si(2x, 1 -
            # α) evaluates to NaN. Si has issues as soon as the
            # argument is a ball containing a multiple of 2π.
            res = Ci(x + π, 2 - α) - Ci(π, 2 - α) + Ci(x, 2 - α) -
                (Ci(2x, 2 - α) + zeta(2 - α))/2 + x*Si(x, 1 - α)
            if π - x < 1e-4
                # When 2x is close to 2π direct evaluation
                # of Si fails. Use that Si(2x, 1 - α) = Si(2x - 2π, 1
                # - α) and the asymptotic expansion.
                y = 2x - 2π
                M = 3
                C, e, P, E = Si_expansion(y, 1 - α, M)
                res -= x/2*(-C*abspow(y, e) + evaluate(P.poly, y) + E*abs(y)^(2M + 1))
            else
                res -= x/2*Si(2x, 1 - α)
            end
            return 2/(π*u0.w(x)*u0(x))*res
        end
    else
        return x -> begin
            a = ArbTools.ubound(x + δ2)

            # Compute with the asymptotic expansion on the whole interval
            res_asymptotic = T021(u0, Ball(), π, x, ϵ = π)

            if π < a || π - x < ϵ
                return res_asymptotic
            end

            part1 = T021(u0, Ball(), a, x; ϵ)

            part2 = T022(u0, Ball(), a, x; rtol, atol, show_trace)

            res = part1 + part2
            if radius(res) < radius(res_asymptotic)
                return res
            else
                return res_asymptotic
            end
        end
    end
end

"""
    T02(u0, Asymptotic())
Returns a function such that `T02(u0, Asymptotic())(x)` computes the
integral T_{0,2} from the paper using an evaluation strategy that
works asymptotically as `x` goes to 0.
"""
function T02(u0::FractionalKdVAnsatz{arb},
             ::Asymptotic;
             N::Integer = 10
             )
    if u0.p != 1
        @warn "T02(u0, Asymptotic()) is not yet rigorous when u0.p != 1"
    end
    Γ = Nemo.gamma
    α = u0.α
    p = u0.p
    π = parent(α)(pi)

    return x -> begin
        if p == 1
            (A, _, P1, E1) = Ci_expansion(x, 2 - α, 2)
            (_, _, _, E1p) = Ci_expansion(2x, 2 - α, 2)
            (B, _, P2, E2) = Si_expansion(x, 1 - α, 1)
            (_, _, _, E2p) = Si_expansion(2x, 1 - α, 1)
            # PROVE: That P3[1] == P3[3] == 0
            (P3, E3) = taylor_with_error(π, setunion(π, π + x), 4) do y
                return Ci(y, 2 - α)
            end

            c_α = A*(1 - parent(α)(2)^(-α)) + B*(1 - parent(α)(2)^(-α - 1))
            K = 2/(π*a0(u0, 0))
            # Ball containing 1 + hat(u0)(x)
            L = ball(parent(α)(1), c(u0, ArbTools.abs_ubound(x))*abspow(x, u0.p0))

            res = (
                K*c_α*L
                + K/2*L*(zeta(-α) - Ci(π, -α))*abspow(x, 1 + α)
                + K*L*abspow(x, α - 1)*x^4*(
                    E3 + E1 - 8*E1p + E2 - 8*E2p
                )
            )

            return res
        end

        S = zero(u0.α)
        for k = reverse(1:N)
            k = parent(α)(k)
            S += k^(α - 1 - p)*(cos(k*x) - 1)*(cosint(p + 1, k*x) - cosintpi(p + 1, k))
        end

        # S is bounded by - later on we want to use this
        S_bound = cosintpi(p + 1, one(p))*(Ci(x, p + 1 - α) - zeta(p + 1 - α)) -
            abspow(x, p)*(Ci(x, 1 - α) - zeta(1 - α))

        L = ball(parent(α)(1), c(u0, ArbTools.abs_ubound(x))*abspow(x, u0.p0))
        return 2L/(π*a0(u0, 0))*abspow(x, α - p)*S, 2L/(π*a0(u0, 0))*abspow(x, α - p)*S2
    end
end

"""
    T021(u0::FractionalKdVAnstaz{arb}, a::arb, x::arb)
Computes the (not yet existing) integral T_{0,2,1} from the paper.

The strategy for evaluation is the same as for T011 except that the
first term is singular and the last two are analytic and their Taylor
expansion is computed at `t = x`.

The integral that needs to be computed in this case is `∫_x^a (y -
x)^s*y^p dy` which is given by `x^(s + p + 1)*(Γ(1 + s)*Γ(-1 - s -
p)/Γ(-p) - B(-1 - s - p, 1 + s; x/a)` where B(a, b; z) is the
incomplete Beta-function. For `p = 1` we instead use the expression
`(a - x)^(1+s)*((a - x)/(2 + s) + x/(1 + s))`.

If `x` is equal or very close to π (determined by `ϵ`) then the Taylor
expansion gives a very poor approximation for `Ci(x + y, -α)`. In this
case we make use of the fact that it's 2π periodic and even, so that
`Ci(x + y, -α) = Ci(x + y - 2π, -α) = Ci(2π - (x + y), -α)`, to be
able to use the asymptotic expansion instead. That gives us the
integral ```∫_x^a (2π - (x + y))^s*t^p dy` which is given by `(2π -
x)^(1 + p + s)*(B(1 + p, 1 + s, a/(2π - x)) - B(1 + p, 1 + s, x/(2π -
x)))```. The value of `x/(2π - x)` will always be less than or equal
to 1 for `x` less than or equal to π, however due to overestimation
the enclosing ball might contain values greater than one, we therefore
have to use `beta_inc_zeroone` to be able to get finite results in
that case.
"""
T021(u0::FractionalKdVAnsatz{arb}, a, x; kwargs...) = T021(u0, Ball(), a, x; kwargs...)

function T021(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ2::arb = parent(u0.α)(1e-4),
              ϵ::arb = parent(u0.α)(1e-1),
              N::Integer = 3,
              )
    return x -> begin
        T021(u0, Ball(), ArbTools.ubound(x + δ2), x, ϵ = ϵ, N = N)
    end
end

function T021(u0::FractionalKdVAnsatz{arb},
              ::Ball,
              a::arb,
              x::arb;
              ϵ::arb = parent(u0.α)(1e-1),
              N::Integer = 3,
              )
    Γ = Nemo.gamma
    α = u0.α
    δ2 = a - x
    π = parent(α)(pi)

    PP = ArbPolyRing(parent(α), :x)

    # Determine if the asymptotic expansion or the Taylor
    # expansion should be used for the second term
    use_asymptotic = π - x < ϵ

    # Analytic terms
    (P, E) = taylor_with_error(x, setunion(x, a), N) do y
        if !use_asymptotic
            return Ci(x + y, -α) - 2Ci(y, -α)
        else
            return -2Ci(y, -α)
        end
    end
    P_restterm = ball(zero(α), E*δ2^N)

    # Singular term
    M = div(N, 2) + 1
    (C, e, P2, P2_E) = Ci_expansion(δ2, -α, M)
    P2_restterm = P2_E*(δ2)^(2M)

    # Compute the integral
    res = zero(α)
    # Integrate the singular term
    # Using ∫_x^a |x - y|^s*y^p dy = ∫_x^a (y - x)^s*y^p dy
    if u0.p == 1
        res += C*abspow(δ2, 1 + e)*(δ2/(2 + e) + x/(1 + e))
    else
        res += C*x^(1 + e + u0.p)*(
            Γ(1 + e)*Γ(-1 - e - u0.p)/Γ(-u0.p)
            - beta_inc_zeroone(-1 -e - u0.p, 1 + e, x/a)
        )
    end

    # Integrate the analytic terms
    full_series = P + P2
    for i = 0:N-1
        if u0.p == 1
            res += full_series[i]*δ2^(1 + i)*(
                δ2/(2 + i) + x/(1 + i)
            )
        else
            res += full_series[i] * x^(1 + i + u0.p) * (
                Γ(parent(α)(1 + i))*Γ(-1 -i - u0.p)/Γ(-u0.p)
                - beta_inc(-1 - i - u0.p, parent(α)(1 + i), x/a)
            )
        end
    end

    # Add the error term
    res += δ2*(P_restterm + P2_restterm)

    if use_asymptotic
        # Handle asymptotic expansion of Ci(x + y, -α)
        # The furthest away from 2π we are is at y = x
        (C, e, P3, P3_E) = Ci_expansion(2π - 2x, -α, M)
        P3_restterm = P2_E*(2π - 2x)^(2M)

        # Add the singular part to the integral
        res += C*(2π - x)^(1 + e + u0.p)*(
            beta_inc_zeroone(1 + u0.p, 1 + e, a/(2π - x))
            - beta_inc_zeroone(1 + u0.p, 1 + e, x/(2π - x))
        )

        for i = 0:2:N-1
            # Only even terms
            res += P3[i]*(2π - x)^(1 + i + u0.p)*(
            beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), a/(2π - x))
            - beta_inc_zeroone(1 + u0.p, parent(α)(1 + i), x/(2π - x))
        )
        end

        # Add error term
        res += δ2*P3_restterm
    end

    # Prove: that the expression inside the absolute value of the
    # integrand is positive
    return res/(parent(u0.α)(π)*u0.w(x)*u0(x))
end

"""
    T022(u0::FractionalKdVAnsatz{arb}, a::arb, x::arb)
Computes the (not yet existing) integral T_{0,2,2} from the paper.
)
"""
T022(u0::FractionalKdVAnsatz{arb}, a, x; kwargs...) = T022(u0, Ball(), a, x; kwargs...)

function T022(u0::FractionalKdVAnsatz{arb},
              ::Ball;
              δ2::arb = parent(u0.α)(1e-4),
              rtol = -1.0,
              atol = -1.0,
              show_trace = false,
              )
    return x -> begin
        T022(u0, Ball(), ArbTools.ubound(x + δ2), x, rtol = rtol, atol = atol, show_trace = show_trace)
    end
end

function T022(u0::FractionalKdVAnsatz{arb},
              ::Ball,
              a::arb,
              x::arb;
              rtol = -1.0,
              atol = -1.0,
              show_trace = false,
              )
    CC = ComplexField(prec(parent(u0.α)))
    mα = CC(-u0.α)
    a = CC(a)
    b = CC(parent(u0.α)(π))

    # PROVE: That there are no branch cuts that interact with the
    # integral
    F(y, analytic = false) = begin
        if isreal(y)
            res = CC(Ci(x - real(y), real(mα)) + Ci(x + real(y), real(mα)) - 2Ci(real(y), real(mα)))
        else
            res = Ci(x - y, mα) + Ci(x + y, mα) - 2Ci(y, mα)
        end
        return ArbTools.real_abs(res, analytic = analytic)*y^u0.p
    end
    res = real(ArbTools.integrate(CC, F, a, b,
                                  rel_tol = rtol,
                                  abs_tol = atol,
                                  eval_limit = 2000,
                                  verbose = Int(show_trace),
                                  ))

    return res/(parent(u0.α)(π)*u0.w(x)*u0(x))

end
