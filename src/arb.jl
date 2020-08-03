function zeta(s; d::Integer = 0)
    if typeof(s) == Int
        s = Float64(s)
    end
    RR = RealField(precision(BigFloat))
    PP = ArbPolyRing(RR, :x)
    s_poly = PP([RR(s), one(RR)])
    a = one(RR)
    res = PP()

    ccall((:arb_poly_zeta_series, Nemo.libarb), Cvoid,
          (Ref{arb_poly}, Ref{arb_poly}, Ref{arb}, Cint, Clong, Clong),
          res, s_poly, a, 0, d + 1, RR.prec)

    convert(typeof(s), coeff(res, d))
end

function stieltjes(type, n::Integer)
    CC = ComplexField(precision(BigFloat))
    res = zero(CC)
    a = one(CC)

    ccall((:acb_dirichlet_stieltjes, Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{fmpz}, Ref{acb}, Clong),
          res, fmpz(n), a, CC.prec)

    convert(type, real(res))
end
function stieltjes(n::Integer)
    stieltjes(Float64, n)
end

# Clausians
# The factor \((-1)^\beta\) used when computing the Clausians comes from
# the fact that we discard the minus sign coming from the derivative
# when defining \(Ci_s^{(1)}(x)\). We could choose to keep this minus
# sign but it should not matter to much in the end.
# \begin{equation}
#   \frac{d}{ds}Ci_{s}(x) = -\sum_{n = 1}^{\infty}\frac{cos(nx)}{n^{s}}\log(n) = -Ci_{s}^{(1)}(x).
# \end{equation}

# Old call method
function Li(z; s = 2, β::Integer = 0)
    CC = ComplexField(precision(BigFloat))
    PP = AcbPolyRing(CC, :x)

    s = PP([CC(s), one(CC)])
    z = CC(z)
    w = PP()

    ccall(("acb_poly_polylog_series", Nemo.libarb), Cvoid,
          (Ref{acb_poly}, Ref{acb_poly}, Ref{acb}, Clong, Clong),
          w, s, z, β + 1, precision(BigFloat))

    coeff(w, β)*factorial(β)
end

function Ci(x; s = 2, β::Integer = 0)
    if x == 0
        return zeta(s, d = β)*(-1)^β
    end

    CC = ComplexField(precision(BigFloat))
    res = real(Li(exp(CC(zero(x), x)), s = s, β = β))

    # We need (-1)^β here to correct for the minus-sign that
    # we get when computing the derivative
    convert(float(typeof(x)), res*(-1)^β)
end

function Si(x; s = 2, β = 0)
    CC = ComplexField(precision(BigFloat))
    res = imag(Li(exp(CC(zero(x), x)), s = s, β = β))

    # We need (-1)^β here to correct for the minus-sign that
    # we get when computing the derivative
    convert(float(typeof(x)), res*(-1)^β)
end

# New call method
function Li(z::acb, s::acb)
    res = parent(z)()
    ccall(("acb_polylog", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Clong), res, s, z, prec(parent(z)))
    res
end

function Li(z::acb, s::Integer)
    res = parent(z)()
    ccall(("acb_polylog_si", Nemo.libarb), Cvoid,
          (Ref{acb}, Clong, Ref{acb}, Clong), res, s, z, prec(parent(z)))
    res
end

function Li(z::T, s) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    res = real(Li(CC(z), CC(s)))
    convert(float(T), res)
end

function Ci(x::acb, s::acb)
    im = x.parent(0, 1)
    0.5(Li(exp(im*x), s) + Li(exp(-im*x), s))
end

function Ci(x::arb, s::arb)
    CC = ComplexField(precision(BigFloat))
    real(Li(exp(CC(zero(x), x)), CC(s)))
end

function Ci(x::arb, s::Integer)
    CC = ComplexField(precision(BigFloat))
    real(Li(exp(CC(zero(x), x)), s))
end

function Ci(x::T, s) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = real(Li(CC(real(z), imag(z)), CC(s)))
    convert(float(T), res)
end

function Ci(x::T, s::Integer) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = real(Li(CC(real(z), imag(z)), s))
    convert(float(T), res)
end

function Si(x::acb, s::acb)
    im = x.parent(0, 1)
    0.5(Li(exp(im*x), s) - Li(exp(-im*x), s))
end

function Si(x::arb, s::arb)
    CC = ComplexField(precision(BigFloat))
    imag(Li(exp(CC(zero(x), x)), CC(s)))
end

function Si(x::arb, s::Integer)
    CC = ComplexField(precision(BigFloat))
    imag(Li(exp(CC(zero(x), x)), s))
end

function Si(x::T, s) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = imag(Li(CC(real(z), imag(z)), CC(s)))
    convert(float(T), res)
end

function Si(x::T, s::Integer) where {T <: Real}
    CC = ComplexField(precision(BigFloat))
    z = exp(im*x)
    res = imag(Li(CC(real(z), imag(z)), s))
    convert(float(T), res)
end

function LCi(x::arb; s = 2, tanh_max = 20)
    s = parent(x)(s) + 1/2
    res = Ci(x, s = s)

    for n in 1:tanh_max
        res += (sqrt(tanh(parent(x)(n))) - 1)*cos(n*x)/(parent(x)(n))^s
    end

    tail_error = (sqrt(tanh(parent(x)(tanh_max))) - 1)*zeta(s)

    res + Nemo.ball(zero(x), one(x))*tail_error
end

function LCi(x; s = 2, tanh_max = 20)
    s = s + 1/2
    res = Ci(x, s = s)

    for n in 1:tanh_max
        res += (sqrt(tanh(n)) - 1)*cos(n*x)/n^s
    end

    res
end

function LSi(x::arb; s = 2, tanh_max = 20)
    s = parent(x)(s) + 1/2
    res = Si(x, s = s)

    for n in 1:tanh_max
        res += (sqrt(tanh(parent(x)(n))) - 1)*sin(n*x)/(parent(x)(n))^s
    end

    tail_error = (sqrt(tanh(parent(x)(tanh_max))) - 1)*zeta(s)

    res + Nemo.ball(zero(x), one(x))*tail_error
end

function LSi(x; s = 2, tanh_max = 20)
    s = s + 1/2
    res = Si(x, s = s)

    for n in 1:tanh_max
        res += (sqrt(tanh(n)) - 1)*sin(n*x)/n^s
    end

    res
end

function hypbeta(a::arb, b::arb, z::arb)
    res = parent(z)()
    ccall(("arb_hypgeom_beta_lower", Nemo.libarb), Cvoid,
          (Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Cint, Clong), res, a, b, z, 0, prec(parent(z)))
    return res
end

function hypbeta(a::acb, b::acb, z::acb)
    res = parent(z)()
    ccall(("acb_hypgeom_beta_lower", Nemo.libarb), Cvoid,
          (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Cint, Clong), res, a, b, z, 0, prec(parent(z)))
    return res
end

function hypbeta(a, b, z::T) where {T}
    CC = ComplexField(precision(BigFloat))
    res = hypbeta(CC(a), CC(b), CC(z))
    accuracy = ArbTools.rel_accuracy_bits(real(res)) < 53 || ArbTools.rel_accuracy_bits(imag(res))
    if accuracy < 53
        @warn "Unsufficient precision $accuracy in hypbeta for z = $z"
        @show res
    end
    return convert(complex(float(T)), res)
end
