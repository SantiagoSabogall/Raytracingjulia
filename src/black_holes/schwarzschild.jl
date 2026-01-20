module Schwarzschild


export Schwarzschild, Omega, metric, inverse_metric, geodesics

mutable struct Schwarzschild


    a::Float64
    EH::Float64
    ISCOco::Float64
    ISCOcounter::Float64

    function Schwarzschild(a::Float64=0.0)
        EH = 2.0
        ISCOco = 6.0
        ISCOcounter = 6.0
        new(a, EH, ISCOco, ISCOcounter)
    end
end


# Angular velocity Ω(r)
function Omega(b::Schwarzschild, r::Real; corotating::Bool=true)
    return 1 / r^(3/2)
end


# Metric components gμν
function metric(b::Schwarzschild, x::Vector{<:Real})
    # x = [t, r, θ, φ]
    r = x[2]
    θ = x[3]

    g_tt = -(1 - 2/r)
    g_rr = -1 / g_tt
    g_thth = r^2
    g_phph = (r * sin(θ))^2
    g_tph = 0.0

    return [g_tt, g_rr, g_thth, g_phph, g_tph]
end


# Inverse metric g^μν
function inverse_metric(b::Schwarzschild, x::Vector{<:Real})
    r = x[2]
    θ = x[3]

    gtt = -1 / (1 - 2/r)
    grr = -1 / gtt
    gthth = 1 / r^2
    gphph = 1 / (r * sin(θ))^2
    gtph = 0.0

    return [gtt, grr, gthth, gphph, gtph]
end


# Geodesic equations dq/dλ = f(q, λ)
function geodesics(b::Schwarzschild, q::Vector{<:Real}, λ::Real)
    # q = [t, r, θ, φ, k_t, k_r, k_th, k_φ]
    t, r, θ, φ = q[1:4]
    k_t, k_r, k_th, k_φ = q[5:8]

    sinθ = sin(θ)
    f = 1 - 2/r

    # Coordinate derivatives
    dtdλ = -k_t / f
    drdλ = f * k_r
    dθdλ = k_th / r^2
    dφdλ = k_φ / (r^2 * sinθ^2)

    # Momentum derivatives
    dk_tdl = 0.0
    dk_rdl = - (k_t / (r - 2))^2 - (k_r / r)^2 + k_th^2 / r^3 + k_φ^2 / (r^3 * sinθ^2)
    dk_thdl = (cos(θ) / sinθ^3) * (k_φ / r)^2
    dk_phdl = 0.0

    return [dtdλ, drdλ, dθdλ, dφdλ,
            dk_tdl, dk_rdl, dk_thdl, dk_phdl]
end



end