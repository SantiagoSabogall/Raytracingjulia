module Kerr  

export KerrBH

import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega

using LinearAlgebra, StaticArrays

mutable struct KerrBH
    a::Float64
    EH::Float64
    ISCOco::Float64
    ISCOcounter::Float64

    function KerrBH(a::Float64)
        EH = 1 + sqrt(1 - a^2)

        Z1 = 1 + (1 - a^2)^(1/3) *
            ((1 + a)^(1/3) + (1 - a)^(1/3))
        Z2 = sqrt(3a^2 + Z1^2)

        ISCOco = 3 + Z2 - sqrt((3 - Z1)*(3 + Z1 + 2Z2))
        ISCOcounter = 3 + Z2 + sqrt((3 - Z1)*(3 + Z1 + 2Z2))

        new(a, EH, ISCOco, ISCOcounter)
    end
end

# ============================
# Angular velocity
# ============================
function Omega(b::KerrBH, r::Real; corotating::Bool=true)
    corotating ? 1 / (r *sqrt(r) + b.a) :
                 -1 / (r * sqrt(r) - b.a)
end

# ============================
# Kerr metric
# ============================
function metric(b::KerrBH, x::SVector{4,Float64})
    r = x[2]
    θ = x[3]

    r2 = r*r
    a2 = b.a*b.a

    sinθ, cosθ = sincos(θ)
    sinθ2 = sinθ*sinθ
    cosθ2 = cosθ*cosθ

    Δ = r2 - 2r + a2
    Σ = r2 + a2*cosθ2

    invΣ = 1.0 / Σ

    g_tt   = -(1 - 2r*invΣ)
    g_rr   = Σ/Δ
    g_θθ   = Σ
    g_φφ   = (r2 + a2 + 2a2*r*sinθ2*invΣ)*sinθ2
    g_tφ   = -2*b.a*r*sinθ2*invΣ

    return g_tt, g_rr, g_θθ, g_φφ, g_tφ
end

# ============================
# Inverse metric
# ============================
function inverse_metric(b::KerrBH, x::AbstractVector)
    r = x[2]
    θ = x[3]

    r2 = r^2
    a2 = b.a^2
    sinθ2 = sin(θ)^2

    Δ = r2 - 2r + a2
    Σ = r2 + a2*cos(θ)^2
    A = (r2 + a2)^2 - Δ*a2*sinθ2

    gtt   = -A/(Δ*Σ)
    grr   = Δ/Σ
    gθθ   = 1/Σ
    gφφ   = (Δ - a2*sinθ2)/(Δ*Σ*sinθ2)
    gtφ   = -2*b.a*r/(Δ*Σ)

    return gtt, grr, gθθ, gφφ, gtφ
end

# ============================
# Photon geodesics (Hamiltonian)
# ============================
using StaticArrays

# 1. Cambiamos la firma para que sea (u, p, t) como espera DiffEq
# u = vector de estado (q), p = parámetros (el objeto KerrBH), t = parámetro afín (λ)
function geodesics(q, b::KerrBH, λ)
    # Desempaquetado con @inbounds para evitar chequeos de límites
    @inbounds begin
        r  = q[2]
        θ  = q[3]
        kt = q[5] 
        kr = q[6]
        kθ = q[7]
        kφ = q[8]
    end

    # Pre-cálculos trigonométricos (usando sincos para velocidad)
    sinθ, cosθ = sincos(θ)
    sinθ2 = sinθ * sinθ
    cosθ2 = cosθ * cosθ
    
    # Manejo de singularidad en los polos (evita división por cero)
    inv_sinθ2 = 1.0 / (sinθ2 + 1e-14)
    inv_sinθ  = 1.0 / (sinθ + 1e-14)

    r2 = r * r
    a  = b.a
    a2 = a * a # Idealmente b.a2 si lo añades al struct
    
    # Geometría base
    Σ = r2 + a2 * cosθ2
    invΣ = 1.0 / Σ
    invΣ2 = invΣ * invΣ
    
    Δ = r2 - 2.0*r + a2
    invΔ = 1.0 / Δ

    # Hamiltoniano y términos auxiliares
    W = -kt * (r2 + a2) - a * kφ
    
    # partΞ optimizado
    partΞ = r2 + (kφ + a * kt)^2 + a2 * (1.0 + kt^2) * cosθ2 + (kφ^2 * cosθ2 * inv_sinθ2)
    Ξ = W^2 - Δ * partΞ

    # Derivadas de Xi
    dΞdE = 2.0 * W * (r2 + a2) + 2.0 * a * Δ * (kφ + a * kt * sinθ2)
    dΞdL = -2.0 * a * W - 2.0 * a * kt * Δ - 2.0 * kφ * Δ * inv_sinθ2
    
    # OJO: Revisa este término dΞdr en tu física, aquí lo mantenemos optimizado
    dΞdr = -4.0 * r * kt * W - 2.0 * (r - 1.0) * partΞ - 2.0 * r * Δ 

    # Coeficientes A, B, C con multiplicaciones en vez de divisiones
    common_term = Ξ * invΔ * invΣ
    dAdr = (r - 1.0) * invΣ - (r * Δ) * invΣ2
    dBdr = -r * invΣ2
    dCdr = (dΞdr * 0.5 * invΔ * invΣ) - (common_term * (r - 1.0) * invΔ) - (r * Ξ * invΔ * invΣ2)

    auxθ = a2 * cosθ * sinθ
    dAdθ = Δ * auxθ * invΣ2
    dBdθ = auxθ * invΣ2
    
    # dCdθ optimizado
    term_trig_C = (kφ^2 * cosθ * inv_sinθ2 * inv_sinθ)
    dCdθ = ((1.0 + kt^2) * auxθ + term_trig_C) * invΣ + (common_term * auxθ * invΣ)

    # Ecuaciones diferenciales (Sistema de 8 ODEs)
    dtdλ = dΞdE * 0.5 * invΔ * invΣ
    drdλ = (Δ * invΣ) * kr
    dθdλ = kθ * invΣ
    dφdλ = -dΞdL * 0.5 * invΔ * invΣ

    dk_t  = 0.0
    dk_r  = -dAdr * kr^2 - dBdr * kθ^2 + dCdr
    dk_th = -dAdθ * kr^2 - dBdθ * kθ^2 + dCdθ
    dk_ph = 0.0

    # 2. RETORNO CRÍTICO: Usamos SVector para evitar alocaciones en el heap
    return SVector{8, Float64}(dtdλ, drdλ, dθdλ, dφdλ, dk_t, dk_r, dk_th, dk_ph)
end

end # module Kerr