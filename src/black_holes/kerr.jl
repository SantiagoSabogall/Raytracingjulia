module Kerr  

export KerrBH

# IMPORTANTE: Usamos 'import' para poder extender las funciones del padre
import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega

using LinearAlgebra

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
    corotating ? 1 / (r^(3/2) + b.a) :
                 -1 / (r^(3/2) - b.a)
end

# ============================
# Kerr metric
# ============================
function metric(b::KerrBH, x::AbstractVector)
    r = x[2]
    θ = x[3]

    r2 = r^2
    a2 = b.a^2
    sinθ2 = sin(θ)^2

    Δ = r2 - 2r + a2
    Σ = r2 + a2*cos(θ)^2

    g_tt   = -(1 - 2r/Σ)
    g_rr   = Σ/Δ
    g_θθ   = Σ
    g_φφ   = (r2 + a2 + 2a2*r*sinθ2/Σ)*sinθ2
    g_tφ   = -2*b.a*r*sinθ2/Σ

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
function geodesics(q, b::KerrBH, λ)
    # q[1]=t, q[2]=r, q[3]=θ, q[4]=φ
    # q[5]=kt, q[6]=kr, q[7]=kθ, q[8]=kφ
    r  = q[2]
    θ  = q[3]
    kt = q[5] 
    kr = q[6]
    kθ = q[7]
    kφ = q[8]

    r2 = r^2
    a2 = b.a^2
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinθ2 = sinθ^2
    cosθ2 = cosθ^2

    Σ = r2 + a2*cosθ2
    Σ2 = Σ^2
    Δ = r2 - 2r + a2

    # W = -kt*(r2 + a2) - a*kφ  (kt es q[5] en Julia, q[4] en Py)
    W = -kt*(r2 + a2) - b.a*kφ
    
    # partXi en Python usa q[7] (kφ) y q[4] (kt)
    partΞ = r2 + (kφ + b.a*kt)^2 + a2*(1 + kt^2)*cosθ2 + (kφ^2 * cosθ2 / sinθ2)
    Ξ = W^2 - Δ*partΞ

    # Derivadas de Xi
    dΞdE = 2*W*(r2 + a2) + 2.0*b.a*Δ*(kφ + b.a*kt*sinθ2)
    dΞdL = -2*b.a*W - 2*b.a*kt*Δ - 2*kφ*Δ/sinθ2
    dΞdr = -4*r*kt*W - 2*(r - 1)*partΞ - 2*r*Δ # Revisa si este 2r*Δ en Py no debería ser Δ * (derivada de partXi)

    # Coeficientes A, B, C
    dAdr = (r - 1)/Σ - (r*Δ)/Σ2
    dBdr = -r/Σ2
    dCdr = dΞdr/(2*Δ*Σ) - (Ξ*(r - 1))/(Σ*Δ^2) - r*Ξ/(Δ*Σ2)

    auxθ = a2*cosθ*sinθ
    dAdθ = Δ*auxθ/Σ2
    dBdθ = auxθ/Σ2
    
    # dCdθ: Ojo al kφ^2 y kt^2
    dCdθ = ((1 + kt^2)*auxθ + (kφ^2 * cosθ / (sinθ2 * sinθ)))/Σ + (Ξ/(Δ*Σ2))*auxθ

    # Ecuaciones diferenciales
    dtdλ  = dΞdE/(2.0*Δ*Σ)
    drdλ  = (Δ/Σ)*kr
    dθdλ  = kθ/Σ
    dφdλ  = -dΞdL/(2.0*Δ*Σ)

    dk_t = 0.0
    dk_r = -dAdr*kr^2 - dBdr*kθ^2 + dCdr
    dk_th = -dAdθ*kr^2 - dBdθ*kθ^2 + dCdθ
    dk_ph = 0.0

    return [dtdλ, drdλ, dθdλ, dφdλ, dk_t, dk_r, dk_th, dk_ph]
end

end # module Kerr