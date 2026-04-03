"""
    KerrPFDM

Provides the metric and Hamiltonian for a Kerr black hole immersed in Perfect Fluid Dark Matter (PFDM).
"""
module KerrPFDM

export KerrPFDMBH


import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega

using LinearAlgebra
using Roots
using StaticArrays

# ============================================================
# Kerr–PFDM Black Hole struct 
# ============================================================
"""
    KerrPFDMBH(a::Float64, k::Float64)

Definition of a Kerr black hole immersed in PFDM.

# Fields
- `a`: Spin parameter, \$a = J/M\$.
- `k`: PFDM intensity parameter.
- `EH`: Radius of the Event horizon, solved numerically where \$\\Delta(r) = 0\$.
- `ISCOco` / `ISCOcounter`: ISCO radii.
"""
mutable struct KerrPFDMBH
    a::Float64
    k::Float64
    EH::Float64
    ISCOco::Union{Float64, Nothing}
    ISCOcounter::Union{Float64, Nothing}

    function KerrPFDMBH(a::Float64, k::Float64)

        # --- Event Horizon Structural constraint (Δ = 0) ---
        fΔ(r) = r^2 + a^2 - 2r*(1 - (k/2)*log(r/abs(k)))
        EH = find_zero(fΔ, 1 + sqrt(1 - a^2))

        # --- ISCO domains: strict analytical existence not guaranteed under PFDM constraints ---
        ISCOco = nothing
        ISCOcounter = nothing

        new(a, k, EH, ISCOco, ISCOcounter)
    end
end
# ============================================================
# Effective PFDM dynamical mass profile
# ============================================================
"""
    m_eff(r, k)

Effective mass profile of the black hole with surrounding PFDM:
\$m_{eff}(r) = 1 - \\frac{k}{2} \\log\\left(\\frac{r}{|k|}\\right)\$.
"""
m_eff(r, k) = 1 - (k/2)*log(max(r, 1e-8)/abs(k))

"""
    dm_dr(r, k)

Radial derivative of the effective mass:
\$\\frac{dm_{eff}}{dr} = -\\frac{k}{2r}\$.
"""
dm_dr(r, k) = -k/(2r)

# ============================
# Angular velocity Ω
# ============================
"""
    Omega(b::KerrPFDMBH, r::Real; corotating::Bool=true)

Calculates the Keplerian angular velocity for a circular orbit in Kerr-PFDM.
"""
function Omega(b::KerrPFDMBH, r::Real; corotating::Bool=true)
    m  = m_eff(r, b.k)
    dm = dm_dr(r, b.k)

    fac = sqrt(m + r*dm)

    corotating ?
        fac / (r^(3/2) + b.a*fac) :
       -fac / (r^(3/2) - b.a*fac)
end

# ============================
# Kerr–PFDM metric
# ============================
"""
    metric(b::KerrPFDMBH, x::AbstractVector)

Evaluates the Kerr-PFDM covariant metric tensor \$g_{\\mu\\nu}\$.
"""
function metric(b::KerrPFDMBH, x::AbstractVector)
    r = x[2]
    θ = x[3]

    r2 = r^2
    a2 = b.a^2
    sinθ2 = sin(θ)^2

    m = m_eff(r, b.k)

    Δ = r2 - 2*m*r + a2
    Σ = r2 + a2*cos(θ)^2

    g_tt = -(1 - 2*m*r/Σ)
    g_rr = Σ/Δ
    g_θθ = Σ
    g_φφ = (r2 + a2 + 2a2*m*r*sinθ2/Σ)*sinθ2
    g_tφ = -2*b.a*m*r*sinθ2/Σ

    return g_tt, g_rr, g_θθ, g_φφ, g_tφ
end

# ============================
# Inverse metric
# ============================
"""
    inverse_metric(b::KerrPFDMBH, x::AbstractVector)

Evaluates the Kerr-PFDM contravariant metric tensor \$g^{\\mu\\nu}\$.
"""
function inverse_metric(b::KerrPFDMBH, x::AbstractVector)
    r = x[2]
    θ = x[3]

    r2 = r^2
    a2 = b.a^2
    sinθ2 = sin(θ)^2

    m = m_eff(r, b.k)

    Δ = r2 - 2*m*r + a2
    Σ = r2 + a2*cos(θ)^2
    A = (r2 + a2)^2 - Δ*a2*sinθ2

    gtt = -A/(Δ*Σ)
    grr = Δ/Σ
    gθθ = 1/Σ
    gφφ = (Δ - a2*sinθ2)/(Δ*Σ*sinθ2)
    gtφ = -2*b.a*m*r/(Δ*Σ)

    return gtt, grr, gθθ, gφφ, gtφ
end

# ============================
# Photon geodesics (Hamiltonian)
# ============================
"""
    geodesics(q, b::KerrPFDMBH, λ)

Evaluates the null geodesics 8-ODE system under the PFDM Hamiltonian constraints.
"""
function geodesics(q, b::KerrPFDMBH, λ)
    # q = [t, r, θ, φ, kt, kr, kθ, kφ]
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

    m  = m_eff(r, b.k)
    dm = dm_dr(r, b.k)

    Σ  = r2 + a2*cosθ2
    Σ2 = Σ^2
    Δ  = r2 - 2*m*r + a2

    # Hamiltonian structure bounds (identical to Kerr metric framework)
    W = -kt*(r2 + a2) - b.a*kφ

    partΞ =
        r2 +
        (kφ + b.a*kt)^2 +
        a2*(1 + kt^2)*cosθ2 +
        (kφ^2 * cosθ2 / sinθ2)

    Ξ = W^2 - Δ*partΞ

    # First-order ODE geometrical derivatives
    dΞdE = 2*W*(r2 + a2) + 2*b.a*Δ*(kφ + b.a*kt*sinθ2)
    dΞdL = -2*b.a*W - 2*b.a*kt*Δ - 2*kφ*Δ/sinθ2
    dΞdr = -4*r*kt*W - 2*(r - m - r*dm)*partΞ - 2*r*Δ

    dAdr = (r - m - r*dm)/Σ - (r*Δ)/Σ2
    dBdr = -r/Σ2
    dCdr = dΞdr/(2*Δ*Σ) -
           Ξ*(r - m - r*dm)/(Σ*Δ^2) -
           r*Ξ/(Δ*Σ2)

    auxθ = a2*cosθ*sinθ
    dAdθ = Δ*auxθ/Σ2
    dBdθ = auxθ/Σ2
    dCdθ =
        ((1 + kt^2)*auxθ +
        (kφ^2*cosθ/(sinθ2*sinθ)))/Σ +
        (Ξ/(Δ*Σ2))*auxθ

    # Evaluated equations of motion
    dtdλ  = dΞdE/(2*Δ*Σ)
    drdλ  = (Δ/Σ)*kr
    dθdλ  = kθ/Σ
    dφdλ  = -dΞdL/(2*Δ*Σ)

    dk_t  = 0.0
    dk_r  = -dAdr*kr^2 - dBdr*kθ^2 + dCdr
    dk_θ  = -dAdθ*kr^2 - dBdθ*kθ^2 + dCdθ
    dk_φ  = 0.0

    return SVector{8, Float64}(dtdλ, drdλ, dθdλ, dφdλ, dk_t, dk_r, dk_θ, dk_φ)
end

end # module KerrPFDM