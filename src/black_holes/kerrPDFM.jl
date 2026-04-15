"""
    KerrPFDM

Provides the analytical definitions of the Kerr-PFDM metric and Hamiltonian
for tracing null geodesics in a rotating black hole spacetime surrounded by
Perfect Fluid Dark Matter (PFDM).

Based on:
  Huang et al. (2026), "Impact of Perfect Fluid Dark Matter on the
  Appearance of Rotating Black Hole", arXiv:2602.00025v1

The key modification with respect to standard Kerr is the replacement of
the constant BH mass M by a radially-dependent effective mass profile:

    m(r) = M - (k/2) * ln(r / |k|)

which modifies Δ and therefore every metric component that depends on it.
The parameter k (PFDM intensity, k > 0 → attractive) is a free parameter.
In the limit k → 0 the spacetime reduces exactly to Kerr.

Units: geometrised (G = c = M = 1).
Coordinates: Boyer-Lindquist (t, r, θ, φ).
"""
module KerrPFDM

export KerrPFDMBH

import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega

using LinearAlgebra, StaticArrays

# ============================================================
# Struct
# ============================================================

"""
    KerrPFDMBH(a, k)

Rotating black hole embedded in Perfect Fluid Dark Matter.

# Fields
- `a`  : Spin parameter (J/M), |a| ≤ 1 for physical BH solutions.
- `k`  : PFDM intensity parameter (k > 0 attractive, k < 0 repulsive).
         The paper restricts to 0 < k ≤ 0.5 for well-defined ISCO.
- `EH` : Outer event horizon radius (largest root of Δ = 0), found numerically.
- `ISCOco`      : Co-rotating ISCO radius (numerical).
- `ISCOcounter` : Counter-rotating ISCO radius (numerical).

# Notes on horizon finding
Δ(r) = r² + a² - 2 m(r) r = r² + a² - 2Mr + k r ln(r/|k|)
has no closed-form roots when k ≠ 0, so we locate them with a simple
bisection search over r ∈ (0, 20].
"""
mutable struct KerrPFDMBH
    a  :: Float64
    k  :: Float64
    EH :: Float64
    ISCOco      :: Float64
    ISCOcounter :: Float64

    function KerrPFDMBH(a::Float64, k::Float64)
        # --- effective mass and its derivatives ---
        m(r)   = 1.0 - (k / 2.0) * log(r / abs(k))   # M = 1
        dm(r)  = -k / (2.0 * r)
        d2m(r) =  k / (2.0 * r^2)

        # --- metric functions ---
        Delta(r) = r^2 + a^2 - 2.0 * m(r) * r
        dDelta(r) = 2.0 * r - 2.0 * m(r) - 2.0 * dm(r) * r   # = 2r - 2m + k

        # ---- event horizon: outermost root of Δ = 0 ----
        EH = _find_horizon(Delta, a)

        # ---- ISCO (numerical) ----
        ISCOco, ISCOcounter = _find_isco(a, k, m, dm, d2m, Delta, dDelta, EH)

        new(a, k, EH, ISCOco, ISCOcounter)
    end
end

# ============================================================
# Internal: horizon finder (bisection)
# ============================================================

"""
    _find_horizon(Delta, a) → r_+

Finds the outermost positive root of Δ(r) = 0 by scanning for a sign
change and then bisecting.  Falls back to the Kerr value if no root is
found (should not happen for physical parameters).
"""
function _find_horizon(Delta::Function, a::Float64)
    # scan coarsely
    r_min = 1e-3
    r_max = 5.0        # r_+ ≤ 2 always for physical BH
    n_scan = 5000
    rs = range(r_min, r_max, length=n_scan)

    r_outer = NaN
    for i in 1:length(rs)-1
        if Delta(rs[i]) * Delta(rs[i+1]) < 0.0
            # bisect
            lo, hi = rs[i], rs[i+1]
            for _ in 1:60
                mid = 0.5*(lo + hi)
                Delta(lo) * Delta(mid) < 0.0 ? (hi = mid) : (lo = mid)
            end
            r_outer = 0.5*(lo + hi)   # keep the outermost crossing
        end
    end

    isnan(r_outer) && (r_outer = 1.0 + sqrt(max(0.0, 1.0 - a^2)))  # fallback Kerr
    return r_outer
end

# ============================================================
# Internal: ISCO finder
# ============================================================

"""
    _find_isco(a, k, m, dm, d2m, Delta, dDelta, EH) → (r_co, r_counter)

Solves the three simultaneous ISCO conditions:
    V(r,E,L) = 0,   ∂_r V = 0,   ∂²_r V = 0
by scanning for the equatorial effective potential  V(r,E,L)
in terms of circular-orbit energy/angular-momentum.

Uses the standard approach: for circular equatorial orbits
    E = E(r),  L = L(r)  from  V = ∂_r V = 0,
then ISCO is where ∂²_r V = 0.  We scan r from rH outward and locate
the first extremum of ∂²_r V.

For robustness, if the numerical search fails we return the Kerr ISCO.
"""
function _find_isco(a, k, m, dm, d2m, Delta, dDelta, EH)
    # We use the determinant-based criterion from the paper (Eq. 25 + marginal stability).
    # Circular equatorial orbit: uᶿ = 0, θ = π/2, gᵣᵣ uʳ = 0
    # Angular velocity Ω = dφ/dt from ∂_r V = 0:
    #
    #   Ω = [-∂_r g_tφ ± sqrt((∂_r g_tφ)² - (∂_r g_tt)(∂_r g_φφ))] / (∂_r g_φφ)
    #
    # Then E, L from normalisation.  ISCO: dE/dr = 0 (or equivalently d²V/dr² = 0).

    function gtt_eq(r)
        mr = m(r)
        Σ  = r^2   # θ = π/2
        return -(1.0 - 2mr*r/Σ)
    end
    function gtphi_eq(r)
        mr = m(r)
        Σ  = r^2
        return -2a*mr*r/Σ
    end
    function gphiphi_eq(r)
        mr  = m(r)
        Σ   = r^2
        return (r^2 + a^2 + 2*a^2*mr*r/Σ)   # sin²θ = 1
    end

    # Keplerian angular velocity (co: +, counter: -)
    function Omega_circ(r, sign)
        # From ∂_r V = 0 at equatorial circular orbit → standard formula
        mr  = m(r)
        dmr = dm(r)
        # ∂_r g_tt, ∂_r g_tφ, ∂_r g_φφ  at θ=π/2
        # Σ = r² → dΣ/dr = 2r
        Σ    = r^2
        invΣ = 1/Σ
        dΣ   = 2r

        dgtt   =  2*(dmr*r + mr)*invΣ - 2*mr*r*(-dΣ*invΣ^2)
        # = 2*(dmr*r + mr)/Σ + 2mr*r*dΣ/Σ²
        dgtt   =  2*(dmr*r + mr)/Σ  -  2*mr*r*dΣ/Σ^2

        dgtphi = -2a*(dmr*r + mr)/Σ + 2a*mr*r*dΣ/Σ^2

        dgphiphi = 2r + 2a^2*(dmr*r + mr)/Σ - 2a^2*mr*r*dΣ/Σ^2

        disc = dgtphi^2 - dgtt*dgphiphi
        disc < 0 && return NaN
        return (-dgtphi + sign*sqrt(disc)) / dgphiphi
    end

    # Specific energy E and angular momentum L for circular orbit
    function EL_circ(r, sign)
        Ω   = Omega_circ(r, sign)
        isnan(Ω) && return NaN, NaN
        gtt   = gtt_eq(r)
        gtphi = gtphi_eq(r)
        gphi  = gphiphi_eq(r)
        denom2 = -(gtt + 2gtphi*Ω + gphi*Ω^2)
        denom2 ≤ 0 && return NaN, NaN
        ut = 1.0 / sqrt(denom2)
        E  = -(gtt + gtphi*Ω) * ut
        L  =  (gphi*Ω + gtphi) * ut
        return E, L
    end

    # We look for dE/dr = 0 for each direction → ISCO
    function find_isco_one(sign)
        r_lo = EH * 1.01
        r_hi = 20.0
        n    = 8000
        rs   = range(r_lo, r_hi, length=n)

        E_prev, _ = EL_circ(rs[1], sign)
        isnan(E_prev) && (E_prev = Inf)

        r_isco = NaN
        for i in 2:n
            Ei, _ = EL_circ(rs[i], sign)
            isnan(Ei) && continue
            if Ei < E_prev
                # E is decreasing here: look for local minimum
                E_prev = Ei
            elseif !isnan(E_prev) && Ei > E_prev
                # passed a minimum → ISCO is at rs[i-1]
                r_isco = rs[i-1]
                break
            end
            E_prev = Ei
        end
        isnan(r_isco) && (r_isco = sign > 0 ? 6.0 : 9.0)   # Kerr fallback
        return r_isco
    end

    r_co      = find_isco_one(+1.0)
    r_counter = find_isco_one(-1.0)
    return r_co, r_counter
end

# ============================================================
# Angular velocity (Keplerian equatorial)
# ============================================================

"""
    Omega(b::KerrPFDMBH, r; corotating=true)

Keplerian angular velocity for equatorial circular orbits in Kerr-PFDM.
Generalises the Kerr formula by using the PFDM effective mass:

    Ω_K = [−∂_r g_tφ ± √((∂_r g_tφ)² − (∂_r g_tt)(∂_r g_φφ))] / ∂_r g_φφ

evaluated at θ = π/2.
"""
function Omega(b::KerrPFDMBH, r::Real; corotating::Bool=true)
    mr  = 1.0 - (b.k / 2.0) * log(r / abs(b.k))
    dmr = -b.k / (2.0 * r)
    a   = b.a

    Σ    = r^2
    dΣ   = 2r
    invΣ = 1.0 / Σ

    # ∂_r metric components at θ = π/2
    dgtt     =  2*(dmr*r + mr)*invΣ  -  2*mr*r*dΣ*invΣ^2
    dgtphi   = -2a*(dmr*r + mr)*invΣ +  2a*mr*r*dΣ*invΣ^2
    dgphiphi =  2r + 2a^2*(dmr*r + mr)*invΣ - 2a^2*mr*r*dΣ*invΣ^2

    disc = dgtphi^2 - dgtt*dgphiphi
    disc = max(disc, 0.0)

    sign = corotating ? +1.0 : -1.0
    return (-dgtphi + sign*sqrt(disc)) / dgphiphi
end

# ============================================================
# Covariant metric  g_μν
# ============================================================

"""
    metric(b::KerrPFDMBH, x::SVector{4,Float64})

Covariant Kerr-PFDM metric in Boyer-Lindquist coordinates.
Returns `(g_tt, g_rr, g_θθ, g_φφ, g_tφ)`.

The only change w.r.t. Kerr is that M is replaced by m(r):
    m(r) = 1 − (k/2) ln(r/|k|)   [M = 1]
"""
function metric(b::KerrPFDMBH, x::SVector{4,Float64})
    r = x[2];  θ = x[3]

    r2   = r * r
    a2   = b.a * b.a
    sinθ, cosθ = sincos(θ)
    sinθ2 = sinθ * sinθ
    cosθ2 = cosθ * cosθ

    # PFDM effective mass
    mr   = 1.0 - (b.k / 2.0) * log(r / abs(b.k))

    Σ    = r2 + a2 * cosθ2
    Δ    = r2 + a2 - 2.0 * mr * r
    invΣ = 1.0 / Σ

    g_tt   = -(1.0 - 2.0*mr*r*invΣ)
    g_rr   =  Σ / Δ
    g_θθ   =  Σ
    g_φφ   = (r2 + a2 + 2.0*a2*mr*r*sinθ2*invΣ) * sinθ2
    g_tφ   = -2.0 * b.a * mr * r * sinθ2 * invΣ

    return g_tt, g_rr, g_θθ, g_φφ, g_tφ
end

# ============================================================
# Contravariant metric  g^μν
# ============================================================

"""
    inverse_metric(b::KerrPFDMBH, x::AbstractVector)

Contravariant Kerr-PFDM metric.
Returns `(g^tt, g^rr, g^θθ, g^φφ, g^tφ)`.

Identical in structure to Kerr but with m(r) in place of M:
    A(r,θ) = (r²+a²)² − Δ a² sin²θ
"""
function inverse_metric(b::KerrPFDMBH, x::AbstractVector)
    r  = x[2];  θ = x[3]

    r2    = r^2
    a2    = b.a^2
    sinθ2 = sin(θ)^2
    cosθ2 = cos(θ)^2

    mr  = 1.0 - (b.k / 2.0) * log(r / abs(b.k))
    Δ   = r2 + a2 - 2.0 * mr * r
    Σ   = r2 + a2 * cosθ2
    A   = (r2 + a2)^2 - Δ * a2 * sinθ2

    invΔΣ = 1.0 / (Δ * Σ)

    gtt   = -A * invΔΣ
    grr   =  Δ / Σ
    gθθ   =  1.0 / Σ
    gφφ   =  (Δ - a2*sinθ2) * invΔΣ / sinθ2
    gtφ   = -2.0 * b.a * mr * r * invΔΣ

    return gtt, grr, gθθ, gφφ, gtφ
end

# ============================================================
# Null geodesics — Hamiltonian formalism
# ============================================================

"""
    geodesics(q, b::KerrPFDMBH, λ)

Right-hand side of the 8-ODE system for null geodesics in Kerr-PFDM.

State vector:  q = (t, r, θ, φ, k_t, k_r, k_θ, k_φ)
where k_μ = p_μ are the covariant momenta (k_t and k_φ are constants).

Hamiltonian:
    H = ½ g^μν k_μ k_ν = 0

Hamilton's equations:
    dxᵘ/dλ =  ∂H/∂k_μ = g^μν k_ν
    dk_μ/dλ = −∂H/∂xᵘ = −½ (∂_μ g^αβ) k_α k_β

---
Key difference from Kerr:
The PFDM mass profile m(r) introduces non-zero derivatives of m w.r.t. r:
    m'(r)  = −k/(2r)
    m''(r) =  k/(2r²)
These appear in ∂_r g^μν and therefore in dk_r/dλ.
∂_θ g^μν is identical in structure to Kerr (Σ still contains the
a²cos²θ term), except evaluated with m(r) instead of M.

---
Returns an `SVector{8,Float64}`.
"""
function geodesics(q, b::KerrPFDMBH, λ)
    @inbounds begin
        r  = q[2];  θ  = q[3]
        kt = q[5];  kr = q[6];  kθ = q[7];  kφ = q[8]
    end

    sinθ, cosθ = sincos(θ)
    sinθ2  = sinθ * sinθ
    cosθ2  = cosθ * cosθ
    inv_sinθ2 = 1.0 / (sinθ2 + 1e-14)

    r2  = r * r
    a   = b.a
    a2  = a * a

    # --- PFDM effective mass and derivatives ---
    mr   = 1.0 - (b.k / 2.0) * log(r / abs(b.k))   # m(r),  M=1
    dmr  = -b.k / (2.0 * r)                          # m'(r)
    d2mr =  b.k / (2.0 * r2)                         # m''(r)

    # --- Kerr-PFDM structure scalars ---
    Σ    = r2 + a2 * cosθ2
    invΣ = 1.0 / Σ
    invΣ2 = invΣ * invΣ

    Δ    = r2 + a2 - 2.0 * mr * r
    invΔ = 1.0 / Δ

    # dΔ/dr = 2r − 2m(r) − 2m'(r)r = 2r − 2m(r) + k
    dΔdr = 2.0*r - 2.0*mr - 2.0*dmr*r

    # dΣ/dr = 2r  (only r-dependence in Σ comes from r²)
    # dΣ/dθ = −2 a² cosθ sinθ
    dΣdr  =  2.0 * r
    dΣdθ  = -2.0 * a2 * cosθ * sinθ

    # --- Contravariant metric components (same expressions as inverse_metric) ---
    A     = (r2 + a2)^2 - Δ * a2 * sinθ2
    invΔΣ = invΔ * invΣ

    gtt   = -A * invΔΣ
    grr   =  Δ * invΣ
    gθθ   =  invΣ
    gφφ   = (Δ - a2*sinθ2) * invΔΣ * inv_sinθ2
    gtφ   = -2.0 * a * mr * r * invΔΣ

    # --- dx^μ/dλ = g^μν k_ν  (contravariant velocities) ---
    # Only gtt, gtφ, gφφ mix t and φ; grr, gθθ are diagonal.
    dtdλ  =  gtt * kt + gtφ * kφ
    drdλ  =  grr * kr
    dθdλ  =  gθθ * kθ
    dφdλ  =  gtφ * kt + gφφ * kφ

    # ---------------------------------------------------------------
    # dk_μ/dλ = −½ (∂_μ g^αβ) k_α k_β
    # Only μ = r and μ = θ are non-trivial (k_t, k_φ = const → 0).
    # ---------------------------------------------------------------

    # ---- Partial derivatives of A w.r.t. r and θ ----
    # A = (r²+a²)² − Δ a² sin²θ
    dAdr = 4.0*r*(r2 + a2) - dΔdr * a2 * sinθ2
    dAdθ = -Δ * a2 * 2.0 * sinθ * cosθ

    # ---- ∂_r g^tt = ∂_r[−A/(Δ Σ)] ----
    # Let f = A, g = ΔΣ  →  ∂_r(f/g) = (∂_r f · g − f · ∂_r g) / g²
    # ∂_r(ΔΣ) = dΔdr·Σ + Δ·dΣdr
    dgtt_dr = -(dAdr*(Δ*Σ) - A*(dΔdr*Σ + Δ*dΣdr)) / (Δ*Σ)^2

    # ---- ∂_r g^rr = ∂_r[Δ/Σ] ----
    dgrr_dr = (dΔdr*Σ - Δ*dΣdr) * invΣ2

    # ---- ∂_r g^θθ = ∂_r[1/Σ] ----
    dgθθ_dr = -dΣdr * invΣ2

    # ---- ∂_r g^φφ = ∂_r[(Δ − a² sin²θ)/(Δ Σ sin²θ)] ----
    # Numerator N = Δ − a² sin²θ,  denominator D = Δ Σ sin²θ
    # dN/dr = dΔdr,   dD/dr = dΔdr·Σ·sin²θ + Δ·dΣdr·sin²θ
    N_phi   = Δ - a2*sinθ2
    D_phi   = Δ * Σ * sinθ2
    dN_dr   = dΔdr
    dD_dr   = dΔdr * Σ * sinθ2 + Δ * dΣdr * sinθ2
    dgφφ_dr = (dN_dr * D_phi - N_phi * dD_dr) / D_phi^2

    # ---- ∂_r g^tφ = ∂_r[−2 a m(r) r / (Δ Σ)] ----
    # Numerator P = −2 a m(r) r,   dP/dr = −2a(m'(r)r + m(r)) = −2a(dmr·r + mr)
    # Denominator Q = Δ Σ,         dQ/dr = dΔdr·Σ + Δ·dΣdr
    P_tphi    = -2.0 * a * mr * r
    dP_dr     = -2.0 * a * (dmr*r + mr)
    Q_tphi    = Δ * Σ
    dQ_dr     = dΔdr * Σ + Δ * dΣdr
    dgtφ_dr   = (dP_dr * Q_tphi - P_tphi * dQ_dr) / Q_tphi^2

    # ---- dk_r/dλ = −½ (∂_r g^μν) k_μ k_ν ----
    dk_r = -0.5 * (dgtt_dr  * kt*kt
                 + 2.0*dgtφ_dr  * kt*kφ
                 + dgrr_dr  * kr*kr
                 + dgθθ_dr  * kθ*kθ
                 + dgφφ_dr  * kφ*kφ)

    # ---------------------------------------------------------------
    # ∂_θ derivatives
    # ---------------------------------------------------------------

    # ∂_θ A = −Δ a² 2 sinθ cosθ  (computed above as dAdθ)
    # ∂_θ(Δ Σ) = Δ · dΣdθ

    dgtt_dθ = -(dAdθ*(Δ*Σ) - A*(Δ*dΣdθ)) / (Δ*Σ)^2

    # ∂_θ g^rr = ∂_θ[Δ/Σ] = −Δ dΣdθ / Σ²
    dgrr_dθ = -Δ * dΣdθ * invΣ2

    # ∂_θ g^θθ = −dΣdθ / Σ²
    dgθθ_dθ = -dΣdθ * invΣ2

    # ∂_θ g^φφ : N = Δ−a²sin²θ, D = Δ Σ sin²θ
    # dN/dθ = −2a² sinθ cosθ
    # dD/dθ = Δ(dΣdθ · sin²θ + Σ · 2sinθcosθ)
    dN_dθ   = -2.0 * a2 * sinθ * cosθ
    dD_dθ   = Δ * (dΣdθ * sinθ2 + Σ * 2.0*sinθ*cosθ)
    dgφφ_dθ = (dN_dθ * D_phi - N_phi * dD_dθ) / D_phi^2

    # ∂_θ g^tφ = ∂_θ[−2am(r)r/(ΔΣ)] = −P_tphi · ∂_θ(ΔΣ) / (ΔΣ)²
    #          = P_tphi · Δ dΣdθ / (ΔΣ)²     (P_tphi already defined)
    dgtφ_dθ = -P_tphi * (Δ * dΣdθ) / Q_tphi^2

    # ---- dk_θ/dλ = −½ (∂_θ g^μν) k_μ k_ν ----
    dk_th = -0.5 * (dgtt_dθ  * kt*kt
                  + 2.0*dgtφ_dθ  * kt*kφ
                  + dgrr_dθ  * kr*kr
                  + dgθθ_dθ  * kθ*kθ
                  + dgφφ_dθ  * kφ*kφ)

    return SVector{8, Float64}(dtdλ, drdλ, dθdλ, dφdλ,
                                0.0,  dk_r, dk_th, 0.0)
end

end  # module KerrPFDM