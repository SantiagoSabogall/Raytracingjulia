"""
    KerrScalaron

Kerr-scalaron (KS) metric for null geodesic ray-tracing in f(R) gravity.

Reference:
  Paul, Bhattacharjee & Kalita (2024), ApJ 964:127
  "Kerr-scalaron Metric and Astronomical Consequences near the
   Galactic Center Black Hole"

─────────────────────────────────────────────────────────────────────────────
Physics summary
─────────────────────────────────────────────────────────────────────────────

Starting from the Schwarzschild-scalaron (SchS) seed metric (Eq. 2):

    ds² = −[1 − 2mα(r)/r] c²dt² + [1 − 2mα(r)/r]⁻¹ dr² + r²dΩ²

with the scalaron factor (Eq. 18, defining α(r)):

    α(r) = 1 + (1/3) exp(−Mψ r)          [units: Mψ in 1/length, r in length]

the Newman–Janis algorithm (Sec. 2.1–2.2) produces the KS line element
(Eq. 22):

    ds² = −[1 − 2m̃(r)r/ρ²] c²dt²
          + (4m̃(r)r sin²θ / ρ²) c dt dφ
          + (ρ²/Δ) dr²
          + ρ² dθ²
          + [r² + a² + 2a²m̃(r)r sin²θ/ρ²] sin²θ dφ²

where (geometrised units, M = G = c = 1):

    m̃(r) ≡ m · α(r) = m [1 + (1/3) exp(−Mψ r)]    effective mass
    ρ²   = r² + a² cos²θ
    Δ    = r² + a² − 2 m̃(r) r

The KS metric has the same tensorial structure as Kerr; only m → m̃(r)
everywhere.  Limits:
  • Mψ → ∞  ⟹  α(r) → 1  ⟹  standard Kerr
  • Mψ → ∞, a → 0  ⟹  Schwarzschild

─────────────────────────────────────────────────────────────────────────────
Derivatives needed for geodesics
─────────────────────────────────────────────────────────────────────────────

    dm̃/dr  = −(m Mψ/3) exp(−Mψ r)        ≡ dtm(r)
    d²m̃/dr² =  (m Mψ²/3) exp(−Mψ r)       ≡ d2tm(r)

    dΔ/dr  = 2r − 2m̃(r) − 2 dtm(r)·r     [differs from Kerr's 2r−2m]

─────────────────────────────────────────────────────────────────────────────
Units
─────────────────────────────────────────────────────────────────────────────
Geometrised: G = c = M = 1.
Coordinates: Boyer–Lindquist (t, r, θ, φ).
Mψ has dimensions of [1/length] = [1/mass] in these units.
"""
module KerrScalaron

export KerrScalaronBH

import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega

using LinearAlgebra, StaticArrays

# ═══════════════════════════════════════════════════════════════════════════
# Struct
# ═══════════════════════════════════════════════════════════════════════════

"""
    KerrScalaronBH(a, Mpsi)

Rotating black hole in f(R) scalaron gravity (Kerr-scalaron metric).

# Fields
- `a`    : Spin parameter J/M (geometrised), 0 ≤ |a| ≤ 1.
- `Mpsi` : Scalaron mass Mψ in units of 1/M (geometrised).
           Physical range from the paper: 10⁻²²–10⁻¹⁶ eV,
           converted to geometrised units for Sgr A* (m ≈ 0.042 au):
           Mψ [au⁻¹] = Mψ [eV] / 8.18×10⁻²⁰.
           In units where M=1: Mψ_geom = Mψ [au⁻¹] × m [au].
           For numerical work with M=1, typical range Mpsi ∈ (0, ~100].
           Mpsi → ∞ recovers Kerr exactly.
- `EH`   : Outer event horizon radius (largest root of Δ=0), numerical.
- `ISCOco`      : Co-rotating ISCO radius (numerical).
- `ISCOcounter` : Counter-rotating ISCO radius (numerical).
"""
mutable struct KerrScalaronBH
    a    :: Float64   # spin parameter
    Mpsi :: Float64   # scalaron mass  (1/length, geometrised)
    EH   :: Float64   # event horizon radius
    ISCOco      :: Float64
    ISCOcounter :: Float64

    function KerrScalaronBH(a::Float64, Mpsi::Float64)
        # define effective mass and its derivatives (M = m = 1)
        mt(r)   =  1.0 + (1.0/3.0) * exp(-Mpsi * r)          # α(r), so m̃ = m·α = α since m=1
        dtm(r)  = -(Mpsi / 3.0)   * exp(-Mpsi * r)            # dm̃/dr
        # d²m̃/dr² not needed in constructor

        Delta(r) = r^2 + a^2 - 2.0 * mt(r) * r

        EH              = _find_horizon(Delta)
        ISCOco, ISCOcounter = _find_isco(a, mt, dtm, Delta, EH)

        new(a, Mpsi, EH, ISCOco, ISCOcounter)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════════

"""
    _find_horizon(Delta) → r_+

Locates the outermost positive root of Δ(r) = 0 via bisection.
Scans r ∈ (ε, 5] with 10 000 steps, then bisects each sign change.
Returns the largest root found (outer horizon).
"""
function _find_horizon(Delta::Function)
    r_min, r_max = 1e-4, 5.0
    n = 10_000
    rs = range(r_min, r_max; length = n)

    r_outer = NaN
    for i in 1:(n - 1)
        if Delta(rs[i]) * Delta(rs[i + 1]) < 0.0
            lo, hi = rs[i], rs[i + 1]
            for _ in 1:64          # 64 bisections → machine precision
                mid = 0.5 * (lo + hi)
                Delta(lo) * Delta(mid) ≤ 0.0 ? (hi = mid) : (lo = mid)
            end
            r_outer = 0.5 * (lo + hi)   # keep updating → last = outermost
        end
    end
    # fallback: Kerr outer horizon
    isnan(r_outer) && (r_outer = 1.0 + sqrt(max(0.0, 1.0 - (r_max / 5)^2)))
    return r_outer
end

"""
    _find_isco(a, mt, dtm, Delta, EH) → (r_co, r_counter)

Locates co- and counter-rotating ISCO radii as minima of the specific
energy E(r) for equatorial circular orbits in the KS metric.

For circular equatorial orbits the angular velocity is:
    Ω = [−∂_r g_tφ ± √((∂_r g_tφ)² − (∂_r g_tt)(∂_r g_φφ))] / ∂_r g_φφ
evaluated at θ = π/2 (sin θ = 1).

The specific energy E(r) has a minimum at the ISCO.
"""
function _find_isco(a::Float64, mt::Function, dtm::Function,
                    Delta::Function, EH::Float64)

    # Equatorial metric derivatives at θ = π/2  (ρ² = r²)
    function dg_equatorial(r)
        mr  = mt(r)
        dmr = dtm(r)
        r2  = r^2
        a2  = a^2
        Σ   = r2                       # ρ² at θ=π/2
        invΣ  = 1.0 / Σ
        invΣ2 = invΣ^2

        # d/dr of metric components at θ = π/2
        # g_tt = −(1 − 2m̃r/ρ²)
        dgtt   =  2.0 * (dmr * r + mr) * invΣ - 2.0 * mr * r * 2r * invΣ2

        # g_tφ = −2a m̃ r sin²θ / ρ²  →  sin²θ=1
        dgtphi = -2.0 * a * ((dmr * r + mr) * invΣ - mr * r * 2r * invΣ2)

        # g_φφ = (r²+a² + 2a²m̃r/ρ²) sin²θ  →  sin²θ=1
        dgphiphi = 2.0 * r + 2.0 * a2 * ((dmr * r + mr) * invΣ - mr * r * 2r * invΣ2)

        return dgtt, dgtphi, dgphiphi
    end

    # Keplerian angular velocity (sign=+1 co-rotating, −1 counter)
    function Omega_circ(r, sign)
        dgtt, dgtphi, dgphiphi = dg_equatorial(r)
        disc = dgtphi^2 - dgtt * dgphiphi
        disc ≤ 0.0 && return NaN
        return (-dgtphi + sign * sqrt(disc)) / dgphiphi
    end

    # Specific energy for circular orbit
    function E_circ(r, sign)
        Ω   = Omega_circ(r, sign)
        isnan(Ω) && return NaN
        mr  = mt(r)
        r2  = r^2
        a2  = a^2
        gtt   = -(1.0 - 2.0 * mr * r / r2)          # θ=π/2
        gtphi = -2.0 * a * mr * r / r2
        gphi  = (r2 + a2 + 2.0 * a2 * mr * r / r2)  # sin²θ=1
        denom2 = -(gtt + 2.0 * gtphi * Ω + gphi * Ω^2)
        denom2 ≤ 0.0 && return NaN
        return -(gtt + gtphi * Ω) / sqrt(denom2)     # E = −(g_tt + g_tφ Ω) uᵗ
    end

    # ISCO = location of minimum specific energy
    function find_one(sign)
        r_lo = EH * 1.001
        r_hi = 30.0
        n    = 12_000
        rs   = range(r_lo, r_hi; length = n)

        E_prev = E_circ(rs[1], sign)
        isnan(E_prev) && (E_prev = Inf)
        r_isco = NaN

        for i in 2:n
            Ei = E_circ(rs[i], sign)
            isnan(Ei) && continue
            if Ei > E_prev           # energy started rising → minimum was previous point
                r_isco = rs[i - 1]
                break
            end
            E_prev = Ei
        end
        # Kerr fallback
        isnan(r_isco) && (r_isco = sign > 0 ? 6.0 : 9.0)
        return r_isco
    end

    return find_one(+1.0), find_one(-1.0)
end

# ═══════════════════════════════════════════════════════════════════════════
# Keplerian angular velocity  Ω_K  (equatorial circular orbits)
# ═══════════════════════════════════════════════════════════════════════════

"""
    Omega(b::KerrScalaronBH, r; corotating=true)

Keplerian angular velocity for equatorial circular orbits in the KS metric.

Generalises the Kerr formula by replacing M with the effective mass m̃(r):

    Ω_K = [−∂_r g_tφ ± √((∂_r g_tφ)² − (∂_r g_tt)(∂_r g_φφ))] / ∂_r g_φφ

evaluated on the equatorial plane θ = π/2.
"""
function Omega(b::KerrScalaronBH, r::Real; corotating::Bool = true)
    mr  = 1.0 + (1.0 / 3.0) * exp(-b.Mpsi * r)
    dmr = -(b.Mpsi / 3.0)   * exp(-b.Mpsi * r)

    r2    = r * r
    a2    = b.a^2
    invr2 = 1.0 / r2

    # ∂_r of equatorial metric components (ρ²=r², dρ²/dr=2r)
    dgtt     =  2.0 * (dmr * r + mr) * invr2 - 4.0 * mr * r * r * invr2^2 / r2
    # simplify: 2(dmr·r + mr)/r² − 4mr/r²
    dgtt     =  2.0 * (dmr * r + mr) / r2 - 4.0 * mr / r2

    dgtphi   = -2.0 * b.a * ((dmr * r + mr) / r2 - 2.0 * mr / r2)
    dgphiphi =  2.0 * r   + 2.0 * a2 * ((dmr * r + mr) / r2 - 2.0 * mr / r2)

    disc  = dgtphi^2 - dgtt * dgphiphi
    disc  = max(disc, 0.0)
    sign  = corotating ? +1.0 : -1.0

    return (-dgtphi + sign * sqrt(disc)) / dgphiphi
end

# ═══════════════════════════════════════════════════════════════════════════
# Covariant metric  g_μν   (Eq. 22 of the paper, in Boyer–Lindquist)
# ═══════════════════════════════════════════════════════════════════════════

"""
    metric(b::KerrScalaronBH, x::SVector{4,Float64})

Covariant KS metric components in Boyer–Lindquist coordinates.
Returns `(g_tt, g_rr, g_θθ, g_φφ, g_tφ)`.

Identical in structure to Kerr but with the scalaron-dressed mass:
    m̃(r) = m [1 + (1/3) exp(−Mψ r)]
"""
function metric(b::KerrScalaronBH, x::SVector{4, Float64})
    r = x[2];  θ = x[3]

    sinθ, cosθ = sincos(θ)
    sinθ2 = sinθ * sinθ
    cosθ2 = cosθ * cosθ
    r2    = r * r
    a2    = b.a * b.a

    mr   = 1.0 + (1.0 / 3.0) * exp(-b.Mpsi * r)   # m̃(r),  m=1

    ρ2   = r2 + a2 * cosθ2
    Δ    = r2 + a2 - 2.0 * mr * r
    invρ2 = 1.0 / ρ2

    g_tt   = -(1.0 - 2.0 * mr * r * invρ2)
    g_rr   =  ρ2 / Δ
    g_θθ   =  ρ2
    g_φφ   = (r2 + a2 + 2.0 * a2 * mr * r * sinθ2 * invρ2) * sinθ2
    g_tφ   = -2.0 * b.a * mr * r * sinθ2 * invρ2

    return g_tt, g_rr, g_θθ, g_φφ, g_tφ
end

# ═══════════════════════════════════════════════════════════════════════════
# Contravariant metric  g^μν
# ═══════════════════════════════════════════════════════════════════════════

"""
    inverse_metric(b::KerrScalaronBH, x::AbstractVector)

Contravariant KS metric components.
Returns `(g^tt, g^rr, g^θθ, g^φφ, g^tφ)`.

Structure identical to Kerr with m → m̃(r):
    A = (r²+a²)² − Δ a² sin²θ
    g^tt = −A/(Δρ²),  g^rr = Δ/ρ²,  g^θθ = 1/ρ²
    g^φφ = (Δ − a²sin²θ)/(Δρ²sin²θ),  g^tφ = −2am̃r/(Δρ²)
"""
function inverse_metric(b::KerrScalaronBH, x::AbstractVector)
    r = x[2];  θ = x[3]

    sinθ2 = sin(θ)^2
    cosθ2 = cos(θ)^2
    r2    = r^2
    a2    = b.a^2

    mr   = 1.0 + (1.0 / 3.0) * exp(-b.Mpsi * r)

    ρ2   = r2 + a2 * cosθ2
    Δ    = r2 + a2 - 2.0 * mr * r
    A    = (r2 + a2)^2 - Δ * a2 * sinθ2

    invΔρ2 = 1.0 / (Δ * ρ2)

    gtt   = -A * invΔρ2
    grr   =  Δ / ρ2
    gθθ   =  1.0 / ρ2
    gφφ   = (Δ - a2 * sinθ2) * invΔρ2 / sinθ2
    gtφ   = -2.0 * b.a * mr * r * invΔρ2

    return gtt, grr, gθθ, gφφ, gtφ
end

# ═══════════════════════════════════════════════════════════════════════════
# Null geodesics — Hamiltonian formalism
# ═══════════════════════════════════════════════════════════════════════════

"""
    geodesics(q, b::KerrScalaronBH, λ)

Right-hand side of the 8-ODE Hamiltonian system for null geodesics.

State vector:
    q = (t, r, θ, φ, k_t, k_r, k_θ, k_φ)
where k_μ = p_μ are covariant momenta.  k_t and k_φ are constants of
motion (metric cyclic in t and φ) and their evolution equations return 0.

Hamiltonian:
    H = ½ g^μν k_μ k_ν = 0   (null condition)

Hamilton's equations:
    ẋ^μ   =  ∂H/∂k_μ  = g^μν k_ν
    k̇_μ   = −∂H/∂x^μ = −½ (∂_μ g^αβ) k_α k_β

─────────────────────────────────────────────────────────────────────────────
Difference from Kerr:
The scalaron factor α(r) = 1 + exp(−Mψr)/3 makes m̃(r) r-dependent,
introducing dm̃/dr and d²m̃/dr² terms in ∂_r g^μν.
∂_θ g^μν has identical structure to Kerr (ρ² still contains a²cos²θ),
evaluated with m̃(r) instead of M.
─────────────────────────────────────────────────────────────────────────────

Returns `SVector{8,Float64}`.
"""
function geodesics(q, b::KerrScalaronBH, λ)
    @inbounds begin
        r  = q[2];  θ  = q[3]
        kt = q[5];  kr = q[6];  kθ = q[7];  kφ = q[8]
    end

    sinθ, cosθ = sincos(θ)
    sinθ2   = sinθ * sinθ
    cosθ2   = cosθ * cosθ
    inv_sinθ2 = 1.0 / (sinθ2 + 1e-14)

    r2 = r * r
    a  = b.a
    a2 = a * a

    # ── KS effective mass and derivatives ───────────────────────────────
    exp_r  = exp(-b.Mpsi * r)
    mr     = 1.0 + (1.0 / 3.0) * exp_r               # m̃(r)
    dtm    = -(b.Mpsi / 3.0) * exp_r                   # dm̃/dr
    d2tm   =  (b.Mpsi^2 / 3.0) * exp_r                 # d²m̃/dr²

    # ── Structural scalars ──────────────────────────────────────────────
    ρ2    = r2 + a2 * cosθ2
    invρ2 = 1.0 / ρ2
    invρ4 = invρ2 * invρ2

    Δ     = r2 + a2 - 2.0 * mr * r
    invΔ  = 1.0 / Δ

    # dΔ/dr = 2r − 2m̃(r) − 2(dm̃/dr)r
    dΔdr  = 2.0 * r - 2.0 * mr - 2.0 * dtm * r

    # dρ²/dr = 2r,   dρ²/dθ = −2a²cosθ sinθ
    dρ2dr = 2.0 * r
    dρ2dθ = -2.0 * a2 * cosθ * sinθ

    # ── Contravariant metric (same as inverse_metric, inlined) ──────────
    A      = (r2 + a2)^2 - Δ * a2 * sinθ2
    invΔρ2 = invΔ * invρ2

    gtt   = -A * invΔρ2
    grr   =  Δ * invρ2
    gθθ   =  invρ2
    gφφ   = (Δ - a2 * sinθ2) * invΔρ2 * inv_sinθ2
    gtφ   = -2.0 * a * mr * r * invΔρ2

    # ── Velocities:  ẋ^μ = g^μν k_ν ────────────────────────────────────
    dtdλ  =  gtt * kt + gtφ * kφ
    drdλ  =  grr * kr
    dθdλ  =  gθθ * kθ
    dφdλ  =  gtφ * kt + gφφ * kφ

    # ════════════════════════════════════════════════════════════════════
    # ∂_r g^μν — needed for dk_r/dλ
    # ════════════════════════════════════════════════════════════════════
    #
    # All r-derivatives come from three sources:
    #   (i)   explicit r in metric expressions
    #   (ii)  Δ(r) through dΔ/dr
    #   (iii) m̃(r) through dm̃/dr
    #
    # ∂_r A = 4r(r²+a²) − (dΔ/dr) a² sin²θ
    dAdr = 4.0 * r * (r2 + a2) - dΔdr * a2 * sinθ2

    # ∂_r [Δ ρ²] = (dΔ/dr) ρ² + Δ (dρ²/dr)
    d_Dro2_dr = dΔdr * ρ2 + Δ * dρ2dr

    # ── ∂_r g^tt = −∂_r[A/(Δρ²)]
    #            = −[dAdr·(Δρ²) − A·∂_r(Δρ²)] / (Δρ²)²
    dgtt_dr = -(dAdr * (Δ * ρ2) - A * d_Dro2_dr) / (Δ * ρ2)^2

    # ── ∂_r g^rr = ∂_r[Δ/ρ²] = (dΔ/dr·ρ² − Δ·dρ²/dr) / ρ⁴
    dgrr_dr = (dΔdr * ρ2 - Δ * dρ2dr) * invρ4

    # ── ∂_r g^θθ = ∂_r[1/ρ²] = −(dρ²/dr) / ρ⁴
    dgθθ_dr = -dρ2dr * invρ4

    # ── ∂_r g^φφ = ∂_r[(Δ−a²sin²θ)/(Δρ²sin²θ)]
    #   N = Δ − a²sin²θ,  dN/dr = dΔ/dr
    #   D = Δ ρ² sin²θ,   dD/dr = (dΔ/dr·ρ² + Δ·dρ²/dr) sin²θ
    N_phi  = Δ - a2 * sinθ2
    dN_dr  = dΔdr
    D_phi  = Δ * ρ2 * sinθ2
    dD_dr  = (dΔdr * ρ2 + Δ * dρ2dr) * sinθ2
    dgφφ_dr = (dN_dr * D_phi - N_phi * dD_dr) / D_phi^2

    # ── ∂_r g^tφ = ∂_r[−2a m̃ r / (Δ ρ²)]
    #   P = −2a m̃ r,  dP/dr = −2a(dm̃/dr·r + m̃)
    #   Q = Δ ρ²,     dQ/dr = d_Dro2_dr
    P_tph  = -2.0 * a * mr * r
    dP_dr  = -2.0 * a * (dtm * r + mr)
    Q_tph  = Δ * ρ2
    dgtφ_dr = (dP_dr * Q_tph - P_tph * d_Dro2_dr) / Q_tph^2

    # ── dk_r/dλ = −½ Σ_{μ,ν} (∂_r g^μν) k_μ k_ν ──────────────────────
    dk_r = -0.5 * (dgtt_dr  * kt * kt
                 + 2.0 * dgtφ_dr * kt * kφ
                 + dgrr_dr  * kr * kr
                 + dgθθ_dr  * kθ * kθ
                 + dgφφ_dr  * kφ * kφ)

    # ════════════════════════════════════════════════════════════════════
    # ∂_θ g^μν — needed for dk_θ/dλ
    # ════════════════════════════════════════════════════════════════════
    #
    # Only ρ²(θ) = r² + a²cos²θ varies with θ.  Δ, m̃, A_part have
    # θ-dependence only through sin²θ and ρ².
    #
    # dρ²/dθ = −2a²cosθ sinθ  (computed above as dρ2dθ)
    # ∂_θ A = −Δ a² · 2 sinθ cosθ
    dAdθ = -Δ * a2 * 2.0 * sinθ * cosθ

    # ∂_θ[Δρ²] = Δ · dρ²/dθ
    d_Dro2_dθ = Δ * dρ2dθ

    dgtt_dθ = -(dAdθ * (Δ * ρ2) - A * d_Dro2_dθ) / (Δ * ρ2)^2

    dgrr_dθ = -Δ * dρ2dθ * invρ4

    dgθθ_dθ = -dρ2dθ * invρ4

    # ∂_θ g^φφ
    #   dN/dθ = −2a² sinθ cosθ
    #   dD/dθ = Δ [dρ²/dθ · sin²θ + ρ² · 2sinθcosθ]
    dN_dθ  = -2.0 * a2 * sinθ * cosθ
    dD_dθ  = Δ * (dρ2dθ * sinθ2 + ρ2 * 2.0 * sinθ * cosθ)
    dgφφ_dθ = (dN_dθ * D_phi - N_phi * dD_dθ) / D_phi^2

    # ∂_θ g^tφ = ∂_θ[−2am̃r/(Δρ²)]  →  only Q=Δρ² depends on θ
    dgtφ_dθ = -P_tph * d_Dro2_dθ / Q_tph^2

    # ── dk_θ/dλ = −½ Σ (∂_θ g^μν) k_μ k_ν ─────────────────────────────
    dk_th = -0.5 * (dgtt_dθ  * kt * kt
                  + 2.0 * dgtφ_dθ * kt * kφ
                  + dgrr_dθ  * kr * kr
                  + dgθθ_dθ  * kθ * kθ
                  + dgφφ_dθ  * kφ * kφ)

    return SVector{8, Float64}(dtdλ, drdλ, dθdλ, dφdλ,
                                0.0,  dk_r, dk_th, 0.0)
end

end  # module KerrScalaron