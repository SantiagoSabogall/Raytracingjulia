"""
    Schwarzschild

Provides the exact analytical implementation of the Schwarzschild 
spacetime metric and its associated null geodesic integrations.
"""
module Schwarzschild

# Enforces generic extension of base Raytracing functions
import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega

using StaticArrays

export SchwarzschildBH

# Immutable struct: prevents heap allocation and boxing overhead during execution.
"""
    SchwarzschildBH(M::Float64=1.0)

Implements the geometry of a spherically symmetric, non-rotating black hole.

# Fields
- `a`: Spin parameter (exactly `0.0` for Schwarzschild).
- `EH`: Radius of the Event Horizon, \$r_+ = 2M\$.
- `ISCOco`: Radius of the Innermost Stable Circular Orbit for co-rotating particles, \$r_{ISCO} = 6M\$.
- `ISCOcounter`: Radius of the ISCO for counter-rotating particles, \$r_{ISCO} = 6M\$.
"""
struct SchwarzschildBH
    a::Float64
    EH::Float64
    ISCOco::Float64
    ISCOcounter::Float64

    function SchwarzschildBH(M::Float64=1.0)
        new(0.0, 2.0*M, 6.0*M, 6.0*M)
    end
end

"""
    Omega(b::SchwarzschildBH, r::Real; corotating::Bool=true)

Calculates the Keplerian angular velocity in the Schwarzschild spacetime:
\$\\Omega_K = \\frac{1}{r^{3/2}}\$ (assuming \$M=1\$).
"""
@inline function Omega(b::SchwarzschildBH, r::Real; corotating::Bool=true)
    return 1.0 / (r * sqrt(r))   # Equivalently r^1.5, evading floating point exponentiation overhead
end

"""
    metric(b::SchwarzschildBH, x::AbstractVector)

Evaluates the covariant Schwarzschild metric \$g_{\\mu\\nu}\$ using spherical 
coordinates \$x^\\mu = (t, r, \\theta, \\phi)\$:
\$ds^2 = -\\left(1 - \\frac{2M}{r}\\right) dt^2 + \\left(1 - \\frac{2M}{r}\\right)^{-1} dr^2 + r^2 d\\theta^2 + r^2 \\sin^2\\theta d\\phi^2\$

Returns `(g_tt, g_rr, g_thth, g_phph, g_tph)`.
"""
@inline function metric(b::SchwarzschildBH, x::AbstractVector)
    r  = x[2]; θ = x[3]
    f  = 1.0 - 2.0/r
    sinθ = sin(θ)
    r2   = r * r
    return -f, 1.0/f, r2, r2 * sinθ * sinθ, 0.0
end

"""
    inverse_metric(b::SchwarzschildBH, x::AbstractVector)

Evaluates the contravariant Schwarzschild metric \$g^{\\mu\\nu}\$.
Returns `(gtt, grr, gthth, gphph, gtph)`.
"""
@inline function inverse_metric(b::SchwarzschildBH, x::AbstractVector)
    r    = x[2]; θ = x[3]
    f    = 1.0 - 2.0/r
    sinθ = sin(θ)
    r2   = r * r
    return -1.0/f, f, 1.0/r2, 1.0/(r2 * sinθ * sinθ), 0.0
end

# Returns an SVector to strictly enforce zero memory allocations on the heap.
"""
    geodesics(q, b::SchwarzschildBH, λ)

Evaluates the right-hand side of the 8-ODE system for null geodesics.

Derived via the Hamiltonian framework:
\$\\mathcal{H} = \\frac{1}{2} \\left[ -\\left(1 - \\frac{2M}{r}\\right)^{-1} p_t^2 + \\left(1 - \\frac{2M}{r}\\right) p_r^2 + \\frac{p_\\theta^2}{r^2} + \\frac{p_\\phi^2}{r^2 \\sin^2\\theta} \\right] = 0\$

Hamilton's equations:
- \$\\dot{x}^\\mu = \\partial \\mathcal{H} / \\partial p_\\mu = g^{\\mu\\nu} p_\\nu\$
- \$\\dot{p}_\\mu = -\\partial \\mathcal{H} / \\partial x^\\mu = -\\frac{1}{2} \\partial_\\mu g^{\\alpha\\beta} p_\\alpha p_\\beta\$

Returns an `SVector{8}` of \$(\\dot{t}, \\dot{r}, \\dot{\\theta}, \\dot{\\phi}, \\dot{p}_t, \\dot{p}_r, \\dot{p}_\\theta, \\dot{p}_\\phi)\$.
"""
@inline function geodesics(q, b::SchwarzschildBH, λ)
    @inbounds begin
        r   = q[2]
        θ   = q[3]
        p_t = q[5]
        p_r = q[6]
        p_θ = q[7]
        p_φ = q[8]
    end

    sinθ, cosθ = sincos(θ)
    sinθ2 = sinθ * sinθ

    r2     = r * r
    f      = 1.0 - 2.0 / r
    inv_f  = 1.0 / f
    inv_r2 = 1.0 / r2
    inv_r3 = inv_r2 / r

    # ── ẋ^μ = ∂H/∂p_μ = g^μν p_ν ──
    dtdλ =  -inv_f * p_t
    drdλ =   f * p_r
    dθdλ =   inv_r2 * p_θ
    dφdλ =   inv_r2 / sinθ2 * p_φ

    # ── ṗ_μ = -∂H/∂x^μ = -½ ∂g^αβ/∂x^μ p_α p_β ──
    # ∂g^tt/∂r =  2/(r²f²),   ∂g^rr/∂r = 2/r²
    # ∂g^θθ/∂r = -2/r³,       ∂g^φφ/∂r = -2/(r³sin²θ)
    dp_t = 0.0
    dp_r = -inv_r2 * (inv_f * inv_f * p_t^2 + p_r^2) +
            inv_r3 * (p_θ^2 + p_φ^2 / sinθ2)
    # ∂g^φφ/∂θ = -2cosθ/(r²sin³θ)
    dp_θ = cosθ / (r2 * sinθ2 * sinθ) * p_φ^2
    dp_φ = 0.0

    return SVector{8, Float64}(dtdλ, drdλ, dθdλ, dφdλ, dp_t, dp_r, dp_θ, dp_φ)
end

end # module Schwarzschild