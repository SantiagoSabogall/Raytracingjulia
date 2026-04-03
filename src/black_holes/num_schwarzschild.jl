"""
    SchwarzschildNumerical

Provides a numerical framework for Schwarzschild spacetime where the lapse function \$N(r)\$
is interpolated from external data constraints.
"""
module SchwarzschildNumerical

using Dierckx
using DelimitedFiles
using StaticArrays

export SchwarzschildNumBH, dr_inverse_metric

import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega



"""
    SchwarzschildNumBH(path_N::String, path_dN::String)

A dynamically loadable numerical Schwarzschild black hole. 
Constructs cubic splines from text files corresponding to \$N(r)\$ and \$\\frac{dN}{dr}\$.
"""
mutable struct SchwarzschildNumBH
    rN::Vector{Float64}
    Nvals::Vector{Float64}
    N::Spline1D

    rdN::Vector{Float64}
    dNvals::Vector{Float64}
    dNdr::Spline1D

    a::Float64
    EH::Float64
    ISCOco::Float64
    ISCOcounter::Float64

    function SchwarzschildNumBH(path_N::String, path_dN::String)
        dataN = readdlm(path_N)
        datadN = readdlm(path_dN)

        rN = Vector{Float64}(dataN[:,1])
        Nvals = Vector{Float64}(dataN[:,2])
        rdN = Vector{Float64}(datadN[:,1])
        dNvals = Vector{Float64}(datadN[:,2])

        N = Spline1D(rN, Nvals)
        dNdr = Spline1D(rdN, dNvals)

        # Default Event Horizon and ISCO values for Schwarzschild (configurable depending on numerical constraints)
        new(rN, Nvals, N, rdN, dNvals, dNdr, 0.0, 2.0, 6.0, 6.0)
    end
end

"""
    Omega(b::SchwarzschildNumBH, r::Real; corotating::Bool=true)

Calculates the Keplerian angular velocity for a numerically defined metric:
\$\\Omega_K = \\frac{1}{r^{3/2}}\$ (assuming mass \$M=1\$ is implicit in the numerical tables).
"""
function Omega(b::SchwarzschildNumBH, r::Real; corotating::Bool=true)
    return 1 / r^(3/2)
end


# ==============================================================================
# Extended Methods
# ==============================================================================

"""
    metric(b::SchwarzschildNumBH, x::AbstractVector)

Evaluates the numerically interpolated Schwarzschild metric:
\$ds^2 = -N(r) dt^2 + \\frac{1}{N(r)} dr^2 + r^2 d\\theta^2 + r^2 \\sin^2\\theta d\\phi^2\$

Returns `(g_tt, g_rr, g_thth, g_phph, g_tph)`.
"""
function metric(b::SchwarzschildNumBH, x::AbstractVector)
    r = x[2]
    θ = x[3]
    Nr = b.N(r)

    g_tt   = -Nr
    g_rr   = 1 / Nr
    g_thth = r^2
    g_phph = (r*sin(θ))^2
    g_tph  = 0.0

    return g_tt, g_rr, g_thth, g_phph, g_tph
end

"""
    inverse_metric(b::SchwarzschildNumBH, x::AbstractVector)

Evaluates the contravariant components of the numerically interpolated metric.
"""
function inverse_metric(b::SchwarzschildNumBH, x::AbstractVector)
    r = x[2]
    θ = x[3]
    Nr = b.N(r)

    gtt   = -1 / Nr
    grr   = Nr
    gthth = 1 / r^2
    gphph = 1 / (r*sin(θ))^2
    gtph  = 0.0

    return gtt, grr, gthth, gphph, gtph
end

# Auxiliary Function (Local to the numerical module computations)
"""
    dr_inverse_metric(b::SchwarzschildNumBH, x::AbstractVector)

Calculates the radial derivative of the contravariant metric components, \$\\partial_r g^{\\mu\\nu}\$, 
using the interpolated derivative \$\\frac{dN}{dr}\$.
"""
function dr_inverse_metric(b::SchwarzschildNumBH, x::AbstractVector)
    r = x[2]
    θ = x[3]
    Nr = b.N(r)
    dNr = b.dNdr(r)

    drgtt   = dNr / Nr^2
    drgrr   = dNr
    drgthth = -2 / r^3
    drgphph = -2 / (r^3*sin(θ)^2)
    drgtph  = 0.0

    return drgtt, drgrr, drgthth, drgphph, drgtph
end

# Corrected signature mapped to (q, b, lambda)
"""
    geodesics(q, b::SchwarzschildNumBH, λ)

Evaluates the right-hand side of the 8-ODE system for null geodesics.

Derived using the fully numeric Hamiltonian formulation:
\$\\dot{p}_\\alpha = -\\frac{1}{2} (\\partial_\\alpha g^{\\mu\\nu}) p_\\mu p_\\nu\$
"""
function geodesics(q, b::SchwarzschildNumBH, λ)
    g = inverse_metric(b, q)
    drg = dr_inverse_metric(b, q)

    # q = [t, r, θ, φ, k_t, k_r, k_th, k_φ]
    dtdλ   = g[1]*q[5]
    drdλ   = g[2]*q[6]
    dθdλ   = g[3]*q[7]
    dφdλ   = g[4]*q[8]

    dk_t   = 0.0
    dk_r   = -0.5*(drg[1]*q[5]^2 + drg[2]*q[6]^2 + drg[3]*q[7]^2 + drg[4]*q[8]^2)
    dk_th  = (cos(q[3]) / sin(q[3])^3) * (q[8] / q[2])^2
    dk_ph  = 0.0

    return SVector{8, Float64}(dtdλ, drdλ, dθdλ, dφdλ, dk_t, dk_r, dk_th, dk_ph)
end

end # module