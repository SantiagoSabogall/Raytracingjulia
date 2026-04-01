module SchwarzschildNumerical

using Dierckx
using DelimitedFiles
using StaticArrays

export SchwarzschildNumBH, dr_inverse_metric

import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega



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

        # EH y ISCO por defecto para Schwarzschild (pueden ajustarse según los datos)
        new(rN, Nvals, N, rdN, dNvals, dNdr, 0.0, 2.0, 6.0, 6.0)
    end
end

function Omega(b::SchwarzschildNumBH, r::Real; corotating::Bool=true)
    return 1 / r^(3/2)
end


# -------------------------
# Métodos Extendidos
# -------------------------

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

# Función auxiliar (no necesita ser extendida si solo se usa aquí)
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

# Firma corregida (q, b, lambda)
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