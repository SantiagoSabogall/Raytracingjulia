module SchwarzschildNumerical

using Dierckx
using DelimitedFiles

export SchwarzschildNumBH, metric, inverse_metric, dr_inverse_metric, geodesics

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

        rN = dataN[:,1]
        Nvals = dataN[:,2]
        rdN = datadN[:,1]
        dNvals = datadN[:,2]

        N = Spline1D(rN, Nvals)
        dNdr = Spline1D(rdN, dNvals)

        new(rN, Nvals, N, rdN, dNvals, dNdr, 0.0, 2.0, 6.0, 6.0)
    end
end

# -------------------------
# Metric
# -------------------------
function metric(b::SchwarzschildNumBH, x::AbstractVector{<:Real})
    r = x[2]
    θ = x[3]

    g_tt   = -b.N(r)
    g_rr   = 1 / b.N(r)
    g_thth = r^2
    g_phph = (r*sin(θ))^2
    g_tph  = 0.0

    return (g_tt, g_rr, g_thth, g_phph, g_tph)
end

function inverse_metric(b::SchwarzschildNumBH, x::AbstractVector{<:Real})
    r = x[2]
    θ = x[3]

    gtt   = -1 / b.N(r)
    grr   = b.N(r)
    gthth = 1 / r^2
    gphph = 1 / (r*sin(θ))^2
    gtph  = 0.0

    return (gtt, grr, gthth, gphph, gtph)
end

function dr_inverse_metric(b::SchwarzschildNumBH, x::AbstractVector{<:Real})
    r = x[2]
    θ = x[3]

    drgtt   = b.dNdr(r) / b.N(r)^2
    drgrr   = b.dNdr(r)
    drgthth = -2 / r^3
    drgphph = -2 / (r^3*sin(θ)^2)
    drgtph  = 0.0

    return (drgtt, drgrr, drgthth, drgphph, drgtph)
end

function geodesics(b::SchwarzschildNumBH, q::AbstractVector{<:Real}, λ::Float64)
    gtt, grr, gthth, gphph, _ = inverse_metric(b, q)
    drgtt, drgrr, drgthth, drgphph, _ = dr_inverse_metric(b, q)

    dtdλ   = gtt*q[5]
    drdλ   = grr*q[6]
    dθdλ   = gthth*q[7]
    dφdλ   = gphph*q[8]

    dk_t   = 0.0
    dk_r   = -0.5*(drgtt*q[5]^2 + drgrr*q[6]^2 + drgthth*q[7]^2 + drgphph*q[8]^2)
    dk_th  = cos(q[3])/(sin(q[3])^3)*(q[8]/q[2])^2
    dk_ph  = 0.0

    return (dtdλ, drdλ, dθdλ, dφdλ, dk_t, dk_r, dk_th, dk_ph)
end

end # module
