using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Revise
using Raytracingjulia
using LinearAlgebra
using NPZ

# ----------------------------
# Black Hole Configuration
# ----------------------------
a = 0.5
blackhole = KerrPDFMBH(a,0.5 )

using StaticArrays

function photon_coords(d::Detector, blackhole, alpha::Float64, beta::Float64;
                       freq::Float64=1.0)

    # Posición inicial (pantalla → BL)
    r = sqrt(alpha*alpha + beta*beta + d.D*d.D)

    arg = (beta * d.sin_iota + d.D * d.cos_iota) / r
    arg = clamp(arg, -1.0, 1.0)

    theta = acos(arg)
    phi   = atan(alpha, (d.D * d.sin_iota - beta * d.cos_iota))

    # 🔥 SVector directamente
    x = @SVector [0.0, r, theta, phi]

    # Métrica
    g_tt, g_rr, g_thth, g_phph, g_tph = metric(blackhole, x)

    # Componentes espaciales
    invD = 1.0 / d.D

    k_th = sqrt(g_thth) * beta * invD
    k_ph = -sqrt(g_phph) * alpha * invD

    # Energía
    denom = sqrt(g_tph*g_tph - g_tt*g_phph)
    k_t = -freq * sqrt(g_phph) / denom +
          alpha * g_tph * invD / sqrt(g_phph)

    # Radial
    term_r = 1.0 - (k_th*k_th / g_thth) - (k_ph*k_ph / g_phph)
    k_r = sqrt(g_rr * max(0.0, term_r))

    # 🔥 Estado completo SIN vcat
    return @SVector [0.0, r, theta, phi, k_t, k_r, k_th, k_ph]
end