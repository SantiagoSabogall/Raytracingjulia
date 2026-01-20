module Detector

export Detector, photon_coords

# ============================================================
# Detector (pantalla del observador)
# ============================================================

mutable struct Detector
    D::Float64
    iota::Float64
    sin_iota::Float64
    cos_iota::Float64
    x_pixels::Int
    y_pixels::Int
    alphaRange::AbstractRange{Float64}
    betaRange::AbstractRange{Float64}

    function Detector(D::Float64, iota::Float64, x_side::Float64;
                      x_pixels::Int=25, ratio::String="16:9")

        # Forzar número par de píxeles
        actual_x_pixels = isodd(x_pixels) ? x_pixels + 1 : x_pixels

        if ratio == "16:9"
            actual_y_pixels = floor(Int, actual_x_pixels * 9 / 16)
            y_side = x_side * 9 / 16
        elseif ratio == "4:3"
            actual_y_pixels = floor(Int, actual_x_pixels * 3 / 4)
            y_side = x_side * 3 / 4
        else
            actual_y_pixels = actual_x_pixels
            y_side = x_side
        end

        alphaRange = range(-x_side, x_side, length=actual_x_pixels)
        betaRange  = range(-y_side,  y_side,  length=actual_y_pixels)

        println("\nScreen size: $actual_x_pixels × $actual_y_pixels pixels")
        println("Total photons: $(actual_x_pixels * actual_y_pixels)")

        new(
            D,
            iota,
            sin(iota),
            cos(iota),
            actual_x_pixels,
            actual_y_pixels,
            alphaRange,
            betaRange
        )
    end
end

# ============================================================
# Condiciones iniciales del fotón
# ============================================================

function photon_coords(d::Detector, blackhole, alpha::Float64, beta::Float64;
                       freq::Float64=1.0)

    # Posición inicial (pantalla → BL)
    r = sqrt(alpha^2 + beta^2 + d.D^2)

    arg = (beta * d.sin_iota + d.D * d.cos_iota) / r
    arg = clamp(arg, -1.0, 1.0)

    theta = acos(arg)
    phi   = atan(alpha, (d.D * d.sin_iota - beta * d.cos_iota))

    x = [0.0, r, theta, phi]

    # Métrica
    g_tt, g_rr, g_thth, g_phph, g_tph = blackhole.metric(x)

    # Componentes espaciales del momento
    k_th = sqrt(g_thth) * beta / d.D
    k_ph = -sqrt(g_phph) * alpha / d.D

    # Energía (normalización nula)
    k_t = -freq * sqrt(g_phph / (g_tph^2 - g_tt * g_phph)) +
          alpha * g_tph / (d.D * sqrt(g_phph))

    term_r = 1.0 - (k_th^2 / g_thth) - (k_ph^2 / g_phph)
    k_r = -sqrt(g_rr * max(0.0, term_r))   # fotón entrante

    k = [k_t, k_r, k_th, k_ph]

    return vcat(x, k)
end

# ============================================================

if abspath(PROGRAM_FILE) == @__FILE__
    println()
    println("THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.")
    println("YOU NEED TO RUN main.jl TO GENERATE THE IMAGE")
    println()
end

end # module Detector
