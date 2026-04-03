"""
    ImagePlane

Defines the observer's screen (detector) setup and the local mapping 
of pixel coordinates to 4-momentum using Bardeen's impact parameters.
"""
module ImagePlane

using StaticArrays
using ..Raytracingjulia: metric


export Detector, photon_coords

# ============================================================
# Detector (pantalla del observador)
# ============================================================

# ============================================================
# Detector (pantalla del observador)
# ============================================================

"""
    Detector

Defines the properties of the observational camera.

# Fields
- `D`: Distance from the black hole to the observer.
- `iota` (\$\\iota\$): Inclination angle of the observer with respect to the black hole's spin axis.
- `sin_iota` / `cos_iota`: Pre-computed trigonometrics.
- `x_pixels` / `y_pixels`: Resolution dimensions.
- `alphaRange` / `betaRange`: Range of pixel impact parameters \$(\\alpha, \\beta)\$.
"""
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

        # Enforce even scaling for pixel spatial resolution
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

# ============================================================
# Condiciones iniciales del fotón
# ============================================================

"""
    photon_coords(d::Detector, blackhole, alpha::Float64, beta::Float64; freq::Float64=1.0)

Calculates the initial 4-position \$x^\\mu = (0, r, \\theta, \\phi)\$ and null 
4-momentum \$k_\\mu = (k_t, k_r, k_\\theta, k_\\phi)\$ for a photon arriving at the 
detector pixel \$(\\alpha, \\beta)\$.

Derived from Bardeen's equations using the observer's local tetrad.
Returns an `SVector{8, Float64}`.
"""
function photon_coords(d::Detector, blackhole, alpha::Float64, beta::Float64;
                       freq::Float64=1.0)

    # Initial spatial projection bounds (Screen → Boyer-Lindquist formulation)
    r = sqrt(alpha*alpha + beta*beta + d.D*d.D)

    arg = (beta * d.sin_iota + d.D * d.cos_iota) / r
    arg = clamp(arg, -1.0, 1.0)

    theta = acos(arg)
    phi   = atan(alpha, (d.D * d.sin_iota - beta * d.cos_iota))

    x = @SVector [0.0, r, theta, phi]

    # Metric covariant array assignments
    g_tt, g_rr, g_thth, g_phph, g_tph = metric(blackhole, x)

    # Spatial null momentum boundary variables
    invD = 1.0 / d.D
    k_th = sqrt(g_thth) * beta * invD
    k_ph = -sqrt(g_phph) * alpha * invD

    # Energy bounds mapping constrained to null geodesic normalization
    k_t = -freq * sqrt(g_phph / (g_tph^2 - g_tt * g_phph)) +
          alpha * g_tph * invD / sqrt(g_phph)

    term_r = 1.0 - (k_th^2 / g_thth) - (k_ph^2 / g_phph)
    k_r = sqrt(g_rr * max(0.0, term_r))   # Inbound trajectory geometry constraints

    # Returned via direct SVector{8} structural typing bypassing dynamic heap allocations
    return SVector{8, Float64}(0.0, r, theta, phi, k_t, k_r, k_th, k_ph)
end



# ============================================================

if abspath(PROGRAM_FILE) == @__FILE__
    println()
    println("THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.")
    println("YOU NEED TO RUN main.jl TO GENERATE THE IMAGE")
    println()
end

end # module Detector
