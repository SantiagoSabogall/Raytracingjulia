"""
    ThinDiskmod

Implements the Novikov-Thorne geometrically thin, optically thick accretion disk model.
"""
module ThinDiskmod

export ThinDisk, intensity,f

using Interpolations: interpolate, Gridded, Linear

# ============================
# Thin disk structure
# ============================
# ============================
# ============================
# Thin disk structure
# ============================
"""
    ThinDisk

Structure representing a Novikov-Thorne accretion disk.

# Fields
- `in_edge`: Inner edge radius (typically ISCO).
- `out_edge`: Outer edge radius.
- `a`: Black hole spin.
- `ISCO`: Innermost Stable Circular Orbit.
- `energy_profile`: Interpolated function for the radial emission profile.
- `intensity`: Final radiative intensity function \$I(r)\$.
"""
mutable struct ThinDisk
    in_edge::Float64
    out_edge::Float64
    a::Float64
    ISCO::Float64
    energy_profile::Function
    intensity::Function
end

# ============================
# Constructor
# ============================
function ThinDisk(
    blackhole;
    corotating::Bool = true,
    R_min::Union{Nothing, Float64} = nothing,
    R_max::Float64 = 20.0
)
    out_edge = float(R_max)
    a = float(blackhole.a)

    ISCO = if hasproperty(blackhole, :ISCOco)
        corotating ? blackhole.ISCOco : blackhole.ISCOcounter
    else
        6.0
    end

    in_edge = R_min === nothing ? ISCO : float(R_min)

    disk = ThinDisk(in_edge, out_edge, a, ISCO, r -> 0.0, r -> 0.0)

    # ============================
    # Precompute Novikov–Thorne profile
    # ============================
    rr = collect(range(in_edge, out_edge; length = 100_000))
    ff = [f(disk, r) for r in rr]
    ff .-= minimum(ff)

    interp = interpolate((rr,), ff, Gridded(Linear()))

    disk.energy_profile = r ->
        (in_edge < r < out_edge) ? interp(r) : 0.0

    disk.intensity = r -> disk.energy_profile(r)

    return disk
end

# ============================
# Novikov–Thorne flux
# ============================
# ============================
# Novikov–Thorne flux
# ============================
"""
    f(d::ThinDisk, r::Real)

Calculates the Novikov-Thorne radiative flux \$\\mathcal{F}(r)\$ for a 
geometrically thin accretion disk according to Page & Thorne (1974).
"""
function f(d::ThinDisk, r::Real)
    a = d.a
    arccos_a = acos(a)

    x0 = sqrt(d.ISCO)
    x  = sqrt(r)

    x1 = 2 * cos((arccos_a - π) / 3)
    x2 = 2 * cos((arccos_a + π) / 3)
    x3 = -2 * cos(arccos_a / 3)

    denom = x^3 - 3x + 2a
    if abs(denom) < 1e-10
        return 0.0
    end

    c = 3 / (2 * x^4 * denom)

    t1 = x - x0 - (3a / 2) * log(x / x0)
    t2 = -((3 * (x1 - a)^2) / (x1 * (x1 - x2) * (x1 - x3))) *
         log((x - x1) / (x0 - x1))
    t3 = -((3 * (x2 - a)^2) / (x2 * (x2 - x1) * (x2 - x3))) *
         log((x - x2) / (x0 - x2))
    t4 = -((3 * (x3 - a)^2) / (x3 * (x3 - x1) * (x3 - x2))) *
         log((x - x3) / (x0 - x3))

    return max(0.0, c * (t1 + t2 + t3 + t4))
end

end # module ThinDisk