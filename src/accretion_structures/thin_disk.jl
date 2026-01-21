module ThinDisk

export ThinDisk, intensity

using Interpolations: interpolate, Gridded, Linear

# ============================
# Thin disk structure
# ============================
mutable struct ThinDisk
    in_edge::Float64
    out_edge::Float64
    a::Float64
    ISCO::Float64
    energy_profile::Function  # Renombrado para evitar conflicto
    intensity::Function       # Campo necesario para run2.jl
end

# ============================
# Constructor
# ============================
function ThinDisk(
    blackhole;
    corotating::Bool=true,
    R_min::Union{Nothing, Float64}=nothing,
    R_max::Float64=20.0
)
    out_edge = float(R_max)
    a = float(blackhole.a)
    
    # Determinar ISCO basado en el tipo de BH (Kerr o Schwarzschild)
    # Si es Schwarzschild, ISCO es usualmente 6.0
    isco_val = hasproperty(blackhole, :ISCOco) ? (corotating ? blackhole.ISCOco : blackhole.ISCOcounter) : 6.0
    
    in_edge = R_min !== nothing ? float(R_min) : isco_val

    # Inicialización temporal
    disk = ThinDisk(in_edge, out_edge, a, isco_val, r -> 0.0, r -> 0.0)

    # Precompute Novikov–Thorne energy profile
    rr = collect(range(disk.in_edge, stop=disk.out_edge, length=10_000))
    ff = [f(disk, r) for r in rr]
    ff .-= minimum(ff)

    interp = interpolate((rr,), ff, Gridded(Linear()))
    
    # Asignamos el perfil físico a energy_profile
    disk.energy_profile = r -> (disk.in_edge <= r <= disk.out_edge) ? interp(r) : 0.0
    
    # Por defecto, la intensidad es el perfil de energía física
    disk.intensity = r -> disk.energy_profile(r)

    return disk
end

# ============================
# Novikov–Thorne flux function
# ============================
function f(d::ThinDisk, r::Real)
    # Evitar singularidades en el borde interno
    if r <= d.ISCO return 0.0 end
    
    a = d.a
    arccos_a = acos(a)
    x0 = sqrt(d.ISCO)
    x = sqrt(r)
    
    # Raíces de la ecuación característica
    x1 = 2 * cos((arccos_a - π) / 3)
    x2 = 2 * cos((arccos_a + π) / 3)
    x3 = -2 * cos(arccos_a / 3)

    c = 3 / (2 * x^4 * (x^3 - 3x + 2a))

    t1 = x - x0 - (3a / 2) * log(x / x0)
    t2 = -((3 * (x1 - a)^2) / (x1 * (x1 - x2) * (x1 - x3))) * log((x - x1) / (x0 - x1))
    t3 = -((3 * (x2 - a)^2) / (x2 * (x2 - x1) * (x2 - x3))) * log((x - x2) / (x0 - x2))
    t4 = -((3 * (x3 - a)^2) / (x3 * (x3 - x1) * (x3 - x2))) * log((x - x3) / (x0 - x3))

    return max(0.0, c * (t1 + t2 + t3 + t4))
end

end # module ThinDisk