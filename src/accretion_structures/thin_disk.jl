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
    energy::Function
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

    if corotating
        ISCO = blackhole.ISCOco
        in_edge = ISCO
    else
        ISCO = blackhole.ISCOcounter
        in_edge = ISCO
    end

    if R_min !== nothing
        in_edge = float(R_min)
    end

    disk = ThinDisk(in_edge, out_edge, a, ISCO, r -> 0.0)

    # Precompute energy profile
    rr = collect(range(disk.in_edge, stop=disk.out_edge, length=100_000))
    ff = f.(Ref(disk), rr)
    ff .-= minimum(ff)

    interp = interpolate((rr,), ff, Gridded(Linear()))
    disk.energy = r -> interp(r)

    return disk
end


# ============================
# Novikov–Thorne flux function
# ============================

function f(d::ThinDisk, r::Real)
    a = d.a
    arccos_a = acos(a)

    x0 = sqrt(d.ISCO)
    x1 = 2 * cos((arccos_a - π) / 3)
    x2 = 2 * cos((arccos_a + π) / 3)
    x3 = -2 * cos(arccos_a / 3)

    x = sqrt(r)

    c = 3 / (2 * x^4 * (x^3 - 3x + 2a))

    t1 = x - x0 - (3a / 2) * log(x / x0)
    t2 = -((3 * (x1 - a)^2) / (x1 * (x1 - x2) * (x1 - x3))) *
          log((x - x1) / (x0 - x1))
    t3 = -((3 * (x2 - a)^2) / (x2 * (x2 - x1) * (x2 - x3))) *
          log((x - x2) / (x0 - x2))
    t4 = -((3 * (x3 - a)^2) / (x3 * (x3 - x1) * (x3 - x2))) *
          log((x - x3) / (x0 - x3))

    return c * (t1 + t2 + t3 + t4)
end


# ============================
# Intensity profile
# ============================

function intensity(d::ThinDisk, r::Real)
    if d.in_edge < r < d.out_edge
        return d.energy(r)
    else
        return 0.0
    end
end


end # module ThinDisk
