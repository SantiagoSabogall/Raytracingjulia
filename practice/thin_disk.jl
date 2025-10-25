using Interpolations: interpolate, Gridded, Linear

mutable struct Structure
    in_edge::Float64
    out_edge::Float64
    a::Float64
    ISCO::Float64
    energy::Any
end

function Structure(blackhole; corotating::Bool=true, R_min::Union{Nothing, Float64}=nothing, R_max::Float64=20.0)

    out_edge_value = float(R_max)
    a = float(blackhole.a)

    if corotating
        ISCO = blackhole.ISCOco
        in_edge_value = blackhole.ISCOco
    else
        ISCO = blackhole.ISCOcounter
        in_edge_value = blackhole.ISCOcounter
    end

    if R_min !== nothing
        in_edge_value = float(R_min)
    end

    s = Structure(in_edge_value, out_edge_value, a, ISCO, nothing)

    rr = collect(range(s.in_edge, stop=s.out_edge, length=100000))
    ff = f.(Ref(s), rr)
    ff .-= minimum(ff)

    energy = interpolate((rr,), ff, Gridded(Linear()))
    s.energy = energy
    return s
end

function f(s::Structure, r::Real)
    a_M = s.a
    arcos_aM = acos(a_M)

    x0 = sqrt(s.ISCO)
    x1 = 2 * cos((arcos_aM - pi) / 3)
    x2 = 2 * cos((arcos_aM + pi) / 3)
    x3 = -2 * cos(arcos_aM / 3)
    x = sqrt(r)

    c = 3 / (2 * (x^4) * (x^3 - 3*x + 2 *a_M))
    
    t1 = x - x0 - (3*s.a/2) * log(x / x0)
    t2 = -((3 * (x1 - a_M)^2) / (x1 * (x1 - x2) * (x1 - x3))) * log((x - x1)/(x0 - x1))
    t3 = -((3 * (x2 - a_M)^2) / (x2 * (x2 - x1) * (x2 - x3))) * log((x - x2)/(x0 - x2))
    t4 = -((3 * (x3 - a_M)^2) / (x3 * (x3 - x1) * (x3 - x2))) * log((x - x3)/(x0 - x3))

    return c * (t1 + t2 + t3 + t4)

end

function intensity(s::Structure, r::Real)
    if r > s.in_edge && r < s.out_edge
        return s.energy(r)
    else
        return 0.0
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    println()
    println("THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.")
    println("YOU NEED TO RUN THE main.jl FILE TO GENERATE THE IMAGE")
    println()
end