module SimpleDisk

export Structure, intensity

mutable struct Structure
    in_edge::Float64
    out_edge::Float64
end

function Structure(blackhole;
                   R_min::Union{Nothing,Float64}=nothing,
                   R_max::Float64=20.0,
                   corotating::Bool=true)

    in_edge_value = if R_min !== nothing
        R_min
    else
        corotating ? blackhole.ISCOco : blackhole.ISCOcounter
    end

    return Structure(float(in_edge_value), float(R_max))
end

function intensity(s::Structure, r::Real)
    if s.in_edge < r < s.out_edge
        return (r - s.out_edge) / (s.in_edge - s.out_edge)
    else
        return 0.0
    end
end

end # module
