module SimpleDisk

export SimpleDisk, intensity

mutable struct SimpleDisk
    in_edge::Float64
    out_edge::Float64
    intensity::Function
end

function default_intensity(r, in_edge, out_edge)
    if in_edge < r < out_edge
        return (r - out_edge) / (in_edge - out_edge)
    else
        return 0.0
    end
end

function SimpleDisk(blackhole;
                   R_min::Union{Nothing,Float64}=nothing,
                   R_max::Float64=20.0,
                   corotating::Bool=true)

    in_edge_value = if R_min !== nothing
        R_min
    else
        isco = corotating ? blackhole.ISCOco : blackhole.ISCOcounter
        isnothing(isco) && error(
            "SimpleDisk: ISCO does not exist for this black hole. Provide R_min explicitly."
        )
        isco
    end

    intensity_func = r ->
        default_intensity(r, float(in_edge_value), float(R_max))

    return SimpleDisk(float(in_edge_value), float(R_max), intensity_func)
end

end # module