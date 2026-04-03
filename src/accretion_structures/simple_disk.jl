"""
    SimpleDiskmod

Provides a simplified toy model for an accretion disk with a linear intensity gradient.
"""
module SimpleDiskmod

export SimpleDisk, intensity

"""
    SimpleDisk

Structure representing a generic localized geometrically flat accretion disk.

# Fields
- `in_edge`: Inner radial edge.
- `out_edge`: Outer radial edge.
- `intensity`: A functional mapping \$r \\mapsto I\$.
"""
mutable struct SimpleDisk
    in_edge::Float64
    out_edge::Float64
    intensity::Function
end

"""
    default_intensity(r, in_edge, out_edge)

A normalized linear gradient intensity decaying from `in_edge` (1.0) towards `out_edge` (0.0).
"""
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