# Acretion structures - simple disk

mutable struct Structure
    in_edge::Float64
    out_edge::Float64
end

function Structure(blackhole; R_min::Union{Nothing, Float64}=nothing, R_max::Float64=20.0, corotating::Bool=true)

    in_edge_value = if R_min !== false
        R_min
    else
        if corotating
            blackhole.ISCOco
        else
            blackhole.ISCOcounter
        end
    end

    return Structure(in_edge_value, R_max)
end


function intensity(s::Structure, r::Real)
    m = (1)/(s.in_edge - s.out_edge)
    I = m * (r - s.out_edge)

    if r > s.in_edge && r < s.out_edge
        return I
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