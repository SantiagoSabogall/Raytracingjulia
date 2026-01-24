module SimpleDisk

export SimpleDisk, intensity

mutable struct SimpleDisk
    in_edge::Float64
    out_edge::Float64
    intensity::Function # Este campo guardará la lógica de brillo
end

# 1. Definimos la función de intensidad FUERA del constructor
# Esta es la lógica que quieres para el disco
function default_intensity(r, in_edge, out_edge)
    if in_edge < r < out_edge
        return (r - out_edge) / (in_edge - out_edge)
    else
        return 0.0
    end
end

# 2. Constructor
function SimpleDisk(blackhole;
                   R_min::Union{Nothing,Float64}=nothing,
                   R_max::Float64=20.0,
                   corotating::Bool=true)

    in_edge_value = if R_min !== nothing
        R_min
    else
        corotating ? blackhole.ISCOco : blackhole.ISCOcounter
    end

    # Creamos una "clausura" (closure) para que la función de intensidad
    # ya sepa cuáles son sus bordes sin tener que pedirlos luego.
    intensity_func = r -> default_intensity(r, float(in_edge_value), float(R_max))

    return SimpleDisk(float(in_edge_value), float(R_max), intensity_func)
end

end # module