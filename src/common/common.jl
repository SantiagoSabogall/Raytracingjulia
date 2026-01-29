module Common

using ..ImagePlane: photon_coords
using ..Raytracingjulia: metric, inverse_metric, Omega, geodesics


using DifferentialEquations
using Base.Threads
using Plots
using JLD2
using Random
using LinearAlgebra
using Base: circshift

export Photon, Image,
       create_photons!,
       create_image!,
       create_image_no_Doppler!,
       create_shadow!,
       save_data,
       plot,
       plot_shadow,
       plot_contours,
       verify_Hamiltonian,
       doppler_shift

# ============================================================
# Photon
# ============================================================

mutable struct Photon
    alpha::Float64
    beta::Float64
    freq::Float64
    i::Union{Int, Nothing}
    j::Union{Int, Nothing}
    iC::Union{Vector{Float64}, Nothing}
    fP::Union{Vector{Float64}, Nothing}

    function Photon(alpha::Float64, beta::Float64, freq::Float64=1.0)
        new(alpha, beta, freq, nothing, nothing, nothing, nothing)
    end
end

# ============================================================
# Geodesic integrators
# ============================================================

function geodesics_integrate(p::Photon, blackhole, acc_structure, detector)

    final_λ = 1.5 * detector.D
    lmbda_range = range(0.0, -final_λ, length=Int(7 * final_λ))
    tspan = (lmbda_range[1], lmbda_range[end])

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    sol = solve(prob, Tsit5(), 
                reltol=1e-8, 
                abstol=1e-8, 
                verbose=false,
                saveat=lmbda_range)
                    


    p.fP = zeros(8)
    I_f = 0.0

    zi  = [cos(u[3]) for u in sol.u]
    zi1 = circshift(zi, -1)
    zi1[end] = 0.0

    indxs = findall(zi .* zi1 .< 0)

    for idx in indxs
        u = sol.u[idx]
        r = u[2]

        if acc_structure.in_edge < r < acc_structure.out_edge
            p.fP = u
            I0 = acc_structure.intensity(r)
            I_f = doppler_shift(p, I0, blackhole)
            break
        end
    end

    return I_f
end

function geodesics_integrate_no_Doppler(p::Photon, blackhole, acc_structure, detector)

    final_λ = 1.5 * detector.D
    lmbda_range = range(0.0, -final_λ, length=Int(7 * final_λ))
    tspan = (lmbda_range[1], lmbda_range[end])

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    sol = solve(prob, Tsit5(), 
                reltol=1e-8, 
                abstol=1e-8, 
                verbose=false,
                saveat=lmbda_range)

    p.fP = zeros(8)

    zi  = [cos(u[3]) for u in sol.u]
    zi1 = circshift(zi, -1)
    zi1[end] = zi[end]

    indxs = findall(zi .* zi1 .< 0)

    for idx in indxs
        u = sol.u[idx]
        r = u[2]

        if acc_structure.in_edge < r < acc_structure.out_edge
            p.fP = u
            return acc_structure.intensity(r)
        end
    end

    return 0.0
end
"""
function shadow_integ(p::Photon, blackhole, detector)

    final_λ = 1.5 * detector.D
    lmbda_range = range(0.0, -final_λ, length=Int(7 * final_λ))
    tspan = (lmbda_range[1], lmbda_range[end])

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    sol = solve(prob, Tsit5(), 
                reltol=1e-8, 
                abstol=1e-8, 
                verbose=false,
                saveat=lmbda_range)

    rs = [u[2] for u in sol.u]
    indxs = findall(r -> r < blackhole.EH + 1e-7, rs)

    return isempty(indxs) ? 100.0 : 0.0
end
"""
function shadow_integ(p::Photon, blackhole, detector)
    final_λ = 1.5 * detector.D
    tspan = (0.0, -final_λ)

    # Reutilizamos el callback de terminación para que sea ultra rápido
    condition(u, t, integrator) = u[2] - (blackhole.EH + 1e-7)
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    
    # IMPORTANTE: save_everystep=false para no llenar la RAM
    sol = solve(prob, Tsit5(), 
                reltol=1e-8, 
                abstol=1e-8, 
                verbose=false,
                callback=cb,
                save_everystep=false)

    # Lógica de la Sombra:
    # Si la integración terminó PORQUE tocó el horizonte, el radio final será ≈ EH
    # Si terminó porque se acabó el tiempo (final_λ), el radio será grande.
    
    r_final = sol.u[end][2]
    
    if r_final < blackhole.EH + 1e-3 # Un margen un poco más generoso
        return 0.0   # Sombra (Negro)
    else
        return 100.0 # Cielo / Fondo (Blanco)
    end
end
# ============================================================
# Physics
# ============================================================

function doppler_shift(p::Photon, I0::Float64, blackhole)

    if p.fP === nothing || all(p.fP .== 0.0)
        return 0.0
    end

    coords = p.fP[1:4]
    g_tt, _, _, g_phph, g_tph = metric(blackhole, coords)

    Ω = Omega(blackhole, p.fP[2])

    k_t   = p.fP[5]
    k_phi = p.fP[8]

    g = sqrt(-g_tt - 2g_tph*Ω - g_phph*Ω^2) /
        (1 + k_phi*Ω / k_t)

    return I0 * g^3
end

function integrate_for_H(p::Photon, blackhole, acc_structure, detector)
    final_lmbda = 1.5 * detector.D
    lmbda_range = range(0.0, -final_lmbda, length = Int(7 * final_lmbda))
    tspan = (lmbda_range[1], lmbda_range[end])

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)

    sol = solve(prob, Tsit5(),
                reltol = 1e-8,
                abstol = 1e-8,
                saveat = lmbda_range)

    # Empezamos con la solución completa
    solution = sol.u
    last_idx = length(sol.u)

    # --- Lógica de corte por Disco de Acreción ---
    zi = [cos(u[3]) for u in sol.u]
    zi1 = circshift(zi, -1)
    zi1[end] = zi[end] # Evitar el 0.0 para que el signo sea consistente

    indxs = findall(zi .* zi1 .< 0)
    for i in indxs
        r = sol.u[i][2]
        if acc_structure.in_edge < r < acc_structure.out_edge
            last_idx = i
            break
        end
    end

    # --- Lógica de corte por Horizonte de Sucesos (EH) ---
    rEH = [u[2] for u in sol.u]
    idxEH = findfirst(r -> r < (blackhole.EH + 0.05), rEH)
    
    if !isnothing(idxEH) && idxEH < last_idx
        last_idx = idxEH
    end

    # Recortamos la solución y los tiempos hasta el punto de impacto/caída
    final_solution = sol.u[1:last_idx]
    final_lambdas = sol.t[1:last_idx]

    # Calculamos el vector de Hamiltonianos
    H = Hamiltonian(final_solution, blackhole)

    return final_lambdas, H
end



function Hamiltonian(sol::Vector{Vector{Float64}}, blackhole)

    H = zeros(length(sol))

    for i in eachindex(sol)
        x = sol[i][1:4]
        p = sol[i][5:8]

        gtt, grr, gθθ, gφφ, gtφ = inverse_metric(blackhole, x)
    H[i] = 0.5 * (
        gtt * p[1]^2 +
        grr * p[2]^2 +
        gθθ * p[3]^2 +
        gφφ * p[4]^2 +
        2 * gtφ * p[1] * p[4]
        )  
    end

    return H
end

# ============================================================
# Image
# ============================================================

mutable struct Image
    blackhole
    acc_structure
    detector
    photon_list::Vector{Photon}
    image_data::Matrix{Float64}

    function Image(blackhole, acc_structure, detector)
        new(blackhole, acc_structure, detector, Photon[], zeros(0,0))
    end
end

function create_photons!(img::Image)

    println("Creating photons...")
    img.photon_list = Photon[]

    for (i, a) in enumerate(img.detector.alphaRange)
        for (j, b) in enumerate(img.detector.betaRange)
            p = Photon(a, b)
            p.iC = photon_coords(img.detector, img.blackhole, a, b)
            p.i = i
            p.j = j
            push!(img.photon_list, p)
        end
    end

    println("Total photons: ", length(img.photon_list))
end

function create_image!(img::Image)

    img.image_data = zeros(img.detector.x_pixels, img.detector.y_pixels)

    println("Integrating with ", nthreads(), " threads...")
    t0 = time()

    @threads for k in eachindex(img.photon_list)
        p = img.photon_list[k]
        img.image_data[p.i, p.j] =
            geodesics_integrate(p, img.blackhole, img.acc_structure, img.detector)
    end

    println("Total time: ", round(time() - t0, digits=2), " s")
end
function create_image_no_Doppler!(img::Image)

    # Inicializar la imagen
    img.image_data = zeros(img.detector.x_pixels, img.detector.y_pixels)

    photon = 1
    println("Integrating trajectories ...")

    start_time = time()

    for p in img.photon_list
        img.image_data[p.i, p.j] =
            geodesics_integrate_no_Doppler(p, img.blackhole, img.acc_structure, img.detector)

        print("\rPhoton # $photon")
        flush(stdout)

        photon += 1
    end

    total_time = time() - start_time

    println("\n\n--- Total time of integration : $total_time seconds ---")
    println(
        "\n--- Time of integration : $(total_time / length(img.photon_list)) seconds/photon ---\n"
    )
end


function create_shadow!(img::Image)
    # 1. Aseguramos que la matriz esté limpia
    img.image_data = zeros(img.detector.x_pixels, img.detector.y_pixels)
    
    n_photons = length(img.photon_list)
    println("Integrating trajectories for shadow on $(nthreads()) threads...")

    start_time = time()

    # 2. El truco mágico: @threads reparte los fotones entre tus núcleos
    @threads for p in img.photon_list
        # Guardamos el resultado directamente en la matriz
        img.image_data[p.i, p.j] = shadow_integ(p, img.blackhole, img.detector)
    end

    total_time = time() - start_time
    
    # 3. Reporte final
    println("\n\nEH radius: $(img.blackhole.EH)")
    println("--- Total time: $(round(total_time, digits=2)) seconds ---")
    println("--- Speed: $(round(total_time/n_photons, digits=6)) seconds/photon ---")
end

function plot(img::Image; save=false, filename=nothing, cmap=:afmhot)
    # 1. Normalizar datos
    maxv = maximum(img.image_data)
    data = maxv > 0 ? img.image_data ./ maxv : img.image_data

    # 2. Crear el heatmap
    plt = heatmap(
        transpose(data),
        aspect_ratio = 1,
        c = cmap,
        axis = false,
        framestyle = :none,
        colorbar = false
    )

    # 3. Lógica de guardado
    if save && !isnothing(filename)
        mkpath("images")
        # Usamos Plots.savefig para evitar que Julia use tu variable 'save'
        Plots.savefig(plt, "images/$filename.png")
        println("Imagen guardada en: images/$filename.png")
    end

    display(plt)
    return plt
end

function plot_shadow(img; save=false, filename=nothing, cmap=:gray)
    # 1. Copiamos los datos para no modificar la imagen original
    data = copy(img.image_data)
    
    # 2. Normalización segura: Evita división por cero si la imagen está vacía
    maxv = maximum(data)
    if maxv > 0
        data .= data ./ maxv
    end

    # 3. Crear el heatmap
    # Usamos clims=(0, 1) para forzar que el negro sea 0 y blanco 1
    plt = heatmap(
        data', # El operador ' es más rápido que transpose() en matrices grandes
        aspect_ratio = :equal,
        c = cmap,
        clims = (0, 1), 
        axis = true,       # Activamos ejes para ver las coordenadas alfa y beta
        framestyle = :box,
        colorbar = false,
        ticks = false      # Quita los números pero deja el marco si lo prefieres
    )

    # Etiquetas de coordenadas de impacto (α, β)
    xlabel!(plt, "α")
    ylabel!(plt, "β")

    # 4. Guardado corregido
    if save && filename !== nothing
        if !isdir("images") mkpath("images") end
        # Aseguramos la extensión .png
        full_path = endswith(filename, ".png") ? "images/$filename" : "images/$filename.png"
        Plots.savefig(plt, full_path)
        println("Sombra guardada en: $full_path")
    end

    display(plt)
    return plt
end


function plot_contours(img; save=false, filename=nothing, cmap=:gray)
    # 1. Preparar los datos (transponer para que X e Y coincidan con el detector)
    data = img.image_data'

    # 2. Crear el gráfico de contornos
    plt = contour(
        data,
        aspect_ratio = :equal,
        c = cmap,
        fill = false,      # Solo líneas, no relleno
        levels = 15,       # Cantidad de líneas de contorno
        lw = 1.2,          # Grosor de línea (linewidth)
        colorbar = false
    )

    # 3. Estética de los ejes (limpio como en tu versión de Python)
    xlabel!(plt, "α")
    ylabel!(plt, "β")
    plot!(plt, grid=true, gridalpha=0.25, ticks=false)

    # 4. Guardado de imagen
    if save && filename !== nothing
        if !isdir("images") mkpath("images") end
        full_path = endswith(filename, ".png") ? "images/$filename" : "images/$filename.png"
        Plots.savefig(plt, full_path)
        println("Contornos guardados en: $full_path")
    end

    display(plt)
    return plt
end

function verify_Hamiltonian(img::Image; n::Int=10)
    println("Verifying Hamiltonian conservation for $n random photons...\n")

    plt = Plots.plot(
        title = "Hamiltonian Conservation",
        xlabel = "λ (Affine Parameter)",
        ylabel = "H",
        legend = :outerright,
        grid = true,
        gridalpha = 0.3
    )

    # Seleccionamos n índices aleatorios sin repetir
    indices = rand(1:length(img.photon_list), n)

    for i in indices
        p = img.photon_list[i]
        
        # Ahora recibe la tupla (lambdas, H_vals)
        lambdas, H_vals = integrate_for_H(p, img.blackhole, img.acc_structure, img.detector)

        Plots.plot!(plt, lambdas, H_vals, label = "Photon #$i")
    end

    # La escala es crucial: queremos ver desviaciones muy pequeñas cerca de 0
    
    display(plt)
end
end # module Common
