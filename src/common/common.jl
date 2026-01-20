module Common

using ..Detector: photon_coords
using ..Kerr: metric, inverse_metric, Omega, geodesics
using ..Schwarzschild: metric, inverse_metric

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

using DifferentialEquations
using Base.Threads
using Plots
using JLD2
using Random
using LinearAlgebra

# ============================================================
# Estructura del Fotón
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
# Integradores de Geodésicas
# ============================================================

function geodesics_integrate(p::Photon, blackhole, acc_structure, detector)
    final_lmbda = 1.5 * detector.D
    lmbda_range = range(0, -final_lmbda, length=Int(7 * final_lmbda))
    tspan = (lmbda_range[1], lmbda_range[end])

    # CORRECCIÓN: Se pasa 'geodesics' como función y 'blackhole' como parámetro de la ODE
    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    sol = solve(prob, Tsit5(), saveat=lmbda_range)

    p.fP = zeros(8)
    I_f = 0.0
    
    # Análisis de cruce por el plano ecuatorial (theta = pi/2 -> cos(theta) = 0)
    zi = [cos(u[3]) for u in sol.u]
    zil = circshift(zi, -1)
    zil[end] = zi[end] # Evitar el 0.0 para no crear un cruce artificial al final

    indxs = findall(zi .* zil .<= 0)

    for idx in indxs
        current_sol = sol.u[idx]
        r = current_sol[2]

        if r < acc_structure.out_edge && r > acc_structure.in_edge
            p.fP = current_sol
            I_0 = acc_structure.intensity(p.fP[2])
            I_f = doppler_shift(p, I_0, blackhole)
            break
        end
    end

    return I_f
end

function geo_inte_no_Doppler(p::Photon, blackhole, acc_structure, detector)
    final_lmbda = 1.5 * detector.D
    lmbda_range = range(0, -final_lmbda, length=Int(7 * final_lmbda))
    tspan = (lmbda_range[1], lmbda_range[end])

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    sol = solve(prob, Tsit5(), saveat=lmbda_range)

    p.fP = zeros(8)
    zi = [cos(u[3]) for u in sol.u]
    zil = circshift(zi, -1)
    zil[end] = zi[end]

    indxs = findall(zi .* zil .<= 0)

    for idx in indxs
        current_sol = sol.u[idx]
        r = current_sol[2]

        if r < acc_structure.out_edge && r > acc_structure.in_edge
            p.fP = current_sol
            break
        end
    end

    if length(p.fP) > 1 && p.fP[2] > 0
        return acc_structure.intensity(p.fP[2])
    else
        return 0.0
    end
end

function shadow_integ(p::Photon, blackhole, detector)
    final_lmbda = 1.5 * detector.D
    lmbda_range = range(0, -final_lmbda, length=Int(7 * final_lmbda))
    tspan = (lmbda_range[1], lmbda_range[end])

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    sol = solve(prob, Tsit5(), saveat=lmbda_range)
    
    r_vals = [u[2] for u in sol.u]
    # Si el radio cae por debajo del Horizonte de Eventos (EH)
    indxs = findall(r_vals .< (blackhole.EH + 1e-5))

    return isempty(indxs) ? 100.0 : 0.0
end

# ============================================================
# Física y Auxiliares
# ============================================================

function doppler_shift(p::Photon, I0::Float64, blackhole)
    if p.fP === nothing || all(p.fP .== 0.0)
        return 0.0
    end
    
    coords = p.fP[1:4]
    g_tt, _, _, g_phph, g_tph = metric(blackhole, coords)
    
    # Aquí usamos la función Omega importada
    omega_disk = Omega(blackhole, p.fP[2]) 
    
    k_t = p.fP[5]
    k_phi = p.fP[8]
    
    # Factor de corrimiento
    g_factor = sqrt(-g_tt - 2*g_tph*omega_disk - g_phph*omega_disk^2) / (1 + k_phi*omega_disk/k_t)
    return I0 * g_factor^3
end

function Hamiltonian(sol_vec::Vector{Vector{Float64}}, blackhole)
    H = zeros(length(sol_vec))

    for i in eachindex(sol_vec)
        x = sol_vec[i][1:4]
        p_mom = sol_vec[i][5:8] 

        # CORRECCIÓN: Llamada funcional a inverse_metric
        gtt, grr, gthth, gphph, gtph = inverse_metric(blackhole, x)

        H[i] = 0.5 * (
            gtt * p_mom[1]^2 +
            grr * p_mom[2]^2 +
            gthth * p_mom[3]^2 +
            gphph * p_mom[4]^2 + 
            2 * gtph * p_mom[1] * p_mom[4]
        )
    end
    return H
end

# ============================================================
# Gestión de Imagen
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
    
    # Acceso a los rangos del detector
    alpha_range = img.detector.alphaRange
    beta_range = img.detector.betaRange

    for (i, a) in enumerate(alpha_range)
        for (j, b) in enumerate(beta_range)
            p = Photon(a, b)
            # CORRECCIÓN: Llamada funcional
            p.iC = photon_coords(img.detector, img.blackhole, a, b)
            p.i = i
            p.j = j
            push!(img.photon_list, p)
        end
    end
    println("Total photons created: ", length(img.photon_list))
end

function create_image!(img::Image)
    img.image_data = zeros(img.detector.x_pixels, img.detector.y_pixels)
    
    println("Integrating trajectories with ", nthreads(), " threads...")
    start_time = time()

    @threads for idx in eachindex(img.photon_list)
        p = img.photon_list[idx]
        img.image_data[p.i, p.j] = geodesics_integrate(p, img.blackhole, img.acc_structure, img.detector)
    end

    total_time = time() - start_time
    println("\n--- Total time: ", round(total_time, digits=2), " seconds ---")
end

# ... (Las demás funciones create_image_no_Doppler!, create_shadow!, etc. 
# seguirían el mismo patrón de corrección que create_image!) ...

# ============================================================
# Visualización y Guardado
# ============================================================

function plot(img::Image; savefig=false, filename=nothing, cmap=:afmhot)
    # Normalización para visualización
    max_val = maximum(img.image_data)
    data_norm = max_val > 0 ? img.image_data ./ max_val : img.image_data
    
    plt = heatmap(
        transpose(data_norm),
        aspect_ratio=1,
        c=cmap,
        axis=false,
        framestyle=:none
    )
    
    if savefig && filename !== nothing
        mkpath("images")
        Plots.savefig(plt, "images/$(filename).png")
        println("Image saved to images/$(filename).png")
    end
    display(plt)
end

end # module Common