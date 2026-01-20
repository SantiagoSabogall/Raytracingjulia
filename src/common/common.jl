module Common

export Photon, Image,
       create_photons!,
       create_image!,
       create_image_no_Doppler!,
       create_shadow!,
       save_data,
       plot,
       plot_shadow,
       plot_contours,
       verify_Hamiltonian


using DifferentialEquations
using Base.Threads
using Plots
using JLD2
using Random
using LinearAlgebra

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

function geodesics_integrate(p::Photon, blackhole, acc_structure, detector)
    final_lmbda = 1.5 * detector.D
    lmbda = range(0, -final_lmbda, length=Int(7 * final_lmbda))
    tspan = (lmbda[1], lmbda[end])

    prob = ODEProblem(blackhole.geodesics, p.iC, tspan)
    sol = solve(prob, Tsit5(), saveat=lmbda)

    p.fP = zeros(8)
    I_f = 0.0
    zi = [cos(u[3]) for u in sol.u]
    zil = circshift(zi, -1)
    zil[end] = 0.0

    indxs = findall(zi .* zil .< 0)

    for i in indxs
        current_sol = sol.u[i]
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
    lmbda = range(0, -final_lmbda, length=Int(7 * final_lmbda))
    tspan = (lmbda[1], lmbda[end])

    prob = ODEProblem(blackhole.geodesics, p.iC, tspan)
    sol = solve(prob, Tsit5(), saveat=lmbda)

    p.fP = zeros(8)
    zi = [cos(u[3]) for u in sol.u]
    zil = circshift(zi, -1)
    zil[end] = 0.0

    indxs = findall(zi .* zil .< 0)

    for i in indxs
        current_sol = sol.u[i]
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
    lmbda = range(0, -final_lmbda, length=Int(7 * final_lmbda))
    tspan = (lmbda[1], lmbda[end])

    prob = ODEProblem(blackhole.geodesics, p.iC, tspan)
    sol = solve(prob, Tsit5(), saveat=lmbda)
    
    # Extraer valores de r
    r_vals = [u[2] for u in sol.u]
    
    indxs = findall(r_vals .< (blackhole.EH + 1e-7))

    if isempty(indxs)
        return 100.0
    else
        return 0.0
    end
end

function doppler_shift(p::Photon, I0::Float64, blackhole)
    if p.fP === nothing || length(p.fP) < 8
        return 0.0
    end
    
    coords = p.fP[1:4]
    g_tt, _, _, g_phph, g_tph = blackhole.metric(coords)
    Omega = blackhole.Omega(p.fP[2])
    
    # Corrección: índices en Julia comienzan en 1, no en 0
    k_t = p.fP[5]
    k_phi = p.fP[8]
    
    g = sqrt(-g_tt - 2*g_tph*Omega - g_phph*Omega^2) / (1 + k_phi*Omega/k_t)
    return I0 * g^3
end

function integrate_for_H(p::Photon, blackhole, acc_structure, detector)
    final_lmbda = 1.5 * detector.D
    lmbda = range(0, -final_lmbda, length=Int(7 * final_lmbda))
    tspan = (lmbda[1], lmbda[end])

    prob = ODEProblem(blackhole.geodesics, p.iC, tspan)
    sol = solve(prob, Tsit5(), saveat=lmbda)

    zi = [cos(u[3]) for u in sol.u]
    zil = circshift(zi, -1)
    zil[end] = 0.0

    indxs = findall(zi .* zil .< 0)
    
    solution = sol.u  # Por defecto, usar toda la solución

    for i in indxs
        current_sol = sol.u[i]
        r = current_sol[2]

        if r < acc_structure.out_edge && r > acc_structure.in_edge
            solution = sol.u[1:i]  # Usar solo hasta este punto
            break
        end
    end

    # Verificar horizonte de eventos
    r_vals = [u[2] for u in sol.u]
    indxsEH = findall(r_vals .< (blackhole.EH + 0.1))
    
    if !isempty(indxsEH)
        i = indxsEH[1]
        solution = sol.u[1:i]
    end

    H = Hamiltonian(solution, blackhole)
    println("Hamiltonian constraint verified = |H_max - H_0| = ", abs(maximum(H) - H[1]))
    return H
end

function Hamiltonian(sol::Vector{Vector{Float64}}, blackhole)
    H = zeros(length(sol))

    for i in eachindex(sol)
        x = sol[i][1:4]
        p_mom = sol[i][5:8]  # Renombrar para evitar conflicto

        gtt, grr, gthth, gphph, gtph = blackhole.inverse_metric(x)

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
    alpha_range = img.detector.alphaRange
    beta_range = img.detector.betaRange

    for (i, a) in enumerate(alpha_range)
        for (j, b) in enumerate(beta_range)
            p = Photon(a, b)
            p.iC = img.detector.photon_coords(img.blackhole, a, b)
            p.i = i
            p.j = j
            push!(img.photon_list, p)
        end
    end
end

function create_image!(img::Image)
    img.image_data = zeros(length(img.detector.alphaRange), length(img.detector.betaRange))
    
    println("Integrating trajectories with ", nthreads(), " threads")
    start_time = time()

    @threads for idx in eachindex(img.photon_list)
        p = img.photon_list[idx]
        img.image_data[p.i, p.j] = geodesics_integrate(p, img.blackhole, img.acc_structure, img.detector)
    end

    total_time = time() - start_time

    println("\n--- Total time of integration : ", total_time, " seconds ---")
    println("--- Time per photon : ", total_time / length(img.photon_list), " seconds/photon ---")
end

function create_image_no_Doppler!(img::Image)
    img.image_data = zeros(length(img.detector.alphaRange), length(img.detector.betaRange))
    
    println("Integrating trajectories with ", nthreads(), " threads")
    start_time = time()

    @threads for idx in eachindex(img.photon_list)
        p = img.photon_list[idx]
        img.image_data[p.i, p.j] = geo_inte_no_Doppler(p, img.blackhole, img.acc_structure, img.detector)
    end

    total_time = time() - start_time

    println("\n--- Total time of integration : ", total_time, " seconds ---")
    println("--- Time per photon : ", total_time / length(img.photon_list), " seconds/photon ---")
end

function create_shadow!(img::Image)
    img.image_data = zeros(length(img.detector.alphaRange), length(img.detector.betaRange))
    
    println("Integrating trajectories with ", nthreads(), " threads")
    start_time = time()

    @threads for idx in eachindex(img.photon_list)
        p = img.photon_list[idx]
        img.image_data[p.i, p.j] = shadow_integ(p, img.blackhole, img.detector)
    end

    total_time = time() - start_time

    println("\n--- Total time of integration : ", total_time, " seconds ---")
    println("--- EH radius: ", img.blackhole.EH)
    println("--- Time per photon : ", total_time / length(img.photon_list), " seconds/photon ---")
end

function save_data(img::Image, filename::String)
    save("$filename.jld2", "image_data", img.image_data)
end

function plot(img::Image; savefig=false, filename=nothing, cmap=:afmhot)
    if maximum(img.image_data) > 0
        image_data = img.image_data ./ maximum(img.image_data)
    else
        image_data = img.image_data
    end
    
    plt = heatmap(
        transpose(image_data),
        aspect_ratio=1,
        c=cmap,
        axis=false,
        framestyle=:none,
        yflip=false
    )
    
    if savefig && filename !== nothing
        mkpath("images")
        savefig(plt, "images/$(filename).png")
    end
    display(plt)
end

function plot_shadow(img::Image; savefig=false, filename=nothing, cmap=:grays)
    if maximum(img.image_data) > 0
        image_data = img.image_data ./ maximum(img.image_data)
    else
        image_data = img.image_data
    end
    
    plt = heatmap(
        transpose(image_data),
        aspect_ratio=1,
        c=cmap,
        axis=false,
        framestyle=:none,
        yflip=false
    )
    
    if savefig && filename !== nothing
        mkpath("images")
        savefig(plt, "images/$(filename).png")
    end
    display(plt)
end

function plot_contours(img::Image; savefig=false, filename=nothing, cmap=:grays)
    plt = contour(
        transpose(img.image_data),
        aspect_ratio=1,
        c=cmap,
        axis=false,
        framestyle=:none,
        yflip=false
    )
    
    if savefig && filename !== nothing
        mkpath("images")
        savefig(plt, "images/$(filename).png")
    end
    display(plt)
end

function verify_Hamiltonian(img::Image, n::Int=10)
    if !hasmethod(img.blackhole.inverse_metric, (Vector{Float64},))
        println("The inverse metric is not defined for this blackhole")
        println("PLEASE CHECK THE BLACKHOLE DEFINITION")
        return
    end

    println("Integrating trajectories...\n")

    plt = plot(
        xlabel="λ",
        ylabel="H",
        ylim=(-2, 2),
        grid=true,
        legend=true
    )

    for photon_id in 1:n
        i = rand(1:length(img.photon_list))
        p = img.photon_list[i]

        H = integrate_for_H(p, img.blackhole, img.acc_structure, img.detector)

        plot!(
            plt,
            H,
            label="Photon #$i"
        )
    end

    display(plt)
    println()
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    println()
    println("THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.")
    println("YOU NEED TO RUN THE main.jl FILE TO GENERATE THE IMAGE")
    println()
end
