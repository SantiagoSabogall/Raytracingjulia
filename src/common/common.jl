"""
    Common

Provides the core integration loop for non-symplectic `Tsit5` ordinary differential 
equation evaluation, observational logic (Doppler shifts, shadowing), and plotting handlers.
"""
module Common

using ..ImagePlane: photon_coords
using ..Raytracingjulia: metric, inverse_metric, Omega, geodesics


using DifferentialEquations
using StaticArrays
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

# ============================================================
# Photon
# ============================================================

"""
    Photon

Mutable container for a bundle of rays. Pre-allocated to prevent heap allocation 
inside the ODE solvers.

# Fields
- `alpha`, `beta`: Impact parameters.
- `freq`: Observed frequency.
- `i`, `j`: Screen pixel indices.
- `iC`: Initial state 8-vector (coordinate and momentum variables).
- `fP`: Final hit-point state on the accretion structure.
"""
mutable struct Photon
    alpha::Float64
    beta::Float64
    freq::Float64
    i::Int
    j::Int
    iC::SVector{8, Float64} # Strictly typed to SVector guaranteeing structural execution stability
    fP::SVector{8, Float64}

    function Photon(alpha::Float64, beta::Float64, freq::Float64=1.0)
        new(alpha, beta, freq, 0, 0, zeros(SVector{8}), zeros(SVector{8}))
    end
end

# ============================================================
# Geodesic integrators
# ============================================================
# ============================================================
# Geodesic integrators
# ============================================================
"""
    geodesics_integrate(p::Photon, blackhole, acc_structure, detector)

Integrates the null geodesics via `DifferentialEquations.jl` using a `Tsit5` solver.
Uses `ContinuousCallback` to accurately detect collisions with the accretion disk or 
the event horizon.

Returns the observed radiative intensity scaled by the Doppler factor \$g^3\$.
"""
function geodesics_integrate(p::Photon, blackhole, acc_structure, detector)

    r_max = 1.2 * detector.D
    tspan = (0.0, -1.5 * detector.D)

    hit_intensity = Ref(0.0)   # 🔥 FIX

    disc_cond(u, t, integrator) = cos(u[3])
    
    function disc_affect!(integrator)
        r = integrator.u[2]
        if acc_structure.in_edge < r < acc_structure.out_edge
            I0 = acc_structure.intensity(r)
            hit_intensity[] = doppler_shift(integrator.u, I0, blackhole)  # 🔥 FIX
            terminate!(integrator)
        end
    end

    horizon_cond(u, t, integrator) = u[2] - (blackhole.EH + 1e-5)
    horizon_affect!(integrator) = terminate!(integrator)

    escape_cond(u, t, integrator) = u[2] - r_max
    escape_affect!(integrator) = terminate!(integrator)

    cb = CallbackSet(
        ContinuousCallback(disc_cond, disc_affect!),
        ContinuousCallback(horizon_cond, horizon_affect!),
        ContinuousCallback(escape_cond, escape_affect!)
    )

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)

    sol = solve(prob, Tsit5(); 
        callback=cb, 
        reltol=1e-6, 
        abstol=1e-6, 
        save_everystep=false,
        verbose=false
    )

    return hit_intensity[]   # 🔥 FIX
end

function geodesics_integrate_no_Doppler(p::Photon, blackhole, acc_structure, detector)

    final_λ = 1.5 * detector.D
    lmbda_range = range(0.0, -final_λ, length=Int(7 * final_λ))
    tspan = (lmbda_range[1], lmbda_range[end])

    horizon_cond(u, t, integrator) = u[2] - (blackhole.EH + 1e-5)
    horizon_affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(horizon_cond, horizon_affect!)

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    sol = solve(prob, Tsit5(), 
                reltol=1e-8, 
                abstol=1e-8, 
                verbose=false,
                saveat=lmbda_range,
                callback=cb)

    p.fP = @SVector zeros(8)

    for i in 1:(length(sol.u) - 1)
        z_curr = cos(sol.u[i][3])
        z_next = cos(sol.u[i+1][3])
        
        if z_curr * z_next < 0
            u = sol.u[i]
            r = u[2]
            if acc_structure.in_edge < r < acc_structure.out_edge
                p.fP = u
                return acc_structure.intensity(r)
            end
        end
    end

    return 0.0
end

function shadow_integ(p::Photon, blackhole, detector)
    final_λ = 1.5 * detector.D
    tspan = (0.0, -final_λ)

    
    condition(u, t, integrator) = u[2] - (blackhole.EH + 1e-7)
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)
    

    sol = solve(prob, Tsit5(), 
                reltol=1e-8, 
                abstol=1e-8, 
                verbose=false,
                callback=cb,
                save_everystep=false)


    r_final = sol.u[end][2]
    
    if r_final < blackhole.EH + 1e-3 
        return 0.0   
    else
        return 100.0 
    end
end
# ============================================================
# Physics
# ============================================================

# ============================================================
# Physics
# ============================================================

"""
    doppler_shift(u_final, I0::Float64, blackhole)

Calculates the kinetic shift of the photon frequency:
\$g = \\frac{\\nu_{obs}}{\\nu_{em}} = \\frac{\\left(k_\\mu u^\\mu\\right)_{obs}}{\\left(k_\\nu u^\\nu\\right)_{em}}\$
Assuming the rest frame of the emitting particle rotates with Keplerian angular velocity \$\\Omega\$.

Returns the invariant Doppler-shifted intensity \$I_{obs} = I_0 g^3\$.
"""
function doppler_shift(u_final, I0::Float64, blackhole)
    coords = @SVector [u_final[1], u_final[2], u_final[3], u_final[4]]
    g_tt, _, _, g_phph, g_tph = metric(blackhole, coords)

    Ω = Omega(blackhole, u_final[2])
    !isfinite(Ω) && return 0.0

    k_t   = u_final[5]
    k_phi = u_final[8]

    # Doppler invariant redshift ratio assignment
    g = sqrt(abs(-g_tt - 2g_tph*Ω - g_phph*Ω^2)) / (1 + k_phi*Ω / k_t)

    return I0 * g^3
end

function integrate_for_H(p::Photon, blackhole, acc_structure, detector)
    final_lmbda = 1.5 * detector.D
    lmbda_range = range(0.0, -final_lmbda, length = Int(7 * final_lmbda))
    tspan = (lmbda_range[1], lmbda_range[end])

    horizon_cond_H(u, t, integrator) = u[2] - (blackhole.EH + 1e-4)
    horizon_affect_H!(integrator) = terminate!(integrator)
    cb_H = ContinuousCallback(horizon_cond_H, horizon_affect_H!)

    prob = ODEProblem(geodesics, p.iC, tspan, blackhole)

    sol = solve(prob, Tsit5(),
                reltol = 1e-8,
                abstol = 1e-8,
                saveat = lmbda_range,
                callback = cb_H)


    last_idx = length(sol.u)

    for i in 1:(length(sol.u) - 1)
        z_curr = cos(sol.u[i][3])
        z_next = cos(sol.u[i+1][3])
        
        if z_curr * z_next < 0
            r = sol.u[i][2]
            if acc_structure.in_edge < r < acc_structure.out_edge
                last_idx = i
                break
            end
        end
    end

    for i in 1:last_idx
        if sol.u[i][2] < (blackhole.EH + 0.05)
            last_idx = i
            break
        end
    end

    final_solution = @view sol.u[1:last_idx]
    final_lambdas = @view sol.t[1:last_idx]


    H = Hamiltonian(final_solution, blackhole)

    return final_lambdas, H
end



"""
    Hamiltonian(sol, blackhole)

Evaluates the Hamiltonian \$\\mathcal{H} = \\frac{1}{2} p_\\mu p^\\mu\$ along a geodesic trajectory.
This serves as a numerical verification constraint: analytic null geodesics require \$\\mathcal{H} \\approx 0\$.
"""
function Hamiltonian(sol, blackhole)

    H = zeros(length(sol))

    for i in eachindex(sol)
        u = sol[i]
        
        gtt, grr, gθθ, gφφ, gtφ = inverse_metric(blackhole, u)
        
        pt = u[5]
        pr = u[6]
        pth = u[7]
        pph = u[8]

        H[i] = 0.5 * (
            gtt * pt^2 +
            grr * pr^2 +
            gθθ * pth^2 +
            gφφ * pph^2 +
            2 * gtφ * pt * pph
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
    n_photons = length(img.photon_list)
    
    println("Integrating trajectories with ", nthreads(), " threads...")
    t0 = time()

    counter = Atomic{Int}(0)
    
    @threads for k in 1:n_photons
        p = img.photon_list[k]
        intensity = geodesics_integrate(p, img.blackhole, img.acc_structure, img.detector)
        img.image_data[p.i, p.j] = intensity
        
        c = atomic_add!(counter, 1) + 1
        if c % 100 == 0 || c == n_photons
            print("\rPhoton # $c / $n_photons")
            flush(stdout)
        end
    end

    total_time = time() - t0
    println("\n\nTotal time: ", round(total_time, digits=2), " s")
    println("Efficiency: ", total_time / n_photons, " s/photon")
end
function create_image_no_Doppler!(img::Image)

   
    img.image_data = zeros(img.detector.x_pixels, img.detector.y_pixels)
    n_photons = length(img.photon_list)

    photon = 1
    println("Integrating trajectories (No Doppler)...")

    start_time = time()

    for p in img.photon_list
        img.image_data[p.i, p.j] =
            geodesics_integrate_no_Doppler(p, img.blackhole, img.acc_structure, img.detector)

        print("\rPhoton # $photon / $n_photons")
        flush(stdout)

        photon += 1
    end

    total_time = time() - start_time

    println("\n\n--- Total time of integration : $total_time seconds ---")
    println(
        "\n--- Time of integration : $(total_time / n_photons) seconds/photon ---\n"
    )
end


function create_shadow!(img::Image)
    
    img.image_data = zeros(img.detector.x_pixels, img.detector.y_pixels)
    
    n_photons = length(img.photon_list)
    println("Integrating trajectories for shadow on $(nthreads()) threads...")

    start_time = time()

    counter = Atomic{Int}(0)
    
    @threads for p in img.photon_list
        
        img.image_data[p.i, p.j] = shadow_integ(p, img.blackhole, img.detector)
        
        c = atomic_add!(counter, 1) + 1
        if c % 100 == 0 || c == n_photons
            print("\rPhoton # $c / $n_photons")
            flush(stdout)
        end
    end

    total_time = time() - start_time
    
    
    println("\n\nEH radius: $(img.blackhole.EH)")
    println("--- Total time: $(round(total_time, digits=2)) seconds ---")
    println("--- Speed: $(round(total_time/n_photons, digits=6)) seconds/photon ---")
end

function plot(img::Image; save=false, filename=nothing, cmap=:afmhot)
    
    maxv = maximum(img.image_data)
    data = maxv > 0 ? img.image_data ./ maxv : img.image_data

    # 2. Render resulting observable spatial heatmap
    plt = heatmap(
        transpose(data),
        aspect_ratio = 1,
        c = cmap,
        axis = false,
        framestyle = :none,
        colorbar = false
    )

   
    if save && !isnothing(filename)
        mkpath("images")
       
        Plots.savefig(plt, "images/$filename.png")
        println("Rendering saved to: images/$filename.png")
    end

    display(plt)
    return plt
end

function plot_shadow(img; save=false, filename=nothing, cmap=:gray)

    data = copy(img.image_data)
    

    maxv = maximum(data)
    if maxv > 0
        data .= data ./ maxv
    end

    plt = heatmap(
        data',
        aspect_ratio = :equal,
        c = cmap,
        clims = (0, 1), 
        axis = true,       
        framestyle = :box,
        colorbar = false,
        ticks = false      
    )


    xlabel!(plt, "α")
    ylabel!(plt, "β")


    if save && filename !== nothing
        if !isdir("images") mkpath("images") end

        full_path = endswith(filename, ".png") ? "images/$filename" : "images/$filename.png"
        Plots.savefig(plt, full_path)
        println("Shadow mapped rendering to: $full_path")
    end

    display(plt)
    return plt
end


function plot_contours(img; save=false, filename=nothing, cmap=:gray)

    data = img.image_data'

    plt = contour(
        data,
        aspect_ratio = :equal,
        c = cmap,
        fill = false,     
        levels = 15,       
        lw = 1.2,          
        colorbar = false
    )


    xlabel!(plt, "α")
    ylabel!(plt, "β")
    plot!(plt, grid=true, gridalpha=0.25, ticks=false)


    if save && filename !== nothing
        if !isdir("images") mkpath("images") end
        full_path = endswith(filename, ".png") ? "images/$filename" : "images/$filename.png"
        Plots.savefig(plt, full_path)
        println("Topological contours cached to: $full_path")
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
        gridalpha = 0.3,
        ylims = (-1, 1)
    )

    indices = rand(1:length(img.photon_list), n)

    for i in indices
        p = img.photon_list[i]

        lambdas, H_vals = integrate_for_H(p, img.blackhole, img.acc_structure, img.detector)

        Plots.plot!(plt, lambdas, H_vals, label = "Photon #$i")

        println("Hamiltonian constraint verified: |H_max - H_0| = ", abs(maximum(H_vals) - H_vals[1]))
    end

    display(plt)
end
end 
