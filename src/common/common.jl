using DifferentialEquations

struct Photon
    alpha::Float64
    beta::Float64
    freq::Float64

    i::Union{Float64, nothing}
    j::Union{Float64, nothing}

    iC::Union{Vector{Float64}, nothing}
    fP::Union{Vector{Float64}, nothing}

    function Photon(alpha::Float64, beta::Float64,freq::Float64=1.0)
    # Pixel coordinates 

    i = nothing
    j = nothing

    # Initial Cartesian Coordinates in the image plane

    iC = nothing

    # Stores the final values of coordinates and Momentum

    fP = nothing

    new(alpha, beta, freq, nothing, nothing, nothing, nothing)

    end
end

function geodesics_integrate(p, blackhole, acc_structure, detector)

    final_lmbda = 1.5 * detector.D
    lmbda = range(0, -final_lmbda, length= Int(7 * final_lmbda))
    tspan = (lmbda[1], lmbda[end])

    prob = ODEProblem(blackhole.geodesics, p.iC, tspan)
    sol = solve(prob, Tsit5(), saveat=lmbda)

    p.fP = zeros(8)
    I_f = 0.0
    zi = [cos(u[3]) for u in sol.u]
    zil = circshift(zi, -1)
    zil[end] = 0.

    indxs = findall(zi .* zil .< 0)

    for i in indxs
        current_sol = sol[i]
        r = current_sol[2]

        if r < acc_structure.out_edge && r > acc_structure.in_edge
            p.fP = current_sol
            I_0 = acc_structure.intensity(p.fP[2])
            I_f = doppler_shift(p, I_0, blackhole)
            break        
    end

    return I_f

end

function geo_inte_no_Doppler(p, blackhole, acc_structure, detector)

    final_lmbda = 1.5 * detector.D
    lmbda = range(0, -final_lmbda, length= Int(7 * final_lmbda))
    tspan = (lmbda[1], lmbda[end])

    prob = ODEProblem(blackhole.geodesics, p.iC, tspan)
    sol = solve(prob, Tsit5(), saveat=lmbda)

    p.fP = zeros(8)
    I_f = 0.0
    zi = [cos(u[3]) for u in sol.u]
    zil = circshift(zi, -1)
    zil[end] = 0.

    indxs = findall(zi .* zil .< 0)

    for i in indxs
        current_sol = sol[i]
        r = current_sol[2]

        if r < acc_structure.out_edge && r > acc_structure.in_edge
            p.fP = current_sol
            break       
        end 
    end

    return acc_structure.intensity(p.fP[2])


end


function shadow_integ(p, blackhole, detector)

    final_lmbda = 1.5 * detector.D
    lmbda = range(0, -final_lmbda, length= Int(7 * final_lmbda))
    tspan = (lmbda[1], lmbda[end])

    prob = ODEProblem(blackhole.geodesics, p.iC, tspan)
    sol = solve(prob, Tsit5(), saveat=lmbda) 
    sol_2 = [u[2] for u in sol.u]

    indxs = findall(sol_2 .< blackhole.EH .+ 1e-7 )

    if length(indxs) == 0
        return 100
    else
        return 0 
    end

end

function doppler_shift(p, I0, blackhole)

    g_tt, _, _, g_phph, g_tph = blackhole.metric(p.fP[:5])
    Omega = blackhole.Omega(p.fP[2])
    g = sqrt(- g_tt - 2*g_tph*Omega - g_phph*Omega^2)/(1 + p.fP[8]*Omega/p.fP[5])

    return I0 * g^3
end


function integrate_for_H(p, blackhole, acc_structure, detector)

    final_lmbda = 1.5 * detector.D
    lmbda = range(0, -final_lmbda, length= Int(7 * final_lmbda))
    tspan = (lmbda[1], lmbda[end])

    prob = ODEProblem(blackhole.geodesics, p.iC, tspan)
    sol = solve(prob, Tsit5(), saveat=lmbda)


    zi = [cos(u[3]) for u in sol.u]
    zil = circshift(zi, -1)
    zil[end] = 0.

    indxs = findall(zi .* zil .< 0)

    for i in indxs
    current_sol = sol[i]
    r = current_sol[2]

        if r < acc_structure.out_edge && r > acc_structure.in_edge
            solution = current_sol
            break       
        end 
    end

    sol_2 = [u[2] for u in sol.u]
    indxsEH = findall(sol_2 .< blackhole.Eh .+0.1)

    for i in indxsEH
        current_sol = sol[i]
        solution = current_sol

    end

    H = Hamilationian(solution, blackhole)
    println("Hamilationian constraint verified = |H_max - H_0| = ", abs(H.max() - H[1]))
    return H
        


function Hamilationian(sol, blackhole):
    H = zeros(length(sol))

    for i in range(length(sol))
        x = 
        p = 
        gtt, grr, gthth, gphph, gtph = blackhole.inverse_metric(x)
        H[i] = 0.5*(gtt*p[1]*p[1] + grr*p[2]*p[2] + gthth*p[3]*p[3] + gphph*p[4]*p[4] + 2*gtph*p[1]*p[4]) 

    end
    return H
end



mutable struct image

    blackhole::
    acc_structure::
    detector::


end