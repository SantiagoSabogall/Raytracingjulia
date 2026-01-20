module Integrado

using DifferentialEquations
using Base: circshift

using ..Common   # ← IMPORTANTE



export geodesics_integrate,
       geodesics_integrate_no_doppler,
       shadow_integrate,
       integrate_for_H

# ============================================================
# Integración general de geodésicas
# ============================================================

function _solve_geodesic(p, blackhole, detector;
                         abstol=1e-9, reltol=1e-9)

    λmax = 1.5 * detector.D
    tspan = (0.0, -λmax)

    # Wrapper compatible con DifferentialEquations.jl
    f(u, λ) = geodesics(blackhole, u, λ)

    prob = ODEProblem(
        f,
        p.iC,
        tspan
    )

    return solve(
        prob,
        Vern9();
        abstol = abstol,
        reltol = reltol,
        save_everystep = true
    )
end

# ============================================================
# Imagen con Doppler
# ============================================================

function geodesics_integrate(p, blackhole, acc_structure, detector)

    sol = _solve_geodesic(p, blackhole, detector)

    zi  = cos.(getindex.(sol.u, 3))
    zil = circshift(zi, -1)
    zil[end] = 0.0

    crossings = findall(zi .* zil .< 0)

    for i in crossings
        u = sol.u[i]
        r = u[2]

        if acc_structure.in_edge < r < acc_structure.out_edge
            p.fP = u
            I0 = acc_structure.intensity(r)
            return doppler_shift(p, I0, blackhole)
        end
    end

    return 0.0
end

# ============================================================
# Imagen sin Doppler
# ============================================================

function geodesics_integrate_no_doppler(p, blackhole, acc_structure, detector)

    sol = _solve_geodesic(p, blackhole, detector)

    zi  = cos.(getindex.(sol.u, 3))
    zil = circshift(zi, -1)
    zil[end] = 0.0

    crossings = findall(zi .* zil .< 0)

    for i in crossings
        u = sol.u[i]
        r = u[2]

        if acc_structure.in_edge < r < acc_structure.out_edge
            p.fP = u
            return acc_structure.intensity(r)
        end
    end

    return 0.0
end

# ============================================================
# Sombra del agujero negro
# ============================================================

function shadow_integrate(p, blackhole, detector)

    sol = _solve_geodesic(p, blackhole, detector)
    r_vals = getindex.(sol.u, 2)

    any(r_vals .< blackhole.EH + 1e-7) && return 0.0
    return 1.0
end

# ============================================================
# Hamiltoniano (chequeo de conservación)
# ============================================================

function integrate_for_H(p, blackhole, acc_structure, detector)

    sol = _solve_geodesic(p, blackhole, detector)
    return Hamiltonian(sol.u, blackhole)
end

end # module Integrado
