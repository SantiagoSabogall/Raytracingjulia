using Plots
using LinearAlgebra

# 1. Definición del perfil de masa PFDM y sus derivadas
# m(r) = M - (k/2) * log(r/|k|)
m(r, k)  = 1.0 - (k / 2.0) * log(r / abs(k))
dm(r, k) = -k / (2.0 * r)

# 2. Función simplificada para encontrar el ISCO (basada en tu código)
function find_isco_value(a, k)
    # Horizonte de eventos aproximado (para empezar la búsqueda)
    # En Kerr-PFDM el horizonte cambia con k, aquí usamos una estimación
    r_start = 1.0 + sqrt(max(0.0, 1.0 - a^2)) + 0.1 
    
    # Función de energía para órbita circular ecuatorial
    function energy_circular(r)
        mr = m(r, k)
        dmr = dm(r, k)
        Σ = r^2
        
        # Derivadas de la métrica (θ = π/2)
        dgtt = 2*(dmr*r + mr)/Σ - 4*mr*r^2/Σ^2
        dgtphi = -2a*(dmr*r + mr)/Σ + 4a*mr*r^2/Σ^2
        dgphiphi = 2r + 2a^2*(dmr*r + mr)/Σ - 4a^2*mr*r^2/Σ^2
        
        # Velocidad angular Ω (prograda)
        disc = dgtphi^2 - dgtt*dgphiphi
        disc < 0 && return NaN
        Ω = (-dgtphi + sqrt(disc)) / dgphiphi
        
        # Componentes de la métrica
        gtt = -(1.0 - 2mr*r/Σ)
        gtphi = -2a*mr*r/Σ
        gphi = (r^2 + a^2 + 2*a^2*mr*r/Σ)
        
        # Energía E
        denom = -(gtt + 2gtphi*Ω + gphi*Ω^2)
        denom <= 0 && return NaN
        ut = 1.0 / sqrt(denom)
        return -(gtt + gtphi*Ω) * ut
    end

    # Búsqueda del mínimo de E (ISCO)
    rs = range(r_start, 100.0, length=2000)
    E_prev = Inf
    for r in rs
        E = energy_circular(r)
        if isnan(E) continue end
        if E > E_prev return r end # El mínimo ocurre justo antes de que suba
        E_prev = E
    end
    return NaN
end

# 3. Generación de la gráfica
ks = range(0.001, 1.2, length=100)
as = [0.3, 0.5, 0.7, 0.8]
colores = [:blue, :orange, :green, :red]

plt = plot(title="The influence of k on the ISCO", xlabel="k", ylabel="ISCO", 
           xlim=(0, 1.4), ylim=(0, 26), legend=:topright)

for (i, a_val) in enumerate(as)
    isco_vals = [find_isco_value(a_val, k_val) for k_val in ks]
    # Filtramos NaNs para una línea limpia
    plot!(plt, ks, isco_vals, label="a=$a_val", color=colores[i], lw=2)
end

display(plt)