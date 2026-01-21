module Schwarzschild

# IMPORTANTE: Usamos 'import' para poder extender las funciones del padre
import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega

export SchwarzschildBH

mutable struct SchwarzschildBH
    a::Float64
    EH::Float64
    ISCOco::Float64
    ISCOcounter::Float64

    function SchwarzschildBH(M::Float64=1.0)
        new(0.0, 2.0*M, 6.0*M, 6.0*M)
    end
end

# Ahora puedes definirlas sin errores de precompilación
function Omega(b::SchwarzschildBH, r::Real; corotating::Bool=true)
    return 1.0 / r^1.5
end

function metric(b::SchwarzschildBH, x::Vector{<:Real})
    r = x[2]; θ = x[3]
    f = 1.0 - 2.0/r
    return [-f, 1.0/f, r^2, (r*sin(θ))^2, 0.0]
end

function inverse_metric(b::SchwarzschildBH, x::Vector{<:Real})
    r = x[2]; θ = x[3]
    f = 1.0 - 2.0/r
    return [-1.0/f, f, 1.0/r^2, 1.0/(r*sin(θ))^2, 0.0]
end

function geodesics(q::Vector{<:Real}, b::SchwarzschildBH, λ::Real)
    # q = [t, r, θ, φ, k_t, k_r, k_th, k_φ]
    r = q[2]; θ = q[3]
    k_t, k_r, k_th, k_φ = q[5:8]
    f = 1.0 - 2.0/r

    # Coordenadas
    dtdλ = -k_t / f
    drdλ = f * k_r
    dθdλ = k_th / r^2
    dφdλ = k_φ / (r*sin(θ))^2

    # Momentos (Ecuaciones simplificadas para Schwarzschild)
    dk_t = 0.0
    dk_r = -(1.0/r^2)*dtdλ^2 + (1.0/(r-2.0)^2)*drdλ^2 + r*dθdλ^2 + r*sin(θ)^2*dφdλ^2
    dk_th = sin(θ)*cos(θ)*dφdλ^2
    dk_φ = 0.0

    return [dtdλ, drdλ, dθdλ, dφdλ, dk_t, dk_r, dk_th, dk_φ]
end

end