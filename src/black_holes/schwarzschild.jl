module Schwarzschild

# IMPORTANTE: Usamos 'import' para poder extender las funciones del padre
import ..Raytracingjulia: metric, inverse_metric, geodesics, Omega

using StaticArrays

export SchwarzschildBH

# Inmutable: los parámetros no cambian tras la construcción.
# Permite al compilador optimizar el acceso a campos sin boxing.
struct SchwarzschildBH
    a::Float64
    EH::Float64
    ISCOco::Float64
    ISCOcounter::Float64

    function SchwarzschildBH(M::Float64=1.0)
        new(0.0, 2.0*M, 6.0*M, 6.0*M)
    end
end

@inline function Omega(b::SchwarzschildBH, r::Real; corotating::Bool=true)
    return 1.0 / (r * sqrt(r))   # r^1.5 sin exponente flotante
end

@inline function metric(b::SchwarzschildBH, x::AbstractVector)
    r  = x[2]; θ = x[3]
    f  = 1.0 - 2.0/r
    sinθ = sin(θ)
    r2   = r * r
    return -f, 1.0/f, r2, r2 * sinθ * sinθ, 0.0
end

@inline function inverse_metric(b::SchwarzschildBH, x::AbstractVector)
    r    = x[2]; θ = x[3]
    f    = 1.0 - 2.0/r
    sinθ = sin(θ)
    r2   = r * r
    return -1.0/f, f, 1.0/r2, 1.0/(r2 * sinθ * sinθ), 0.0
end

# Out-of-place: retorna SVector para cero alocaciones en el heap
@inline function geodesics(q, b::SchwarzschildBH, λ)
    @inbounds begin
        r    = q[2]
        θ    = q[3]
        k_t  = q[5]
        k_r  = q[6]
        k_th = q[7]
        k_φ  = q[8]
    end

    sinθ, cosθ = sincos(θ)
    f     = 1.0 - 2.0/r
    r2    = r * r
    sinθ2 = sinθ * sinθ
    inv_r2    = 1.0 / r2
    inv_rm2sq = 1.0 / ((r - 2.0) * (r - 2.0))

    dtdλ = -k_t / f
    drdλ =  f * k_r
    dθdλ =  k_th * inv_r2
    dφdλ =  k_φ  * inv_r2 / sinθ2

    dk_t  = 0.0
    dk_r  = -inv_r2 * dtdλ^2 + inv_rm2sq * drdλ^2 + r * dθdλ^2 + r * sinθ2 * dφdλ^2
    dk_th = sinθ * cosθ * dφdλ^2
    dk_φ  = 0.0

    return SVector{8, Float64}(dtdλ, drdλ, dθdλ, dφdλ, dk_t, dk_r, dk_th, dk_φ)
end

end # module Schwarzschild