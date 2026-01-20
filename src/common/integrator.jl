module Integrators

export rk45, rk45_optimized

using LinearAlgebra

"""
    rk45(f, t0, y0, t1; kwargs...)

Integrate y'(t) = f(t, y) from `t0` to `t1` using adaptive Dormand–Prince RK45.

Returns `(T, Y)` where:
- `T` is a vector of times
- `Y` is a vector of state vectors
"""
function rk45(f, t0::Float64, y0::Vector{Float64}, t1::Float64;
              atol::Float64=1e-9, rtol::Float64=1e-9, h0::Float64=1e-2,
              h_min::Float64=1e-12, h_max::Float64=1.0,
              max_steps::Int=1_000_000,
              t_eval::Union{Nothing, Vector{Float64}}=nothing)

    y = copy(y0)
    n = length(y)
    t = t0
    forward = t1 ≥ t
    sgn = forward ? 1.0 : -1.0

    using_eval = t_eval !== nothing
    if using_eval
        T_eval = copy(t_eval)
        Ti = Float64[]
        Yi = Vector{Float64}[]
    end

    C = (0.0, 1//5, 3//10, 4//5, 8//9, 1.0, 1.0)
    A = (
        (),
        (1//5,),
        (3//40, 9//40),
        (44//45, -56//15, 32//9),
        (19372//6561, -25360//2187, 64448//6561, -212//729),
        (9017//3168, -355//33, 46732//5247, 49//176, -5103//18656),
        (35//384, 0, 500//1113, 125//192, -2187//6784, 11//84),
    )

    B5 = (35//384, 0, 500//1113, 125//192, -2187//6784, 11//84, 0)
    B4 = (5179//57600, 0, 7571//16695, 393//640, -92097//339200, 187//2100, 1//40)

    h = clamp(h0, h_min, h_max) * sgn
    T = Float64[t]
    Y = Vector{Float64}[copy(y)]

    accepted = 0

    while (forward ? t < t1 : t > t1) && accepted < max_steps
        h = clamp(abs(h), h_min, h_max) * sgn
        if forward && t + h > t1 || !forward && t + h < t1
            h = t1 - t
        end

        k = Vector{Vector{Float64}}(undef, 7)
        k[1] = f(t, y)

        for i in 2:7
            yi = copy(y)
            for j in 1:length(A[i])
                yi .+= h * float(A[i][j]) .* k[j]
            end
            k[i] = f(t + float(C[i]) * h, yi)
        end

        y5 = copy(y)
        y4 = copy(y)
        for i in 1:7
            y5 .+= h * float(B5[i]) .* k[i]
            y4 .+= h * float(B4[i]) .* k[i]
        end

        err = y5 .- y4
        en = sqrt(sum((err[i] / (atol + rtol * max(abs(y[i]), abs(y5[i]))))^2 for i in 1:n) / n)

        if en ≤ 1.0
            t += h
            y = y5
            accepted += 1
            push!(T, t)
            push!(Y, copy(y))
            h *= clamp(0.9 * en^(-0.2), 0.2, 5.0)
        else
            h *= clamp(0.9 * en^(-0.25), 0.2, 1.0)
        end
    end

    return T, Y
end


"""
    rk45_optimized(f!, t0, y0, t1; kwargs...)

In-place optimized RK45 integrator.
"""
function rk45_optimized(f!, t0::Float64, y0::Vector{Float64}, t1::Float64;
                        atol::Float64=1e-9, rtol::Float64=1e-9, h0::Float64=1e-2,
                        h_min::Float64=1e-12, h_max::Float64=1.0,
                        max_steps::Int=1_000_000)

    n = length(y0)
    y = copy(y0)
    t = t0
    forward = t1 ≥ t
    sgn = forward ? 1.0 : -1.0

    k = [zeros(n) for _ in 1:7]
    ytmp = zeros(n)
    y5 = zeros(n)
    y4 = zeros(n)
    err = zeros(n)

    h = clamp(h0, h_min, h_max) * sgn
    T = Float64[t]
    Y = Vector{Float64}[copy(y)]

    accepted = 0

    while (forward ? t < t1 : t > t1) && accepted < max_steps
        h = clamp(abs(h), h_min, h_max) * sgn
        if forward && t + h > t1 || !forward && t + h < t1
            h = t1 - t
        end

        f!(k[1], t, y)

        for i in 2:7
            copy!(ytmp, y)
            for j in 1:length(A[i])
                @. ytmp += h * float(A[i][j]) * k[j]
            end
            f!(k[i], t + float(C[i]) * h, ytmp)
        end

        copy!(y5, y)
        copy!(y4, y)
        for i in 1:7
            @. y5 += h * float(B5[i]) * k[i]
            @. y4 += h * float(B4[i]) * k[i]
        end

        @. err = y5 - y4
        en = sqrt(sum((err[i] / (atol + rtol * max(abs(y[i]), abs(y5[i]))))^2 for i in 1:n) / n)

        if en ≤ 1.0
            t += h
            copy!(y, y5)
            accepted += 1
            push!(T, t)
            push!(Y, copy(y))
            h *= clamp(0.9 * en^(-0.2), 0.2, 5.0)
        else
            h *= clamp(0.9 * en^(-0.25), 0.2, 1.0)
        end
    end

    return T, Y
end


if abspath(PROGRAM_FILE) == @__FILE__
    println("This file defines numerical integrators.")
    println("Run main.jl to generate physical results.")
end

end # module Integrators
