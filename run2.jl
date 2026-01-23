# 1. Cargar tus módulos (ajusta las rutas si es necesario)
include("src/black_holes/kerr.jl")
include("src/detectors/image_plane.jl")
include("src/common/common.jl")
include("src/accretion_structures/thin_disk.jl")

using .Kerr, .Detector, .Common, .ThinDisk
using DifferentialEquations

# 2. Configuración mínima igual a la de tu script de Python
a = 0.5
bh = Kerr.BlackHole(a)
dist = 100.0
inc = (pi/180) * 85.0
det = Detector.detector(D=dist, iota=inc, x_pixels=400, x_side=25)

# 3. EL TEST (Copia y pega esto)
test_alpha = 5.0 
test_beta = 2.0  # Un poco de altura para que cruce el plano

println("--- INICIANDO DEBUG DE FOTÓN ---")
p = Common.Photon(test_alpha, test_beta)
# Generar condiciones iniciales
p.iC = Detector.photon_coords(det, bh, test_alpha, test_beta)

println("Condiciones Iniciales (iC): ", p.iC)

final_λ = 1.5 * det.D
tspan = (0.0, -final_λ)

# Usamos la función de geodésicas que definiste en Kerr
prob = ODEProblem(Kerr.geodesics, p.iC, tspan, bh)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# 4. Verificación de cruce de plano (Theta)
thetas = [u[3] for u in sol.u]
println("Rango de theta: $(minimum(thetas)) a $(maximum(thetas))")
# El plano ecuatorial es pi/2 ≈ 1.5708
cruzó = any(thetas .> 1.5708) && any(thetas .< 1.5708)
println("¿Cruzó el plano ecuatorial (pi/2)? ", cruzó)

# 5. Verificación de radios
radios = [u[2] for u in sol.u]
println("Radio inicial: ", radios[1])
println("Radio mínimo alcanzado: ", minimum(radios))

# 6. Verificación de impacto con el disco
# Intentamos correr la lógica de integración completa para este fotón
acc = ThinDisk.ThinDisk(bh)
I_f = Common.geodesics_integrate(p, bh, acc, det)
println("Intensidad final detectada: ", I_f)