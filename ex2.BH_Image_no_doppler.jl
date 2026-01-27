############################################################
# run_schwarzschild_num.jl
############################################################

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# ----------------------------
# Revise para recargar cambios sin reiniciar Julia
# ----------------------------
using Revise
using Raytracingjulia
using LinearAlgebra
using NPZ

############################################################
# Configuración del agujero negro numérico (Schwarzschild)
############################################################

# Aquí pasamos las rutas de tus datos numéricos
path_N   ="src/black_holes/numerical_data/schwarzschild_data/N.txt"
path_dN  = "src/black_holes/numerical_data/schwarzschild_data/derN.txt"

blackhole = SchwarzschildNumBH(path_N, path_dN)

############################################################
# Detector
############################################################
D = 100.0
iota = 85 * pi / 180
x_side = 25.0
x_pixels = 80

detector = Detector(
    D,
    iota,
    x_side;
    x_pixels = x_pixels,
    ratio = "16:9"
)

############################################################
# Estructura de acreción
############################################################
acc_structure = ThinDisk.ThinDisk(blackhole)

############################################################
# Imagen
############################################################
image = Image(blackhole, acc_structure, detector)

############################################################
# Crear fotones y trazar rayos
############################################################
println("Generando fotones...")
create_photons!(image)

println("Trazando rayos...")
create_image!(image)

############################################################
# Guardar imagen
############################################################
filename = "num_schwarzschild_80x45_thin_disk.png"
println("Mostrando imagen de prueba y guardando como $filename ...")
plot(image, save=true, filename=filename)
