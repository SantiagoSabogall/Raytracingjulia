using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Revise
using Raytracingjulia
using LinearAlgebra
using NPZ

# Configuración del agujero negro
# ----------------------------
a = 0.5
Mpsi = 0.5 # Parámetro M de scalaron (Float64)
blackhole = KerrScalaronBH(a, Mpsi)

# ----------------------------
# Detector
# ----------------------------
D = 100.0
iota = 85 * pi / 180
x_side = 25.0
x_pixels = 300

detector = Detector(
    D,
    iota,
    x_side;
    x_pixels = x_pixels,
    ratio = "16:9"
)

# ----------------------------
# Estructura de acreción
# ----------------------------
#  Estructura de acrecion para KerrPFDM
acc_structure = SimpleDisk(
    blackhole;
    R_min = 3.0
)

# acc_structure = ThinDisk(
#     blackhole
# )
# ----------------------------
# Imagen
# ----------------------------
image = Image(blackhole, acc_structure, detector)


create_photons!(image)


create_image!(image)

filename = "KerrScalaron_a_0.5_1920"

println("Mostrando imagen de prueba...")
plot(image, save=true, filename=filename)
