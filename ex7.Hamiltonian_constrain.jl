using Pkg
Pkg.activate(@__DIR__)
#Pkg.instantiate()

using Revise
using Raytracingjulia


# ----------------------------
# Configuración del agujero negro
# ----------------------------
a = 0.5
blackhole = KerrBH(a)

# ----------------------------
# Detector
# ----------------------------
D = 100.0
iota = 85 * pi / 180
x_side = 25.0
x_pixels = 90

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
acc_structure = ThinDisk.ThinDisk(blackhole)

# ----------------------------
# Imagen
# ----------------------------
image = Image(blackhole, acc_structure, detector)


create_photons!(image)



verify_Hamiltonian(image, n=10)
