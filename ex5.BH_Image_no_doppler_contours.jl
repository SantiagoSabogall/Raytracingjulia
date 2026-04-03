using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Revise
using Raytracingjulia
using LinearAlgebra
using NPZ

# ----------------------------
# Black Hole Configuration
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
# Image Processing
# ----------------------------
image = Image(blackhole, acc_structure, detector)


create_photons!(image)


create_image_no_Doppler!(image)

filename = "Kerr_a_0.5_900_thin_disk_test.png"

println("Rendering test image...")
plot_contours(image, save=true, filename=filename)
