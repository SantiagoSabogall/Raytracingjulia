using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Raytracingjulia
using LinearAlgebra
using NPZ

# =========================
# BLACK HOLE DEFINITION
# =========================
a = 0.5
blackhole = KerrBH(a)

# =========================
# DETECTOR PARAMETERS
# =========================
D = 100.0
iota = 85 * pi / 180
x_side = 25.0
x_pixels = 1920

detector = Detector(
    D,
    iota,
    x_side;
    x_pixels = x_pixels,
    ratio = "16:9"
)

# =========================
# ACCRETION STRUCTURE
# =========================
acc_structure = ThinDisk(blackhole)

# =========================
# IMAGE GENERATION
# =========================
image = Image(blackhole, acc_structure, detector)

image.create_photons()
image.create_image()

# =========================
# SAVE DATA
# =========================
filename = "Kerr_a_0.5_1920x1080"

mkpath("images_data")
NPZ.npy_save("images_data/$filename.npy", image.image_data)

image.plot(savefig=true, filename=filename)
