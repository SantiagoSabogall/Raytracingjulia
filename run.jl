using Pkg
Pkg.activate(@__DIR__)
# Pkg.instantiate() # Solo necesitas correr esto la primera vez o si cambias dependencias

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
x_pixels = 900

# CORRECCIÓN 1: Usamos el nombre del Struct (probablemente ImagePlane), no del Módulo
detector = Detector.Detector( 
    D,
    iota,
    x_side;
    x_pixels = x_pixels,
    ratio = "16:9"
)

# =========================
# ACCRETION STRUCTURE
# =========================
acc_structure = ThinDisk.ThinDisk(blackhole)

# =========================
# IMAGE GENERATION
# =========================
image = Image(blackhole, acc_structure, detector)

# CORRECCIÓN 2: Sintaxis de Julia (Funciones con signo de exclamación)
println("Generando fotones...")
create_photons!(image)

println("Trazando rayos...")
create_image!(image)

# =========================
# SAVE DATA
# =========================
#filename = "Kerr_a_0.5_1920x1080"

#mkpath("images_data")
# Accedemos al campo .image_data del struct image
#NPZ.npy_save("images_data/$filename.npy", image.image_data)

# CORRECCIÓN 3: Plotting
#println("Guardando imagen...")
#plot(image, savefig=true, filename=filename)

println("Mostrando imagen de prueba...")

# savefig=false evita que se cree el archivo .png
# filename no es necesario si savefig es false
plot(image, savefig=false)