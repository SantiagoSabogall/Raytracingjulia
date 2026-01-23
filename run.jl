using Pkg
Pkg.activate(@__DIR__)
# Pkg.instantiate() # Solo necesitas correr esto la primera vez o si cambias dependencias
using Raytracingjulia
using LinearAlgebra
using NPZ




a = 0.5
blackhole = KerrBH(a)


D =100.0
iota = 85 * pi / 180
x_side = 25.0
x_pixels = 100


detector = Detector.Detector( 
    D,
    iota,
    x_side;
    x_pixels = x_pixels,
    ratio = "16:9"
)                                                                                                                   


acc_structure = ThinDisk.ThinDisk(blackhole)


image = Image(blackhole, acc_structure, detector)

println("Generando fotones...")
create_photons!(image)

println("Trazando rayos...")
create_image!(image)

filename = "Kerr_a_0.5_100"

#mkpath("images_data")
# Accedemos al campo .image_data del struct image
#NPZ.npy_save("images_data/$filename.npy", image.image_data)

# CORRECCIÓN 3: Plotting
#println("Guardando imagen...")
#plot(image, savefig=true, filename=filename)

println("Mostrando imagen de prueba...")

plot(image, filename=filename)

