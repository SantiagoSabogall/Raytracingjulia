using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Raytracingjulia

# La primera ejecución incluye compilación; repite para comparar tiempos.
blackhole = KerrPFDMBH(0.5, 0.5)
detector = Detector(100.0, 85 * pi / 180, 25.0;
                    x_pixels=90, ratio="16:9")
disk = SimpleDisk(blackhole; R_min=3.0)
image = Image(blackhole, disk, detector)

create_photons!(image)
elapsed = @elapsed create_image!(image)
println("Tiempo de create_image!: $(round(elapsed, digits=3)) s")
println("Tiempo por rayo: $(round(elapsed / length(image.photon_list), sigdigits=4)) s")
