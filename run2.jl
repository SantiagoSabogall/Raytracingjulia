using Pkg
Pkg.activate(@__DIR__)

using Raytracingjulia
using LinearAlgebra
using NPZ
using Plots

# ============================================================
# 1. OTRO AGUJERO NEGRO: Schwarzschild (a = 0)
# ============================================================
# Si no tienes un struct SchwarzschildBH, puedes usar Kerr con a = 0.0001
bh = SchwarzschildBH() 

# =========================
# 2. CONFIGURACIÓN DEL DETECTOR
# =========================
D = 40.0                # Más cerca para reducir errores de propagación
iota = 45 * pi / 180    # Inclinación moderada (45°) para ver el disco desde arriba
x_side = 40.0           # Campo de visión amplio para no "perder" el disco
x_pixels = 400

detector = Detector.Detector(D, iota, x_side; x_pixels = x_pixels, ratio = "1:1")

# =========================
# 3. OTRO DISCO: Disco Simple (Intensidad Constante)
# =========================
# Usaremos un disco que emite luz uniforme para verificar geometría
acc_structure = ThinDisk.ThinDisk(bh)
acc_structure.in_edge = 6.0    # ISCO de Schwarzschild es 6.0
acc_structure.out_edge = 30.0  # Un disco muy grande
# Forzamos intensidad constante de 1.0 para pruebas
acc_structure.intensity = (r) -> 1.0 

# =========================
# 4. EJECUCIÓN
# =========================
img = Image(bh, acc_structure, detector)

println("--- Simulando Schwarzschild con Disco Simple ---")
create_photons!(img)
create_image!(img)

# =========================
# 5. DIAGNÓSTICO Y VISUALIZACIÓN
# =========================
max_int = maximum(img.image_data)
println("Intensidad máxima: $max_int")

if max_int > 0
    # Usamos un mapa de colores diferente (Inferno) para variar
    plot(img, savefig=true, filename="Schwarzschild_Test", cmap=:inferno)
else
    println("⚠️ Sigue saliendo negro. Probemos con iota = 0 (vista desde el polo).")
end

println("Presiona ENTER para salir..."); readline()