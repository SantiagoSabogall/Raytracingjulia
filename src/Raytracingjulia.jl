module Raytracingjulia


include("black_holes/schwarzschild.jl")
include("black_holes/kerr.jl")
include("black_holes/num_schwarzschild.jl")
#include("black_holes/scalar_hair_BH.jl")

using .Schwarzschild
using .Kerr
using .SchwarzschildNumerical

# 1. PRIMERO: Common (Estructuras de datos base, tipos de fotones, etc.)
# Es vital porque los integradores y los agujeros negros suelen usar tipos definidos aquí.
include("detectors/image_plane.jl")
using .Detector



# 2. SEGUNDO: Integradores
# Los integradores necesitan a Common, pero no necesariamente a los agujeros negros todavía.

include("common/common.jl")
using .Common

include("common/integrator.jl")
include("common/integrator2.jl")
using .Integrado
using .Integrators

# 3. TERCERO: Black Holes
# Las métricas usan los tipos de Common y a veces funciones del Integrador.


# 4. CUARTO: Estructuras de acreción y detectores
include("accretion_structures/simple_disk.jl")
include("accretion_structures/thin_disk.jl")
using .SimpleDisk
using .ThinDisk



# ============================================================
# Detector
# ============================================================

include("detectors/image_plane.jl")
using .Detector

# ============================================================


# ============================================================
# Public API
# ============================================================

export
    # --- Black holes ---
    Schwarzschild,
    KerrBH,
    SchwarzschildNumBH,

    # --- Accretion ---
    SimpleDisk,
    ThinDisk,

    # --- Detector ---
    Detector,
    photon_coords,

    # --- Image / Ray tracing ---
    Image,
    create_photons!,
    create_image!,
    create_image_no_Doppler!,
    create_shadow!,
    plot,
    plot_shadow,
    plot_contours,
    save_data,
    verify_Hamiltonian

end # module
