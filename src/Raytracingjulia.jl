module Raytracingjulia

# ============================================================
# Black holes
# ============================================================

include("black_holes/schwarzschild.jl")
include("black_holes/kerr.jl")
include("black_holes/num_schwarzschild.jl")
#include("black_holes/scalar_hair_BH.jl")

using .Schwarzschild
using .Kerr
using .SchwarzschildNumerical
#using .ScalarHairBH

# ============================================================
# Accretion structures
# ============================================================

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
# Integrators
# ============================================================

include("common/integrator.jl")
include("common/integrator2.jl")
using .Integrado

# ============================================================
# Common (Image, Photon, plotting, etc.)
# ============================================================

include("common/common.jl")
using .Common

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
