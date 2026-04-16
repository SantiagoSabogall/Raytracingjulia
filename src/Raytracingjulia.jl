"""
    Raytracingjulia

A high-performance black hole raytracer framework formulated in Julia.

Uses Hamiltonian mechanics to trace null geodesics backwards from a
detector to an accretion structure, evaluating the radiative transfer under 
the geometric optics approximation. 

# Conventions
- Metric signature is `(-, +, +, +)`.
- Units are natural: `G = c = M = 1`.
"""
module Raytracingjulia

using LinearAlgebra
using Plots
using DifferentialEquations

# ==============================================================================
# 1. Base Function Declarations
# Ensures consistent method extensions across all derived submodules.
# ==============================================================================

"""
    metric(b, x::AbstractVector)

Evaluates the covariant metric tensor elements \$g_{\\mu\\nu}\$ at 
the 4-position `x` for a given black hole `b`.
"""
function metric end

"""
    inverse_metric(b, x::AbstractVector)

Evaluates the contravariant metric tensor elements \$g^{\\mu\\nu}\$ at 
the 4-position `x` for a given black hole `b`.
"""
function inverse_metric end

"""
    geodesics(q, b, λ)

Evaluates the right-hand side of the geodesic equations derived from the null Hamiltonian:
\$\\mathcal{H} = \\frac{1}{2} g^{\\mu\\nu} p_\\mu p_\\nu = 0\$

Returns an `SVector` with terms \$dx^\\mu/d\\lambda\$ and \$dp_\\mu/d\\lambda\$.
"""
function geodesics end

"""
    Omega(b, r::Real; corotating::Bool=true)

Calculates the Keplerian angular velocity \$\\Omega_K = \\frac{d\\phi}{dt}\$ 
for a test particle in a circular equatorial orbit at radius `r`.
"""
function Omega end

export metric, inverse_metric, geodesics, Omega

# ==============================================================================
# 2. Submodule Inclusions
# ==============================================================================
include("black_holes/schwarzschild.jl")
include("black_holes/kerr.jl")
include("black_holes/num_schwarzschild.jl") 
include("black_holes/kerrPDFM.jl")
include("black_holes/kerrscalaron.jl")

include("detectors/image_plane.jl")
include("accretion_structures/thin_disk.jl")
include("accretion_structures/simple_disk.jl")
include("common/common.jl")

# ==============================================================================
# 3. User-Facing Exports
# ==============================================================================
using .Schwarzschild: SchwarzschildBH
using .Kerr: KerrBH
using .KerrPFDM: KerrPFDMBH
using .KerrScalaron: KerrScalaronBH
using .SchwarzschildNumerical: SchwarzschildNumBH, dr_inverse_metric
using .ImagePlane: Detector
using .ThinDiskmod: ThinDisk
using .SimpleDiskmod: SimpleDisk
using .Common: Image, create_photons!, create_image!,create_image_no_Doppler!,create_shadow!, plot_shadow, plot_contours, verify_Hamiltonian, doppler_shift, plot

export SchwarzschildBH, KerrBH, KerrPFDMBH, KerrScalaronBH, SchwarzschildNumBH, dr_inverse_metric
export Detector, ThinDisk, SimpleDisk, Image
export create_photons!, create_image!,create_image_no_Doppler!,create_shadow!, plot_shadow, plot_contours, verify_Hamiltonian, doppler_shift, plot

end # module