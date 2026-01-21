module Raytracingjulia

using LinearAlgebra
using Plots
using DifferentialEquations

# --- 1. DECLARACIÓN DE FUNCIONES BASE ---
# Esto permite que todos los módulos extiendan la MISMA función
function metric end
function inverse_metric end
function geodesics end
function Omega end

export metric, inverse_metric, geodesics, Omega

# --- 2. INCLUSIÓN DE SUBMÓDULOS ---
include("black_holes/schwarzschild.jl")
include("black_holes/kerr.jl")
include("black_holes/num_schwarzschild.jl") # El que acabamos de agregar
include("detectors/image_plane.jl")
include("accretion_structures/thin_disk.jl")
include("common/common.jl")

# --- 3. EXPORTACIÓN DE TIPOS Y MÉTODOS PARA EL USUARIO ---
using .Schwarzschild: SchwarzschildBH
using .Kerr: KerrBH
using .SchwarzschildNumerical: SchwarzschildNumBH
using .Detector: Detector
using .ThinDisk: ThinDisk
using .Common: Image, create_photons!, create_image!, plot

export SchwarzschildBH, KerrBH, SchwarzschildNumBH
export Detector, ThinDisk, Image
export create_photons!, create_image!, plot

end # module