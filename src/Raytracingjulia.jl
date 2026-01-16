module Raytracingjulia

# --- Common utilities ---
include("common/common.jl")
include("common/integrator.jl")

# --- Accretion structures ---
include("accretion_structures/simple_disk.jl")
include("accretion_structures/thin_disk.jl")

# --- Black holes ---
include("black_holes/schwarzschild.jl")
include("black_holes/kerr.jl")
include("black_holes/scalar_hair_BH.jl")
include("black_holes/num_schwarzschild.jl")

# --- Detectors ---
include("detectors/image_plane.jl")

end
