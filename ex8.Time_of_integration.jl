using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Revise
using Raytracingjulia
using LinearAlgebra
using NPZ

# ----------------------------
# Configuración del agujero negro
# ----------------------------
a = 0.5
blackhole = KerrPDFMBH(a,0.5 )
