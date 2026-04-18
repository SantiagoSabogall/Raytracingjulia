# Raytracingjulia

**A High-Performance General Relativistic Ray Tracer in Julia**

`Raytracingjulia` is a rigorous, high-performance computational astrophysics library designed for mapping null geodesics in curved spacetimes under the strict framework of Hamiltonian mechanics. Designed with zero-allocation integration loops and type stability in mind, it provides the tools necessary to simulate observations of black hole accretion flows, shadows, and photon rings suitable for academic research and publication.

## Core Principles

1. **Scientific Rigor**: Trace geometries using exact analytic derivations of the Christoffel symbols (or direct Cartesian Hamiltonian bounds) preventing unphysical gauge artifacts. The integration strictly maps 4-momentum states using $\mathcal{H} = \frac{1}{2}g^{\mu\nu}p_\mu p_\nu = 0$.
2. **Computational Excellence**: Leverages `StaticArrays.jl` to ensure the 8-ODE photon geodesic systems allocate entirely on the stack (zero allocations in the integration loop). Threads over observer screen (`Detector`) pixels to achieve massive computational speedups.
3. **Astrophysical Fidelity**: Accretion maps, including the standard Page & Thorne (1974) / Novikov-Thorne radiative flux solutions for geometrically thin disks, coupled with general relativistic Doppler shift evaluations.

## Implemented Spacetimes

- **Schwarzschild**: Exact analytic geometry of a spherically symmetric, non-rotating local topology.
- **Kerr**: Axially symmetric, rotating spacetime derived in standard Boyer-Lindquist coordinates.
- **Kerr-PFDM**: Extensions of the Kerr black hole immersed in Perfect Fluid Dark Matter (PFDM), modifying the effective mass term $m_{eff}(r)$.
- **Numerical Schwarzschild**: Interpolated configurations for numerically constrained lapse combinations using cubic splines—ideal for evaluating beyond-GR / phenomenological metrics.

## Methodology & Conventions

- **Units**: We adopt the standard natural geometric unit convention where gravitational constant $G=1$, speed of light $c=1$, and total ADM mass $M=1$.
- **Metric Signature**: The metric signature conventionally utilized throughout the code is $(-,+,+,+)$.
- **Observer Paradigm**: Imaging conforms to Bardeen's formulation of the ZAMO local tetrad for projecting observational impact parameters $(\alpha, \beta)$ to global 4-momentum null coordinates. 

## Project Structure

```
src/
├── Raytracingjulia.jl          # Core base definitions 
├── accretion_structures/       # Disk emission geometries 
│   ├── simple_disk.jl          # Linear intensity mappings
│   └── thin_disk.jl            # Novikov-Thorne Page-Thorne flux models
├── black_holes/                # Curvature constraints
│   ├── kerr.jl                 # Spin metrics
│   ├── kerrPDFM.jl             # PFDM extensions
│   ├── num_schwarzschild.jl    # Numerical splines models
│   └── schwarzschild.jl        # Standard M=1 spherical metrics
├── detectors/
│   └── image_plane.jl          # Bardeen observer tetrad definitions
└── common/
    └── common.jl               # ODE Integration loops, Doppler shifting, Hamiltonians
```

## Reproducibility and Data Availability 

This code emphasizes transparent research. All implementations natively integrate comprehensive Julia docstrings including mathematical definitions in LaTeX for rapid in-terminal verification (`?KerrBH`, `?geodesics`, etc.). The code handles deterministic photon ray bundles via `DifferentialEquations.jl` employing adaptive TSit5 ODE solvers.

## Acknowledgements

`Raytracingjulia` is a rigorous Julia adaptation and extension inspired by the Python codebase `andromeda`, originally developed by Eduard Larrañaga. We sincerely acknowledge his foundational work on the geometric and numerical formulation of null geodesics mapping.

If you use or reference this codebase, please ensure appropriate attribution to the original author:

```bibtex
@software{larranaga2025andromeda,
  author = {Larrañaga, Eduard},
  title = {andromeda:},
  url = {https://github.com/ashcat2005/andromeda},
  version = {main branch},
  year = {2025}
}
```