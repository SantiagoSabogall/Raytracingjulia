# Raytracingjulia

Raytracingjulia integra geodésicas nulas desde un detector hacia un agujero negro para producir imágenes de discos de acreción y sombras.

## Inicio rápido

Desde la raíz del repositorio:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=. ex3.BH_Shadow.jl
```

El ejemplo guarda la sombra en `images/`. El `README.md` de la raíz enumera los demás ejemplos y sus resultados.

## Flujo de cálculo

1. Construye un espacio-tiempo (`SchwarzschildBH`, `KerrBH`, `KerrPFDMBH`, `KerrScalaronBH` o `SchwarzschildNumBH`).
2. Define un `Detector` y, para imágenes de emisión, un `SimpleDisk` o `ThinDisk`.
3. Crea `Image(blackhole, disk, detector)` y llama a `create_photons!`.
4. Usa `create_image!`, `create_image_no_Doppler!` o `create_shadow!`, y representa los datos con `plot`, `plot_contours` o `plot_shadow`.

El disco se detecta por un cruce ecuatorial continuo. Los rayos que llegan al horizonte o escapan sin cruzar el disco tienen intensidad cero. La sombra clasifica captura y escape sin usar la emisión del disco.

## Construir esta documentación

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Abre `docs/build/index.html` después de generar el sitio. La [referencia de API](api.md) contiene los tipos y funciones principales.
