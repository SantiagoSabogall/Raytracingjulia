# Raytracingjulia

Trazador de rayos nulos para visualizar sombras de agujeros negros y discos de acreción en Julia. Integra geodésicas mediante el Hamiltoniano nulo con `DifferentialEquations.jl` y representa las imágenes con `Plots.jl`.

## Requisitos e instalación

- Julia 1.10 o 1.11 según `Project.toml`. El desarrollo local también se ha probado con Julia 1.12, que aún no figura en `[compat]`.
- Una copia local de este repositorio.

Desde la raíz del proyecto:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Si tu terminal tiene un alias antiguo de `julia`, compruébalo con `type julia` y usa la ruta del ejecutable instalado o elimina ese alias.

## Primer resultado

Para trazar una sombra de Kerr:

```bash
julia -t auto --project=. ex3.BH_Shadow.jl
```

El archivo se guarda en `images/Kerr_a_0.5_90_shadow_test.png`. La sombra aparece negra y el fondo claro. `-t auto` permite usar los hilos disponibles; los scripts activan este proyecto también cuando se ejecutan sin `--project=.`.

Para una imagen del disco con Doppler:

```bash
julia -t auto --project=. ex1.BH_Image.jl
```

`ex1` usa una resolución de 1920 píxeles de ancho y puede tardar. Reduce `x_pixels` dentro del script para una prueba rápida. Las figuras se escriben en `images/` respecto al directorio desde el que lanzas Julia.

| Script | Resultado |
| --- | --- |
| `ex1.BH_Image.jl` | Imagen de disco Kerr-PFDM con Doppler. |
| `ex2.BH_Image_no_doppler.jl` | Imagen de disco Kerr sin Doppler. |
| `ex3.BH_Shadow.jl` | Sombra binaria de Kerr. |
| `ex4.BH_Image_contours.jl` | Contornos de intensidad con Doppler. |
| `ex5.BH_Image_no_doppler_contours.jl` | Contornos de intensidad sin Doppler. |
| `ex6.read_plot.jl` | Archivo vacío; aún no implementa una operación. |
| `ex7.Hamiltonian_constrain.jl` | Gráfico de conservación del Hamiltoniano. |
| `ex8.Time_of_integration.jl` | Tiempo de trazado de una imagen pequeña Kerr-PFDM. |

## Modelos y convenciones

Se usa la firma métrica `(-,+,+,+)` y unidades geométricas `G=c=1`. Los modelos Kerr, Kerr-PFDM y KerrScalaron fijan `M=1`; `SchwarzschildBH(M)` permite otra masa y la utiliza en la métrica, las geodésicas y la velocidad orbital.

| Espacio-tiempo | Constructor |
| --- | --- |
| Schwarzschild | `SchwarzschildBH()` o `SchwarzschildBH(2.0)` |
| Kerr | `KerrBH(a)` |
| Kerr con materia oscura fluida | `KerrPFDMBH(a, k)` |
| Kerr Scalaron | `KerrScalaronBH(a, Mpsi)` |
| Métrica estática numérica | `SchwarzschildNumBH(path_N, path_dN; EH=2.0, ISCOco=6.0, ISCOcounter=ISCOco)` |

Las tablas numéricas son archivos de texto con dos columnas, radio y valor: `N(r)` en `path_N` y `dN/dr` en `path_dN`. Los radios deben estar ordenados para la interpolación. La velocidad orbital se calcula a partir de `dN/dr`; para una métrica distinta de Schwarzschild con `M=1`, proporciona sus valores físicos de `EH` e ISCO en lugar de usar los predeterminados.

`Detector(D, iota, x_side; x_pixels=..., ratio="16:9")` define la distancia, inclinación en radianes y extensión horizontal del plano de imagen. `SimpleDisk` usa un perfil lineal de intensidad; `ThinDisk` usa un perfil de Novikov-Thorne. Ambos aceptan `R_min`, `R_max` y `corotating` (predeterminado `true`). El sentido de rotación se usa en el corrimiento Doppler.

`create_photons!` prepara los rayos. `create_image!` asigna intensidad solo a los rayos que cruzan el disco y aplica el corrimiento Doppler; `create_image_no_Doppler!` usa el mismo evento de cruce sin ese factor. `create_shadow!` clasifica rayos capturados o escapados. Las pruebas incluyen dualidad métrica, conservación del Hamiltoniano, límites analíticos y casos de impacto/escape.

## Documentación ampliada

La guía y referencia de API se generan con [Documenter.jl](https://documenter.juliadocs.org/stable/). Desde la raíz:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Abre `docs/build/index.html` en un navegador. `docs/Project.toml` mantiene Documenter separado de las dependencias de la biblioteca. Consulta [la guía](docs/src/index.md) y [la API](docs/src/api.md) para ampliar la documentación.

## Origen

Proyecto inspirado en [`andromeda`](https://github.com/ashcat2005/andromeda), de Eduard Larrañaga. Consulta [LICENSE](LICENSE) para la licencia de este repositorio.
