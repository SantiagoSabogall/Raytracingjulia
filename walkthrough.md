# Notas históricas de rendimiento

Este documento conserva el contexto de una refactorización previa del trazado con Doppler. Describe decisiones de implementación, pero **no es un benchmark actual**: los tiempos dependen de Julia, hardware, métrica, resolución y tolerancias.

## Diseño vigente

`create_image!` prepara un `ODEProblem`, callbacks y un conjunto de integradores `Tsit5` antes de recorrer los píxeles. Para cada rayo, reutiliza un integrador mediante `reinit!` y `solve!`. El callback de cruce ecuatorial marca una intersección válida solo si el radio está dentro del disco. La intensidad y el corrimiento Doppler se calculan después de resolver el rayo; los eventos de horizonte y escape no emiten.

Esta estructura evita crear un integrador nuevo por cada píxel en el trazado con Doppler. No implica que la ejecución completa tenga cero asignaciones: las estructuras, callbacks, gráficos y otras operaciones pueden reservar memoria. El trazado sin Doppler usa actualmente un `solve` por rayo y debe evaluarse por separado.

## Medición reproducible

`ex8.Time_of_integration.jl` mide el tiempo de `create_image!` con un escenario Kerr-PFDM pequeño:

```bash
julia -t auto --project=. ex8.Time_of_integration.jl
```

Ejecuta el script al menos dos veces si quieres separar el costo de compilación inicial. Para comparar optimizaciones, fija la versión de Julia, el número de hilos, los parámetros del detector, la métrica y las tolerancias; registra también el error numérico o diferencias en la imagen.
