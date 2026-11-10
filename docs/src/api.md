# Referencia de API

## Geometría y detector

```@docs
Raytracingjulia.SchwarzschildBH
Raytracingjulia.KerrBH
Raytracingjulia.KerrPFDMBH
Raytracingjulia.KerrScalaronBH
Raytracingjulia.SchwarzschildNumBH
Raytracingjulia.Detector
```

## Métricas y geodésicas

```@docs
Raytracingjulia.metric
Raytracingjulia.inverse_metric
Raytracingjulia.geodesics
Raytracingjulia.Omega
```

## Discos

```@docs
Raytracingjulia.SimpleDisk
Raytracingjulia.ThinDisk
```

## Trazado e imágenes

`Image(blackhole, disk, detector)` contiene los rayos y la matriz `image_data`. Primero ejecuta `create_photons!(image)`; después elige una de las funciones de integración:

| Función | Resultado |
| --- | --- |
| `create_image!(image)` | Intensidad con Doppler. |
| `create_image_no_Doppler!(image)` | Intensidad sin Doppler. |
| `create_shadow!(image)` | Clasificación de captura y escape. |

`plot(image)`, `plot_contours(image)` y `plot_shadow(image)` generan los gráficos correspondientes. Aceptan `save=true` y `filename="nombre"`; se guardan en `images/`.
