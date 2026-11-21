# K-Means: de las estrellas a los municipios

Presentación didáctica e interactiva sobre K-Means (Estadística Multivariada) y su simulación con **municipios sintéticos** de Colombia.

> **Los datos son inventados.** No son datos reales de Colombia ni deben leerse como hallazgos empíricos. El contorno del mapa es una silueta esquemática, no cartografía oficial.

## Contenido
| Ruta | Qué es |
|---|---|
| `docs/index.html` | La presentación (un solo archivo, funciona sin internet). Con GitHub Pages queda publicada. |
| `R/simulacion_kmeans.R` | Las 7 etapas de la simulación en R: exploración, estandarización, K-Means, elección de K (coeficiente de fusión, codo, silueta), perfiles y K-medias difuso. |
| `data/municipios_sinteticos.csv` | 1.100 municipios sintéticos: `id`, 7 variables, `lon`, `lat`. |
| `data/poly_colombia.csv` | Polígono del contorno usado en los mapas. |
| `python/generar_datos.py` | Cómo se generaron los datos (semilla 2026). |
| `figuras/` | Salidas de ejemplo del script de R. |

## Variables
`pct_va_agro`, `pct_area_agro`, `va_agro_pc`, `pob_rural`, `rendimiento`, `conectividad`, `nbi` (nombres ilustrativos; valores simulados).

## Cómo reproducir
En R (probado con R 4.3.3), desde la raíz del repositorio:
```r
source("R/simulacion_kmeans.R")
```
Solo usa el paquete `cluster` (incluido con R). Las figuras se guardan en `figuras/`.

Para regenerar los datos (Python con numpy y scipy), desde `python/`:
```
python generar_datos.py
```

## Diferencias esperadas con la presentación
- El orden/etiqueta de los grupos y pequeños detalles pueden variar: los centros iniciales se sortean distinto.
- El coeficiente de fusión en R usa Ward sobre los 1.100 municipios; la presentación usa una submuestra de 440 para correr rápido en el navegador.

## Fuentes y créditos
Basada en la clase de K-Means y en el capítulo 11 del libro *Análisis Estadístico de Datos Multivariados*. Las citas de la diapositiva 4 están marcadas como recurso pedagógico con verificación pendiente.