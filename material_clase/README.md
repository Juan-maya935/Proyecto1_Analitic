# Material de clase

Scripts originales del profesor Johann A. Ospina, traídos del respaldo del curso.
Sirven para una sola cosa: **verificar que solo usamos métodos vistos en clase.**

El enunciado exige "únicamente los scripts, procedimientos y metodologías trabajados
en clase", y librerías extra solo para gráficos, mapas o Excel. Esa lista no está
escrita en ninguna parte — hay que sacarla de estos archivos.

## Qué hay aquí

| Archivo | Qué aporta |
|---|---|
| `Practice Geostatistics/Ejemplo4_Geoestadística.R` | **El molde del proyecto.** Mismo departamento, mismos rasters, mismo flujo |
| `Practice Geostatistics/Ejemplo1–3_Geoestadística.R` | Kriging a mano: `optim()`, `solve()`, `dist()` |
| `Practic_Geostatistics.R` | Camino alterno con `gstat` sobre el dataset `meuse` |
| `carpeta/Script_Spatial_Analytics.R` | Moran's I y Geary's C calculados a mano, matriz `W` por umbral |
| `carpeta/Script_Class5.R` | Manejo de rasters con `terra` |
| `Script_AFM*.R` | Análisis factorial múltiple (otra unidad) |

## Librerías con precedente de clase

✅ `terra`, `sp`, `gstat`, `geodata`, `chirps`, `ggplot2`, `gridExtra`, `geostan`,
`openxlsx`, `dplyr` — y todo lo de base R (`lm`, `optim`, `solve`, `dist`, `tapply`).

❌ Sin precedente: `automap`, `spdep`, `sf`, `caret`, `randomForest`, `mgcv`, `tmap`,
`leaflet`, el paquete viejo `raster` (en clase es `terra`), y `terra::interpIDW`.

## Lo que NO está en el repo

Los libros y papers del curso (109 MB) se dejaron fuera **a propósito**: este repo es
público y son obras con ISBN. Están en el respaldo del curso, no aquí.

## Verificación de los datos del proyecto

Comprobado que `data/datos_proyecto_1.zip` es el original del profesor, sin tocar:

```
md5  d1bd258c058004a7877617817027b60a   zip del profesor == zip del repo
md5  813557728ea4ceee9b5977e37bc68789   enunciado PDF idéntico
CRC  229 de 229 archivos .tif extraídos coinciden con el zip
```

Para regenerar los `.tif` (no se versionan, son derivados):

```bash
unzip data/datos_proyecto_1.zip -d data/
```
