# Criterio metodológico — qué cuenta como "visto en clase"

Este documento existe porque el enunciado impone una restricción que no está escrita en
ninguna parte salvo en los propios scripts del profesor:

> Utilizar **únicamente** los scripts, procedimientos y metodologías trabajados durante
> las clases. Solo se permite el uso de librerías adicionales para la elaboración de
> gráficos y para la carga de mapas o archivos de Excel. Para los demás procedimientos,
> deberán utilizarse las herramientas y funciones trabajadas en clase.
>
> — `docs/Proyecto1_AnalíticaDeDatos_2026-2S.pdf`, condiciones de entrega

"Lo visto en clase" no es una opinión: es un conjunto finito de archivos. Abajo está el
inventario, qué contiene cada uno, y qué se sigue de ahí para este proyecto.

---

## 1. El material de referencia

Todo vive en `ANALÍTICA DE DATOS - 1/Material/` del vault del curso. Los que aplican a
este proyecto:

| Archivo | Qué contiene | Librerías |
|---|---|---|
| `Practic_Geostatistics.R` | Kriging ordinario sobre el dataset `meuse` | `sp`, `gstat` |
| `Practice Geostatistics/Ejemplo1_Geoestadística.R` | Geoestadística a mano, sin paquetes | base R |
| `Practice Geostatistics/Ejemplo2_Geoestadística.R` | Ídem | base R |
| `Practice Geostatistics/Ejemplo3_Geoestadística.R` | Descarga de datos CHIRPS | `geodata`, `chirps`, `terra` |
| `Practice Geostatistics/Ejemplo4_Geoestadística.R` | **Temperatura máxima del Valle del Cauca desde rasters** | `terra` |
| `carpeta/Script_Spatial_Analytics.R` | Moran's I y Geary's C calculados a mano | `ggplot2`, `geostan`, `gridExtra`, `openxlsx`, `dplyr` |
| `carpeta/Script_Class5.R` | Manejo de vector/raster, CRS, autocorrelación | `terra` |

Documentos de apoyo: `Lecture_SpatialStatistics.pdf`, `Class_5_Handle_Spatial_Data.pdf`,
`Patrones_Puntuales.pdf`, `Presentacion_Modelling_ppp.pdf`,
`Geoestadística_Datos_Espaciales_Espacio-Temporales_y_Funcionales.pdf`.

---

## 2. Hay dos caminos, y no son equivalentes

### Camino A — con `gstat`, sobre `meuse`

`Practic_Geostatistics.R` resuelve kriging ordinario apoyándose en el paquete:

```r
ve <- variogram(logZn ~ 1, meuse, cutoff = 1300, width = 90)   # línea 55
vt <- vgm(psill = 0.12, model = "Sph", range = 850, nugget = 0.01)  # línea 67
va <- fit.variogram(ve, vt)                                     # línea 72
ok <- krige(logZn ~ 1, locations = meuse, newdata = meuse.grid, model = va)  # línea 86
cv <- krige.cv(logZn ~ 1, locations = meuse, model = va)        # línea 111
rmse <- sqrt(mean(cv$residual^2))                               # línea 115
```

Es geoestadística clásica, pero sobre un dataset de juguete (`meuse`, contaminación por
Zn en suelo) y con objetos `sp`, no rasters.

### Camino B — a mano, sobre el Valle del Cauca

`Ejemplo4_Geoestadística.R` no carga `gstat` en ningún momento. Su única librería es
`terra` (línea 4). Todo el aparato geoestadístico está implementado explícitamente:

```r
D <- as.matrix(dist(datos[, c("x_km", "y_km")]))                # línea 110
pares <- which(upper.tri(D), arr.ind = TRUE)                    # línea 114
gamma_ij <- (residuales[pares[,1]] - residuales[pares[,2]])^2 / 2  # línea 116
gamma_emp <- tapply(gamma_ij, bins, mean)                       # línea 124

sce_exponencial <- function(theta, h, gamma_emp) { ... }        # línea 138
ajuste <- optim(par = c(c0, c1, phi), fn = sce_exponencial, method = "L-BFGS-B")  # 146

C_hat <- function(h) sill_hat * exp(-h / phi_hat)               # línea 175
kriging_universal <- function(s0, coords_xy, resid, Sigma, C_hat) {
  A <- rbind(cbind(Sigma, rep(1, n_loc)), c(rep(1, n_loc), 0))  # línea 184
  sol <- solve(A, b)                                            # línea 186
}
```

### Por qué B es el molde de este proyecto

No es preferencia estética. El `Ejemplo4` es, literalmente, este mismo problema:

- **Mismo territorio.** Su encabezado dice *"Temperatura máxima en el valle del cauca"*.
- **Mismo tipo de insumo.** Parte de un directorio de `.tif` apilados con `rast()`
  (líneas 10-15), exactamente como nuestros `chirps_semanal_valle_<año>.tif`.
- **Mismo problema de borde.** Recupera el límite departamental desde la máscara de NA
  con `as.polygons()` (líneas 21-23) — nosotros tenemos el GADM, pero el paso es el mismo.
- **Mismo flujo completo**, y en el mismo orden que exige el enunciado (EDA → Modelación
  → Validación → Predicción):

  ```
  #Cargar imagenes → #Borde del valle del cauca → #Puntos de muestreo → #Proyección
  #EDA espacial → #Extraer tendencia → #Distancias → #Semivariagrama nube y muestral
  #Ajuste de modelo → #Predicción → #Mapa de interpolación → #Validación cruzada
  ```

El profesor ya resolvió el ejercicio con estos datos y por este camino. Seguirlo no es
solo lo más seguro frente a la restricción: es lo que hace la sustentación defendible,
porque cada paso tiene un precedente citable en su propio material.

---

## 3. Whitelist operativa

**Librerías con precedente en clase:** `terra`, `sp`, `gstat`, `ggplot2`, `gridExtra`,
`openxlsx`, `dplyr`, `geostan`, `geodata`, `chirps`.

**Funciones con precedente**, por bloque:

| Bloque | Funciones | Fuente |
|---|---|---|
| Carga y manejo raster | `rast`, `vect`, `values`, `crs`, `res`, `ncell`, `mask`, `as.polygons`, `spatSample`, `extract`, `crds`, `writeRaster` | `Script_Class5.R`, `Ejemplo4` |
| EDA | `summary`, `hist`, `plot`, `boxplot`, `par(mfrow=)`, `sd`, `quantile` | `Practic_Geostatistics.R` 22-28, `Ejemplo4` 77-100 |
| Autocorrelación | matriz `W` por umbral de distancia, Moran e I de Geary **por fórmula explícita**; `acf` para la temporal | `Script_Spatial_Analytics.R` 53-128, `Script_Class5.R` |
| Tendencia | `lm`, `predict`, `residuals`, `coef` | `Ejemplo4` 104-106 |
| Variograma | `dist`, `upper.tri`, `cut`, `tapply`, `optim` | `Ejemplo4` 110-153 |
| Kriging | `solve`, `rbind`/`cbind` para el sistema, función de covarianza propia | `Ejemplo4` 175-194 |
| Validación | LOOCV manual con bucle, o `krige.cv`; `rmse <- sqrt(mean(e^2))` | `Ejemplo4` 231-245, `Practic_Geostatistics.R` 111-115 |

**Sin precedente — no usar:** `automap` (y en particular `autofitVariogram`), `spdep`
(incluido `moran.test`), `sf`, `caret`, `randomForest`, `mgcv`/`gam`, `tmap`, `leaflet`,
y el paquete `raster` (en clase se usa `terra`, que es su reemplazo).

`lm()`, `optim()`, `solve()`, `dist()` y demás son base R: siempre disponibles.

---

## 4. Estado del script actual

`Proyecto_1_AnaliticaDeDatos.R` sigue el Camino B de forma fiel. La correspondencia con
`Ejemplo4` es sección por sección:

| Script del proyecto | `Ejemplo4` |
|---|---|
| §5 Modelación de tendencia y residuales | `#Extraer tendencia` (104-106) |
| §6 Semivariograma experimental y ajuste | `#Semivariagrama` + `#Ajuste de modelo` (113-161) |
| §7 Formulación del estimador kriging universal | `#Predicción` (175-194) |
| §8 Validación cruzada LOOCV | `#Validación cruzada` (231-245) |
| §9 Predicción sobre la grilla departamental | `#Mapa de interpolación` (199-228) |

Usa la misma parametrización (`c0_hat`, `c1_hat`, `phi_hat`, `sill_hat`,
`rango_practico <- 3 * phi_hat`), la misma nube de pares con `upper.tri`, el mismo modelo
exponencial y la misma función de covarianza. Metodológicamente está alineado.

### Observación pendiente

La línea 21 declara seis paquetes:

```r
paquetes_requeridos <- c("terra", "sp", "gstat", "ggplot2", "gridExtra", "corrplot")
```

De esos, **solo `terra` se usa**. No hay una sola llamada a `variogram()`, `vgm()`,
`krige()`, `coordinates()`, `spplot()`, `ggplot()`, `grid.arrange()` ni `corrplot()` en
todo el archivo. Las otras cinco se cargan y nunca se invocan.

Esto importa por dos razones:

1. **`corrplot` no tiene precedente en clase.** Es el único paquete del script que no
   aparece en ningún material del curso. Se ampararía en la excepción de "elaboración de
   gráficos", pero es una discusión que no hace falta dar si no se está usando.
2. **Declarar `gstat` sugiere un camino que el script no toma.** Un lector que revise la
   cabecera esperará kriging con `gstat` y encontrará kriging a mano. Vale la pena que la
   lista de paquetes refleje lo que el código hace.

Sugerencia: dejar `paquetes_requeridos <- c("terra")`, y volver a añadir lo que se
necesite cuando efectivamente se use.

---

## 5. Regla práctica

Antes de incorporar una función al entregable, buscarla en los archivos de la sección 1.
Si no aparece, hay tres salidas: reemplazarla por el equivalente que sí está,
implementarla a mano como hace el profesor, o —si es estrictamente de graficación—
justificarla por la excepción del enunciado.
