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

## 4. Estado del script entregable

El entregable de código es **`R/proyecto1.R`**, archivo único y reproducible desde un
clon limpio. Declara **una sola librería**: `library(terra)`. Todo lo demás es base R.

Sigue el Camino B (geoestadística a mano), con correspondencia sección por sección
contra `Ejemplo4_Geoestadística.R`:

| `R/proyecto1.R` | `Ejemplo4` |
|---|---|
| Paso 2 — Puntos de muestreo | `#Borde del valle` + `#Puntos de muestreo` (21-40) |
| Paso 4 — Modelo de tendencia (deriva) | `#Extraer tendencia` (103-106) |
| Paso 5 — Autocorrelación (Moran y Geary) | `Script_Spatial_Analytics.R` (matriz `W` a mano) |
| Paso 6 — Semivariograma y ajuste | `#Semivariagrama` + `#Ajuste de modelo` (113-161) |
| Paso 7 — Kriging universal + LOOCV | `#Predicción` (175-194), `#Validación cruzada` (231-246) |
| Paso 8 — Predicción sobre la grilla | `#Mapa de interpolación` (199-228) |

Misma parametrización (`c0`, `c1`, `phi`, `sill`, `rango_practico <- 3 * phi`), misma
nube de pares con `upper.tri`, misma función de covarianza y mismo sistema kriging
resuelto con `solve()`.

---

## 5. Objeciones frecuentes, y dónde se responden

Preguntas que surgieron en la discusión del grupo, con la línea exacta del material que
las contesta. Sirven para la sustentación.

### «LOOCV no lo vimos en clase»

Sí se vio, en cuatro scripts:

| Script | Línea |
|---|---|
| `Practic_Geostatistics.R` | 108-111 — `krige.cv()`, comentario literal *"Leave-one-out cross-validation"* |
| `Ejemplo1_Geoestadística.R` | 205, 223 |
| `Ejemplo3_Geoestadística.R` | 275, 296 |
| `Ejemplo4_Geoestadística.R` | 231-246 — el bucle manual |

`Ejemplo4` 231-246 es exactamente el procedimiento del entregable: saca el punto *i*,
reajusta el `lm`, krigea los residuos sin *i*, suma la tendencia, y cierra con
`rmse_loo <- sqrt(mean(errores_loo^2))`.

### «Es solo kriging, y el proyecto es multivariado»

El modelo no es kriging ordinario simple: es **kriging universal**. La deriva se estima
con `lm(precip ~ x_km + y_km + altitud)` y el kriging ordinario va sobre sus residuos;
la predicción final vuelve a sumar la tendencia. La covariable entra por la deriva, y el
modelo tiene tres predictores. El molde es `Ejemplo4` 103-106 (`m_trend`, `residuales`)
más 179-209 (`kriging_universal`, `tendencia_k + krg_k["prediccion"]`).

### «¿Por qué temperatura y radiación quedaron fuera?»

Por cobertura, no por preferencia. Medido dentro del área de estudio (688 celdas
válidas), en la semana de trabajo:

| Variable | Celdas con dato | % del área | Valores distintos |
|---|---|---|---|
| Precipitación | 688 | 100 % | 688 |
| Altitud | 688 | 100 % | 688 |
| Temperatura | 123 | 17,9 % | 56 |
| Radiación | 28 | 4,1 % | 3 |
| **Temperatura Y radiación a la vez** | **2** | **0,3 %** | — |

Causa: CHIRPS es nativo 0,05° y NASA POWER es 0,5°, diez veces más grueso. El remuestreo
a la grilla fina reproyectó la rejilla pero no rellenó los huecos, así que dentro del
departamento solo sobrevivieron los centros originales. La misma cobertura aparece en
los archivos anuales (`power_semanal_valle_<año>.tif`, las 104 bandas), así que no es un
defecto de la climatología.

Consecuencias concretas:

1. Un modelo con ambas covariables se ajustaría con **2 observaciones**.
2. La radiación tiene **3 valores distintos** en todo el Valle: es un escalón
   costa/valle/cordillera, o sea posición geográfica, que ya entra por `x_km`. Su
   correlación aparente con la lluvia (r = −0,974 sobre 28 celdas) no es interpretable.
3. La temperatura sí aporta señal real (`cor(temp, altitud) = −0,361`, o sea no es un
   sustituto de la altitud; test F p = 0,030), pero sube el R² de 0,9539 a 0,9557:
   **+0,0018**. A cambio, dejaría sin covariable al 82 % del departamento, y la
   predicción sobre la grilla —que el enunciado exige— quedaría con un hueco.

Por eso ambas se documentan en el EDA y se excluyen de la modelación, con el criterio
explícito de `R/proyecto1.R` 209-227: cobertura ≥ 90 % del área **y** más de 10 valores
distintos.

### «Metamos kNN»

`kNN` **no aparece en ningún script del curso**. La única mención de vecindad es
`Script_Spatial_Analytics.R` 56 — *"Definir vecinos basados en distancia euclidiana"*—
y es para construir la matriz `W` de Moran y Geary, no un predictor. Incorporarlo sería
salirse del enunciado.

### «Y el AFM»

`Script_AFM.R` y `Script_AFM_Ejemplo2.R` sí son material del curso, implementados a mano
con descomposición en valores singulares (`acp_svd`), sin `FactoMineR`. Son legítimos
para análisis exploratorio multivariado, pero no producen predicción espacial sobre una
grilla ni validación geoestadística, que es lo que piden los puntos 3 y 4 del alcance.

---

## 6. Regla práctica

Antes de incorporar una función al entregable, buscarla en los archivos de la sección 1.
Si no aparece, hay tres salidas: reemplazarla por el equivalente que sí está,
implementarla a mano como hace el profesor, o —si es estrictamente de graficación—
justificarla por la excepción del enunciado.
