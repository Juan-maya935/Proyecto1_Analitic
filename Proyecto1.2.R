# ==========================================================
# PROYECTO ANALÍTICA DE DATOS
# ==========================================================

library(terra)


ruta <- "C:/Descargas3/Analitica de datos/proyecto 1AnaliticaDeDatos/datos_proyecto_1/datos_proyecto_1/imagenes_semanales"


# ==========================================================
# 1. DATASET DE PRECIPITACIÓN - CHIRPS
# ==========================================================

archivos_chirps <- list.files(
  path = ruta,
  pattern = "^chirps_semanal_valle_[0-9]{4}\\.tif$",
  full.names = TRUE
)

archivos_chirps


chirps <- rast(archivos_chirps)

chirps


# ----------------------------------------------------------
# Renombrar capas
# ----------------------------------------------------------

nombres_chirps <- c()

for (anio in 2010:2025) {
  
  nombres_chirps <- c(
    nombres_chirps,
    paste0(
      "prec_",
      anio,
      "_sem_",
      sprintf("%02d", 1:52)
    )
  )
}

names(chirps) <- nombres_chirps


# ----------------------------------------------------------
# Revisión estructural
# ----------------------------------------------------------

nlyr(chirps)

dim(chirps)

res(chirps)

ext(chirps)

crs(chirps)

names(chirps)[1:10]


# ----------------------------------------------------------
# Limpieza
# ----------------------------------------------------------

chirps[chirps < 0] <- NA


# ----------------------------------------------------------
# Métricas generales por capa
# ----------------------------------------------------------

min_chirps <- global(
  chirps,
  "min",
  na.rm = TRUE
)

max_chirps <- global(
  chirps,
  "max",
  na.rm = TRUE
)

media_chirps <- global(
  chirps,
  "mean",
  na.rm = TRUE
)

sd_chirps <- global(
  chirps,
  "sd",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Cantidad de valores válidos y NA por capa
# ----------------------------------------------------------

validos_chirps <- global(
  !is.na(chirps),
  "sum",
  na.rm = TRUE
)

na_chirps <- global(
  is.na(chirps),
  "sum",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Resumen de todas las capas
# ----------------------------------------------------------

summary(min_chirps)

summary(max_chirps)

summary(media_chirps)

summary(sd_chirps)

summary(validos_chirps)

summary(na_chirps)


# ----------------------------------------------------------
# Revisar primera semana
# ----------------------------------------------------------

summary(
  values(chirps[[1]])
)

plot(
  chirps[[1]],
  main = "Precipitación - 2010 Semana 01"
)


# ==========================================================
# 2. DATASET POWER
# TEMPERATURA Y RADIACIÓN
# ==========================================================

archivos_power <- list.files(
  path = ruta,
  pattern = "^power_semanal_valle_[0-9]{4}\\.tif$",
  full.names = TRUE
)

archivos_power


# ==========================================================
# 4. SEPARAR TEMPERATURA Y RADIACIÓN
# ==========================================================

lista_temperatura <- list()

lista_radiacion <- list()


for (i in 1:length(archivos_power)) {
  
  power_anio <- rast(
    archivos_power[i]
  )
  
  
  anio <- sub(
    ".*power_semanal_valle_([0-9]{4})\\.tif",
    "\\1",
    archivos_power[i]
  )
  
  
  # --------------------------------------------------------
  # Capas 1 a 52 = temperatura
  # --------------------------------------------------------
  
  temp_anio <- power_anio[[1:52]]
  
  
  names(temp_anio) <- paste0(
    "temp_",
    anio,
    "_sem_",
    sprintf("%02d", 1:52)
  )
  
  
  lista_temperatura[[i]] <- temp_anio
  
  
  # --------------------------------------------------------
  # Capas 53 a 104 = radiación
  # --------------------------------------------------------
  
  rad_anio <- power_anio[[53:104]]
  
  
  names(rad_anio) <- paste0(
    "rad_",
    anio,
    "_sem_",
    sprintf("%02d", 1:52)
  )
  
  
  lista_radiacion[[i]] <- rad_anio
}


# ==========================================================
# 5. CREAR RASTER COMPLETO DE TEMPERATURA
# ==========================================================

temperatura <- rast(
  lista_temperatura
)

temperatura


# ----------------------------------------------------------
# Revisar estructura
# ----------------------------------------------------------

nlyr(temperatura)

dim(temperatura)

res(temperatura)

ext(temperatura)

crs(temperatura)

names(temperatura)[1:10]


# ----------------------------------------------------------
# Métricas de temperatura
# ----------------------------------------------------------

min_temperatura <- global(
  temperatura,
  "min",
  na.rm = TRUE
)

max_temperatura <- global(
  temperatura,
  "max",
  na.rm = TRUE
)

media_temperatura <- global(
  temperatura,
  "mean",
  na.rm = TRUE
)

sd_temperatura <- global(
  temperatura,
  "sd",
  na.rm = TRUE
)


validos_temperatura <- global(
  !is.na(temperatura),
  "sum",
  na.rm = TRUE
)

na_temperatura <- global(
  is.na(temperatura),
  "sum",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Resumen de las métricas
# ----------------------------------------------------------

summary(min_temperatura)

summary(max_temperatura)

summary(media_temperatura)

summary(sd_temperatura)

summary(validos_temperatura)

summary(na_temperatura)


# ----------------------------------------------------------
# Primera semana
# ----------------------------------------------------------

summary(
  values(temperatura[[1]])
)

plot(
  temperatura[[1]],
  main = "Temperatura - 2010 Semana 01"
)


# ==========================================================
# 6. CREAR RASTER COMPLETO DE RADIACIÓN
# ==========================================================

radiacion <- rast(
  lista_radiacion
)

radiacion


# ----------------------------------------------------------
# Revisar estructura
# ----------------------------------------------------------

nlyr(radiacion)

dim(radiacion)

res(radiacion)

ext(radiacion)

crs(radiacion)

names(radiacion)[1:10]


# ----------------------------------------------------------
# Métricas de radiación
# ----------------------------------------------------------

min_radiacion <- global(
  radiacion,
  "min",
  na.rm = TRUE
)

max_radiacion <- global(
  radiacion,
  "max",
  na.rm = TRUE
)

media_radiacion <- global(
  radiacion,
  "mean",
  na.rm = TRUE
)

sd_radiacion <- global(
  radiacion,
  "sd",
  na.rm = TRUE
)


validos_radiacion <- global(
  !is.na(radiacion),
  "sum",
  na.rm = TRUE
)

na_radiacion <- global(
  is.na(radiacion),
  "sum",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Resumen de las métricas
# ----------------------------------------------------------

summary(min_radiacion)

summary(max_radiacion)

summary(media_radiacion)

summary(sd_radiacion)

summary(validos_radiacion)

summary(na_radiacion)


# ----------------------------------------------------------
# Primera semana
# ----------------------------------------------------------

summary(
  values(radiacion[[1]])
)

plot(
  radiacion[[1]],
  main = "Radiación - 2010 Semana 01"
)

# ==========================================================
# VISUALIZACIÓN DE TEMPERATURA
# ==========================================================

# Primera semana: 2010 - semana 01
plot(
  temperatura[[1]],
  main = "Temperatura - 2010 Semana 01"
)


# ----------------------------------------------------------
# Ejemplo: semana 30 del año 2010
# ----------------------------------------------------------

plot(
  temperatura[[30]],
  main = "Temperatura - 2010 Semana 30"
)


# ----------------------------------------------------------
# Comparar varias semanas del mismo año
# ----------------------------------------------------------

par(mfrow = c(2, 2))

plot(
  temperatura[[1]],
  main = "2010 - Semana 01"
)

plot(
  temperatura[[13]],
  main = "2010 - Semana 13"
)

plot(
  temperatura[[26]],
  main = "2010 - Semana 26"
)

plot(
  temperatura[[39]],
  main = "2010 - Semana 39"
)

par(mfrow = c(1, 1))

# Cantidad de semanas con datos válidos por celda
conteo_temp <- app(
  temperatura,
  fun = function(x) sum(!is.na(x))
)

plot(
  conteo_temp,
  main = "Disponibilidad de datos de temperatura"
)




# ==========================================================
# 4. DATASET DE ALTITUD
# ==========================================================

altitud <- rast(
  file.path(
    ruta,
    "altitud_valle.tif"
  )
)

names(altitud) <- "altitud"

altitud


# ----------------------------------------------------------
# Revisar estructura
# ----------------------------------------------------------

nlyr(altitud)

dim(altitud)

res(altitud)

ext(altitud)

crs(altitud)


# ----------------------------------------------------------
# Métricas de altitud
# ----------------------------------------------------------

min_altitud <- global(
  altitud,
  "min",
  na.rm = TRUE
)

max_altitud <- global(
  altitud,
  "max",
  na.rm = TRUE
)

media_altitud <- global(
  altitud,
  "mean",
  na.rm = TRUE
)

sd_altitud <- global(
  altitud,
  "sd",
  na.rm = TRUE
)

validos_altitud <- global(
  !is.na(altitud),
  "sum",
  na.rm = TRUE
)

na_altitud <- global(
  is.na(altitud),
  "sum",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Mostrar métricas
# ----------------------------------------------------------

min_altitud

max_altitud

media_altitud

sd_altitud

validos_altitud

na_altitud


# ----------------------------------------------------------
# Resumen estadístico
# ----------------------------------------------------------

summary(
  values(altitud)
)


# ----------------------------------------------------------
# Visualización
# ----------------------------------------------------------

plot(
  altitud,
  main = "Altitud - Valle del Cauca"
)


# ==========================================================
# 4. CLIMATOLOGÍA DE TEMPERATURA
# ==========================================================

temp_climatologia <- rast(
  file.path(
    ruta,
    "power_temp_climatologia_semanal_valle.tif"
  )
)

temp_climatologia


# ----------------------------------------------------------
# Revisar estructura
# ----------------------------------------------------------

nlyr(temp_climatologia)

dim(temp_climatologia)

res(temp_climatologia)

ext(temp_climatologia)

crs(temp_climatologia)


# ----------------------------------------------------------
# Métricas
# ----------------------------------------------------------

global(
  temp_climatologia,
  "min",
  na.rm = TRUE
)

global(
  temp_climatologia,
  "max",
  na.rm = TRUE
)

global(
  temp_climatologia,
  "mean",
  na.rm = TRUE
)

global(
  temp_climatologia,
  "sd",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Valores válidos
# ----------------------------------------------------------

global(
  !is.na(temp_climatologia),
  "sum",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Valores NA
# ----------------------------------------------------------

global(
  is.na(temp_climatologia),
  "sum",
  na.rm = TRUE
)


plot(
  temp_climatologia[[1]],
  main = "Climatología de temperatura - Semana 01"
)


# ==========================================================
# 9. CLIMATOLOGÍA DE RADIACIÓN
# ==========================================================

rad_climatologia <- rast(
  file.path(
    ruta,
    "power_radiacion_climatologia_semanal_valle.tif"
  )
)

rad_climatologia


# ----------------------------------------------------------
# Revisar estructura
# ----------------------------------------------------------

nlyr(rad_climatologia)

dim(rad_climatologia)

res(rad_climatologia)

ext(rad_climatologia)

crs(rad_climatologia)


# ----------------------------------------------------------
# Métricas
# ----------------------------------------------------------

global(
  rad_climatologia,
  "min",
  na.rm = TRUE
)

global(
  rad_climatologia,
  "max",
  na.rm = TRUE
)

global(
  rad_climatologia,
  "mean",
  na.rm = TRUE
)

global(
  rad_climatologia,
  "sd",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Valores válidos
# ----------------------------------------------------------

global(
  !is.na(rad_climatologia),
  "sum",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Valores NA
# ----------------------------------------------------------

global(
  is.na(rad_climatologia),
  "sum",
  na.rm = TRUE
)


plot(
  rad_climatologia[[1]],
  main = "Climatología de radiación - Semana 01"
)


# ==========================================================
# 10. CLIMATOLOGÍA DE PRECIPITACIÓN
# ==========================================================

prec_climatologia <- rast(
  file.path(
    ruta,
    "chirps_climatologia_semanal_valle.tif"
  )
)

prec_climatologia


# ----------------------------------------------------------
# Eliminar posibles valores negativos
# ----------------------------------------------------------

prec_climatologia[
  prec_climatologia < 0
] <- NA


# ----------------------------------------------------------
# Revisar estructura
# ----------------------------------------------------------

nlyr(prec_climatologia)

dim(prec_climatologia)

res(prec_climatologia)

ext(prec_climatologia)

crs(prec_climatologia)


# ----------------------------------------------------------
# Métricas
# ----------------------------------------------------------

global(
  prec_climatologia,
  "min",
  na.rm = TRUE
)

global(
  prec_climatologia,
  "max",
  na.rm = TRUE
)

global(
  prec_climatologia,
  "mean",
  na.rm = TRUE
)

global(
  prec_climatologia,
  "sd",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Valores válidos
# ----------------------------------------------------------

global(
  !is.na(prec_climatologia),
  "sum",
  na.rm = TRUE
)


# ----------------------------------------------------------
# Valores NA
# ----------------------------------------------------------

global(
  is.na(prec_climatologia),
  "sum",
  na.rm = TRUE
)


plot(
  prec_climatologia[[1]],
  main = "Climatología precipitación - Semana 01"
)


# ==========================================================
# 11. COMPARACIÓN ESPACIAL INICIAL
# ==========================================================

par(
  mfrow = c(2, 2)
)

plot(
  chirps[[1]],
  main = "Precipitación"
)

plot(
  temperatura[[1]],
  main = "Temperatura"
)

plot(
  radiacion[[1]],
  main = "Radiación"
)

plot(
  altitud,
  main = "Altitud"
)

par(
  mfrow = c(1, 1)
)

# ==========================================================
# BLOQUE 2 - TEMPERATURA
# ANÁLISIS DE DEPENDENCIA ESPACIAL
# ==========================================================

# ==========================================================
# 1. TEMPERATURA MEDIA Y DISPONIBILIDAD DE INFORMACIÓN
# ==========================================================

# Temperatura media de las 832 semanas disponibles
temperatura_media <- mean(
  temperatura,
  na.rm = TRUE
)

plot(
  temperatura_media,
  main = "Temperatura media 2010-2025"
)


# Número de observaciones válidas por celda
n_obs_temp <- app(
  temperatura,
  fun = function(x) sum(!is.na(x))
)

plot(
  n_obs_temp,
  main = "Número de observaciones de temperatura por celda"
)

global(
  n_obs_temp,
  c("min", "max", "mean"),
  na.rm = TRUE
)


# ==========================================================
# 2. CONVERSIÓN DE CELDAS CON TEMPERATURA A PUNTOS
# ==========================================================

temp_puntos <- as.points(
  temperatura_media,
  values = TRUE,
  na.rm = TRUE
)

temp_puntos


# Convertimos a data.frame
coords_temp <- crds(temp_puntos)

datos_temp <- data.frame(
  lon = coords_temp[, 1],
  lat = coords_temp[, 2],
  temperatura = values(temp_puntos)[, 1]
)


# ==========================================================
# 3. INCORPORACIÓN DE LA ALTITUD
# ==========================================================

alt_temp <- extract(
  altitud,
  temp_puntos
)

datos_temp$altitud <- alt_temp[, 2]


# Revisión antes de eliminar faltantes
summary(datos_temp)

colSums(
  is.na(datos_temp)
)


# Conservamos únicamente observaciones completas
datos_temp <- na.omit(datos_temp)

nrow(datos_temp)


# ==========================================================
# 4. RELACIÓN TEMPERATURA - ALTITUD
# ==========================================================

cor_temp_alt <- cor(
  datos_temp$altitud,
  datos_temp$temperatura,
  use = "complete.obs"
)

cor_temp_alt


modelo_altitud <- lm(
  temperatura ~ altitud,
  data = datos_temp
)

summary(modelo_altitud)


plot(
  datos_temp$altitud,
  datos_temp$temperatura,
  pch = 19,
  xlab = "Altitud (m)",
  ylab = "Temperatura (°C)",
  main = "Relación temperatura - altitud"
)

abline(
  modelo_altitud,
  col = "red",
  lwd = 2
)


# ==========================================================
# 5. CONVERSIÓN DE COORDENADAS GEOGRÁFICAS A KM
# ==========================================================

lat0 <- mean(datos_temp$lat)

lon0 <- mean(datos_temp$lon)


datos_temp$x_km <- (
  datos_temp$lon - lon0
) * 111.320 * cos(lat0 * pi / 180)


datos_temp$y_km <- (
  datos_temp$lat - lat0
) * 110.574


head(datos_temp)


# ==========================================================
# 6. MODELO DE TENDENCIA
#    TEMPERATURA ~ ALTITUD + POSICIÓN
# ==========================================================

modelo_tendencia <- lm(
  temperatura ~ altitud + x_km + y_km,
  data = datos_temp
)

summary(modelo_tendencia)


# Los residuos representan la variación que no fue
# explicada por altitud ni por la tendencia espacial.
residuales_temp <- residuals(
  modelo_tendencia
)


# ==========================================================
# 7. MATRIZ DE DISTANCIAS
# ==========================================================

D_temp <- as.matrix(
  dist(
    datos_temp[, c("x_km", "y_km")]
  )
)

round(
  D_temp[1:5, 1:5],
  2
)


# ==========================================================
# 8. NUBE DEL SEMIVARIOGRAMA
# ==========================================================

pares_temp <- which(
  upper.tri(D_temp),
  arr.ind = TRUE
)


h_ij_temp <- D_temp[
  upper.tri(D_temp)
]


gamma_ij_temp <- (
  residuales_temp[pares_temp[, 1]] -
    residuales_temp[pares_temp[, 2]]
)^2 / 2


plot(
  h_ij_temp,
  gamma_ij_temp,
  pch = 20,
  col = "grey50",
  xlab = "Distancia (km)",
  ylab = "Semivarianza",
  main = "Nube de variograma - Temperatura"
)


# ==========================================================
# 9. SEMIVARIOGRAMA EXPERIMENTAL
# ==========================================================

n_bins <- 10

bins_temp <- cut(
  h_ij_temp,
  breaks = n_bins
)


gamma_emp_temp <- tapply(
  gamma_ij_temp,
  bins_temp,
  mean
)


N_h_temp <- tapply(
  gamma_ij_temp,
  bins_temp,
  length
)


h_medio_temp <- tapply(
  h_ij_temp,
  bins_temp,
  mean
)


semivariograma_temp <- na.omit(
  data.frame(
    h = h_medio_temp,
    gamma = gamma_emp_temp,
    N = N_h_temp
  )
)

semivariograma_temp


plot(
  semivariograma_temp$h,
  semivariograma_temp$gamma,
  pch = 19,
  type = "b",
  xlab = "Distancia (km)",
  ylab = "Semivarianza empírica",
  main = "Semivariograma experimental - Temperatura"
)


# ==========================================================
# 10. SEMIVARIOGRAMA PARA DISTANCIAS CORTAS
# ==========================================================

# Se utiliza como distancia máxima la mitad de la
# distancia máxima observada entre puntos.

distancia_maxima <- max(
  h_ij_temp,
  na.rm = TRUE
)

distancia_corte <- distancia_maxima / 2

distancia_maxima
distancia_corte


seleccion <- h_ij_temp <= distancia_corte

h_corto <- h_ij_temp[seleccion]

gamma_corto <- gamma_ij_temp[seleccion]


bins_corto <- cut(
  h_corto,
  breaks = n_bins
)


gamma_emp_corto <- tapply(
  gamma_corto,
  bins_corto,
  mean
)


N_h_corto <- tapply(
  gamma_corto,
  bins_corto,
  length
)


h_medio_corto <- tapply(
  h_corto,
  bins_corto,
  mean
)


semivariograma_corto <- na.omit(
  data.frame(
    h = h_medio_corto,
    gamma = gamma_emp_corto,
    N = N_h_corto
  )
)

semivariograma_corto


plot(
  semivariograma_corto$h,
  semivariograma_corto$gamma,
  pch = 19,
  type = "b",
  xlab = "Distancia (km)",
  ylab = "Semivarianza empírica",
  main = "Semivariograma - Temperatura (distancias cortas)"
)


# ==========================================================
# 11. VALIDACIÓN CRUZADA DEL MODELO DE TENDENCIA
#     Leave-One-Out Cross Validation
# ==========================================================

n <- nrow(datos_temp)

errores_loo <- numeric(n)

predicciones_loo <- numeric(n)


for (i in 1:n) {
  
  # Retiramos una observación
  datos_entrenamiento <- datos_temp[-i, ]
  
  
  # Ajustamos el modelo con las observaciones restantes
  modelo_sin_i <- lm(
    temperatura ~ altitud + x_km + y_km,
    data = datos_entrenamiento
  )
  
  
  # Predecimos la observación retirada
  pred_i <- predict(
    modelo_sin_i,
    newdata = datos_temp[i, , drop = FALSE]
  )
  
  
  predicciones_loo[i] <- pred_i
  
  
  # Error = observado - predicho
  errores_loo[i] <-
    datos_temp$temperatura[i] -
    pred_i
}


validacion <- data.frame(
  observado = datos_temp$temperatura,
  predicho = predicciones_loo,
  error = errores_loo
)


# Métricas de validación
rmse_loo <- sqrt(
  mean(
    errores_loo^2,
    na.rm = TRUE
  )
)


mae_loo <- mean(
  abs(errores_loo),
  na.rm = TRUE
)


error_medio <- mean(
  errores_loo,
  na.rm = TRUE
)


cor_obs_pred <- cor(
  validacion$observado,
  validacion$predicho,
  use = "complete.obs"
)


r2_cv <- 1 -
  sum(
    errores_loo^2,
    na.rm = TRUE
  ) /
  sum(
    (
      validacion$observado -
        mean(validacion$observado, na.rm = TRUE)
    )^2,
    na.rm = TRUE
  )


# Resultados de la validación
metricas_validacion <- data.frame(
  RMSE = rmse_loo,
  MAE = mae_loo,
  Error_medio = error_medio,
  Correlacion = cor_obs_pred,
  R2 = r2_cv
)

metricas_validacion


# Observado vs predicho
plot(
  validacion$observado,
  validacion$predicho,
  pch = 19,
  xlab = "Temperatura observada (°C)",
  ylab = "Temperatura predicha (°C)",
  main = "Validación cruzada - Temperatura"
)

abline(
  0,
  1,
  col = "red",
  lwd = 2
)


# Errores
plot(
  validacion$observado,
  validacion$error,
  pch = 19,
  xlab = "Temperatura observada (°C)",
  ylab = "Error (observado - predicho)",
  main = "Errores de validación cruzada"
)

abline(
  h = 0,
  col = "red",
  lty = 2
)


# ==========================================================
# 12. AJUSTE DEL MODELO EXPONENCIAL AL SEMIVARIOGRAMA
# ==========================================================

# Función de suma de cuadrados del error
sce_exponencial <- function(theta, h, gamma_emp) {
  
  c0 <- theta[1]
  c1 <- theta[2]
  phi <- theta[3]
  
  gamma_teo <- c0 +
    c1 * (1 - exp(-h / phi))
  
  sum(
    (gamma_emp - gamma_teo)^2
  )
}


# Diferentes valores iniciales para phi
phis_iniciales <- quantile(
  semivariograma_corto$h,
  probs = c(
    0.1,
    0.25,
    0.5,
    0.75,
    1
  )
)


ajustes_candidatos <- lapply(
  phis_iniciales,
  function(phi0) {
    
    optim(
      par = c(
        c0 = 0,
        c1 = var(residuales_temp),
        phi = phi0
      ),
      
      fn = sce_exponencial,
      
      h = semivariograma_corto$h,
      
      gamma_emp = semivariograma_corto$gamma,
      
      method = "L-BFGS-B",
      
      lower = c(
        0,
        0,
        1
      ),
      
      upper = c(
        max(semivariograma_corto$gamma),
        5 * max(semivariograma_corto$gamma),
        3 * max(semivariograma_corto$h)
      )
    )
  }
)


# Selección del ajuste con menor error
sce_por_arranque <- sapply(
  ajustes_candidatos,
  function(a) a$value
)

sce_por_arranque


ajuste_temp <-
  ajustes_candidatos[[which.min(sce_por_arranque)]]


# Parámetros estimados
c0_hat <- unname(
  ajuste_temp$par["c0"]
)

c1_hat <- unname(
  ajuste_temp$par["c1"]
)

phi_hat <- unname(
  ajuste_temp$par["phi"]
)


# Sill y rango práctico
sill_hat <- c0_hat + c1_hat

rango_practico <- 3 * phi_hat


parametros_variograma <- data.frame(
  nugget = c0_hat,
  sill = sill_hat,
  phi = phi_hat,
  rango_practico_km = rango_practico
)

parametros_variograma


# Curva del modelo exponencial
curva_h <- seq(
  0,
  max(semivariograma_corto$h),
  length.out = 100
)


curva_gamma <- c0_hat +
  c1_hat * (
    1 - exp(-curva_h / phi_hat)
  )


plot(
  semivariograma_corto$h,
  semivariograma_corto$gamma,
  pch = 19,
  xlab = "Distancia (km)",
  ylab = "Semivarianza",
  main = "Semivariograma vs modelo exponencial"
)

lines(
  curva_h,
  curva_gamma,
  col = "red",
  lwd = 2
)


# ==========================================================
# 13. SELECCIÓN DE LA DISTANCIA DE VECINDAD
# ==========================================================

# Se evalúan radios cercanos debido a que la exploración
# inicial mostró mejor desempeño en distancias cortas.

distancias_prueba <- c(
  8, 10, 12, 14,
  16, 18, 20
)


resultados_distancia <- data.frame(
  distancia = distancias_prueba,
  RMSE = NA,
  MAE = NA,
  n_predicciones = NA
)


for (j in 1:length(distancias_prueba)) {
  
  distancia_max <- distancias_prueba[j]
  
  errores <- rep(
    NA,
    nrow(datos_temp)
  )
  
  
  for (i in 1:nrow(datos_temp)) {
    
    # Retirar el punto que se desea predecir
    entrenamiento <- datos_temp[-i, ]
    
    
    # Modelo de tendencia
    modelo_i <- lm(
      temperatura ~ altitud + x_km + y_km,
      data = entrenamiento
    )
    
    
    # Residuos del modelo
    residuos_i <- residuals(
      modelo_i
    )
    
    
    # Distancias entre el punto retirado
    # y las observaciones restantes
    dist_i <- sqrt(
      (
        entrenamiento$x_km -
          datos_temp$x_km[i]
      )^2 +
        (
          entrenamiento$y_km -
            datos_temp$y_km[i]
        )^2
    )
    
    
    # Vecinos dentro del radio evaluado
    vecinos <- which(
      dist_i <= distancia_max
    )
    
    
    # Se exige un mínimo de tres vecinos
    if (length(vecinos) >= 3) {
      
      tendencia_i <- predict(
        modelo_i,
        newdata = datos_temp[i, , drop = FALSE]
      )
      
      
      # Mayor peso para observaciones cercanas
      pesos <- 1 / dist_i[vecinos]
      
      
      residuo_predicho <- weighted.mean(
        residuos_i[vecinos],
        pesos
      )
      
      
      pred_i <-
        tendencia_i +
        residuo_predicho
      
      
      errores[i] <-
        datos_temp$temperatura[i] -
        pred_i
    }
  }
  
  
  # Métricas para cada radio
  resultados_distancia$RMSE[j] <- sqrt(
    mean(
      errores^2,
      na.rm = TRUE
    )
  )
  
  
  resultados_distancia$MAE[j] <- mean(
    abs(errores),
    na.rm = TRUE
  )
  
  
  resultados_distancia$n_predicciones[j] <-
    sum(!is.na(errores))
}


resultados_distancia


# ==========================================================
# 14. GRÁFICOS PARA SELECCIONAR LA VECINDAD
# ==========================================================

plot(
  resultados_distancia$distancia,
  resultados_distancia$RMSE,
  type = "b",
  pch = 19,
  xlab = "Distancia máxima de vecindad (km)",
  ylab = "RMSE (°C)",
  main = "Distancia de vecindad vs error de predicción"
)


plot(
  resultados_distancia$distancia,
  resultados_distancia$n_predicciones,
  type = "b",
  pch = 19,
  xlab = "Distancia máxima de vecindad (km)",
  ylab = "Número de puntos predichos",
  main = "Distancia de vecindad vs cobertura"
)


# ==========================================================
# 15. DISTANCIA SELECCIONADA
# ==========================================================

# La validación mostró que:
# 12 km presentó el menor RMSE, pero no alcanzó cobertura total.
# 14 km fue la menor distancia con cobertura de 123/123 puntos
# y mantuvo un error bajo.

distancia_optima <- 14

distancia_optima


# ==========================================================
# IMPUTACIÓN DE TEMPERATURA CAPA POR CAPA
# ==========================================================


# ==========================================================
# 1. PARÁMETROS DEFINIDOS EN LA VALIDACIÓN PREVIA
# ==========================================================

distancia_optima <- 14      # km
min_vecinos <- 3


# ==========================================================
# 2. PREPARAR TODAS LAS CELDAS DEL ÁREA
# ==========================================================

# Utilizamos la altitud como referencia espacial.
# Solo se consideran celdas donde existe información
# de altitud.

celdas_area <- as.points(
  altitud,
  values = TRUE,
  na.rm = TRUE
)


# Coordenadas geográficas
coords_area <- crds(celdas_area)


# Dataset espacial general
base_area <- data.frame(
  lon = coords_area[, 1],
  lat = coords_area[, 2],
  altitud = values(celdas_area)[, 1]
)


# ==========================================================
# 3. CONVERTIR COORDENADAS A KILÓMETROS
# ==========================================================

lat0 <- mean(
  base_area$lat,
  na.rm = TRUE
)

lon0 <- mean(
  base_area$lon,
  na.rm = TRUE
)


base_area$x_km <- (
  base_area$lon - lon0
) * 111.320 * cos(lat0 * pi / 180)


base_area$y_km <- (
  base_area$lat - lat0
) * 110.574


head(base_area)


# ==========================================================
# 4. EXTRAER TEMPERATURA EN TODAS LAS CELDAS DEL ÁREA
# ==========================================================

temp_area <- extract(
  temperatura,
  celdas_area
)


# La primera columna corresponde al ID
temp_valores <- temp_area[, -1, drop = FALSE]


dim(temp_valores)


# ==========================================================
# 5. CREAR OBJETOS PARA GUARDAR RESULTADOS
# ==========================================================

# Copia de los valores originales
temp_imputada <- temp_valores


# Máscara:
# 0 = dato original
# 1 = dato imputado
# NA = no fue posible estimarlo

mascara_imputacion <- matrix(
  NA,
  nrow = nrow(temp_valores),
  ncol = ncol(temp_valores)
)


# Donde originalmente existe temperatura
mascara_imputacion[
  !is.na(as.matrix(temp_valores))
] <- 0


# ==========================================================
# 6. IMPUTACIÓN CAPA POR CAPA
# ==========================================================

for (capa in 1:ncol(temp_valores)) {
  
  cat(
    "Procesando capa",
    capa,
    "de",
    ncol(temp_valores),
    "\n"
  )
  
  
  # --------------------------------------------------------
  # Valores de temperatura de la capa actual
  # --------------------------------------------------------
  
  temp_capa <- temp_valores[, capa]
  
  
  # --------------------------------------------------------
  # Identificar datos observados y faltantes
  # --------------------------------------------------------
  
  idx_observados <- which(
    !is.na(temp_capa)
  )
  
  idx_faltantes <- which(
    is.na(temp_capa)
  )
  
  
  # --------------------------------------------------------
  # Verificar que existan suficientes observaciones
  # para ajustar el modelo
  # --------------------------------------------------------
  
  if (
    length(idx_observados) < 4 ||
    length(idx_faltantes) == 0
  ) {
    
    next
  }
  
  
  # ========================================================
  # 6.1 DATOS OBSERVADOS PARA ESTA SEMANA
  # ========================================================
  
  datos_obs <- data.frame(
    temperatura = temp_capa[idx_observados],
    altitud = base_area$altitud[idx_observados],
    x_km = base_area$x_km[idx_observados],
    y_km = base_area$y_km[idx_observados]
  )
  
  
  # Eliminar cualquier registro incompleto
  datos_obs <- na.omit(
    datos_obs
  )
  
  
  if (nrow(datos_obs) < 4) {
    next
  }
  
  
  # ========================================================
  # 6.2 MODELO DE TENDENCIA
  # ========================================================
  
  modelo_capa <- tryCatch(
    
    lm(
      temperatura ~ altitud + x_km + y_km,
      data = datos_obs
    ),
    
    error = function(e) NULL
  )
  
  
  # Si no fue posible ajustar el modelo,
  # pasar a la siguiente capa
  if (is.null(modelo_capa)) {
    next
  }
  
  
  # ========================================================
  # 6.3 RESIDUOS DE LOS DATOS OBSERVADOS
  # ========================================================
  
  residuos_obs <- residuals(
    modelo_capa
  )
  
  
  # Coordenadas de los puntos observados
  x_obs <- datos_obs$x_km
  
  y_obs <- datos_obs$y_km
  
  
  # ========================================================
  # 6.4 RECORRER CELDAS FALTANTES
  # ========================================================
  
  for (idx in idx_faltantes) {
    
    
    # ------------------------------------------------------
    # La celda debe tener altitud válida
    # ------------------------------------------------------
    
    if (
      is.na(base_area$altitud[idx])
    ) {
      next
    }
    
    
    # ------------------------------------------------------
    # Distancia de la celda faltante a todas
    # las observaciones existentes
    # ------------------------------------------------------
    
    distancias <- sqrt(
      (x_obs - base_area$x_km[idx])^2 +
        (y_obs - base_area$y_km[idx])^2
    )
    
    
    # ------------------------------------------------------
    # Seleccionar vecinos dentro de 14 km
    # ------------------------------------------------------
    
    vecinos <- which(
      distancias <= distancia_optima
    )
    
    
    # ------------------------------------------------------
    # Solo predecir si existen al menos 3 vecinos
    # ------------------------------------------------------
    
    if (
      length(vecinos) >= min_vecinos
    ) {
      
      
      # ====================================================
      # A. PREDICCIÓN DE LA TENDENCIA
      # ====================================================
      
      nueva_celda <- data.frame(
        altitud = base_area$altitud[idx],
        x_km = base_area$x_km[idx],
        y_km = base_area$y_km[idx]
      )
      
      
      tendencia_pred <- predict(
        modelo_capa,
        newdata = nueva_celda
      )
      
      
      # ====================================================
      # B. COMPONENTE ESPACIAL
      # ====================================================
      
      # Los vecinos más cercanos reciben mayor peso
      pesos <- 1 / distancias[vecinos]
      
      
      residuo_pred <- weighted.mean(
        residuos_obs[vecinos],
        pesos,
        na.rm = TRUE
      )
      
      
      # ====================================================
      # C. TEMPERATURA FINAL ESTIMADA
      # ====================================================
      
      temp_predicha <-
        tendencia_pred +
        residuo_pred
      
      
      # Guardar predicción
      temp_imputada[idx, capa] <-
        temp_predicha
      
      
      # Marcar que este valor fue imputado
      mascara_imputacion[idx, capa] <- 1
    }
  }
}


# ==========================================================
# 7. CONVERTIR RESULTADOS NUEVAMENTE A RASTER
# ==========================================================

# Crear lista para almacenar las capas reconstruidas

lista_temp_completa <- vector(
  "list",
  ncol(temp_imputada)
)


lista_mascara <- vector(
  "list",
  ncol(temp_imputada)
)


for (capa in 1:ncol(temp_imputada)) {
  
  
  # --------------------------------------------------------
  # TEMPERATURA
  # --------------------------------------------------------
  
  datos_xyz <- data.frame(
    x = base_area$lon,
    y = base_area$lat,
    temperatura = temp_imputada[, capa]
  )
  
  
  raster_capa <- rast(
    datos_xyz,
    type = "xyz",
    crs = crs(temperatura)
  )
  
  
  names(raster_capa) <-
    names(temperatura)[capa]
  
  
  lista_temp_completa[[capa]] <-
    raster_capa
  
  
  # --------------------------------------------------------
  # MÁSCARA DE IMPUTACIÓN
  # --------------------------------------------------------
  
  mascara_xyz <- data.frame(
    x = base_area$lon,
    y = base_area$lat,
    tipo = mascara_imputacion[, capa]
  )
  
  
  mascara_capa <- rast(
    mascara_xyz,
    type = "xyz",
    crs = crs(temperatura)
  )
  
  
  names(mascara_capa) <-
    names(temperatura)[capa]
  
  
  lista_mascara[[capa]] <-
    mascara_capa
}


# ==========================================================
# 8. UNIR TODAS LAS CAPAS
# ==========================================================

temperatura_completa <- rast(
  lista_temp_completa
)


mascara_temp_imputada <- rast(
  lista_mascara
)


# ==========================================================
# 9. REVISAR RESULTADO
# ==========================================================

temperatura_completa

nlyr(
  temperatura_completa
)

names(
  temperatura_completa
)[1:10]


# ==========================================================
# 10. CONTAR CUÁNTOS DATOS FUERON IMPUTADOS
# ==========================================================

n_imputados <- global(
  mascara_temp_imputada == 1,
  "sum",
  na.rm = TRUE
)

n_imputados


# ==========================================================
# 11. COMPARAR NA ANTES Y DESPUÉS
# ==========================================================

na_antes <- global(
  is.na(temperatura),
  "sum",
  na.rm = TRUE
)


na_despues <- global(
  is.na(temperatura_completa),
  "sum",
  na.rm = TRUE
)


head(
  data.frame(
    capa = names(temperatura),
    NA_antes = na_antes[, 1],
    NA_despues = na_despues[, 1]
  )
)


# ==========================================================
# 12. COMPARACIÓN VISUAL DE UNA SEMANA
# ==========================================================

par(
  mfrow = c(1, 3)
)


plot(
  temperatura[[1]],
  main = "Temperatura original"
)


plot(
  temperatura_completa[[1]],
  main = "Temperatura con imputación"
)


plot(
  mascara_temp_imputada[[1]],
  main = "Origen del dato\n0 = original | 1 = imputado"
)


par(
  mfrow = c(1, 1)
)




