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



