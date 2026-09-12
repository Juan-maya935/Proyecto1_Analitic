################################################################################
# UNIVERSIDAD AUTONOMA DE OCCIDENTE
# Facultad de Ingenieria - Ingenieria de Datos e Inteligencia Artificial
# Analitica de Datos (2026-2S) - Prof. Johann A. Ospina
#
# PROYECTO 1
# Modelamiento de la precipitacion semanal en el Valle del Cauca
# y evaluacion del papel de la correlacion espacial.
#
# Autor: Cesar Armando Reyes Oliveros - 2236379
#
# Script unico y reproducible. Solo `terra` (libreria de clase, Script_Class5.R)
# y funciones base de R.
################################################################################

rm(list = ls())
set.seed(2026)
options(scipen = 999)

library(terra)

# Rutas relativas a la raiz del proyecto -> el script corre con
#   Rscript R/proyecto1.R   desde /Proyecto1_Analitic
DIR_DATOS <- "data/datos_proyecto_1/imagenes_semanales"
DIR_SALIDA <- "resultados"
if (!dir.exists(DIR_SALIDA)) dir.create(DIR_SALIDA, recursive = TRUE)

stopifnot(dir.exists(DIR_DATOS))


# ---- Helper de figuras ------------------------------------------------------
# Problema: si solo se llama plot(), la figura vive en el panel de RStudio y se
# pierde. Si solo se llama png()/dev.off(), va al archivo pero no se ve.
# Esta funcion hace las dos cosas: evalua el mismo bloque de dibujo dos veces,
# una contra la pantalla y otra contra el PNG.
#
# Uso:
#   figura("fig02_borde.png", {
#     plot(precip)
#     plot(borde, add = TRUE)
#   })

figura <- function(nombre, expr, ancho = 2000, alto = 1800, res = 300) {
  bloque <- substitute(expr)          # captura el codigo SIN ejecutarlo
  entorno <- parent.frame()           # ...para correrlo donde fue escrito

  eval(bloque, entorno)               # 1) a la pantalla

  ruta <- file.path(DIR_SALIDA, nombre)
  png(ruta, width = ancho, height = alto, res = res)
  eval(bloque, entorno)               # 2) al archivo
  dev.off()

  cat(sprintf("  [figura] %s\n", ruta))
  invisible(ruta)
}


################################################################################
# PASO 1. AUDITORIA DE LOS DATOS
#
# Antes de modelar hay que saber que hay. No se asume nada del enunciado:
# se verifica la geometria, la cobertura real y los valores centinela.
################################################################################

cat("\n================ PASO 1: AUDITORIA DE LOS DATOS ================\n\n")

# ---- 1.1 Carga de las capas -------------------------------------------------
# Climatologia = promedio historico 2010-2025 por semana ISO (52 bandas).
# Se usa la climatologia y no un anio suelto para que la senal espacial no
# quede dominada por el ruido de un evento particular.

r_precip <- rast(file.path(DIR_DATOS, "chirps_climatologia_semanal_valle.tif"))
r_temp   <- rast(file.path(DIR_DATOS, "power_temp_climatologia_semanal_valle.tif"))
r_rad    <- rast(file.path(DIR_DATOS, "power_radiacion_climatologia_semanal_valle.tif"))
r_alt    <- rast(file.path(DIR_DATOS, "altitud_valle.tif"))


# ---- 1.2 Verificacion de la geometria comun ---------------------------------
# El enunciado afirma que todas las capas comparten grilla. Se comprueba:
# si no fuera cierto, no se podrian cruzar pixel a pixel.

geometria <- function(r) {
  e <- as.vector(ext(r))
  c(filas = nrow(r), cols = ncol(r), bandas = nlyr(r),
    res = unname(res(r)[1]),
    xmin = unname(e[1]), xmax = unname(e[2]),
    ymin = unname(e[3]), ymax = unname(e[4]))
}

tabla_geom <- rbind(
  Precipitacion = geometria(r_precip),
  Temperatura   = geometria(r_temp),
  Radiacion     = geometria(r_rad),
  Altitud       = geometria(r_alt)
)

cat("--- 1.2 Geometria de las capas ---\n")
print(tabla_geom)

# Prueba formal: todas las filas de geometria (menos el numero de bandas)
# deben ser identicas a la primera.
cols_comparables <- setdiff(colnames(tabla_geom), "bandas")
grilla_comun <- all(apply(tabla_geom[, cols_comparables], 2,
                          function(x) all(abs(x - x[1]) < 1e-9)))
cat(sprintf("\nGrilla comun a las cuatro capas: %s\n", ifelse(grilla_comun, "SI", "NO")))
cat(sprintf("Pixeles por banda: %d x %d = %d\n",
            nrow(r_precip), ncol(r_precip), ncell(r_precip)))
cat(sprintf("Tamano de pixel: %.3f grados (~%.1f km)\n",
            res(r_precip)[1], res(r_precip)[1] * 111.32))


# ---- 1.3 Valores centinela --------------------------------------------------
# CHIRPS codifica "sin dato" como -9999 / -69993 en vez de NA. Si no se
# limpian, cualquier media queda destruida.

v_precip_crudo <- values(r_precip)
v_precip_crudo <- v_precip_crudo[!is.na(v_precip_crudo)]

cat("\n--- 1.3 Valores centinela en CHIRPS ---\n")
cat(sprintf("Minimo crudo observado en todas las bandas: %.1f\n", min(v_precip_crudo)))
cat(sprintf("Celdas con valor negativo (imposible en lluvia): %d\n",
            sum(v_precip_crudo < 0)))

# Regla de limpieza: la precipitacion es una cantidad no negativa.
# Todo valor < 0 es centinela, no dato.
r_precip[r_precip < 0] <- NA

cat(sprintf("Minimo despues de la limpieza: %.2f mm\n",
            min(values(r_precip), na.rm = TRUE)))


# Semana usada solo para los diagnosticos de cobertura de este paso.
# La semana definitiva de trabajo se elige con evidencia en 1.6.
SEMANA_DIAG <- 42

# ---- 1.4 Mascara del area de estudio ----------------------------------------
# El rectangulo de 1443 pixeles no es el departamento: incluye oceano Pacifico
# y territorio vecino. El area de estudio son las celdas con altitud valida
# (SRTM solo trae tierra firme) Y con precipitacion valida.
# Este es el denominador honesto para medir cobertura.

mascara <- !is.na(r_alt) & !is.na(r_precip[[SEMANA_DIAG]])
N_VALLE <- sum(values(mascara), na.rm = TRUE)

cat("\n--- 1.4 Mascara del area de estudio ---\n")
cat(sprintf("Celdas del rectangulo:        %d\n", ncell(r_alt)))
cat(sprintf("Celdas del area de estudio:   %d (%.1f%% del rectangulo)\n",
            N_VALLE, 100 * N_VALLE / ncell(r_alt)))


# ---- 1.5 Cobertura efectiva de cada variable --------------------------------
# "Misma grilla" no implica "misma cantidad de dato". NASA POWER es nativo
# 0.5 grados y CHIRPS 0.05 grados: al llevar POWER a la grilla fina, solo los
# centros originales quedaron con valor. Se mide DENTRO del area de estudio.

cobertura <- function(r, banda = 1) {
  v <- values(mask(r[[banda]], mascara, maskvalues = c(FALSE, NA)))
  v <- v[!is.na(v)]
  c(celdas_con_dato = length(v),
    pct_area = 100 * length(v) / N_VALLE,
    valores_unicos = length(unique(round(v, 6))))
}

tabla_cob <- rbind(
  Precipitacion = cobertura(r_precip, SEMANA_DIAG),
  Altitud       = cobertura(r_alt,    1),
  Temperatura   = cobertura(r_temp,   SEMANA_DIAG),
  Radiacion     = cobertura(r_rad,    SEMANA_DIAG)
)

cat(sprintf("\n--- 1.5 Cobertura efectiva dentro del area (semana %d) ---\n", SEMANA_DIAG))
print(round(tabla_cob, 2))

write.csv(data.frame(Variable = rownames(tabla_cob), round(tabla_cob, 2)),
          file.path(DIR_SALIDA, "tabla_cobertura_variables.csv"), row.names = FALSE)


# ---- 1.5b Criterio de seleccion de covariables ------------------------------
# Una covariable entra al modelo solo si describe el area completa y tiene
# variabilidad real. Dos filtros, ambos necesarios:
#   (a) cobertura >= 90% del area de estudio  -> no hay que inventar dato
#   (b) mas de 10 valores distintos           -> es un gradiente, no un escalon

UMBRAL_PCT <- 90
UMBRAL_UNICOS <- 10

decision <- ifelse(tabla_cob[, "pct_area"] >= UMBRAL_PCT &
                   tabla_cob[, "valores_unicos"] > UMBRAL_UNICOS,
                   "SE USA", "SE DESCARTA")

cat("\n--- 1.5b Decision sobre covariables ---\n")
print(data.frame(Variable = rownames(tabla_cob),
                 Pct_area = round(tabla_cob[, "pct_area"], 1),
                 Valores_unicos = tabla_cob[, "valores_unicos"],
                 Decision = decision, row.names = NULL))


# ---- 1.6 Eleccion de la semana de estudio -----------------------------------
# No se asume cual es la semana lluviosa: se calcula el promedio espacial de
# cada una de las 52 semanas y se toma el maximo.

media_semanal <- sapply(1:nlyr(r_precip),
                        function(k) mean(values(r_precip[[k]]), na.rm = TRUE))

SEMANA <- which.max(media_semanal)

cat("\n--- 1.6 Ciclo anual de precipitacion ---\n")
cat(sprintf("Semana mas lluviosa (climatologia 2010-2025): %d (%.1f mm)\n",
            SEMANA, media_semanal[SEMANA]))
cat(sprintf("Semana mas seca: %d (%.1f mm)\n",
            which.min(media_semanal), min(media_semanal)))

figura("fig01_ciclo_anual.png", {
  par(mar = c(4.5, 4.5, 3, 1))
  plot(1:52, media_semanal, type = "o", pch = 19, cex = 0.6, col = "#08519c",
       xlab = "Semana ISO", ylab = "Precipitacion media del departamento (mm)",
       main = "Ciclo anual de precipitacion - Valle del Cauca (CHIRPS 2010-2025)",
       las = 1)
  abline(v = SEMANA, col = "red", lty = 2, lwd = 2)
  text(SEMANA, media_semanal[SEMANA], labels = paste0("  Semana ", SEMANA),
       pos = 4, col = "red", cex = 0.8)
  grid(col = "grey85")
}, ancho = 2000, alto = 1100)

cat(sprintf("\n[OK] Paso 1 terminado. Semana de trabajo fijada en %d.\n", SEMANA))


################################################################################
# PASO 2. PUNTOS DE MUESTREO
#
# El molde del curso (Ejemplo4_Geoestadistica.R, lineas 21-35) construye el
# borde real del area a partir de la mascara y sortea puntos dentro de el.
################################################################################

cat("\n================ PASO 2: PUNTOS DE MUESTREO ================\n\n")

# ---- 2.1 Borde real del departamento ----------------------------------------
# as.polygons() convierte los pixeles TRUE de la mascara en un poligono.
# Devuelve DOS poligonos: el de los TRUE (atributo 1) y el de los FALSE (0).
# Hay que quedarse con el 1; si no, el "borde" incluye el oceano Pacifico.

# CHIRPS tiene 787 celdas con dato y SRTM 688: CHIRPS desborda el departamento
# y cubre mar abierto. Sin recortar, el mapa pinta lluvia sobre el Pacifico.
# maskvalues = c(FALSE, NA) deja fuera tanto los FALSE como los NA de la mascara.
precip <- mask(r_precip[[SEMANA]], mascara, maskvalues = c(FALSE, NA))

cat(sprintf("Celdas de CHIRPS antes del recorte: %d\n",
            sum(!is.na(values(r_precip[[SEMANA]])))))
cat(sprintf("Celdas tras recortar al area de estudio: %d\n",
            sum(!is.na(values(precip)))))

borde <- as.polygons(mascara, dissolve = TRUE)
borde <- borde[borde[[1]] == 1, ]

cat(sprintf("Poligonos devueltos por as.polygons(): %d\n", nrow(as.polygons(mascara, dissolve = TRUE))))
cat(sprintf("Area del borde conservado: %.0f km2\n", expanse(borde, unit = "km")))

figura("fig02_borde_valle.png", {
  par(mar = c(4, 4, 3, 4))
  plot(precip, main = "Precipitacion semana 44 (mm) - Valle del Cauca",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde, add = TRUE, border = "black", lwd = 1.5)
}, ancho = 1800, alto = 1900)


# ---- 2.2 Sorteo de puntos de observacion ------------------------------------
# Se simula una red de estaciones sorteando celdas dentro del borde.
# Mismo procedimiento de Ejemplo4_Geoestadistica.R (lineas 29-35).
#
# Por que no usar las 688 celdas:
#   (a) el sistema kriging invierte una matriz (n+1)x(n+1) por cada punto
#       predicho y por cada iteracion del LOOCV; con n=688 es inviable;
#   (b) un semivariograma modela un proceso muestreado en estaciones. Con las
#       688 celdas ya se tiene el mapa y no hay nada que interpolar.

set.seed(2026)
n_puntos <- 70

puntos <- spatSample(borde, size = n_puntos, method = "random")

datos <- data.frame(
  lon     = crds(puntos)[, 1],
  lat     = crds(puntos)[, 2],
  precip  = extract(precip, puntos)[, 2],   # [,2] porque extract() antepone ID
  altitud = extract(r_alt,  puntos)[, 2]
)
datos <- na.omit(datos)

# Proyeccion a km. El variograma mide semivarianza contra DISTANCIA, y un grado
# de longitud no vale lo mismo que uno de latitud:
#   1 grado de latitud  = 110.574 km (constante)
#   1 grado de longitud = 111.320 km * cos(latitud)
lat0 <- mean(datos$lat)
lon0 <- mean(datos$lon)
datos$x_km <- (datos$lon - lon0) * 111.320 * cos(lat0 * pi / 180)
datos$y_km <- (datos$lat - lat0) * 110.574

n <- nrow(datos)
cat(sprintf("\nEstaciones simuladas: %d\n", n))
cat(sprintf("Pares para el semivariograma: %d\n", n * (n - 1) / 2))


# ---- 2.3 Estadistica descriptiva --------------------------------------------

media_p  <- mean(datos$precip)
mediana_p<- median(datos$precip)
sd_p     <- sd(datos$precip)
cv_p     <- 100 * sd_p / media_p
sw       <- shapiro.test(datos$precip)

cat("\n--- Precipitacion en las estaciones (mm/semana) ---\n")
print(summary(datos$precip))
cat(sprintf("Desviacion estandar : %.2f mm\n", sd_p))
cat(sprintf("Coef. de variacion  : %.1f %%\n", cv_p))
cat(sprintf("Media / mediana     : %.2f  (>1 = cola a la derecha)\n", media_p / mediana_p))
cat(sprintf("Shapiro-Wilk        : W = %.4f, p = %.3e  -> %s\n",
            sw$statistic, sw$p.value,
            ifelse(sw$p.value < 0.05, "SE RECHAZA normalidad", "no se rechaza")))

cat("\n--- Correlacion de Pearson con la precipitacion ---\n")
print(round(cor(datos[, c("precip", "altitud", "x_km", "y_km")])[1, ], 3))


################################################################################
# PASO 3. ANALISIS EXPLORATORIO ESPACIAL (EDA)
################################################################################

cat("\n================ PASO 3: EDA ESPACIAL ================\n")

# ---- 3.1 Distribucion univariada --------------------------------------------
figura("fig03_distribucion.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  hist(datos$precip, breaks = 12, col = "#9ecae1", border = "white",
       main = "(A) Distribucion de la precipitacion",
       xlab = "Precipitacion (mm/semana)", ylab = "Frecuencia", las = 1)
  abline(v = media_p,   col = "red",   lwd = 2, lty = 2)
  abline(v = mediana_p, col = "black", lwd = 2, lty = 3)
  legend("topright", bty = "n", cex = 0.75,
         legend = c(sprintf("Media = %.1f", media_p),
                    sprintf("Mediana = %.1f", mediana_p)),
         col = c("red", "black"), lty = c(2, 3), lwd = 2)
  boxplot(datos$precip, col = "#9ecae1", las = 1,
          main = "(B) Diagrama de caja",
          ylab = "Precipitacion (mm/semana)")
  par(mfrow = c(1, 1))
}, ancho = 2400, alto = 1200)


# ---- 3.2 Mapa de posting ----------------------------------------------------
# Tamano del simbolo proporcional al valor: muestra DONDE esta cada magnitud.
figura("fig04_posting.png", {
  par(mar = c(4.5, 4.5, 3, 1))
  tam <- 0.6 + (datos$precip - min(datos$precip)) / diff(range(datos$precip)) * 2.4
  plot(borde, border = "grey40", las = 1,
       main = "Mapa de posting: estaciones simuladas (semana 44)",
       xlab = "Longitud", ylab = "Latitud")
  points(datos$lon, datos$lat, pch = 21, bg = "#3182bd", cex = tam)
  legend("topright", bty = "o", bg = "white", cex = 0.75, pt.cex = c(0.8, 1.8, 3.0),
         pch = 21, pt.bg = "#3182bd",
         legend = c("~35 mm", "~110 mm", "~190 mm"))
}, ancho = 1800, alto = 1900)


# ---- 3.3 Relacion con las covariables ---------------------------------------
figura("fig05_covariables.png", {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  for (v in c("x_km", "y_km", "altitud")) {
    etiqueta <- switch(v,
      x_km = "Coordenada Este (km)",
      y_km = "Coordenada Norte (km)",
      altitud = "Altitud SRTM (m)")
    plot(datos[[v]], datos$precip, pch = 21, bg = "#31a354", cex = 1.2, las = 1,
         xlab = etiqueta, ylab = "Precipitacion (mm/semana)",
         main = sprintf("r = %.3f", cor(datos[[v]], datos$precip)))
    abline(lm(datos$precip ~ datos[[v]]), col = "red", lwd = 2)
  }
  par(mfrow = c(1, 1))
}, ancho = 2800, alto = 1100)


################################################################################
# PASO 4. MODELO DE TENDENCIA (DERIVA DE GRAN ESCALA)
#
# Kriging universal = tendencia deterministica + residuo espacialmente
# correlacionado. Aqui se estima la primera parte con lm(), igual que
# Ejemplo4_Geoestadistica.R (seccion "Extraer tendencia").
################################################################################

cat("\n================ PASO 4: TENDENCIA ================\n")

# ---- 4.1 Confusion entre altitud y longitud ---------------------------------
# La correlacion simple altitud-precipitacion es NEGATIVA (-0.651), pero la
# altitud esta fuertemente correlacionada con la longitud (0.755): la costa es
# baja Y occidental, la cordillera es alta Y oriental.
# Al controlar por posicion, el signo del efecto de la altitud se invierte.

cat(sprintf("\ncor(altitud, precip) = %.3f   (simple)\n", cor(datos$altitud, datos$precip)))
cat(sprintf("cor(altitud, x_km)   = %.3f   (colinealidad)\n", cor(datos$altitud, datos$x_km)))

mod_A <- lm(precip ~ x_km + y_km, data = datos)
mod_B <- lm(precip ~ x_km + y_km + altitud, data = datos)

cat("\n--- Modelo A: solo coordenadas ---\n")
print(round(summary(mod_A)$coefficients, 4))
cat(sprintf("R2 = %.4f   R2 ajustado = %.4f\n",
            summary(mod_A)$r.squared, summary(mod_A)$adj.r.squared))

cat("\n--- Modelo B: coordenadas + altitud ---\n")
print(round(summary(mod_B)$coefficients, 4))
cat(sprintf("R2 = %.4f   R2 ajustado = %.4f\n",
            summary(mod_B)$r.squared, summary(mod_B)$adj.r.squared))

cat("\n--- Comparacion de modelos anidados (test F) ---\n")
print(anova(mod_A, mod_B))

# Se conserva el modelo B: la altitud aporta informacion significativa (p<0.05)
# y su signo positivo es fisicamente coherente (ascenso orografico: a igual
# longitud, mas altura implica mas condensacion).
modelo_tendencia <- mod_B
datos$tendencia  <- fitted(modelo_tendencia)
datos$residual   <- residuals(modelo_tendencia)

# ---- 4.2 Diagnostico de los residuales --------------------------------------
sw_res <- shapiro.test(datos$residual)
cat(sprintf("\nResiduales: media = %.2e, sd = %.2f mm\n",
            mean(datos$residual), sd(datos$residual)))
cat(sprintf("Shapiro-Wilk sobre residuales: W = %.4f, p = %.4f -> %s\n",
            sw_res$statistic, sw_res$p.value,
            ifelse(sw_res$p.value < 0.05, "SE RECHAZA normalidad",
                   "NO se rechaza normalidad")))

figura("fig06_residuales.png", {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  hist(datos$residual, breaks = 12, col = "#fdae6b", border = "white", las = 1,
       main = "(A) Residuales del modelo", xlab = "Residual (mm)")
  qqnorm(datos$residual, pch = 21, bg = "#fdae6b", las = 1,
         main = "(B) Normal Q-Q")
  qqline(datos$residual, col = "red", lwd = 2)
  plot(datos$tendencia, datos$residual, pch = 21, bg = "#fdae6b", cex = 1.2, las = 1,
       xlab = "Valor ajustado (mm)", ylab = "Residual (mm)",
       main = "(C) Residual vs ajustado")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  par(mfrow = c(1, 1))
}, ancho = 2800, alto = 1100)
