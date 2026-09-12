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

# ---- Helper de mapas --------------------------------------------------------
# plot() sobre un SpatVector ajusta la region de dibujo al poligono exacto y no
# deja margen: cualquier legend() se sale del area y se corta. Este helper abre
# un lienzo vacio con espacio extra abajo y dibuja el borde encima.

mapa_base <- function(borde, titulo, ...) {
  plot(borde, border = "grey40", las = 1, main = titulo,
       xlab = "Longitud", ylab = "Latitud", ...)
}

# terra::plot deja el area de dibujo pegada al borde del poligono, asi que
# legend("topleft") se recorta por arriba. Esta funcion devuelve un punto
# seguro dentro del hueco del noroeste del departamento.
# OJO: as.vector(ext()) devuelve un vector CON NOMBRES (xmin, xmax, ymin, ymax).
# Sin unname(), c(x = e[1], ...) produce el nombre "x.xmin" y lg["x"] da NA;
# legend(NA, NA, ...) no dibuja nada y no lanza ningun error.
esquina_leyenda <- function(borde, baja = 0.08) {
  e <- unname(as.vector(ext(borde)))
  c(x = e[1], y = e[4] - baja * (e[4] - e[3]))
}


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
  mapa_base(borde, "Mapa de posting: estaciones simuladas (semana 44)")
  points(datos$lon, datos$lat, pch = 21, bg = "#3182bd", cex = tam)
  lg <- esquina_leyenda(borde)
  legend(lg["x"], lg["y"], bty = "o", bg = "white", cex = 0.8, pt.cex = c(0.8, 1.8, 3.0),
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


# ---- 4.3 Representacion visual de la confusion ------------------------------
# Cuatro paneles que muestran por que la correlacion simple de la altitud
# enganya, y cual es su efecto real.

# Bandas de longitud (terciles) para el analisis estratificado
cortes_x <- quantile(datos$x_km, c(0, 1/3, 2/3, 1))
datos$banda <- cut(datos$x_km, breaks = cortes_x, include.lowest = TRUE,
                   labels = c("Oeste (pacifico)", "Centro", "Este (valle)"))

col_banda <- c("#1a9850", "#fdae61", "#d73027")
pal_lon <- colorRampPalette(c("#1a9850", "#fdae61", "#d73027"))(100)
idx_lon <- round((datos$x_km - min(datos$x_km)) / diff(range(datos$x_km)) * 99) + 1

cat("\n--- Pendiente precip~altitud DENTRO de cada banda de longitud ---\n")
for (b in levels(datos$banda)) {
  sub <- datos[datos$banda == b, ]
  pend <- coef(lm(precip ~ altitud, data = sub))[2]
  cat(sprintf("  %-18s n=%2d   pendiente = %+.4f mm/m\n", b, nrow(sub), pend))
}
cat(sprintf("  %-18s n=%2d   pendiente = %+.4f mm/m  <- agrupado\n",
            "TODOS JUNTOS", nrow(datos), coef(lm(precip ~ altitud, data = datos))[2]))

# --- 4.3a La colinealidad: altitud contra longitud ---
figura("fig07a_colinealidad.png", {
  par(mar = c(4.5, 4.5, 3.5, 1))
  plot(datos$x_km, datos$altitud, pch = 21, bg = pal_lon[idx_lon], cex = 1.5, las = 1,
       xlab = "Coordenada Este (km)", ylab = "Altitud SRTM (m)",
       main = sprintf("Altitud y longitud van juntas (r = %.3f)",
                      cor(datos$x_km, datos$altitud)))
  abline(lm(altitud ~ x_km, data = datos), col = "black", lwd = 2)
  legend("topleft", bty = "n", cex = 0.8, pch = 21, pt.bg = col_banda,
         legend = levels(datos$banda))
}, ancho = 1800, alto = 1500)

# --- 4.3b La correlacion simple, enganyosa ---
figura("fig07b_correlacion_simple.png", {
  par(mar = c(4.5, 4.5, 3.5, 1))
  plot(datos$altitud, datos$precip, pch = 21, bg = pal_lon[idx_lon], cex = 1.5, las = 1,
       xlab = "Altitud SRTM (m)", ylab = "Precipitacion (mm/semana)",
       main = sprintf("Correlacion simple r = %.3f: 'mas alto = mas seco' (FALSO)",
                      cor(datos$altitud, datos$precip)))
  abline(lm(precip ~ altitud, data = datos), col = "red", lwd = 2.5)
  legend("topright", bty = "n", cex = 0.8, pch = 21, pt.bg = col_banda,
         legend = levels(datos$banda), title = "Banda de longitud")
}, ancho = 1800, alto = 1500)

# --- 4.3c Estratificado por banda de longitud ---
figura("fig07c_estratificado.png", {
  par(mar = c(4.5, 4.5, 3.5, 1))
  plot(datos$altitud, datos$precip, type = "n", las = 1,
       xlab = "Altitud SRTM (m)", ylab = "Precipitacion (mm/semana)",
       main = "Estratificado por longitud: el efecto de la altitud NO es homogeneo")
  for (k in seq_along(levels(datos$banda))) {
    sub <- datos[datos$banda == levels(datos$banda)[k], ]
    points(sub$altitud, sub$precip, pch = 21, bg = col_banda[k], cex = 1.5)
    aj <- lm(precip ~ altitud, data = sub); xs <- range(sub$altitud)
    lines(xs, predict(aj, data.frame(altitud = xs)), col = col_banda[k], lwd = 3)
  }
  legend("topright", bty = "n", cex = 0.85, lwd = 3, col = col_banda,
         legend = sprintf("%s: %+.4f mm/m", levels(datos$banda),
                          sapply(levels(datos$banda), function(b)
                            coef(lm(precip ~ altitud, data = datos[datos$banda == b, ]))[2])))
}, ancho = 1800, alto = 1500)

# --- 4.3d Grafico de variable anyadida: el efecto parcial real ---
figura("fig07d_variable_anyadida.png", {
  par(mar = c(4.5, 4.5, 3.5, 1))
  res_precip  <- residuals(lm(precip  ~ x_km + y_km, data = datos))
  res_altitud <- residuals(lm(altitud ~ x_km + y_km, data = datos))
  plot(res_altitud, res_precip, pch = 21, bg = "#4575b4", cex = 1.5, las = 1,
       xlab = "Altitud | quitada la posicion (m)",
       ylab = "Precipitacion | quitada la posicion (mm)",
       main = sprintf("Efecto parcial de la altitud = %+.4f mm/m",
                      coef(lm(res_precip ~ res_altitud))[2]))
  abline(lm(res_precip ~ res_altitud), col = "blue", lwd = 2.5)
  abline(h = 0, v = 0, col = "grey70", lty = 3)
}, ancho = 1800, alto = 1500)


# ---- 4.4 Las mismas relaciones, sobre el mapa -------------------------------
# Los graficos anteriores son en el espacio de las variables. Estos son en el
# espacio geografico: muestran DONDE ocurre cada cosa.

# --- 4.4a Las tres bandas de longitud sobre el territorio ---
figura("fig08a_mapa_bandas.png", {
  par(mar = c(4.5, 4.5, 3.5, 1))
  mapa_base(borde, "Bandas de longitud del analisis estratificado")
  for (k in seq_along(levels(datos$banda))) {
    sub <- datos[datos$banda == levels(datos$banda)[k], ]
    points(sub$lon, sub$lat, pch = 21, bg = col_banda[k], cex = 1.6)
  }
  abline(v = lon0 + cortes_x[2:3] / (111.320 * cos(lat0 * pi / 180)),
         col = "grey30", lty = 2, lwd = 2)
  lg <- esquina_leyenda(borde)
  legend(lg["x"], lg["y"], bty = "o", bg = "white", cex = 0.8, pch = 21, pt.cex = 1.4,
         pt.bg = col_banda, legend = levels(datos$banda))
}, ancho = 1700, alto = 1900)

# --- 4.4b Altitud del terreno con las estaciones encima ---
figura("fig08b_mapa_altitud.png", {
  par(mar = c(4, 4, 3.5, 4))
  plot(mask(r_alt, mascara, maskvalues = c(FALSE, NA)),
       col = terrain.colors(50), las = 1,
       main = "Altitud SRTM (m) y estaciones simuladas",
       xlab = "Longitud", ylab = "Latitud")
  plot(borde, add = TRUE, border = "black", lwd = 1.2)
  points(datos$lon, datos$lat, pch = 3, col = "black", cex = 0.9, lwd = 1.4)
}, ancho = 1700, alto = 1900)

# --- 4.4c Residuales sobre el mapa: azul negativo, rojo positivo ---
# Si los colores aparecen agrupados en manchas, hay autocorrelacion espacial
# remanente y el kriging tiene trabajo que hacer.
figura("fig08c_mapa_residuales.png", {
  par(mar = c(4.5, 4.5, 3.5, 1))
  tam_res <- 0.7 + abs(datos$residual) / max(abs(datos$residual)) * 2.3
  col_res <- ifelse(datos$residual >= 0, "#d73027", "#4575b4")
  mapa_base(borde, "Residuales de la tendencia (rojo +, azul -)")
  points(datos$lon, datos$lat, pch = 21, bg = col_res, cex = tam_res)
  lg <- esquina_leyenda(borde)
  legend(lg["x"], lg["y"], bty = "o", bg = "white", cex = 0.78, pch = 21,
         pt.bg = c("#d73027", "#4575b4"), pt.cex = 1.6,
         legend = c("Subestima el modelo (+)", "Sobreestima el modelo (-)"))
}, ancho = 1700, alto = 1900)

# --- 4.4d Observado contra ajustado, lado a lado en el mapa ---
figura("fig08d_mapa_obs_vs_tendencia.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1))
  rango <- range(c(datos$precip, datos$tendencia))
  escala <- function(v) 0.6 + (v - rango[1]) / diff(rango) * 2.4
  mapa_base(borde, "(A) Precipitacion observada")
  points(datos$lon, datos$lat, pch = 21, bg = "#3182bd", cex = escala(datos$precip))
  mapa_base(borde, "(B) Tendencia ajustada por el modelo")
  points(datos$lon, datos$lat, pch = 21, bg = "#31a354", cex = escala(datos$tendencia))
  par(mfrow = c(1, 1))
}, ancho = 3000, alto = 1700)


################################################################################
# PASO 5. AUTOCORRELACION ESPACIAL (MORAN Y GEARY)
#
# Es el corazon del enunciado: "evaluar el papel de la correlacion espacial".
# Indices calculados a mano segun carpeta/Script_Spatial_Analytics.R.
################################################################################

cat("\n================ PASO 5: AUTOCORRELACION ESPACIAL ================\n")

coords <- as.matrix(datos[, c("x_km", "y_km")])
D <- as.matrix(dist(coords))       # matriz de distancias, idioma de Ejemplo4

# ---- 5.1 Matriz de pesos e indices ------------------------------------------
# OJO: el script de clase estandariza las filas con
#   W <- matrix(ifelse(row_sums == 0, 0, W / row_sums), nrow = n, ncol = n)
# ifelse() devuelve un objeto del largo del TEST (n), no del resultado (n^2):
# se queda con la primera columna y matrix() la recicla, asi que la
# estandarizacion NO ocurre y las filas no suman 1. Aqui se usa sweep().

pesos_espaciales <- function(D, q) {
  n <- nrow(D)
  umbral <- quantile(D[upper.tri(D)], q)
  W <- matrix(0, n, n)
  W[D <= umbral] <- 1
  diag(W) <- 0
  filas <- rowSums(W)
  list(W = sweep(W, 1, ifelse(filas == 0, 1, filas), "/"),
       umbral = unname(umbral), vecinos = mean(filas))
}

moran <- function(z, W) {
  n <- length(z); zc <- z - mean(z)
  (n / sum(W)) * (sum(W * outer(zc, zc)) / sum(zc^2))
}

geary <- function(z, W) {
  n <- length(z); zc <- z - mean(z)
  ((n - 1) / (2 * sum(W))) * (sum(W * outer(z, z, function(a, b) (a - b)^2)) / sum(zc^2))
}

# Significancia por permutacion: se re-etiquetan los valores al azar sobre las
# mismas posiciones. Si la estructura fuera casual, el I observado caeria dentro
# de la nube de I permutados.
p_permutacion <- function(z, W, n_perm = 999) {
  obs <- moran(z, W)
  set.seed(2026)
  nulos <- replicate(n_perm, moran(sample(z), W))
  (1 + sum(nulos >= obs)) / (n_perm + 1)
}

# ---- 5.2 Sensibilidad al umbral de vecindad ---------------------------------
# El script de clase usa el percentil 5 de las distancias. La eleccion cambia el
# resultado, asi que se reporta un rango en vez de un solo numero.

tabla_ac <- do.call(rbind, lapply(c(0.05, 0.10, 0.25), function(q) {
  pw <- pesos_espaciales(D, q)
  data.frame(
    percentil    = q * 100,
    umbral_km    = pw$umbral,
    vecinos_prom = pw$vecinos,
    I_precip     = moran(datos$precip,   pw$W),
    C_precip     = geary(datos$precip,   pw$W),
    I_residual   = moran(datos$residual, pw$W),
    C_residual   = geary(datos$residual, pw$W),
    p_residual   = p_permutacion(datos$residual, pw$W))
}))

cat("\n--- Indices por umbral de vecindad ---\n")
print(format(tabla_ac, digits = 4))
cat(sprintf("\nValor esperado bajo aleatoriedad: E[I] = %.4f   |   E[C] = 1\n", -1 / (n - 1)))

write.csv(tabla_ac, file.path(DIR_SALIDA, "tabla_autocorrelacion.csv"), row.names = FALSE)

# ---- 5.3 Diagrama de dispersion de Moran ------------------------------------
# Eje x: valor estandarizado. Eje y: promedio de sus vecinos (rezago espacial).
# La pendiente de la recta ES el indice de Moran.

pw_ref <- pesos_espaciales(D, 0.10)
z_std  <- as.numeric(scale(datos$residual))
lag_z  <- as.numeric(pw_ref$W %*% z_std)
I_ref  <- moran(datos$residual, pw_ref$W)

figura("fig09_moran_scatter.png", {
  par(mar = c(4.5, 4.5, 3.5, 1))
  col_cuad <- ifelse(z_std >= 0 & lag_z >= 0, "#d73027",
              ifelse(z_std <  0 & lag_z <  0, "#4575b4", "grey65"))
  plot(z_std, lag_z, pch = 21, bg = col_cuad, cex = 1.4, las = 1,
       xlab = "Residual estandarizado", ylab = "Promedio de los vecinos (rezago espacial)",
       main = sprintf("Diagrama de Moran de los residuales (I = %.3f)", I_ref))
  abline(h = 0, v = 0, col = "grey60", lty = 3)
  abline(lm(lag_z ~ z_std), col = "black", lwd = 2.5)
  legend("topleft", bty = "n", cex = 0.75, pch = 21,
         pt.bg = c("#d73027", "#4575b4", "grey65"),
         legend = c("Alto rodeado de alto", "Bajo rodeado de bajo", "Discordante"))
}, ancho = 1800, alto = 1600)

# ---- 5.4 Correlograma: Moran por anillo de distancia ------------------------
# Muestra hasta donde llega la dependencia espacial. El corte por cero anticipa
# el rango practico que estimara el semivariograma.

anillos <- seq(0, max(D) * 0.7, length.out = 9)
correlograma <- do.call(rbind, lapply(seq_len(length(anillos) - 1), function(k) {
  W <- matrix(0, n, n)
  W[D > anillos[k] & D <= anillos[k + 1]] <- 1
  diag(W) <- 0
  filas <- rowSums(W)
  if (sum(filas) == 0) return(NULL)
  W <- sweep(W, 1, ifelse(filas == 0, 1, filas), "/")
  data.frame(h_medio = mean(anillos[k:(k + 1)]),
             pares = sum(W > 0),
             I = moran(datos$residual, W))
}))

cat("\n--- Correlograma de los residuales ---\n")
print(format(correlograma, digits = 4))

figura("fig10_correlograma.png", {
  par(mar = c(4.5, 4.5, 3.5, 1))
  plot(correlograma$h_medio, correlograma$I, type = "o", pch = 19, cex = 1.3,
       col = "#08519c", lwd = 2, las = 1,
       xlab = "Distancia h (km)", ylab = "Indice de Moran I",
       main = "Correlograma de los residuales de la tendencia")
  abline(h = -1 / (n - 1), col = "red", lty = 2, lwd = 2)
  text(max(correlograma$h_medio) * 0.75, -1 / (n - 1), pos = 3, cex = 0.75,
       col = "red", labels = "E[I] bajo aleatoriedad")
  grid(col = "grey85")
}, ancho = 1900, alto = 1400)


################################################################################
# PASO 6. SEMIVARIOGRAMA EMPIRICO Y AJUSTE DE MODELOS TEORICOS
#
# Respaldo del material:
#
# - Lecture_SpatialStatistics.pdf, ec. (11):
#     gamma(h) = 1/(2 N(h)) * SUM (z_i - z_j)^2
#   y tres modelos validos, ec. (12) esferico, (13) exponencial, (14) gaussiano,
#   parametrizados por pepita c0, meseta c0+c1 y rango a.
#
# - Geoestadistica (libro del curso):
#     "En caso de encontrarse tendencia en los datos, esta tendencia se modela
#      con modelos de regresion polinomicos [...] y el semivariograma se
#      construye con los residuales obtenidos"
#     "en general, se usan los rezagos espaciales hasta la mitad de la maxima
#      distancia"
#     "Es importante elegir una distancia maxima [...] en caso de que se observe
#      un comportamiento erratico a distancias mayores"
#
# Se llevan DOS tendencias en paralelo (lineal y cuadratica) y se comparan.
################################################################################

cat("\n================ PASO 6: SEMIVARIOGRAMA ================\n")

# ---- 6.1 Las dos tendencias en competencia ----------------------------------
formulas_tendencia <- list(
  lineal     = precip ~ x_km + y_km + altitud,
  cuadratica = precip ~ x_km + I(x_km^2) + y_km + altitud
)

modelos_tend <- lapply(formulas_tendencia, function(f) lm(f, data = datos))
residuales_tend <- lapply(modelos_tend, residuals)

cat("\n--- Comparacion de tendencias ---\n")
print(do.call(rbind, lapply(names(modelos_tend), function(nm) {
  m <- modelos_tend[[nm]]; r <- residuals(m)
  pw <- pesos_espaciales(D, 0.10)
  data.frame(tendencia = nm,
             R2_ajustado = round(summary(m)$adj.r.squared, 4),
             sd_residual = round(sd(r), 2),
             shapiro_p   = round(shapiro.test(r)$p.value, 4),
             moran_residual = round(moran(r, pw$W), 3))
})))

# ---- 6.2 Semivariograma empirico --------------------------------------------
# Ecuacion (11) de la Lecture, agrupando los pares en intervalos de distancia.

semivariograma_empirico <- function(resid, D, cutoff, n_bins = 12) {
  pares <- which(upper.tri(D), arr.ind = TRUE)
  h     <- D[upper.tri(D)]
  gamma <- (resid[pares[, 1]] - resid[pares[, 2]])^2 / 2
  en    <- h <= cutoff
  bins  <- cut(h[en], breaks = n_bins)
  na.omit(data.frame(h = tapply(h[en], bins, mean),
                     gamma = tapply(gamma[en], bins, mean),
                     N = tapply(gamma[en], bins, length)))
}

D_MAX <- max(D)
CUTOFF_LIBRO <- D_MAX / 2        # regla general del libro
CUTOFF_CORR  <- 50               # donde el correlograma cruza cero (Paso 5)

cat(sprintf("\nDistancia maxima entre estaciones: %.1f km\n", D_MAX))
cat(sprintf("Cutoff regla del libro (mitad):    %.1f km\n", CUTOFF_LIBRO))
cat(sprintf("Cutoff por comportamiento erratico: %.1f km  (correlograma, Paso 5)\n", CUTOFF_CORR))

# ---- 6.3 Los tres modelos validos de la Lecture -----------------------------
mod_esferico <- function(h, c0, c1, a)
  ifelse(h <= a, c0 + c1 * (1.5 * (h / a) - 0.5 * (h / a)^3), c0 + c1)
mod_exponencial <- function(h, c0, c1, a) c0 + c1 * (1 - exp(-h / a))
mod_gaussiano   <- function(h, c0, c1, a) c0 + c1 * (1 - exp(-(h / a)^2))

MODELOS <- list(esferico = mod_esferico, exponencial = mod_exponencial,
                gaussiano = mod_gaussiano)

# Ajuste por minimos cuadrados ordinarios sobre el semivariograma empirico
# (libro, seccion 3.1), con optim L-BFGS-B y varios arranques porque los
# modelos no son lineales en los parametros.
ajustar_modelo <- function(sv, fn) {
  sce <- function(th) sum((sv$gamma - fn(sv$h, th[1], th[2], th[3]))^2)
  arranques <- quantile(sv$h, c(0.1, 0.25, 0.5, 0.75))
  cand <- lapply(arranques, function(a0)
    optim(c(c0 = 0, c1 = max(sv$gamma), a = a0), sce, method = "L-BFGS-B",
          lower = c(0, 0.01, 1),
          upper = c(max(sv$gamma), 5 * max(sv$gamma), 3 * max(sv$h))))
  mejor <- cand[[which.min(sapply(cand, function(x) x$value))]]
  c(c0 = unname(mejor$par[1]), c1 = unname(mejor$par[2]),
    a = unname(mejor$par[3]), SCE = mejor$value)
}

# ---- 6.4 Rejilla de ajustes: 2 tendencias x 2 cutoffs x 3 modelos -----------
rejilla <- expand.grid(tendencia = names(residuales_tend),
                       cutoff = c(CUTOFF_LIBRO, CUTOFF_CORR),
                       modelo = names(MODELOS), stringsAsFactors = FALSE)

resultados_sv <- do.call(rbind, lapply(seq_len(nrow(rejilla)), function(i) {
  fila <- rejilla[i, ]
  sv <- semivariograma_empirico(residuales_tend[[fila$tendencia]], D, fila$cutoff)
  aj <- ajustar_modelo(sv, MODELOS[[fila$modelo]])
  # Rango practico: distancia a la que gamma alcanza el 95% de la meseta
  rp <- switch(fila$modelo, esferico = aj["a"],
               exponencial = 3 * aj["a"], gaussiano = sqrt(3) * aj["a"])
  data.frame(tendencia = fila$tendencia, cutoff_km = round(fila$cutoff, 1),
             modelo = fila$modelo, pepita = round(aj["c0"], 1),
             meseta = round(aj["c0"] + aj["c1"], 1),
             rango_practico_km = round(rp, 1), SCE = round(aj["SCE"], 1),
             row.names = NULL)
}))

cat("\n--- Ajustes: 2 tendencias x 2 cutoffs x 3 modelos ---\n")
print(resultados_sv[order(resultados_sv$tendencia, resultados_sv$cutoff_km,
                          resultados_sv$SCE), ])
write.csv(resultados_sv, file.path(DIR_SALIDA, "tabla_semivariogramas.csv"),
          row.names = FALSE)


# ---- 6.5 Diagnostico: meseta contra varianza muestral -----------------------
# Bajo estacionariedad de segundo orden (libro, Definicion 4) la meseta del
# semivariograma debe aproximar la varianza del proceso. Si la meseta la supera
# claramente, el semivariograma sigue creciendo: queda tendencia sin modelar y
# el supuesto no se cumple.

varianzas <- sapply(residuales_tend, var)
cat("\n--- Meseta estimada contra varianza de los residuales ---\n")
diag_meseta <- do.call(rbind, lapply(seq_len(nrow(resultados_sv)), function(i) {
  f <- resultados_sv[i, ]
  data.frame(tendencia = f$tendencia, cutoff = f$cutoff_km, modelo = f$modelo,
             meseta = f$meseta,
             varianza = round(varianzas[[f$tendencia]], 1),
             razon = round(f$meseta / varianzas[[f$tendencia]], 2),
             row.names = NULL)
}))
diag_meseta$veredicto <- ifelse(abs(diag_meseta$razon - 1) <= 0.25, "coherente",
                         ifelse(diag_meseta$razon > 1.25, "meseta > varianza: NO estacionario",
                                "meseta < varianza"))
print(diag_meseta[order(diag_meseta$tendencia, diag_meseta$cutoff), ])
write.csv(diag_meseta, file.path(DIR_SALIDA, "tabla_diagnostico_meseta.csv"), row.names = FALSE)


# ---- 6.6 Figuras del semivariograma -----------------------------------------
col_mod <- c(esferico = "#d73027", exponencial = "#1a9850", gaussiano = "#4575b4")

figura_variograma <- function(nombre_tend, cutoff, titulo, archivo) {
  resid <- residuales_tend[[nombre_tend]]
  sv <- semivariograma_empirico(resid, D, cutoff)
  pares_idx <- which(upper.tri(D), arr.ind = TRUE)
  h_all <- D[upper.tri(D)]
  g_all <- (resid[pares_idx[, 1]] - resid[pares_idx[, 2]])^2 / 2
  en <- h_all <= cutoff

  figura(archivo, {
    par(mar = c(4.5, 4.5, 3.5, 1))
    plot(h_all[en], g_all[en], pch = 20, col = "grey80", cex = 0.6, las = 1,
         xlab = "Distancia h (km)", ylab = expression(gamma(h)),
         main = titulo, ylim = c(0, quantile(sv$gamma, 1) * 1.9))
    abline(h = var(resid), col = "grey30", lty = 3, lwd = 2)
    hs <- seq(0, cutoff, length.out = 250)
    for (m in names(MODELOS)) {
      aj <- ajustar_modelo(sv, MODELOS[[m]])
      lines(hs, MODELOS[[m]](hs, aj["c0"], aj["c1"], aj["a"]),
            col = col_mod[m], lwd = 2.5)
    }
    points(sv$h, sv$gamma, pch = 19, col = "black", cex = 1.5)
    legend("bottomright", bty = "o", bg = "white", cex = 0.75,
           legend = c("nube de pares", "semivariograma empirico",
                      names(MODELOS), "varianza muestral"),
           col = c("grey80", "black", col_mod, "grey30"),
           pch = c(20, 19, NA, NA, NA, NA),
           lty = c(0, 0, 1, 1, 1, 3), lwd = c(0, 0, 2.5, 2.5, 2.5, 2))
  }, ancho = 2000, alto = 1500)
}

figura_variograma("lineal", CUTOFF_CORR,
  sprintf("Tendencia LINEAL - cutoff %g km", CUTOFF_CORR),
  "fig11a_variograma_lineal.png")
figura_variograma("cuadratica", CUTOFF_CORR,
  sprintf("Tendencia CUADRATICA - cutoff %g km", CUTOFF_CORR),
  "fig11b_variograma_cuadratica.png")
figura_variograma("cuadratica", CUTOFF_LIBRO,
  sprintf("Tendencia CUADRATICA - cutoff %.0f km (regla del libro)", CUTOFF_LIBRO),
  "fig11c_variograma_cuadratica_libro.png")
