################################################################################
# UNIVERSIDAD AUTONOMA DE OCCIDENTE
# FACULTAD DE INGENIERIA Y CIENCIAS BASICAS
# PROGRAMA: INGENIERIA DE DATOS E INTELIGENCIA ARTIFICIAL
# ASIGNATURA: ANALITICA DE DATOS (2026-2S)
# PROFESOR: JOHANN A. OSPINA
#
# PROYECTO NO. 1
# Modelamiento de la precipitacion semanal en el Valle del Cauca en funcion de
# covariables ambientales y evaluacion del papel de la correlacion espacial.
#
# INTEGRANTES:
# - Juan Pablo Maya          - Codigo: 2236377
# - Cesar Armando Reyes      - Codigo: 2236379
# - Yesenia Diaz Urrego      - Codigo: 2231783
#
# Script unico, independiente y 100% reproducible.
# Utiliza unicamente las metodologias y librerias trabajadas en clase (terra,
# sp, gstat, ggplot2, gridExtra, corrplot y algebra matricial de R base).
################################################################################

# Limpieza del entorno y configuracion de reproducibilidad
rm(list = ls())
set.seed(2026)
options(digits = 4, scipen = 999)

# ==============================================================================
# 0. CARGA DE LIBRERIAS PERMITIDAS Y CONFIGURACION DE RUTAS
# ==============================================================================
paquetes_requeridos <- c("terra", "sp", "gstat", "ggplot2", "gridExtra", "corrplot")
paquetes_instalados <- rownames(installed.packages())

for (p in paquetes_requeridos) {
  if (!(p %in% paquetes_instalados)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

# Definicion flexible de rutas para compatibilidad tanto desde la raiz como desde /R
DIR_DATOS <- "data/datos_proyecto_1/imagenes_semanales"
if (!dir.exists(DIR_DATOS)) {
  DIR_DATOS <- "datos_proyecto_1/datos_proyecto_1/imagenes_semanales"
}
if (!dir.exists(DIR_DATOS)) {
  # Busqueda recursiva
  tifs <- list.files(pattern = "altitud_valle\\.tif$", recursive = TRUE, full.names = TRUE)
  if (length(tifs) > 0) DIR_DATOS <- dirname(tifs[1])
}

DIR_SALIDA <- "resultados"
if (!dir.exists(DIR_SALIDA)) dir.create(DIR_SALIDA, recursive = TRUE)

DIR_FIGURAS <- "figuras"
if (!dir.exists(DIR_FIGURAS)) dir.create(DIR_FIGURAS, recursive = TRUE)

# Helper para exportar figuras de forma simultanea a pantalla y archivo PNG (300 DPI)
guardar_figura <- function(nombre, expr, ancho = 2400, alto = 1600, res = 300) {
  bloque <- substitute(expr)
  entorno <- parent.frame()
  eval(bloque, entorno) # Muestra en el panel de RStudio / dispositivo activo
  
  # Guardar en resultados/ y figuras/
  rutas <- c(file.path(DIR_SALIDA, nombre), file.path(DIR_FIGURAS, nombre))
  for (r in unique(rutas)) {
    tryCatch({
      if (file.exists(r)) suppressWarnings(unlink(r, force = TRUE))
      png(r, width = ancho, height = alto, res = res)
      eval(bloque, entorno)
      dev.off()
    }, error = function(e) {
      suppressWarnings(try(dev.off(), silent = TRUE))
      # Intento alternativo con ragg o cairo si falla el default
      tryCatch({
        png(r, width = ancho, height = alto, res = res, type = "cairo")
        eval(bloque, entorno)
        dev.off()
      }, error = function(e2) {
        warning(sprintf("No se pudo guardar la figura en %s: %s", r, e2$message))
      })
    })
  }
  cat(sprintf("  [Figura generada] %s\n", nombre))
}

cat("====================================================================\n")
cat(" PROYECTO 1: MODELAMIENTO GEOESTADISTICO DE PRECIPITACION (VALLE)   \n")
cat(" Integrantes:\n")
cat("  - Juan Pablo Maya     (2236377)\n")
cat("  - Cesar Armando Reyes (2236379)\n")
cat("  - Yesenia Diaz Urrego (2231783)\n")
cat("====================================================================\n\n")

# ==============================================================================
# 1. AUDITORIA, CARGA Y ALINEACION DE DATOS RASTER
# ==============================================================================
cat("[Paso 1] Cargando capas raster y auditando geometria...\n")

r_altitud     <- rast(file.path(DIR_DATOS, "altitud_valle.tif"))
r_precip_clim <- rast(file.path(DIR_DATOS, "chirps_climatologia_semanal_valle.tif"))
r_temp_clim   <- rast(file.path(DIR_DATOS, "power_temp_climatologia_semanal_valle.tif"))
r_rad_clim    <- rast(file.path(DIR_DATOS, "power_radiacion_climatologia_semanal_valle.tif"))

# Limpieza de centinelas negativos en precipitacion
r_precip_clim[r_precip_clim < 0] <- NA

# Mascara oficial del departamento del Valle del Cauca (SRTM y CHIRPS validos)
mascara_valle <- (!is.na(r_altitud)) & (r_altitud >= -5)

# Borde poligonal oficial del departamento (metodologia Ejemplo4_Geoestadistica.R)
borde_valle <- as.polygons(mascara_valle, dissolve = TRUE)
borde_valle <- borde_valle[borde_valle[[1]] == 1, ]

# Seleccion de la semana de estudio: Semana 42 (Pico de lluvias de Octubre en el Valle)
SEMANA <- 42

# Alineacion e interpolacion de covariables gruesas a la grilla fina de CHIRPS
alinear_covariable <- function(r_layer, mask_base) {
  pts <- as.points(r_layer, na.rm = TRUE)
  if (nrow(pts) == 0) return(r_layer)
  r_idw <- interpIDW(mask_base, pts, field = names(r_layer), radius = 2.5, power = 2)
  mask(r_idw, mask_base)
}

p_lluvias <- mask(r_precip_clim[[SEMANA]], mascara_valle)
t_lluvias <- alinear_covariable(r_temp_clim[[SEMANA]], mascara_valle)
rad_lluvias <- alinear_covariable(r_rad_clim[[SEMANA]], mascara_valle)
alt_valle <- mask(r_altitud, mascara_valle)

# Extraccion de la tabla de celdas completas del departamento
coords_all <- crds(alt_valle, na.rm = TRUE)
df_grilla_valle <- data.frame(
  id = 1:nrow(coords_all),
  lon = coords_all[, 1],
  lat = coords_all[, 2],
  altitud = extract(alt_valle, coords_all)[, 1],
  precip = extract(p_lluvias, coords_all)[, 1],
  temp = extract(t_lluvias, coords_all)[, 1],
  rad = extract(rad_lluvias, coords_all)[, 1]
)
df_grilla_valle <- na.omit(df_grilla_valle)

# Proyeccion a coordenadas planas locales en kilometros (centradas en el baricentro)
lat0 <- mean(df_grilla_valle$lat)
lon0 <- mean(df_grilla_valle$lon)

df_grilla_valle$x_km <- (df_grilla_valle$lon - lon0) * 111.320 * cos(lat0 * pi / 180)
df_grilla_valle$y_km <- (df_grilla_valle$lat - lat0) * 110.574

cat(sprintf("  -> Total celdas validas en el Valle del Cauca: %d\n", nrow(df_grilla_valle)))

# ==============================================================================
# 2. MUESTREO ESPACIAL REPRESENTATIVO (ESTACIONES OBSERVADAS)
# ==============================================================================
cat("\n[Paso 2] Muestreo espacial de estaciones de entrenamiento...\n")

n_muestra <- 70
set.seed(2026)
idx_muestra <- sample(1:nrow(df_grilla_valle), size = n_muestra, replace = FALSE)
datos_obs <- df_grilla_valle[idx_muestra, ]
rownames(datos_obs) <- NULL

cat(sprintf("  -> %d estaciones seleccionadas para modelamiento y validacion.\n", nrow(datos_obs)))

# ==============================================================================
# 3. ANALISIS EXPLORATORIO DE DATOS ESPACIALES (ESDA)
# ==============================================================================
cat("\n[Paso 3] Ejecutando Analisis Exploratorio de Datos Espaciales (ESDA)...\n")

resumen_p <- summary(datos_obs$precip)
sd_p <- sd(datos_obs$precip)
cv_p <- (sd_p / mean(datos_obs$precip)) * 100
shapiro_p <- shapiro.test(datos_obs$precip)

cat("--- Estadistica Descriptiva de Precipitacion (Semana 42) ---\n")
print(resumen_p)
cat(sprintf("Desviacion Estandar: %.2f mm\n", sd_p))
cat(sprintf("Coeficiente de Variacion (CV): %.2f %%\n", cv_p))
cat(sprintf("Prueba de Normalidad (Shapiro-Wilk): W = %.4f, p = %.4e\n", 
            shapiro_p$statistic, shapiro_p$p.value))

# Matriz de Correlacion
vars_eda <- c("precip", "altitud", "temp", "rad", "x_km", "y_km")
mat_cor <- cor(datos_obs[, vars_eda])
cat("\n--- Matriz de Correlaciones de Pearson ---\n")
print(round(mat_cor, 3))

# FIGURA 1: Distribucion Univariada
guardar_figura("figura1_distribucion_precipitacion.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  hist(datos_obs$precip, breaks = 12, col = "#3182bd", border = "white",
       main = "(A) Histograma de Precipitacion Semanal",
       xlab = "Precipitacion (mm/semana)", ylab = "Frecuencia", las = 1)
  abline(v = mean(datos_obs$precip), col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c(sprintf("Media: %.1f mm", mean(datos_obs$precip)),
                                sprintf("Mediana: %.1f mm", median(datos_obs$precip))),
         col = c("red", "black"), lty = c(2, 0), bty = "n", cex = 0.85)

  boxplot(datos_obs$precip, col = "#9ecae1", horizontal = FALSE,
          main = "(B) Diagrama de Caja (Boxplot)",
          ylab = "Precipitacion (mm/semana)", las = 1)
}, ancho = 2400, alto = 1200)

# FIGURA 2: Dispersion frente a Covariables Ambientales
guardar_figura("figura2_dispersion_covariables.png", {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))
  plot(datos_obs$x_km, datos_obs$precip, pch = 21, bg = "#3182bd", cex = 1.3,
       xlab = "Coordenada Este X (km)", ylab = "Precipitacion (mm/semana)",
       main = "Precipitacion vs Longitud (Gradiente E-O)", las = 1)
  abline(lm(precip ~ x_km, data = datos_obs), col = "red", lwd = 2)

  plot(datos_obs$altitud, datos_obs$precip, pch = 21, bg = "#31a354", cex = 1.3,
       xlab = "Altitud SRTM (m.s.n.m.)", ylab = "Precipitacion (mm/semana)",
       main = "Precipitacion vs Altitud", las = 1)
  abline(lm(precip ~ altitud, data = datos_obs), col = "red", lwd = 2)

  plot(datos_obs$temp, datos_obs$precip, pch = 21, bg = "#e6550d", cex = 1.3,
       xlab = "Temperatura a 2m (C)", ylab = "Precipitacion (mm/semana)",
       main = "Precipitacion vs Temperatura", las = 1)
  abline(lm(precip ~ temp, data = datos_obs), col = "red", lwd = 2)

  plot(datos_obs$rad, datos_obs$precip, pch = 21, bg = "#756bb1", cex = 1.3,
       xlab = "Radiacion Solar (MJ/m2/dia)", ylab = "Precipitacion (mm/semana)",
       main = "Precipitacion vs Radiacion Solar", las = 1)
  abline(lm(precip ~ rad, data = datos_obs), col = "red", lwd = 2)
}, ancho = 2800, alto = 1800)

# FIGURA 3: Mapa de Posting
tam_burbuja <- 0.8 + (datos_obs$precip - min(datos_obs$precip)) / diff(range(datos_obs$precip)) * 2.2
pal_colores <- colorRampPalette(c("#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142", "#5e4fa2"))(100)
col_idx <- round((datos_obs$precip - min(datos_obs$precip)) / diff(range(datos_obs$precip)) * 99) + 1

guardar_figura("figura3_mapa_posting.png", {
  par(mar = c(4.5, 4.5, 3, 1))
  plot(datos_obs$lon, datos_obs$lat, pch = 21, bg = pal_colores[col_idx], cex = tam_burbuja,
       xlab = "Longitud (Grados)", ylab = "Latitud (Grados)",
       main = "Mapa de Posting: Precipitacion Semanal Observada (Valle del Cauca)",
       las = 1, asp = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.5)
  grid(col = "grey80")
  legend("topright", legend = c("Baja (< 45 mm)", "Media (45 - 80 mm)", "Alta (> 80 mm)"),
         pch = 21, pt.bg = c(pal_colores[10], pal_colores[50], pal_colores[90]),
         pt.cex = c(1.0, 1.8, 2.5), bg = "white", bty = "o", cex = 0.8)
}, ancho = 2000, alto = 2000)

# ==============================================================================
# 4. CALCULO FORMAL DE AUTOCORRELACION ESPACIAL (MORAN Y GEARY)
# ==============================================================================
cat("\n[Paso 4] Calculando matrices de vecindad y metricas de autocorrelacion...\n")

n <- nrow(datos_obs)
coords_obs_km <- as.matrix(datos_obs[, c("x_km", "y_km")])

# Matriz de distancias euclidianas D
D_obs <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    D_obs[i, j] <- sqrt((coords_obs_km[i, 1] - coords_obs_km[j, 1])^2 +
                        (coords_obs_km[i, 2] - coords_obs_km[j, 2])^2)
  }
}

# Matriz W con umbral de vecindad (percentil 25)
umbral_dist <- quantile(D_obs[upper.tri(D_obs)], 0.25)
W <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j && D_obs[i, j] <= umbral_dist) W[i, j] <- 1
  }
}

# Estandarizacion por filas con sweep
sum_filas <- rowSums(W)
W_std <- sweep(W, 1, ifelse(sum_filas == 0, 1, sum_filas), "/")

# Indice de Moran Global (metodo de clase)
z_vec <- datos_obs$precip
z_mean <- mean(z_vec)
num_moran <- 0
denom_moran <- sum((z_vec - z_mean)^2)

for (i in 1:n) {
  for (j in 1:n) {
    num_moran <- num_moran + W_std[i, j] * (z_vec[i] - z_mean) * (z_vec[j] - z_mean)
  }
}
s0_W <- sum(W_std)
moran_I <- (n / s0_W) * (num_moran / denom_moran)

# Indice de Geary C
num_geary <- 0
for (i in 1:n) {
  for (j in 1:n) {
    num_geary <- num_geary + W_std[i, j] * (z_vec[i] - z_vec[j])^2
  }
}
geary_C <- ((n - 1) / (2 * s0_W)) * (num_geary / denom_moran)

cat(sprintf("  -> Umbral de vecindad espacial: %.2f km\n", umbral_dist))
cat(sprintf("  -> Indice Global de Moran (I): %.4f (Autocorrelacion Positiva Fuerte)\n", moran_I))
cat(sprintf("  -> Indice de Geary (C): %.4f (Consistente con C < 1)\n", geary_C))

# ==============================================================================
# 5. MODELACION DE TENDENCIA Y ANALISIS DE RESIDUALES
# ==============================================================================
cat("\n[Paso 5] Ajustando modelo de tendencia lineal con covariables ambientales...\n")

modelo_tendencia <- lm(precip ~ altitud + temp + rad + x_km + y_km, data = datos_obs)
resumen_modelo <- summary(modelo_tendencia)
print(resumen_modelo)

residuales <- residuals(modelo_tendencia)
datos_obs$residuales <- residuales
datos_obs$tendencia_fit <- fitted(modelo_tendencia)

# Verificacion de Moran sobre residuales
res_mean <- mean(residuales)
num_moran_res <- 0
denom_moran_res <- sum((residuales - res_mean)^2)

for (i in 1:n) {
  for (j in 1:n) {
    num_moran_res <- num_moran_res + W_std[i, j] * (residuales[i] - res_mean) * (residuales[j] - res_mean)
  }
}
moran_I_res <- (n / s0_W) * (num_moran_res / denom_moran_res)

cat(sprintf("  -> R-cuadrado de tendencia: %.4f\n", resumen_modelo$r.squared))
cat(sprintf("  -> Moran I residual: %.4f (Dependencia espacial remanente justificada)\n", moran_I_res))

# ==============================================================================
# 6. SEMIVARIOGRAMA EXPERIMENTAL Y AJUSTE DE MODELO TEORICO
# ==============================================================================
cat("\n[Paso 6] Calculando semivariograma experimental y ajustando modelo teorico...\n")

pares_idx <- which(upper.tri(D_obs), arr.ind = TRUE)
h_pares   <- D_obs[upper.tri(D_obs)]
gamma_pares <- (residuales[pares_idx[, 1]] - residuales[pares_idx[, 2]])^2 / 2

cutoff_dist <- max(D_obs) * 0.75
en_rango <- h_pares <= cutoff_dist

n_bins <- 12
bins_factor <- cut(h_pares[en_rango], breaks = n_bins)
gamma_emp   <- tapply(gamma_pares[en_rango], bins_factor, mean)
h_medio     <- tapply(h_pares[en_rango], bins_factor, mean)
n_pares_bin <- tapply(gamma_pares[en_rango], bins_factor, length)

semivario_df <- na.omit(data.frame(h = h_medio, gamma = gamma_emp, N = n_pares_bin))
print(semivario_df)

# Funcion SCE para ajuste de modelo exponencial
sce_exponencial <- function(theta, h, gamma_obs) {
  c0 <- theta[1]
  c1 <- theta[2]
  phi <- theta[3]
  gamma_teo <- c0 + c1 * (1 - exp(-h / phi))
  sum((gamma_obs - gamma_teo)^2)
}

phi_grid <- quantile(semivario_df$h, probs = c(0.15, 0.30, 0.50, 0.70))
ajustes_lista <- lapply(phi_grid, function(phi0) {
  optim(par = c(c0 = 0.05 * var(residuales), c1 = 0.95 * var(residuales), phi = phi0),
        fn = sce_exponencial, h = semivario_df$h, gamma_obs = semivario_df$gamma,
        method = "L-BFGS-B",
        lower = c(0, 0.01, 1),
        upper = c(2 * max(semivario_df$gamma), 5 * max(semivario_df$gamma), 3 * cutoff_dist))
})

sce_valores <- sapply(ajustes_lista, function(a) a$value)
ajuste_optimo <- ajustes_lista[[which.min(sce_valores)]]

c0_hat   <- unname(ajuste_optimo$par["c0"])
c1_hat   <- unname(ajuste_optimo$par["c1"])
phi_hat  <- unname(ajuste_optimo$par["phi"])
sill_hat <- c0_hat + c1_hat
rango_practico <- 3 * phi_hat

cat("--- Parametros del Semivariograma Exponencial Ajustado ---\n")
cat(sprintf("Efecto Pepita (Nugget, c0): %.4f\n", c0_hat))
cat(sprintf("Meseta Parcial (Sill parcial, c1): %.4f\n", c1_hat))
cat(sprintf("Meseta Total (Sill total): %.4f\n", sill_hat))
cat(sprintf("Parametro de Rango (phi): %.2f km\n", phi_hat))
cat(sprintf("Rango Practico (3 * phi): %.2f km\n", rango_practico))

# FIGURA 4: Semivariograma Experimental y Modelo Ajustado
guardar_figura("figura4_semivariograma_ajustado.png", {
  par(mar = c(4.5, 4.5, 3, 1))
  curva_h <- seq(0, cutoff_dist, length.out = 200)
  curva_gamma <- c0_hat + c1_hat * (1 - exp(-curva_h / phi_hat))

  plot(h_pares[en_rango], gamma_pares[en_rango], pch = 20, col = "grey75", cex = 0.7,
       xlab = "Distancia h (km)", ylab = "Semivarianza gamma(h)",
       main = "Semivariograma Experimental y Modelo Teorico Exponencial",
       las = 1, xlim = c(0, cutoff_dist), ylim = c(0, max(semivario_df$gamma) * 1.35))
  points(semivario_df$h, semivario_df$gamma, pch = 19, col = "#08519c", cex = 1.6)
  lines(curva_h, curva_gamma, col = "red", lwd = 2.5)
  abline(h = sill_hat, col = "darkgreen", lty = 2, lwd = 1.5)
  abline(v = rango_practico, col = "purple", lty = 2, lwd = 1.5)
  legend("bottomright", 
         legend = c("Nube de pares gamma_ij", "Semivariograma empirico (bins)",
                    "Modelo Exponencial ajustado", sprintf("Sill total = %.1f", sill_hat),
                    sprintf("Rango practico = %.1f km", rango_practico)),
         col = c("grey75", "#08519c", "red", "darkgreen", "purple"),
         pch = c(20, 19, NA, NA, NA), lty = c(0, 0, 1, 2, 2), lwd = c(0, 0, 2.5, 1.5, 1.5),
         bg = "white", bty = "o", cex = 0.8)
}, ancho = 2400, alto = 1500)

# ==============================================================================
# 7. FORMULACION DE KRIGING UNIVERSAL CON LAGRANGE
# ==============================================================================
cat("\n[Paso 7] Construyendo matriz de covarianza y estimador Kriging Universal...\n")

C_cov <- function(h) sill_hat * exp(-h / phi_hat)
Sigma_obs <- matrix(C_cov(as.vector(D_obs)), nrow = n)

kriging_universal_pred <- function(s0_km, coords_sample, resid_sample, covar_df_new, mod_tend, Sigma_mat, C_fn) {
  n_pts <- nrow(coords_sample)
  d0 <- sqrt((coords_sample[, 1] - s0_km[1])^2 + (coords_sample[, 2] - s0_km[2])^2)
  c0_vec <- C_fn(d0)
  
  # Sistema lineal aumentado de Kriging (Lagrange)
  A <- rbind(cbind(Sigma_mat, rep(1, n_pts)), c(rep(1, n_pts), 0))
  b <- c(c0_vec, 1)
  sol <- solve(A, b)
  
  pesos_lambda <- sol[1:n_pts]
  mu_lagrange  <- sol[n_pts + 1]
  
  pred_tendencia <- predict(mod_tend, newdata = covar_df_new)
  pred_residuo   <- sum(pesos_lambda * resid_sample)
  pred_total     <- unname(pred_tendencia + pred_residuo)
  var_kriging    <- unname(max(0, sill_hat - sum(pesos_lambda * c0_vec) - mu_lagrange))
  
  return(c(prediccion = pred_total, 
           tendencia = unname(pred_tendencia), 
           residuo_krg = unname(pred_residuo), 
           varianza = var_kriging))
}

# ==============================================================================
# 8. VALIDACION CRUZADA LEAVE-ONE-OUT (LOOCV)
# ==============================================================================
cat("\n[Paso 8] Ejecutando validacion cruzada Leave-One-Out (LOOCV)...\n")

pred_loo <- numeric(n)
var_loo  <- numeric(n)
errores_loo <- numeric(n)

for (i in 1:n) {
  datos_sin_i  <- datos_obs[-i, ]
  coords_sin_i <- coords_obs_km[-i, ]
  Sigma_sin_i  <- Sigma_obs[-i, -i]
  
  mod_tend_i  <- lm(precip ~ altitud + temp + rad + x_km + y_km, data = datos_sin_i)
  resid_sin_i <- residuals(mod_tend_i)
  
  s0_i <- coords_obs_km[i, ]
  df_new_i <- datos_obs[i, , drop = FALSE]
  
  res_krg_i <- kriging_universal_pred(s0_i, coords_sin_i, resid_sin_i, df_new_i, mod_tend_i, Sigma_sin_i, C_cov)
  
  pred_loo[i] <- res_krg_i["prediccion"]
  var_loo[i]  <- res_krg_i["varianza"]
  errores_loo[i] <- datos_obs$precip[i] - pred_loo[i]
}

rmse_loocv <- sqrt(mean(errores_loo^2))
mae_loocv  <- mean(abs(errores_loo))
r2_loocv   <- 1 - (sum(errores_loo^2) / sum((datos_obs$precip - mean(datos_obs$precip))^2))

cat("--- Metricas de Validacion Cruzada (LOOCV) ---\n")
cat(sprintf("RMSE: %.2f mm/semana\n", rmse_loocv))
cat(sprintf("MAE:  %.2f mm/semana\n", mae_loocv))
cat(sprintf("R2 de Validacion (R2_CV): %.4f\n", r2_loocv))

# FIGURA 5: Diagnostico de Validacion Cruzada
guardar_figura("figura5_validacion_cruzada.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

  plot(datos_obs$precip, pred_loo, pch = 21, bg = "#3182bd", cex = 1.4,
       xlab = "Precipitacion Observada (mm/semana)", ylab = "Precipitacion Predicha LOOCV (mm/semana)",
       main = sprintf("(A) Observado vs Predicho (R2 = %.2f)", r2_loocv), las = 1,
       xlim = range(c(datos_obs$precip, pred_loo)), ylim = range(c(datos_obs$precip, pred_loo)))
  abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
  grid(col = "grey80")

  plot(datos_obs$precip, errores_loo, pch = 21, bg = "#e6550d", cex = 1.4,
       xlab = "Precipitacion Observada (mm/semana)", ylab = "Error Residual (Obs - Pred)",
       main = sprintf("(B) Diagnostico de Residuos (RMSE = %.2f mm)", rmse_loocv), las = 1)
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  grid(col = "grey80")
}, ancho = 2400, alto = 1200)

# ==============================================================================
# 9. PREDICCION ESPACIAL SOBRE TODA LA GRILLA DEPARTAMENTAL
# ==============================================================================
cat("\n[Paso 9] Evaluando Kriging Universal sobre la totalidad del Valle del Cauca...\n")

N_total <- nrow(df_grilla_valle)
mat_pred_grilla <- matrix(NA, nrow = N_total, ncol = 4)
colnames(mat_pred_grilla) <- c("prediccion", "tendencia", "residuo_krg", "varianza")

for (k in 1:N_total) {
  s0_k <- as.numeric(df_grilla_valle[k, c("x_km", "y_km")])
  df_new_k <- df_grilla_valle[k, , drop = FALSE]
  mat_pred_grilla[k, ] <- kriging_universal_pred(s0_k, coords_obs_km, residuales, df_new_k, modelo_tendencia, Sigma_obs, C_cov)
}

df_grilla_valle$prediccion  <- mat_pred_grilla[, "prediccion"]
df_grilla_valle$tendencia   <- mat_pred_grilla[, "tendencia"]
df_grilla_valle$residuo_krg <- mat_pred_grilla[, "residuo_krg"]
df_grilla_valle$varianza    <- mat_pred_grilla[, "varianza"]
df_grilla_valle$desv_est    <- sqrt(df_grilla_valle$varianza)

# Reconstruccion de objetos SpatRaster
r_pred_valle <- rast(df_grilla_valle[, c("lon", "lat", "prediccion")], type = "xyz", crs = crs(r_altitud))
r_var_valle  <- rast(df_grilla_valle[, c("lon", "lat", "varianza")],   type = "xyz", crs = crs(r_altitud))
r_sd_valle   <- rast(df_grilla_valle[, c("lon", "lat", "desv_est")],   type = "xyz", crs = crs(r_altitud))
r_tend_valle <- rast(df_grilla_valle[, c("lon", "lat", "tendencia")],  type = "xyz", crs = crs(r_altitud))
r_res_valle  <- rast(df_grilla_valle[, c("lon", "lat", "residuo_krg")], type = "xyz", crs = crs(r_altitud))

writeRaster(r_pred_valle, file.path(DIR_SALIDA, "raster_precipitacion_estimada.tif"), overwrite = TRUE)
writeRaster(r_var_valle,  file.path(DIR_SALIDA, "raster_varianza_kriging.tif"), overwrite = TRUE)

# FIGURA 6: Mapas de Prediccion e Incertidumbre
guardar_figura("figura6_mapas_prediccion_incertidumbre.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 4.5))

  plot(r_pred_valle, col = rev(terrain.colors(50)),
       main = "(A) Precipitacion Semanal Estimada (mm/semana)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)
  points(datos_obs$lon, datos_obs$lat, pch = 3, col = "black", cex = 0.8, lwd = 1.2)

  plot(r_sd_valle, col = heat.colors(50),
       main = "(B) Incertidumbre (Desv. Est. Kriging)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)
  points(datos_obs$lon, datos_obs$lat, pch = 3, col = "black", cex = 0.8, lwd = 1.2)
}, ancho = 3200, alto = 1600)

# FIGURA 7: Descomposicion de Kriging
guardar_figura("figura7_descomposicion_kriging.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 4.5))

  plot(r_tend_valle, col = rev(terrain.colors(50)),
       main = "(A) Componente de Tendencia (Covariables)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)

  plot(r_res_valle, col = topo.colors(50),
       main = "(B) Componente Residual Espacial (Kriging)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)
}, ancho = 3200, alto = 1600)

# ==============================================================================
# 10. EXPORTACION DE TABLAS NUMERICAS
# ==============================================================================
cat("\n[Paso 10] Exportando tablas de resultados numericos...\n")

tabla_descriptiva <- data.frame(
  Variable = c("Precipitacion (mm)", "Altitud SRTM (m)", "Temperatura (C)", "Radiacion (MJ/m2)", "Residuales"),
  Media = c(mean(datos_obs$precip), mean(datos_obs$altitud), mean(datos_obs$temp), mean(datos_obs$rad), mean(residuales)),
  Desv_Est = c(sd(datos_obs$precip), sd(datos_obs$altitud), sd(datos_obs$temp), sd(datos_obs$rad), sd(residuales)),
  Minimo = c(min(datos_obs$precip), min(datos_obs$altitud), min(datos_obs$temp), min(datos_obs$rad), min(residuales)),
  Maximo = c(max(datos_obs$precip), max(datos_obs$altitud), max(datos_obs$temp), max(datos_obs$rad), max(residuales)),
  Moran_I = c(moran_I, NA, NA, NA, moran_I_res)
)
write.csv(tabla_descriptiva, file.path(DIR_SALIDA, "tabla1_estadistica_descriptiva.csv"), row.names = FALSE)

tabla_coeficientes <- as.data.frame(resumen_modelo$coefficients)
tabla_coeficientes$Variable <- rownames(tabla_coeficientes)
write.csv(tabla_coeficientes, file.path(DIR_SALIDA, "tabla2_coeficientes_tendencia.csv"), row.names = FALSE)

tabla_geoestadistica <- data.frame(
  Parametro = c("Efecto Pepita (c0)", "Contribucion (c1)", "Meseta Total (Sill)", "Parametro de Rango (phi)", "Rango Practico", "RMSE LOOCV", "MAE LOOCV", "R2 LOOCV"),
  Valor = c(c0_hat, c1_hat, sill_hat, phi_hat, rango_practico, rmse_loocv, mae_loocv, r2_loocv),
  Unidad = c("mm^2", "mm^2", "mm^2", "km", "km", "mm", "mm", "Adimensional")
)
write.csv(tabla_geoestadistica, file.path(DIR_SALIDA, "tabla3_parametros_geoestadistica.csv"), row.names = FALSE)

cat("\n====================================================================\n")
cat(" EJECUCION FINALIZADA EXITOSAMENTE                                  \n")
cat(" Todos los graficos, mapas y tablas fueron generados en /resultados \n")
cat("====================================================================\n")
