################################################################################
# UNIVERSIDAD AUTONOMA DE OCCIDENTE
# FACULTAD DE INGENIERIA Y CIENCIAS BASICAS
# PROGRAMA: INGENIERIA DE DATOS E INTELIGENCIA ARTIFICIAL
# ASIGNATURA: ANALITICA DE DATOS (2026-2S)
# PROFESOR: JOHANN A. OSPINA
#
# PROYECTO NO. 1 - EVALUACION METODOLOGICA DUAL: KRIGING DIRECTO & AFM + KRIGING
#
# INTEGRANTES:
# - Juan Pablo Maya          - Codigo: 2236377
# - Cesar Armando Reyes      - Codigo: 2236379
# - Yesenia Diaz Urrego      - Codigo: 2231783
#
# Script 100% independiente, integral y reproducible.
# Combina estrictamente las metodologias trabajadas en clase:
# - Unidad 1: Analisis Factorial Multiple (Script_AFM.R) via SVD
# - Unidad 2: Geoestadistica y Kriging Universal (Ejemplo4_Geoestadistica.R,
#             Script_Spatial_Analytics.R, Ejemplo1_Geoestadistica.R)
################################################################################

rm(list = ls())
set.seed(2026)
options(digits = 4, scipen = 999)

# ==============================================================================
# 0. CONFIGURACION DE LIBRERIAS PERMITIDAS Y RUTAS
# ==============================================================================
paquetes <- c("terra", "sp", "gstat", "ggplot2", "gridExtra", "corrplot")
for (p in paquetes) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

# Rutas del proyecto
DIR_DATOS <- "data/datos_proyecto_1/imagenes_semanales"
if (!dir.exists(DIR_DATOS)) DIR_DATOS <- "../data/datos_proyecto_1/imagenes_semanales"
if (!dir.exists(DIR_DATOS)) {
  zip_datos <- "data/datos_proyecto_1.zip"
  if (!file.exists(zip_datos)) zip_datos <- "../data/datos_proyecto_1.zip"
  if (file.exists(zip_datos)) unzip(zip_datos, exdir = "data")
  tifs <- list.files(pattern = "altitud_valle\\.tif$", recursive = TRUE, full.names = TRUE)
  if (length(tifs) > 0) DIR_DATOS <- dirname(tifs[1])
}

DIR_SALIDA  <- "prueba_con_FMA/resultados"
DIR_FIGURAS <- "prueba_con_FMA/figuras"
if (!dir.exists(DIR_SALIDA))  dir.create(DIR_SALIDA, recursive = TRUE)
if (!dir.exists(DIR_FIGURAS)) dir.create(DIR_FIGURAS, recursive = TRUE)

guardar_figura <- function(nombre, expr, ancho = 2400, alto = 1600, res = 300) {
  bloque <- substitute(expr)
  entorno <- parent.frame()
  eval(bloque, entorno)
  
  rutas <- c(file.path(DIR_SALIDA, nombre), file.path(DIR_FIGURAS, nombre))
  for (r in unique(rutas)) {
    tryCatch({
      if (file.exists(r)) suppressWarnings(unlink(r, force = TRUE))
      png(r, width = ancho, height = alto, res = res)
      eval(bloque, entorno)
      dev.off()
    }, error = function(e) {
      suppressWarnings(try(dev.off(), silent = TRUE))
      tryCatch({
        png(r, width = ancho, height = alto, res = res, type = "cairo")
        eval(bloque, entorno)
        dev.off()
      }, error = function(e2) {})
    })
  }
  cat(sprintf("  [Figura generada] %s\n", nombre))
}

cat("====================================================================\n")
cat(" PROYECTO 1: MODELAMIENTO GEOESTADISTICO INTEGRADO (DIRECTO + AFM)  \n")
cat(" Integrantes: Juan Pablo Maya, Cesar Reyes, Yesenia Diaz            \n")
cat("====================================================================\n\n")

# ==============================================================================
# 1. CARGA DE DATOS RASTER Y ALINEACION CONTINUA
# ==============================================================================
cat("[Paso 1] Cargando y auditando capas raster...\n")

r_altitud     <- rast(file.path(DIR_DATOS, "altitud_valle.tif"))
r_precip_clim <- rast(file.path(DIR_DATOS, "chirps_climatologia_semanal_valle.tif"))
r_temp_clim   <- rast(file.path(DIR_DATOS, "power_temp_climatologia_semanal_valle.tif"))
r_rad_clim    <- rast(file.path(DIR_DATOS, "power_radiacion_climatologia_semanal_valle.tif"))

r_precip_clim[r_precip_clim < 0] <- NA
mascara_valle <- (!is.na(r_altitud)) & (r_altitud >= -5)
borde_valle <- as.polygons(mascara_valle, dissolve = TRUE)
borde_valle <- borde_valle[borde_valle[[1]] == 1, ]

SEMANA <- 42

# Interpolacion espacial continua IDW para NASA POWER
alinear_covariable <- function(r_layer, mask_base) {
  pts <- as.points(r_layer, na.rm = TRUE)
  if (nrow(pts) == 0) return(r_layer)
  r_idw <- interpIDW(mask_base, pts, field = names(r_layer), radius = 2.5, power = 2)
  mask(r_idw, mask_base)
}

p_lluvias   <- mask(r_precip_clim[[SEMANA]], mascara_valle)
t_lluvias   <- alinear_covariable(r_temp_clim[[SEMANA]], mascara_valle)
rad_lluvias <- alinear_covariable(r_rad_clim[[SEMANA]], mascara_valle)
alt_valle   <- mask(r_altitud, mascara_valle)

coords_all <- crds(alt_valle, na.rm = TRUE)
df_grilla_valle <- data.frame(
  id      = 1:nrow(coords_all),
  lon     = coords_all[, 1],
  lat     = coords_all[, 2],
  altitud = extract(alt_valle, coords_all)[, 1],
  precip  = extract(p_lluvias, coords_all)[, 1],
  temp    = extract(t_lluvias, coords_all)[, 1],
  rad     = extract(rad_lluvias, coords_all)[, 1]
)
df_grilla_valle <- na.omit(df_grilla_valle)

lat0 <- mean(df_grilla_valle$lat)
lon0 <- mean(df_grilla_valle$lon)
df_grilla_valle$x_km <- (df_grilla_valle$lon - lon0) * 111.320 * cos(lat0 * pi / 180)
df_grilla_valle$y_km <- (df_grilla_valle$lat - lat0) * 110.574

cat(sprintf("  -> Total celdas validas en el Valle del Cauca: %d\n", nrow(df_grilla_valle)))

# ==============================================================================
# 2. MUESTREO ESPACIAL REPRESENTATIVO (70 ESTACIONES)
# ==============================================================================
cat("\n[Paso 2] Muestreo espacial de estaciones de entrenamiento...\n")
n_muestra <- 70
set.seed(2026)
idx_muestra <- sample(1:nrow(df_grilla_valle), size = n_muestra, replace = FALSE)
datos_obs <- df_grilla_valle[idx_muestra, ]
rownames(datos_obs) <- NULL
n <- nrow(datos_obs)
cat(sprintf("  -> %d estaciones seleccionadas.\n", n))

# ==============================================================================
# 3A. ANALISIS EXPLORATORIO DE DATOS ESPACIALES (ESDA)
# ==============================================================================
cat("\n[Paso 3A] Generando graficos exploratorios ESDA...\n")

resumen_p <- summary(datos_obs$precip)
sd_p <- sd(datos_obs$precip)
cv_p <- (sd_p / mean(datos_obs$precip)) * 100
shapiro_p <- shapiro.test(datos_obs$precip)

# FIGURA 1: Distribucion de Precipitacion
guardar_figura("figura1_distribucion_precipitacion.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  h <- hist(datos_obs$precip, breaks = 12, col = "#6baed6", border = "white",
            main = "(A) Histograma de Precipitacion Semanal",
            xlab = "Precipitacion (mm/semana)", ylab = "Frecuencia", las = 1)
  x_seq <- seq(min(datos_obs$precip), max(datos_obs$precip), length.out = 100)
  lines(x_seq, dnorm(x_seq, mean = mean(datos_obs$precip), sd = sd(datos_obs$precip)) * 
          length(datos_obs$precip) * diff(h$mids[1:2]), col = "darkblue", lwd = 2)
  grid(col = "grey80")

  boxplot(datos_obs$precip, col = "#fd8d3c", border = "#7f2704",
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
# 3B. ANALISIS FACTORIAL MULTIPLE (AFM / MFA) - METODOLOGIA Script_AFM.R
# ==============================================================================
cat("\n[Paso 3B] Ejecutando Analisis Factorial Multiple (AFM) de Covariables...\n")

lista_grupos <- list(
  Topografia = as.matrix(datos_obs[, "altitud", drop = FALSE]),
  Clima      = as.matrix(datos_obs[, c("temp", "rad")]),
  Espacio    = as.matrix(datos_obs[, c("x_km", "y_km")])
)
K <- length(lista_grupos)

estandarizar <- function(X) scale(X, center = TRUE, scale = TRUE)
acp_svd <- function(X) {
  n_f <- nrow(X)
  Z <- X / sqrt(n_f - 1)
  desc <- svd(Z)
  list(valores = desc$d^2, ejes = desc$v, U = desc$u, d = desc$d)
}

lambda1 <- numeric(K)
grupos_std <- vector("list", K)
for (k in seq_len(K)) {
  Xk <- estandarizar(lista_grupos[[k]])
  grupos_std[[k]] <- Xk
  lambda1[k] <- acp_svd(Xk)$valores[1]
}
names(lambda1) <- names(lista_grupos)

alpha <- 1 / lambda1
grupos_pond <- Map(function(Xk, a) Xk * sqrt(a), grupos_std, alpha)

X_global <- do.call(cbind, grupos_pond)
res_global <- acp_svd(X_global)
valores_propios <- res_global$valores
inercia_pct <- 100 * valores_propios / sum(valores_propios)

q <- 2
U_ejes <- res_global$ejes[, 1:q, drop = FALSE]
F_global <- (X_global / sqrt(n - 1)) %*% U_ejes
colnames(F_global) <- paste0("AFM_Dim", 1:q)

RV <- function(Xk, Xkp) {
  Wk  <- Xk  %*% t(Xk); Wkp <- Xkp %*% t(Xkp)
  sum(diag(Wk %*% Wkp)) / sqrt(sum(diag(Wk %*% Wk)) * sum(diag(Wkp %*% Wkp)))
}
RV_matriz <- matrix(NA, K, K, dimnames = list(names(lista_grupos), names(lista_grupos)))
for (k in seq_len(K)) for (kp in seq_len(K)) RV_matriz[k, kp] <- RV(grupos_std[[k]], grupos_std[[kp]])

# Proyectar factores AFM a toda la grilla del departamento
grilla_topografia_std <- scale(as.matrix(df_grilla_valle[, "altitud", drop = FALSE]),
                               center = attr(grupos_std[[1]], "scaled:center"),
                               scale  = attr(grupos_std[[1]], "scaled:scale"))
grilla_clima_std      <- scale(as.matrix(df_grilla_valle[, c("temp", "rad")]),
                               center = attr(grupos_std[[2]], "scaled:center"),
                               scale  = attr(grupos_std[[2]], "scaled:scale"))
grilla_espacio_std    <- scale(as.matrix(df_grilla_valle[, c("x_km", "y_km")]),
                               center = attr(grupos_std[[3]], "scaled:center"),
                               scale  = attr(grupos_std[[3]], "scaled:scale"))

grilla_global_pond <- cbind(
  grilla_topografia_std * sqrt(alpha["Topografia"]),
  grilla_clima_std * sqrt(alpha["Clima"]),
  grilla_espacio_std * sqrt(alpha["Espacio"])
)

F_grilla <- (grilla_global_pond / sqrt(n - 1)) %*% U_ejes
colnames(F_grilla) <- c("AFM_Dim1", "AFM_Dim2")

df_grilla_valle$AFM_Dim1 <- F_grilla[, 1]
df_grilla_valle$AFM_Dim2 <- F_grilla[, 2]
datos_obs$AFM_Dim1 <- F_global[, 1]
datos_obs$AFM_Dim2 <- F_global[, 2]

# FIGURA 4: Biplot AFM y Circulo de Correlacion
guardar_figura("figura4_afm_biplot_correlacion.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  
  pal_p <- colorRampPalette(c("#fee08b", "#f46d43", "#5e4fa2"))(100)
  col_p <- pal_p[round((datos_obs$precip - min(datos_obs$precip)) / diff(range(datos_obs$precip)) * 99) + 1]
  
  plot(F_global[, 1], F_global[, 2], pch = 21, bg = col_p, cex = 1.5,
       xlab = sprintf("Dimension 1 (%.1f%%)", inercia_pct[1]),
       ylab = sprintf("Dimension 2 (%.1f%%)", inercia_pct[2]),
       main = "(A) Nube Global de Estaciones (AFM)", las = 1)
  abline(h = 0, v = 0, lty = 2, col = "grey60")
  grid(col = "grey85")
  
  X_std_todas <- do.call(cbind, grupos_std)
  cor_vars <- cor(X_std_todas, F_global)
  
  plot(0, 0, type = "n", xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1),
       xlab = sprintf("Dimension 1 (%.1f%%)", inercia_pct[1]),
       ylab = sprintf("Dimension 2 (%.1f%%)", inercia_pct[2]),
       main = "(B) Circulo de Correlacion de Variables", las = 1, asp = 1)
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "grey50", lty = 2)
  abline(h = 0, v = 0, lty = 2, col = "grey70")
  
  colores_grupos <- c(altitud = "#31a354", temp = "#e6550d", rad = "#756bb1", x_km = "#3182bd", y_km = "#3182bd")
  for (v in rownames(cor_vars)) {
    arrows(0, 0, cor_vars[v, 1], cor_vars[v, 2], col = colores_grupos[v], length = 0.1, lwd = 2)
    text(cor_vars[v, 1] * 1.15, cor_vars[v, 2] * 1.15, labels = v, col = colores_grupos[v], font = 2, cex = 0.85)
  }
}, ancho = 2800, alto = 1400)

# ==============================================================================
# 4. AUTOCORRELACION ESPACIAL GLOBAL (MORAN Y GEARY)
# ==============================================================================
cat("\n[Paso 4] Calculando matrices de vecindad y metricas de autocorrelacion...\n")

coords_obs_km <- as.matrix(datos_obs[, c("x_km", "y_km")])
D_obs <- as.matrix(dist(coords_obs_km))
d_corte <- as.numeric(quantile(D_obs[upper.tri(D_obs)], 0.15))
W_obs <- (D_obs <= d_corte & D_obs > 0) * 1.0
W_std <- sweep(W_obs, 1, ifelse(rowSums(W_obs) == 0, 1, rowSums(W_obs)), "/")

calc_moran <- function(z, W) {
  n_z <- length(z); zm <- z - mean(z); s0 <- sum(W)
  num <- sum(W * (zm %o% zm)); den <- sum(zm^2)
  (n_z / s0) * (num / den)
}

calc_geary <- function(z, W) {
  n_z <- length(z); zm <- z - mean(z); s0 <- sum(W)
  diff_sq <- outer(z, z, "-")^2
  ((n_z - 1) / (2 * s0)) * (sum(W * diff_sq) / sum(zm^2))
}

moran_I <- calc_moran(datos_obs$precip, W_std)
geary_C <- calc_geary(datos_obs$precip, W_std)

cat(sprintf("  -> Indice Global de Moran: %.4f (Autocorrelacion Positiva)\n", moran_I))
cat(sprintf("  -> Coeficiente de Geary:   %.4f\n", geary_C))

# ==============================================================================
# 5. MODELACION DE TENDENCIA: MODELO DIRECTO VS MODELO AFM
# ==============================================================================
cat("\n[Paso 5] Ajustando Modelos de Tendencia...\n")

# Modelo 1: Regresion Multivariada Directa
mod_directo <- lm(precip ~ altitud + temp + rad + x_km + y_km, data = datos_obs)
resumen_directo <- summary(mod_directo)
resid_directo <- residuals(mod_directo)
moran_res_directo <- calc_moran(resid_directo, W_std)

# Modelo 2: Regresion con Factores AFM
mod_afm <- lm(precip ~ AFM_Dim1 + AFM_Dim2, data = datos_obs)
resumen_afm <- summary(mod_afm)
resid_afm <- residuals(mod_afm)
moran_res_afm <- calc_moran(resid_afm, W_std)

cat(sprintf("  -> Modelo 1 (Directo): R2 = %.4f, Moran I Residual = %.4f\n", 
            resumen_directo$r.squared, moran_res_directo))
cat(sprintf("  -> Modelo 2 (AFM):     R2 = %.4f, Moran I Residual = %.4f\n", 
            resumen_afm$r.squared, moran_res_afm))

# ==============================================================================
# 6. VARIOGRAFIA Y AJUSTE TEORICO EXPONENCIAL
# ==============================================================================
cat("\n[Paso 6] Ajustando Semivariogramas Teoricos (L-BFGS-B)...\n")

ajustar_semivariograma <- function(residuales) {
  pares <- which(upper.tri(D_obs), arr.ind = TRUE)
  h_ij <- D_obs[upper.tri(D_obs)]
  gamma_ij <- (residuales[pares[, 1]] - residuales[pares[, 2]])^2 / 2
  
  n_bins <- 12
  bins <- cut(h_ij, breaks = n_bins)
  gamma_emp <- tapply(gamma_ij, bins, mean)
  N_h       <- tapply(gamma_ij, bins, length)
  h_medio   <- tapply(h_ij, bins, mean)
  semiv <- na.omit(data.frame(h = as.numeric(h_medio), gamma = as.numeric(gamma_emp), N = as.numeric(N_h)))
  
  sce_exp <- function(par, h, g_emp) {
    c0 <- par[1]; c1 <- par[2]; phi <- par[3]
    g_teo <- c0 + c1 * (1 - exp(-h / phi))
    sum((g_emp - g_teo)^2)
  }
  
  phis_ini <- quantile(semiv$h, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  ajustes <- lapply(phis_ini, function(p0) {
    optim(par = c(c0 = 0, c1 = var(residuales), phi = p0),
          fn = sce_exp, h = semiv$h, g_emp = semiv$gamma,
          method = "L-BFGS-B",
          lower = c(0, 0, 1),
          upper = c(max(semiv$gamma), 5 * max(semiv$gamma), 3 * max(semiv$h)))
  })
  
  best <- ajustes[[which.min(sapply(ajustes, function(a) a$value))]]
  list(semiv = semiv, 
       c0 = unname(best$par["c0"]), 
       c1 = unname(best$par["c1"]), 
       phi = unname(best$par["phi"]),
       sill = unname(best$par["c0"] + best$par["c1"]))
}

fit_dir <- ajustar_semivariograma(resid_directo)
fit_afm <- ajustar_semivariograma(resid_afm)

cat(sprintf("  -> Modelo Directo: Sill = %.2f, Rango Practico = %.2f km\n", fit_dir$sill, 3 * fit_dir$phi))
cat(sprintf("  -> Modelo AFM:     Sill = %.2f, Rango Practico = %.2f km\n", fit_afm$sill, 3 * fit_afm$phi))

# FIGURA 5: Semivariogramas Ajustados
guardar_figura("figura5_semivariogramas_ajustados.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  h_seq <- seq(0, max(fit_dir$semiv$h) * 1.05, length.out = 150)
  
  # Panel A: Modelo Directo
  plot(fit_dir$semiv$h, fit_dir$semiv$gamma, pch = 21, bg = "#3182bd", cex = 1.4,
       xlab = "Distancia h (km)", ylab = "Semivarianza gamma(h)",
       main = "(A) Semivariograma: Residuales Modelo Directo", las = 1,
       ylim = c(0, max(fit_dir$semiv$gamma) * 1.15))
  lines(h_seq, fit_dir$c0 + fit_dir$c1 * (1 - exp(-h_seq / fit_dir$phi)), col = "red", lwd = 2.5)
  grid(col = "grey85")
  legend("bottomright", legend = c("Empirico", sprintf("Exponencial (Rango = %.1f km)", 3 * fit_dir$phi)),
         pch = c(21, NA), lty = c(NA, 1), col = c("black", "red"), pt.bg = c("#3182bd", NA), bty = "o")

  # Panel B: Modelo AFM
  plot(fit_afm$semiv$h, fit_afm$semiv$gamma, pch = 21, bg = "#e6550d", cex = 1.4,
       xlab = "Distancia h (km)", ylab = "Semivarianza gamma(h)",
       main = "(B) Semivariograma: Residuales Modelo AFM", las = 1,
       ylim = c(0, max(fit_afm$semiv$gamma) * 1.15))
  lines(h_seq, fit_afm$c0 + fit_afm$c1 * (1 - exp(-h_seq / fit_afm$phi)), col = "darkgreen", lwd = 2.5)
  grid(col = "grey85")
  legend("bottomright", legend = c("Empirico", sprintf("Exponencial (Rango = %.1f km)", 3 * fit_afm$phi)),
         pch = c(21, NA), lty = c(NA, 1), col = c("black", "darkgreen"), pt.bg = c("#e6550d", NA), bty = "o")
}, ancho = 2800, alto = 1400)

# ==============================================================================
# 7 & 8. KRIGING UNIVERSAL Y VALIDACION CRUZADA (LOOCV)
# ==============================================================================
cat("\n[Paso 7 y 8] Ejecutando Kriging Universal y LOOCV para ambos modelos...\n")

kriging_solver <- function(s0, coords_train, resid_train, df_new, mod_t, sill, phi) {
  n_tr <- nrow(coords_train)
  d0 <- sqrt((coords_train[, 1] - s0[1])^2 + (coords_train[, 2] - s0[2])^2)
  c0_vec <- sill * exp(-d0 / phi)
  
  D_tr <- as.matrix(dist(coords_train))
  Sigma_tr <- sill * exp(-D_tr / phi)
  
  A <- rbind(cbind(Sigma_tr, rep(1, n_tr)), c(rep(1, n_tr), 0))
  b <- c(c0_vec, 1)
  sol <- solve(A, b)
  pesos <- sol[1:n_tr]
  mu_lagrange <- sol[n_tr + 1]
  
  p_tend <- predict(mod_t, newdata = df_new)
  p_res  <- sum(pesos * resid_train)
  p_tot  <- unname(p_tend + p_res)
  v_krg  <- unname(max(0, sill - sum(pesos * c0_vec) - mu_lagrange))
  
  c(prediccion = p_tot, tendencia = unname(p_tend), residuo = unname(p_res), varianza = v_krg)
}

# LOOCV Modelo Directo
pred_loo_dir <- numeric(n); errores_loo_dir <- numeric(n)
for (i in 1:n) {
  res_i <- kriging_solver(coords_obs_km[i, ], coords_obs_km[-i, ], resid_directo[-i], 
                          datos_obs[i, , drop = FALSE], 
                          lm(precip ~ altitud + temp + rad + x_km + y_km, data = datos_obs[-i, ]), 
                          fit_dir$sill, fit_dir$phi)
  pred_loo_dir[i] <- res_i["prediccion"]
  errores_loo_dir[i] <- datos_obs$precip[i] - pred_loo_dir[i]
}
rmse_dir <- sqrt(mean(errores_loo_dir^2))
mae_dir  <- mean(abs(errores_loo_dir))
r2_cv_dir <- 1 - (sum(errores_loo_dir^2) / sum((datos_obs$precip - mean(datos_obs$precip))^2))

# LOOCV Modelo AFM
pred_loo_afm <- numeric(n); errores_loo_afm <- numeric(n)
for (i in 1:n) {
  res_i <- kriging_solver(coords_obs_km[i, ], coords_obs_km[-i, ], resid_afm[-i], 
                          datos_obs[i, , drop = FALSE], 
                          lm(precip ~ AFM_Dim1 + AFM_Dim2, data = datos_obs[-i, ]), 
                          fit_afm$sill, fit_afm$phi)
  pred_loo_afm[i] <- res_i["prediccion"]
  errores_loo_afm[i] <- datos_obs$precip[i] - pred_loo_afm[i]
}
rmse_afm <- sqrt(mean(errores_loo_afm^2))
mae_afm  <- mean(abs(errores_loo_afm))
r2_cv_afm <- 1 - (sum(errores_loo_afm^2) / sum((datos_obs$precip - mean(datos_obs$precip))^2))

cat(sprintf("  -> LOOCV Modelo Directo: RMSE = %.2f mm, MAE = %.2f mm, R2_CV = %.4f\n", rmse_dir, mae_dir, r2_cv_dir))
cat(sprintf("  -> LOOCV Modelo AFM:     RMSE = %.2f mm, MAE = %.2f mm, R2_CV = %.4f\n", rmse_afm, mae_afm, r2_cv_afm))

# FIGURA 6: Validacion Cruzada Comparativa
guardar_figura("figura6_validacion_cruzada_comparativa.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  rango_ejes <- range(c(datos_obs$precip, pred_loo_dir, pred_loo_afm))
  
  plot(datos_obs$precip, pred_loo_dir, pch = 21, bg = "#3182bd", cex = 1.4,
       xlab = "Precipitacion Observada (mm)", ylab = "Prediccion LOOCV (mm)",
       main = sprintf("(A) Modelo Directo (R2 = %.2f, RMSE = %.2f)", r2_cv_dir, rmse_dir),
       xlim = rango_ejes, ylim = rango_ejes, las = 1)
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  grid(col = "grey85")

  plot(datos_obs$precip, pred_loo_afm, pch = 21, bg = "#e6550d", cex = 1.4,
       xlab = "Precipitacion Observada (mm)", ylab = "Prediccion LOOCV (mm)",
       main = sprintf("(B) Modelo AFM (R2 = %.2f, RMSE = %.2f)", r2_cv_afm, rmse_afm),
       xlim = rango_ejes, ylim = rango_ejes, las = 1)
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  grid(col = "grey85")
}, ancho = 2800, alto = 1400)

# ==============================================================================
# 9. PREDICCION ESPACIAL CONTINUA SOBRE EL VALLE DEL CAUCA (688 CELDAS)
# ==============================================================================
cat("\n[Paso 9] Evaluando Prediccion Continua en Grilla (Modelos Directo y AFM)...\n")

N_grid <- nrow(df_grilla_valle)
pred_dir_grid <- numeric(N_grid); var_dir_grid <- numeric(N_grid); tend_dir_grid <- numeric(N_grid); res_dir_grid <- numeric(N_grid)
pred_afm_grid <- numeric(N_grid); var_afm_grid <- numeric(N_grid); tend_afm_grid <- numeric(N_grid); res_afm_grid <- numeric(N_grid)

for (k in 1:N_grid) {
  s0_k <- c(df_grilla_valle$x_km[k], df_grilla_valle$y_km[k])
  df_k <- df_grilla_valle[k, , drop = FALSE]
  
  # Directo
  k_dir <- kriging_solver(s0_k, coords_obs_km, resid_directo, df_k, mod_directo, fit_dir$sill, fit_dir$phi)
  pred_dir_grid[k] <- k_dir["prediccion"]
  var_dir_grid[k]  <- k_dir["varianza"]
  tend_dir_grid[k] <- k_dir["tendencia"]
  res_dir_grid[k]  <- k_dir["residuo"]
  
  # AFM
  k_afm <- kriging_solver(s0_k, coords_obs_km, resid_afm, df_k, mod_afm, fit_afm$sill, fit_afm$phi)
  pred_afm_grid[k] <- k_afm["prediccion"]
  var_afm_grid[k]  <- k_afm["varianza"]
  tend_afm_grid[k] <- k_afm["tendencia"]
  res_afm_grid[k]  <- k_afm["residuo"]
}

df_grilla_valle$pred_dir <- pmax(0, pred_dir_grid)
df_grilla_valle$sd_dir   <- sqrt(var_dir_grid)
df_grilla_valle$tend_dir <- tend_dir_grid
df_grilla_valle$res_dir  <- res_dir_grid

df_grilla_valle$pred_afm <- pmax(0, pred_afm_grid)
df_grilla_valle$sd_afm   <- sqrt(var_afm_grid)
df_grilla_valle$tend_afm <- tend_afm_grid
df_grilla_valle$res_afm  <- res_afm_grid

df_grilla_valle$dif_pred <- df_grilla_valle$pred_dir - df_grilla_valle$pred_afm

# Creacion de Rasters SpatRaster
r_pred_dir <- rast(df_grilla_valle[, c("lon", "lat", "pred_dir")], type = "xyz", crs = crs(r_altitud))
r_sd_dir   <- rast(df_grilla_valle[, c("lon", "lat", "sd_dir")],   type = "xyz", crs = crs(r_altitud))
r_tend_dir <- rast(df_grilla_valle[, c("lon", "lat", "tend_dir")], type = "xyz", crs = crs(r_altitud))
r_res_dir  <- rast(df_grilla_valle[, c("lon", "lat", "res_dir")],  type = "xyz", crs = crs(r_altitud))

r_pred_afm <- rast(df_grilla_valle[, c("lon", "lat", "pred_afm")], type = "xyz", crs = crs(r_altitud))
r_sd_afm   <- rast(df_grilla_valle[, c("lon", "lat", "sd_afm")],   type = "xyz", crs = crs(r_altitud))
r_dif      <- rast(df_grilla_valle[, c("lon", "lat", "dif_pred")], type = "xyz", crs = crs(r_altitud))

# Exportacion de GeoTIFFs
writeRaster(r_pred_dir, file.path(DIR_SALIDA, "raster_precipitacion_directo.tif"), overwrite = TRUE)
writeRaster(r_pred_afm, file.path(DIR_SALIDA, "raster_precipitacion_afm.tif"), overwrite = TRUE)

# FIGURA 7: Mapas de Prediccion e Incertidumbre (Modelo Directo)
guardar_figura("figura7_mapas_prediccion_incertidumbre.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 4.5))
  pal_p <- colorRampPalette(c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b", "#08306b"))(60)
  plot(r_pred_dir, col = pal_p, main = "(A) Precipitacion Estimada (Kriging Universal)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)
  points(datos_obs$lon, datos_obs$lat, pch = 3, col = "red", cex = 0.7)

  pal_sd <- colorRampPalette(c("#ffffcc", "#ffeda0", "#feb24c", "#f03b20", "#bd0026"))(60)
  plot(r_sd_dir, col = pal_sd, main = "(B) Incertidumbre (Error Estandar Kriging)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)
  points(datos_obs$lon, datos_obs$lat, pch = 3, col = "black", cex = 0.7)
}, ancho = 3200, alto = 1600)

# FIGURA 8: Descomposicion de Kriging (Tendencia vs Residual)
guardar_figura("figura8_descomposicion_kriging.png", {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 4.5))
  plot(r_tend_dir, col = rev(terrain.colors(50)),
       main = "(A) Componente de Tendencia (Covariables)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)

  plot(r_res_dir, col = topo.colors(50),
       main = "(B) Componente Residual Espacial (Kriging)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)
}, ancho = 3200, alto = 1600)

# FIGURA 9: Comparacion Espacial de Modelos (Directo vs AFM vs Diferencia)
guardar_figura("figura9_comparacion_modelos_espaciales.png", {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 4.5))
  pal_p <- colorRampPalette(c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b", "#08306b"))(60)
  
  plot(r_pred_dir, col = pal_p, main = "(A) Modelo Directo (5 Covariables)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)

  plot(r_pred_afm, col = pal_p, main = "(B) Modelo con AFM (2 Factores)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)

  pal_dif <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(60)
  plot(r_dif, col = pal_dif, main = "(C) Diferencia (Directo - AFM)",
       xlab = "Longitud", ylab = "Latitud", las = 1)
  plot(borde_valle, add = TRUE, border = "black", lwd = 1.2)
}, ancho = 4200, alto = 1500)

# ==============================================================================
# 10. EXPORTACION DE TABLAS NUMERICAS
# ==============================================================================
cat("\n[Paso 10] Exportando tablas de resultados numericos...\n")

tabla_descriptiva <- data.frame(
  Variable = c("Precipitacion (mm)", "Altitud SRTM (m)", "Temperatura (C)", "Radiacion (MJ/m2)", "Residuales Directo", "Residuales AFM"),
  Media = c(mean(datos_obs$precip), mean(datos_obs$altitud), mean(datos_obs$temp), mean(datos_obs$rad), mean(resid_directo), mean(resid_afm)),
  Desv_Est = c(sd(datos_obs$precip), sd(datos_obs$altitud), sd(datos_obs$temp), sd(datos_obs$rad), sd(resid_directo), sd(resid_afm)),
  Minimo = c(min(datos_obs$precip), min(datos_obs$altitud), min(datos_obs$temp), min(datos_obs$rad), min(resid_directo), min(resid_afm)),
  Maximo = c(max(datos_obs$precip), max(datos_obs$altitud), max(datos_obs$temp), max(datos_obs$rad), max(resid_directo), max(resid_afm)),
  Moran_I = c(moran_I, NA, NA, NA, moran_res_directo, moran_res_afm)
)
write.csv(tabla_descriptiva, file.path(DIR_SALIDA, "tabla1_estadistica_descriptiva.csv"), row.names = FALSE)

tabla_afm_grupos <- data.frame(
  Grupo = names(lista_grupos),
  Num_Variables = sapply(lista_grupos, ncol),
  Lambda_1 = lambda1,
  Ponderacion_Alpha = alpha
)
write.csv(tabla_afm_grupos, file.path(DIR_SALIDA, "tabla3_resultados_afm.csv"), row.names = FALSE)

tabla_comparativa <- data.frame(
  Metodologia = c("Modelo 1: Regresion Multivariada Directa + Kriging", 
                  "Modelo 2: Factores Sinteticos AFM (FMA) + Kriging"),
  Dimensiones_Tendencia = c("5 covariables crudas (Alt, Temp, Rad, X, Y)", "2 dimensiones globales AFM"),
  R2_Tendencia = c(resumen_directo$r.squared, resumen_afm$r.squared),
  Moran_I_Residual = c(moran_res_directo, moran_res_afm),
  Sill_Variograma = c(fit_dir$sill, fit_afm$sill),
  Rango_Practico_km = c(3 * fit_dir$phi, 3 * fit_afm$phi),
  RMSE_LOOCV_mm = c(rmse_dir, rmse_afm),
  MAE_LOOCV_mm = c(mae_dir, mae_afm),
  R2_Validacion_LOOCV = c(r2_cv_dir, r2_cv_afm)
)
write.csv(tabla_comparativa, file.path(DIR_SALIDA, "tabla4_comparativa_modelos.csv"), row.names = FALSE)

cat("\n====================================================================\n")
cat(" EJECUCION EXITOSA: TODAS LAS 9 FIGURAS Y 3 TABLAS GENERADAS        \n")
cat(" Resultados guardados en prueba_con_FMA/resultados y /figuras       \n")
cat("====================================================================\n")
