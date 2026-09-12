###################################################################
# Ejemplo 3: EDA espacial básico y Geoestadística
# Analítica de Datos
# Johann A. Ospina
# Universidad Autónoma de Occidente
# 26 de agosto de 2026
# Nota: las funciones están con la teoría vista en clase son correctas, sin embargo,
# no quiere decir que los resultados sean apropiados. Usted debe leer y ajustar
# lo necesario para que los resultados sean apropiados
###################################################################


#Instalar librerias
install.packages("geodata")
install.packages("chirps")
install.packages("terra")

#Cargar librerias
library(geodata)
library(chirps)
library(terra)


#Fijar periodo
fecha_fin <- Sys.Date()
fecha_inicio_5anios <- fecha_fin - 5 * 365          
n_semanas_demo <- 8                                 # semanas de prueba
fecha_inicio_demo <- fecha_fin - n_semanas_demo * 7

print(c(inicio_5_anios = as.character(fecha_inicio_5anios), fin = as.character(fecha_fin)))
print(c(inicio_demo = as.character(fecha_inicio_demo), fin = as.character(fecha_fin)))


#Borde del valle del cauca

borde_ok <- FALSE
if (paquetes_ok) {
  borde_ok <- tryCatch({
    library(geodata)
    departamentos_col <- gadm(country = "COL", level = 1, path = tempdir())  # departamentos de Colombia
    valle <<- departamentos_col[departamentos_col$NAME_1 == "Valle del Cauca", ]  # se fija el borde real
    TRUE
  }, error = function(e) FALSE)
}
print(borde_ok)

if (!borde_ok) {
  # Borde aproximado de respaldo (poligono simplificado del Valle del Cauca)
  vertices_valle <- matrix(c(
    -77.15, 4.85, -76.60, 4.85, -76.10, 4.30, -75.90, 3.75,
    -76.10, 3.05, -76.75, 3.05, -77.15, 3.90, -77.15, 4.85
  ), ncol = 2, byrow = TRUE)
  valle <- vect(vertices_valle, type = "polygons", crs = "EPSG:4326")
}
plot(valle, main = "Borde del Valle del Cauca")


#Descargar precipitación diaria de satélite (CHIRPS)

precip_ok <- FALSE
if (paquetes_ok) {
  precip_ok <- tryCatch({
    library(chirps)
    precip_diaria <- get_chirps(valle, dates = c(as.character(fecha_inicio_demo), as.character(fecha_fin)),
                                server = "CHC", as.raster = TRUE)     # SpatRaster: 1 capa por dia
    
    fechas_dias <- seq(fecha_inicio_demo, fecha_fin, by = "day")[1:nlyr(precip_diaria)]
    indice_semana <- (as.numeric(fechas_dias - fecha_inicio_demo) %/% 7) + 1
    precip_semanal <<- tapp(precip_diaria, index = indice_semana, fun = sum)  # suma de 7 dias por capa
    names(precip_semanal) <<- paste0("semana_", 1:nlyr(precip_semanal))
    TRUE
  }, error = function(e) FALSE)
}
print(precip_ok)

if (!precip_ok) {
  # Simulacion de respaldo: pila de imagenes semanales con dependencia espacial,
  # sobre la misma extension del borde fijado (misma estructura que CHIRPS real)
  n_semanas_sim <- n_semanas_demo
  extension_valle <- ext(valle)
  base_grid <- rast(extension_valle, nrows = 40, ncols = 40, crs = "EPSG:4326")
  
  set.seed(2024)
  celdas <- xyFromCell(base_grid, 1:ncell(base_grid))
  D_celdas <- as.matrix(dist(celdas))
  Sigma_celdas <- 25 * exp(-D_celdas / 0.4)
  L_celdas <- chol(Sigma_celdas + diag(1e-6, nrow(celdas)))
  
  capas <- vector("list", n_semanas_sim)
  for (s in 1:n_semanas_sim) {
    tendencia <- 90 - 15 * (celdas[, 2] - mean(celdas[, 2]))          # mas lluvia hacia el sur, ejemplo
    valores <- pmax(tendencia + as.vector(t(L_celdas) %*% rnorm(nrow(celdas))), 0)
    capa_s <- base_grid
    values(capa_s) <- valores
    capas[[s]] <- capa_s
  }
  precip_semanal <- rast(capas)
  names(precip_semanal) <- paste0("semana_", 1:n_semanas_sim)
}

precip_semanal <- mask(precip_semanal, valle)       # se recorta exactamente al borde del departamento
plot(precip_semanal[[1]], main = "Precipitacion semanal (semana 1) sobre el Valle del Cauca")


# Muestreo de puntos aleatorios

set.seed(2026)
puntos_valle <- spatSample(valle, size = 40, method = "random")   # 40 puntos dentro del poligono real
extraido <- extract(precip_semanal, puntos_valle)                
coords_xy_grados <- crds(puntos_valle)                           

datos_precip <- data.frame(
  ID = extraido$ID,
  lon = coords_xy_grados[, 1],
  lat = coords_xy_grados[, 2],
  precip_media = rowMeans(extraido[, -1, drop = FALSE], na.rm = TRUE)  # promedio semanal en la ventana
)
print(head(datos_precip))


#Proyección en km

lat0 <- mean(datos_precip$lat)
datos_precip$x_km <- (datos_precip$lon - mean(datos_precip$lon)) * 111.320 * cos(lat0 * pi / 180)
datos_precip$y_km <- (datos_precip$lat - lat0) * 110.574
print(head(datos_precip[, c("lat", "lon", "x_km", "y_km")]))


#EDA espacial
print(summary(datos_precip$precip_media))
print(sd(datos_precip$precip_media) / mean(datos_precip$precip_media))  # coeficiente de variacion

hist(datos_precip$precip_media, breaks = 10, col = "steelblue",
     main = "Histograma de precipitacion semanal promedio", xlab = "Precipitacion (mm/semana)")

tam <- (datos_precip$precip_media - min(datos_precip$precip_media)) /
  diff(range(datos_precip$precip_media)) * 2 + 0.5

plot(datos_precip$lon, datos_precip$lat, cex = tam, pch = 21, bg = "dodgerblue",
     xlab = "Longitud", ylab = "Latitud",
     main = "Mapa de posting: precipitacion semanal promedio")

par(mfrow = c(1, 2))
plot(datos_precip$lat, datos_precip$precip_media, xlab = "Latitud", ylab = "Precipitacion",
     main = "Precipitacion vs Latitud")
plot(datos_precip$lon, datos_precip$precip_media, xlab = "Longitud", ylab = "Precipitacion",
     main = "Precipitacion vs Longitud")



#Tendencia
m_trend <- lm(precip_media ~ x_km + y_km, data = datos_precip)
residuales <- residuals(m_trend)


#Matrix D
n <- nrow(datos_precip)
D <- as.matrix(dist(datos_precip[, c("x_km", "y_km")]))

#Semivariograma nube
pares <- which(upper.tri(D), arr.ind = TRUE)
h_ij <- D[upper.tri(D)]
gamma_ij <- (residuales[pares[, 1]] - residuales[pares[, 2]])^2 / 2

plot(h_ij, gamma_ij, pch = 20, col = "grey50",
     xlab = "Distancia (km)", ylab = "Semivarianza",
     main = "Nube de variograma (residuales de precipitacion)")

n_bins <- 10 
bins <- cut(h_ij, breaks = n_bins)
gamma_emp <- tapply(gamma_ij, bins, mean)
N_h <- tapply(gamma_ij, bins, length)
h_medio <- tapply(h_ij, bins, mean)

semivariograma <- na.omit(data.frame(h = h_medio, gamma = gamma_emp, N = N_h))
print(semivariograma)

plot(semivariograma$h, semivariograma$gamma, pch = 19, col = "navy",
     type = "b",
     xlab = "Distancia (km)", ylab = "Semivarianza empirica",
     main = "Semivariograma experimental (precipitacion)")


#Ajuste del modelo

sce_exponencial <- function(theta, h, gamma_emp) {
  c0 <- theta[1]; c1 <- theta[2]; phi <- theta[3]
  gamma_teo <- c0 + c1 * (1 - exp(-h / phi))
  sum((gamma_emp - gamma_teo)^2)
}

ajuste <- optim(
  par = c(c0 = 0.01, c1 = var(residuales), phi = median(semivariograma$h)),
  fn = sce_exponencial, h = semivariograma$h, gamma_emp = semivariograma$gamma,
  method = "L-BFGS-B", lower = c(0, 0, 1)
)


c0_hat <- unname(ajuste$par["c0"])
c1_hat <- unname(ajuste$par["c1"])
phi_hat <- unname(ajuste$par["phi"])
sill_hat <- c0_hat + c1_hat
rango_practico <- 3 * phi_hat


curva_h <- seq(0, max(semivariograma$h), length.out = 100)
curva_gamma <- c0_hat + c1_hat * (1 - exp(-curva_h / phi_hat))

plot(semivariograma$h, semivariograma$gamma, pch = 19, col = "navy",
     xlab = "Distancia (km)", ylab = "Semivarianza",
     main = "Semivariograma vs modelo exponencial ajustado")
lines(curva_h, curva_gamma, col = "red", lwd = 2)


#Kriging Universal
C_hat <- function(h) sill_hat * exp(-h / phi_hat)
Sigma <- matrix(C_hat(as.vector(D)), nrow = n)

kriging_universal <- function(s0, coords_xy, resid, Sigma, C_hat) {
  n_loc <- nrow(coords_xy)
  d0 <- sqrt((coords_xy[, 1] - s0[1])^2 + (coords_xy[, 2] - s0[2])^2)
  c0_vec <- C_hat(d0)
  
  A <- rbind(cbind(Sigma, rep(1, n_loc)), c(rep(1, n_loc), 0))
  b <- c(c0_vec, 1)
  sol <- solve(A, b)
  
  pesos <- sol[1:n_loc]
  lambda <- sol[n_loc + 1]
  prediccion_residual <- unname(sum(pesos * resid))
  varianza_krg <- unname(sill_hat - sum(pesos * c0_vec) - lambda)
  
  c(prediccion = prediccion_residual, varianza = varianza_krg)
}

coords_xy <- as.matrix(datos_precip[, c("x_km", "y_km")])


#Mapa de predicción
grid_x <- seq(min(datos_precip$x_km), max(datos_precip$x_km), length.out = 30)
grid_y <- seq(min(datos_precip$y_km), max(datos_precip$y_km), length.out = 30)
malla <- expand.grid(x_km = grid_x, y_km = grid_y)

pred_malla <- matrix(NA, nrow = nrow(malla), ncol = 2)

for (k in 1:nrow(malla)) {
  s0_k <- c(malla$x_km[k], malla$y_km[k])
  tendencia_k <- predict(m_trend, newdata = malla[k, , drop = FALSE])
  krg_k <- kriging_universal(s0_k, coords_xy, residuales, Sigma, C_hat)
  pred_malla[k, ] <- c(tendencia_k + krg_k["prediccion"], krg_k["varianza"])
}
malla$prediccion <- pred_malla[, 1]
malla$varianza <- pred_malla[, 2]

z_matriz <- matrix(malla$prediccion, nrow = length(grid_x), ncol = length(grid_y))
contour(grid_x, grid_y, z_matriz, xlab = "x (km)", ylab = "y (km)",
        main = "Superficie interpolada de precipitacion (kriging universal)")
points(datos_precip$x_km, datos_precip$y_km, pch = 19, col = "red", cex = 0.6)

# Recorte del mapa final al borde real del departamento (terra, solo para el mapa)
malla_grados <- data.frame(
  lon = malla$x_km / (111.320 * cos(lat0 * pi / 180)) + mean(datos_precip$lon),
  lat = malla$y_km / 110.574 + lat0,
  prediccion = malla$prediccion
)

#Volver a raster
raster_prediccion <- rast(malla_grados, type = "xyz", crs = "EPSG:4326")
raster_prediccion_recortado <- mask(raster_prediccion, valle)

plot(raster_prediccion_recortado, main = "Precipitacion estimada, recortada al Valle del Cauca")
plot(valle, add = TRUE, border = "black")


#Validación cruzada

errores_loo <- numeric(n)

for (i in 1:n) {
  datos_sin_i <- datos_precip[-i, ]
  modelo_sin_i <- lm(precip_media ~ x_km + y_km, data = datos_sin_i)
  resid_sin_i <- residuals(modelo_sin_i)
  Sigma_sin_i <- Sigma[-i, -i]
  
  tendencia_i <- predict(modelo_sin_i, newdata = datos_precip[i, , drop = FALSE])
  krg_i <- kriging_universal(coords_xy[i, ], coords_xy[-i, ], resid_sin_i, Sigma_sin_i, C_hat)
  pred_i <- tendencia_i + krg_i["prediccion"]
  errores_loo[i] <- datos_precip$precip_media[i] - pred_i
}

rmse_loo <- sqrt(mean(errores_loo^2))
print(rmse_loo)

plot(datos_precip$precip_media, errores_loo, pch = 19, col = "darkgreen",
     xlab = "Precipitacion observada/extraida", ylab = "Error (obs - predicho)",
     main = "Diagnostico de validacion cruzada (leave-one-out)")
abline(h = 0, col = "red", lty = 2)
