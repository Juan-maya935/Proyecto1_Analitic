# Ejemplo 4: EDA espacial y geoestadística
# Temperatura máxima en el valle del cauca

library(terra)


#Cargar imaenes
carpeta_imagenes <- "imagenes_tmax_valle_cauca"     # ajuste esta ruta si es necesario

archivos_tif <- list.files(carpeta_imagenes, pattern = "\\.tif$", full.names = TRUE)
if (length(archivos_tif) == 0) {
}


pila_tmax <- rast(archivos_tif)              # se apilan todas las imagenes en un solo objeto
names(pila_tmax) <- tools::file_path_sans_ext(basename(archivos_tif))   # nombre = fecha del archivo
pila_tmax


#Borde del valle del cauca
mascara_valida <- !is.na(pila_tmax[[1]])                 # TRUE donde hay dato, FALSE donde es NA
borde_recuperado <- as.polygons(mascara_valida, dissolve = TRUE)
borde_recuperado <- borde_recuperado[borde_recuperado[[1]] == 1, ]  # solo el area con dato valido

plot(pila_tmax[[1]], main = "Primera imagen cargada, con el borde real del Valle del Cauca")
plot(borde_recuperado, add = TRUE, border = "black", lwd = 1.5)


#Puntos de muestreo
set.seed(2026)
n_puntos_finales <- 40
n_puntos_muestreados <- 80

puntos <- spatSample(borde_recuperado, size = n_puntos_muestreados, method = "random")
serie_extraida <- extract(pila_tmax, puntos)   # ID + una columna por mes disponible
coords_xy_grados <- crds(puntos)               # longitud y latitud de cada punto

n <- nrow(serie_extraida)
n


#Temperatura promedio
tmax_promedio <- rowMeans(serie_extraida[, -1, drop = FALSE], na.rm = TRUE)

tiempo_en_meses <- 1:(ncol(serie_extraida) - 1)
tendencia_decada <- sapply(1:n, function(i) {
  valores_i <- as.numeric(serie_extraida[i, -1])
  tryCatch({
    unname(coef(lm(valores_i ~ tiempo_en_meses))[2]) * 12 * 10
  }, error = function(e) NA_real_)      # por si algun punto quedo con muy pocos datos validos
})

datos <- data.frame(
  ID = serie_extraida$ID,
  lon = coords_xy_grados[, 1],
  lat = coords_xy_grados[, 2],
  tmax_promedio = tmax_promedio,
  tendencia_decada = tendencia_decada
)

datos <- na.omit(datos)     # se descartan puntos sin datos suficientes en su serie
datos <- datos[1:min(n_puntos_finales, nrow(datos)), ]   # se conservan como maximo los pedidos
datos

n <- nrow(datos)            # se actualiza n al numero final de puntos validos

head(datos)


#Proyección
lat0 <- mean(datos$lat)
datos$x_km <- (datos$lon - mean(datos$lon)) * 111.320 * cos(lat0 * pi / 180)
datos$y_km <- (datos$lat - lat0) * 110.574
head(datos)


#EDA espacial
summary(datos$tmax_promedio)
summary(datos$tendencia_decada)
sd(datos$tmax_promedio) / mean(datos$tmax_promedio)

par(mfrow = c(1, 2))
hist(datos$tmax_promedio, breaks = 10, col = "tomato",
     main = "Temperatura maxima promedio", xlab = "Grados C")
hist(datos$tendencia_decada, breaks = 10, col = "orange",
     main = "Tendencia de calentamiento", xlab = "Grados C por decada")


tam <- (datos$tmax_promedio - min(datos$tmax_promedio)) /
  diff(range(datos$tmax_promedio)) * 2 + 0.5
plot(datos$lon, datos$lat, cex = tam, pch = 21, bg = "red",
     xlab = "Longitud", ylab = "Latitud",
     main = "Mapa de posting: temperatura maxima promedio")

par(mfrow = c(1, 2))
plot(datos$lat, datos$tmax_promedio, xlab = "Latitud", ylab = "Temperatura",
     main = "Temperatura vs Latitud")
plot(datos$lon, datos$tmax_promedio, xlab = "Longitud", ylab = "Temperatura",
     main = "Temperatura vs Longitud")
par(mfrow = c(1, 1))


#Extraer tendencia
m_trend <- lm(tendencia_decada ~ x_km + y_km, data = datos)
print(summary(m_trend))
residuales <- residuals(m_trend)


#Distancias
D <- as.matrix(dist(datos[, c("x_km", "y_km")]))


#Semivariagrama nube y muestral
pares <- which(upper.tri(D), arr.ind = TRUE)      # pares unicos i<j
h_ij <- D[upper.tri(D)]                           # distancia de cada par
gamma_ij <- (residuales[pares[, 1]] - residuales[pares[, 2]])^2 / 2   # semivarianza cruda

plot(h_ij, gamma_ij, pch = 20, col = "grey50",
     xlab = "Distancia (km)", ylab = "Semivarianza",
     main = "Nube de variograma (residuales)")

n_bins <- 10
bins <- cut(h_ij, breaks = n_bins)
gamma_emp <- tapply(gamma_ij, bins, mean)         # semivarianza media por bin
N_h <- tapply(gamma_ij, bins, length)             # numero de pares por bin
h_medio <- tapply(h_ij, bins, mean)               # distancia media por bin

semivariograma <- na.omit(data.frame(h = h_medio, gamma = gamma_emp, N = N_h))
print(semivariograma)

plot(semivariograma$h, semivariograma$gamma, pch = 19, col = "navy", type = "b",
     xlab = "Distancia (km)", ylab = "Semivarianza empirica",
     main = "Semivariograma experimental")


#Ajuste de modelo

sce_exponencial <- function(theta, h, gamma_emp) {
  c0 <- theta[1]; c1 <- theta[2]; phi <- theta[3]
  gamma_teo <- c0 + c1 * (1 - exp(-h / phi))
  sum((gamma_emp - gamma_teo)^2)
}

phis_iniciales <- quantile(semivariograma$h, probs = c(0.1, 0.25, 0.5, 0.75, 1))
ajustes_candidatos <- lapply(phis_iniciales, function(phi0) {
  optim(par = c(c0 = 0, c1 = var(residuales), phi = phi0),
        fn = sce_exponencial, h = semivariograma$h, gamma_emp = semivariograma$gamma,
        method = "L-BFGS-B",
        lower = c(0, 0, 1), upper = c(max(semivariograma$gamma), 5 * max(semivariograma$gamma), 3 * max(semivariograma$h)))
})
sce_por_arranque <- sapply(ajustes_candidatos, function(a) a$value)
print(sce_por_arranque)                 # suma de cuadrados de cada arranque (diagnostico)
ajuste <- ajustes_candidatos[[which.min(sce_por_arranque)]]   # se conserva el mejor de todos
print(ajuste$par)

c0_hat <- unname(ajuste$par["c0"])       # nugget
c1_hat <- unname(ajuste$par["c1"])       # contribucion (sill - nugget)
phi_hat <- unname(ajuste$par["phi"])     # parametro de rango
sill_hat <- c0_hat + c1_hat              # umbral (sill)
rango_practico <- 3 * phi_hat            # distancia de independencia aproximada (km)
print(c(nugget = c0_hat, sill = sill_hat, rango_practico_km = rango_practico))


curva_h <- seq(0, max(semivariograma$h), length.out = 100)
curva_gamma <- c0_hat + c1_hat * (1 - exp(-curva_h / phi_hat))

plot(semivariograma$h, semivariograma$gamma, pch = 19, col = "navy",
     xlab = "Distancia (km)", ylab = "Semivarianza",
     main = "Semivariograma vs modelo exponencial ajustado")
lines(curva_h, curva_gamma, col = "red", lwd = 2)


#Predicción

C_hat <- function(h) sill_hat * exp(-h / phi_hat)
dim(C_hat)
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

coords_xy <- as.matrix(datos[, c("x_km", "y_km")])


#Mapa de interpolación
grid_x <- seq(min(datos$x_km), max(datos$x_km), length.out = 30)
grid_y <- seq(min(datos$y_km), max(datos$y_km), length.out = 30)
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
        main = "Tendencia de calentamiento interpolada (kriging universal)")
points(datos$x_km, datos$y_km, pch = 19, col = "red", cex = 0.6)

#Grilla en terminos de coordenadas geográficas
malla_grados <- data.frame(
  lon = malla$x_km / (111.320 * cos(lat0 * pi / 180)) + mean(datos$lon),
  lat = malla$y_km / 110.574 + lat0,
  prediccion = malla$prediccion
)
raster_prediccion <- rast(malla_grados, type = "xyz", crs = crs(pila_tmax))
raster_prediccion_recortado <- mask(raster_prediccion, borde_recuperado)
plot(raster_prediccion_recortado, main = "")
plot(borde_recuperado, add = TRUE, border = "black")


#Validación cruxadas
errores_loo <- numeric(n)
for (i in 1:n) {
  datos_sin_i <- datos[-i, ]
  modelo_sin_i <- lm(tendencia_decada ~ x_km + y_km, data = datos_sin_i)
  resid_sin_i <- residuals(modelo_sin_i)
  Sigma_sin_i <- Sigma[-i, -i]

  tendencia_i <- predict(modelo_sin_i, newdata = datos[i, , drop = FALSE])
  krg_i <- kriging_universal(coords_xy[i, ], coords_xy[-i, ], resid_sin_i, Sigma_sin_i, C_hat)
  pred_i <- tendencia_i + krg_i["prediccion"]
  errores_loo[i] <- datos$tendencia_decada[i] - pred_i
}

rmse_loo <- sqrt(mean(errores_loo^2))


plot(datos$tendencia_decada, errores_loo, pch = 19, col = "darkgreen",
     xlab = "Tendencia observada/extraida", ylab = "Error (obs - predicho)",
     main = "Diagnostico de validacion cruzada (leave-one-out)")
abline(h = 0, col = "red", lty = 2)


