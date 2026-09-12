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



set.seed(2026)                                  # reproducibilidad del muestreo
idx <- sample(1:nrow(quakes), 60)               # submuestra de 60 sismos
datos <- quakes[idx, c("long", "lat", "depth")] # longitud, latitud y profundidad
rownames(datos) <- NULL
head(datos)


#############################################################
# EDA espacial
#############################################################

resumen_depth <- summary(datos$depth)
resumen_depth

cv_depth <- sd(datos$depth) / mean(datos$depth)
cv_depth

# Histograma de la variable regionalizada
hist(datos$depth, breaks = 12, col = "lightblue",
     main = "Histograma de profundidad (km)", xlab = "Profundidad")

# Mapa de "posting": tamano del punto proporcional al valor de depth
tam <- (datos$depth - min(datos$depth)) / diff(range(datos$depth)) * 2 + 0.5
plot(datos$long, datos$lat, cex = tam, pch = 21, bg = "steelblue",
     xlab = "Longitud", ylab = "Latitud",
     main = "Mapa de posting: tamano ~ profundidad")

# Dispersion de la variable contra cada coordenada (deteccion de tendencia)
par(mfrow = c(1, 2))
plot(datos$long, datos$depth, xlab = "Longitud", ylab = "Profundidad",
     main = "Depth vs Longitud")
plot(datos$lat, datos$depth, xlab = "Latitud", ylab = "Profundidad",
     main = "Depth vs Latitud")


# Tendencia lineal en longitud y latitud
tendencia <- lm(depth ~ long + lat, data = datos)
tendencia

residuales <- residuals(tendencia)
datos$residuales <- residuales

#Matrid D

n <- nrow(datos)
D <- matrix(0, n, n)                        

for (i in 1:n) {
  for (j in 1:n) {
    D[i, j] <- sqrt((datos$long[i] - datos$long[j])^2 +
                      (datos$lat[i]  - datos$lat[j])^2) 
  }
}

dim(D)


############################################################
#Semivariograma experimental
############################################################

pares <- which(upper.tri(D), arr.ind = TRUE)    
h_ij   <- D[upper.tri(D)]                       # distancias de esos pares
gamma_ij <- (residuales[pares[, 1]] - residuales[pares[, 2]])^2 / 2 

# Nube de semivariograma nube
plot(h_ij, gamma_ij, pch = 20, col = "grey50",
     xlab = "Distancia h", ylab = "Semivarianza",
     main = "Semivariogram cloud")

# Agrupacion en intervalos de distancia (bins)
n_bins <- 10
bins   <- cut(h_ij, breaks = n_bins)
gamma_emp <- tapply(gamma_ij, bins, mean)       
N_h       <- tapply(gamma_ij, bins, length) 
h_medio   <- tapply(h_ij, bins, mean)       

semivariograma <- data.frame(h = h_medio, gamma = gamma_emp, N = N_h)
semivariograma <- na.omit(semivariograma)

plot(semivariograma$h, semivariograma$gamma, pch = 19, col = "darkblue",
     type = "b",
     xlab = "Distancia h", ylab = "Semivarianza empirica",
     main = "Semivariograma experimental")

###########################################################
#Ajuste de modelo teórico
###########################################################

# Suma de cuadrados entre el semivariograma empirico y el modelo teorico
sce_exponencial <- function(theta, h, gamma_emp) {
  c0 <- theta[1]; c1 <- theta[2]; phi <- theta[3]
  gamma_teo <- c0 + c1 * (1 - exp(-h / phi))
  sum((gamma_emp - gamma_teo)^2)
}

ajuste <- optim(
  par = c(c0 = 0, c1 = var(residuales), phi = median(semivariograma$h)),
  fn  = sce_exponencial, h = semivariograma$h, gamma_emp = semivariograma$gamma,
  method = "L-BFGS-B", lower = c(0, 0, 0.01)
)

#Limited-memory Broyden-Fletcher-Goldfarb-Shanno with Box (L-BFGS-B)
#Optimización numérica quasi-Newton Rapson


c0_hat  <- unname(ajuste$par["c0"])             # nugget
c1_hat  <- unname(ajuste$par["c1"])             # contribucion (sill - nugget)
phi_hat <- unname(ajuste$par["phi"])            # parametro de rango
sill_hat  <- c0_hat + c1_hat                    # umbral (sill)
rango_practico <- 3 * phi_hat                   # rango practico (95% del sill)


#curva ajustada sobre el semivariograma empirico
curva_h <- seq(0, max(semivariograma$h), length.out = 100)
curva_gamma <- c0_hat + c1_hat * (1 - exp(-curva_h / phi_hat))

plot(semivariograma$h, semivariograma$gamma, pch = 19, col = "darkblue",
     xlab = "Distancia h", ylab = "Semivarianza",
     main = "Semivariograma empirico vs modelo exponencial")
lines(curva_h, curva_gamma, col = "red", lwd = 2)

#################################################################
#Función de covarianza
#################################################################


# C(h) = sill - gamma(h), evaluada sobre toda la matriz de distancias D
C_hat   <- function(h) sill_hat * exp(-h / phi_hat)
Sigma   <- matrix(C_hat(as.vector(D)), nrow = n)
dim(Sigma)

####################################################################
#Predicción con Kriging Ordinario
####################################################################

kriging_ordinario <- function(s0, coords, z, Sigma, C_hat) {
  n_loc <- nrow(coords)
  d0 <- sqrt((coords[, 1] - s0[1])^2 + (coords[, 2] - s0[2])^2)  
  c0_vec <- C_hat(d0)                                            
  
  # Sistema aumentado con la restriccion de pesos que suman 1 (Lagrange)
  A <- rbind(cbind(Sigma, rep(1, n_loc)), c(rep(1, n_loc), 0))
  b <- c(c0_vec, 1)
  sol <- solve(A, b)                                             # pesos + multiplicador
  
  pesos   <- sol[1:n_loc]
  lambda  <- sol[n_loc + 1]
  prediccion   <- unname(sum(pesos * z))
  varianza_krg <- unname(sill_hat - sum(pesos * c0_vec) - lambda)
  
  c(prediccion = prediccion, varianza = varianza_krg)
}

# Prediccion en un punto no observado
s0_ejemplo <- c(mean(datos$long), mean(datos$lat))
coords_mat <- as.matrix(datos[, c("long", "lat")])
resultado_krg <- kriging_ordinario(s0_ejemplo, coords_mat, datos$depth, Sigma, C_hat)



#Mapa de predicción

grid_long <- seq(min(datos$long), max(datos$long), length.out = 15)
grid_lat  <- seq(min(datos$lat),  max(datos$lat),  length.out = 15)
malla     <- expand.grid(long = grid_long, lat = grid_lat)

pred_malla <- matrix(NA, nrow = nrow(malla), ncol = 2)

for (k in 1:nrow(malla)) {
  s0_k <- c(malla$long[k], malla$lat[k])
  pred_malla[k, ] <- kriging_ordinario(s0_k, coords_mat, datos$depth, Sigma, C_hat)
}
malla$prediccion <- pred_malla[, 1]
malla$varianza   <- pred_malla[, 2]

# Curvas de nivel
z_matriz <- matrix(malla$prediccion, nrow = length(grid_long), ncol = length(grid_lat))
contour(grid_long, grid_lat, z_matriz,
        xlab = "Longitud", ylab = "Latitud",
        main = "Superficie interpolada (kriging ordinario)")
points(datos$long, datos$lat, pch = 19, col = "red", cex = 0.6)

# Mapa de la varianza de predicción
v_matriz <- matrix(malla$varianza, nrow = length(grid_long), ncol = length(grid_lat))
contour(grid_long, grid_lat, v_matriz,
        xlab = "Longitud", ylab = "Latitud",
        main = "Varianza de kriging")
points(datos$long, datos$lat, pch = 19, col = "red", cex = 0.6)


#Validación cruzada

errores_loo <- numeric(n)

for (i in 1:n) {
  coords_sin_i <- coords_mat[-i, ]
  z_sin_i      <- datos$depth[-i]
  Sigma_sin_i  <- Sigma[-i, -i]
  pred_i <- kriging_ordinario(coords_mat[i, ], coords_sin_i, z_sin_i, Sigma_sin_i, C_hat)
  errores_loo[i] <- datos$depth[i] - pred_i["prediccion"]
}

rmse_loo <- sqrt(mean(errores_loo^2))            # raiz del error cuadratico medio


# Grafico de valores observados vs error de prediccion
plot(datos$depth, errores_loo, pch = 19, col = "darkgreen",
     xlab = "Profundidad observada", ylab = "Error (obs - predicho)",
     main = "Validacion cruzada (leave-one-out)")
abline(h = 0, col = "red", lty = 2)
