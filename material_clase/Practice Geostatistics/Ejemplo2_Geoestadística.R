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


#Generación de datos
set.seed(2019)                                   # reproducibilidad
n_total <- 150                                   # numero de estaciones
x <- runif(n_total, 300000, 700000)              # coordenada este (m)
y <- runif(n_total, 5000000, 5400000)            # coordenada norte (m)

D_total <- as.matrix(dist(cbind(x, y)))          # distancias para simular dependencia
Sigma_sim <- 0.09 * exp(-D_total / 60000)        # covarianza exponencial "verdadera"
L_sim <- chol(Sigma_sim + diag(1e-8, n_total))   # factorizacion de Cholesky


# El intercepto se centra segun x,y para que el pH quede en un rango realista (~4-9)
x_coef <- -3.2e-6                                # pendiente en x (gradiente NE-SO)
y_coef <- -4.0e-6                                # pendiente en y (gradiente NE-SO)
intercepto <- 8.5 - x_coef * mean(x) - y_coef * mean(y)
tendencia_real <- intercepto + x_coef * x + y_coef * y
pH <- tendencia_real + as.vector(t(L_sim) %*% rnorm(n_total)) + rnorm(n_total, sd = 0.15)

ca_geo <- data.frame(ID = 1:n_total, x = x, y = y, pH = pH)

miss_idx <- sample(1:n_total, 15)               
ca_geo$pH[miss_idx] <- NA


########################################################
#EDA espacial
########################################################

print(str(ca_geo))                               # estructura del data frame
print(names(ca_geo))                             # variables disponibles
print(summary(ca_geo$pH))                        # resumen numerico del pH

hist(ca_geo$pH, breaks = 12, col = "lightblue",
     main = "Histograma de pH", xlab = "pH")


#Datos faltantes
miss <- is.na(ca_geo$pH)                         # vector logico de faltantes
print(table(miss))


#Mapa de puntos

rango_pH <- range(ca_geo$pH[!miss])
colores <- colorRampPalette(c("red", "yellow", "darkgreen"))(100)
paleta_idx <- round((ca_geo$pH[!miss] - rango_pH[1]) / diff(rango_pH) * 99) + 1

plot(ca_geo$x[!miss], ca_geo$y[!miss], pch = 19, col = colores[paleta_idx],
     xlab = "x (m)", ylab = "y (m)", main = "Mapa de pH observado")
points(ca_geo$x[miss], ca_geo$y[miss], pch = 4, col = "black")  # estaciones sin dato


#Verificación de tendencia
m_trend <- lm(pH ~ x + y, data = ca_geo)         # tendencia lineal en las coordenadas
print(summary(m_trend))


#Predicción en estaciones faltantes

ca_geo_miss <- ca_geo[miss, ]                    # estaciones con pH faltante
predicciones_tendencia <- predict(m_trend, newdata = ca_geo_miss, se.fit = TRUE)
print(predicciones_tendencia$fit)

# Probabilidad de que el suelo sea alcalino (pH > 7) segun la tendencia
pAlkaline <- 1 - pnorm(7, mean = predicciones_tendencia$fit, sd = predicciones_tendencia$se.fit)
print(pAlkaline)

hist(pAlkaline, breaks = 8, col = "khaki",
     main = "Probabilidad de suelo alcalino (solo tendencia)", xlab = "P(pH > 7)")


#Semivariograma nube

datos_obs <- ca_geo[!miss, ]                     # solo estaciones con dato observado
n <- nrow(datos_obs)

D <- as.matrix(dist(datos_obs[, c("x", "y")]))   

pares <- which(upper.tri(D), arr.ind = TRUE)     # pares unicos i<j
h_ij <- D[upper.tri(D)]                          
gamma_ij <- (datos_obs$pH[pares[, 1]] - datos_obs$pH[pares[, 2]])^2 / 2 

cutoff <- 300000                                 # distancia maxima de analisis (m)
en_cutoff <- h_ij <= cutoff

plot(h_ij[en_cutoff], gamma_ij[en_cutoff], pch = 20, col = "grey50",
     xlab = "Distancia (m)", ylab = "Semivarianza",
     main = "Nube de variograma (pH crudo)")


#Semivariograma empírico

n_bins <- 12
bins <- cut(h_ij[en_cutoff], breaks = n_bins)
gamma_emp <- tapply(gamma_ij[en_cutoff], bins, mean)   # semivarianza media por bin
N_h <- tapply(gamma_ij[en_cutoff], bins, length)       # numero de pares por bin
h_medio <- tapply(h_ij[en_cutoff], bins, mean)         # distancia media por bin

variograma_crudo <- na.omit(data.frame(h = h_medio, gamma = gamma_emp, N = N_h))
print(variograma_crudo)

plot(variograma_crudo$h, variograma_crudo$gamma, pch = 19, col = "darkblue",
     type = "b",
     xlab = "Distancia (m)", ylab = "Semivarianza",
     main = "Variograma empirico (pH crudo, sin tendencia)")



residuales <- residuals(lm(pH ~ x + y, data = datos_obs)) 
gamma_res <- (residuales[pares[, 1]] - residuales[pares[, 2]])^2 / 2

gamma_emp_res <- tapply(gamma_res[en_cutoff], bins, mean)
variograma_residual <- na.omit(data.frame(h = h_medio, gamma = gamma_emp_res, N = N_h))
print(variograma_residual)

plot(variograma_residual$h, variograma_residual$gamma, pch = 19, col = "forestgreen",
     type = "b",
     xlab = "Distancia (m)", ylab = "Semivarianza",
     main = "Variograma de los residuales (tendencia removida)")



#Ajuste por MCO no lineales

sce_exponencial <- function(theta, h, gamma_emp) {
  c0 <- theta[1]; c1 <- theta[2]; phi <- theta[3]  # nugget, contribucion, rango
  gamma_teo <- c0 + c1 * (1 - exp(-h / phi))
  sum((gamma_emp - gamma_teo)^2)
}

ajuste <- optim(
  par = c(c0 = 0.02, c1 = var(residuales), phi = median(variograma_residual$h)),
  fn = sce_exponencial, h = variograma_residual$h, gamma_emp = variograma_residual$gamma,
  method = "L-BFGS-B", lower = c(0, 0, 1000)
)


c0_hat <- unname(ajuste$par["c0"])               # nugget
c1_hat <- unname(ajuste$par["c1"])               # contribucion (sill - nugget)
phi_hat <- unname(ajuste$par["phi"])             # parametro de rango
sill_hat <- c0_hat + c1_hat                      # umbral (sill)
rango_practico <- 3 * phi_hat                    # distancia de independencia aproximada


curva_h <- seq(0, cutoff, length.out = 100)
curva_gamma <- c0_hat + c1_hat * (1 - exp(-curva_h / phi_hat))

plot(variograma_residual$h, variograma_residual$gamma, pch = 19, col = "forestgreen",
     xlab = "Distancia (m)", ylab = "Semivarianza",
     main = "Variograma de residuales vs modelo exponencial ajustado")
lines(curva_h, curva_gamma, col = "red", lwd = 2)


#Predicción con kriging universal

C_hat <- function(h) sill_hat * exp(-h / phi_hat)   # funcion de covarianza ajustada
Sigma_res <- matrix(C_hat(as.vector(D)), nrow = n)  # covarianza de los residuales observados

kriging_universal <- function(s0, coords, resid, Sigma, C_hat) {
  n_loc <- nrow(coords)
  d0 <- sqrt((coords[, 1] - s0[1])^2 + (coords[, 2] - s0[2])^2)  # distancias al punto nuevo
  c0_vec <- C_hat(d0)
  
  A <- rbind(cbind(Sigma, rep(1, n_loc)), c(rep(1, n_loc), 0))   # sistema con Lagrange
  b <- c(c0_vec, 1)
  sol <- solve(A, b)
  
  pesos <- sol[1:n_loc]
  lambda <- sol[n_loc + 1]
  prediccion_residual <- unname(sum(pesos * resid))
  varianza_krg <- unname(sill_hat - sum(pesos * c0_vec) - lambda)
  
  c(prediccion = prediccion_residual, varianza = varianza_krg)
}

coords_obs <- as.matrix(datos_obs[, c("x", "y")])
pred_krg_pH <- numeric(nrow(ca_geo_miss))
var_krg_pH <- numeric(nrow(ca_geo_miss))

for (k in 1:nrow(ca_geo_miss)) {
  s0_k <- c(ca_geo_miss$x[k], ca_geo_miss$y[k])
  tendencia_k <- predict(m_trend, newdata = ca_geo_miss[k, ])          # componente de tendencia
  krg_k <- kriging_universal(s0_k, coords_obs, residuales, Sigma_res, C_hat)
  pred_krg_pH[k] <- tendencia_k + krg_k["prediccion"]                  # tendencia + residual interpolado
  var_krg_pH[k] <- krg_k["varianza"]
}

comparacion <- data.frame(
  ID = ca_geo_miss$ID,
  pred_tendencia = predicciones_tendencia$fit,
  pred_kriging = pred_krg_pH,
  varianza_kriging = var_krg_pH
)


# Probabilidad de exceder pH = 7 usando la prediccion completa de kriging
pAlkaline_krg <- 1 - pnorm(7, mean = pred_krg_pH, sd = sqrt(var_krg_pH))
print(pAlkaline_krg)


#Mapa de superficie

grid_x <- seq(min(ca_geo$x), max(ca_geo$x), length.out = 15)
grid_y <- seq(min(ca_geo$y), max(ca_geo$y), length.out = 15)
malla <- expand.grid(x = grid_x, y = grid_y)

pred_malla <- matrix(NA, nrow = nrow(malla), ncol = 2)
for (k in 1:nrow(malla)) {
  s0_k <- c(malla$x[k], malla$y[k])
  tendencia_k <- predict(m_trend, newdata = malla[k, , drop = FALSE])
  krg_k <- kriging_universal(s0_k, coords_obs, residuales, Sigma_res, C_hat)
  pred_malla[k, ] <- c(tendencia_k + krg_k["prediccion"], krg_k["varianza"])
}
malla$prediccion <- pred_malla[, 1]
malla$varianza <- pred_malla[, 2]

z_matriz <- matrix(malla$prediccion, nrow = length(grid_x), ncol = length(grid_y))
contour(grid_x, grid_y, z_matriz, xlab = "x (m)", ylab = "y (m)",
        main = "Superficie interpolada de pH (kriging universal)")
points(datos_obs$x, datos_obs$y, pch = 19, col = "red", cex = 0.6)

v_matriz <- matrix(malla$varianza, nrow = length(grid_x), ncol = length(grid_y))
contour(grid_x, grid_y, v_matriz, xlab = "x (m)", ylab = "y (m)",
        main = "Varianza de kriging")
points(datos_obs$x, datos_obs$y, pch = 19, col = "red", cex = 0.6)

