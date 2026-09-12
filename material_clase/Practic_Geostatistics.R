## ============================================================
## PRÁCTICA DE GEOESTADÍSTICA - DATASET MEUSE
## Kriging ordinario sobre concentración de Zn en el suelo
## ============================================================



paquetes <- c("sp", "gstat")
instalados <- rownames(installed.packages())
for (p in paquetes) {
  if (!(p %in% instalados)) install.packages(p)
}

library(sp)
library(gstat)


# Datos "meuse": polución del suelo con concentraciones de metales pesados
data(meuse)
str(meuse)  # algunas variables son continuas, otras categóricas (factor)

##Análisis exploratorio (EDA) 
summary(meuse$zinc)

hist(meuse$zinc, main = "Distribución de Zn", xlab = "Zn (ppm)")
plot(meuse$x, meuse$zinc, xlab = "Coordenada X", ylab = "Zn (ppm)")
plot(meuse$y, meuse$zinc, xlab = "Coordenada Y", ylab = "Zn (ppm)")
boxplot(meuse$zinc, main = "Boxplot de Zn")

##La distribución de Zn está sesgada a la derecha; se aplica log10
# para aproximarla a una distribución normal, requisito habitual
# para el ajuste del variograma y el kriging ordinario.

meuse$logZn <- log10(meuse$zinc)
hist(meuse$logZn, breaks = 16, main = "Distribución de log10(Zn)",
     xlab = "log10(Zn)")

##Conversión a objeto espacial -------------------------------------
coordinates(meuse) <- c("x", "y")
class(meuse)  # SpatialPointsDataFrame

# Visualización de puntos muestreados
plot(meuse, asp = 1, pch = 1)  # asp = 1: mismas escalas en X e Y

# Margen del río Meuse, como referencia espacial
data(meuse.riv)
lines(meuse.riv)

# Tamaño del punto proporcional a la concentración de Zn
plot(meuse, asp = 1, cex = 4 * meuse$zinc / max(meuse$zinc), pch = 1,
     main = "Concentración de Zn (tamaño proporcional)")

## Variograma empírico
# Distancia entre saltos (width) de 90 m, hasta un cutoff de 1300 m
ve <- variogram(logZn ~ 1, meuse, cutoff = 1300, width = 90)
plot(ve, plot.numbers = TRUE, asp = 1, pch = 19, type = "l")

##Ajuste del modelo de variograma ----------------------------------
# Parámetros estimados visualmente a partir del variograma empírico:
#   Rango (range)          850 m  -> distancia a partir de la cual
#                                      no hay dependencia espacial
#   Pepita (nugget)        0.01   -> semivarianza a distancia 0
#   Meseta total (sill)    0.13   -> semivarianza en el rango
#   Meseta parcial         0.12   -> sill - nugget


vt <- vgm(psill = 0.12, model = "Sph", range = 850, nugget = 0.01)
vt
plot(ve, pl = TRUE, model = vt)

# Ajuste automático por mínimos cuadrados ponderados sobre el modelo inicial
va <- fit.variogram(ve, vt)
va
plot(ve, pl = TRUE, model = va)

# Suma de cuadrados del error del ajuste (menor es mejor)
attr(va, "SSErr")

##Grilla de predicción ----------------------------------------------
data(meuse.grid)  # malla regular de 40 m x 40 m incluida en el dataset
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE  # se define como objeto raster
plot(meuse.grid, main = "Grilla de predicción")

#Kriging ordinario --------------------------------------------------
ok <- krige(logZn ~ 1, locations = meuse, newdata = meuse.grid, model = va)

# Retransformación a la escala original (ppm)
ok$pred <- 10^(ok$var1.pred)
str(ok)

# Capa de puntos muestreados para superponer en el mapa de predicción
pts.s <- list("sp.points", meuse, pch = 3, col = "black")

spplot(ok,
       "var1.pred",                       # predicción en escala log10
       col.regions = rev(heat.colors(50)), # paleta de colores
       main = "Predicción OK, log10-ppm Zn",
       sp.layout = list(pts.s))            # puntos de muestreo superpuestos

# Mapa de la varianza del error de predicción (incertidumbre)
spplot(ok,
       "var1.var",
       col.regions = rev(heat.colors(50)),
       main = "Varianza de predicción OK",
       sp.layout = list(pts.s))

#Validación cruzada del modelo --------------------------------------
# Leave-one-out cross-validation: evalúa qué tan bien predice el modelo
# ajustado en ubicaciones donde sí se conoce el valor real.
cv <- krige.cv(logZn ~ 1, locations = meuse, model = va)
summary(cv)

# RMSE en la escala transformada (log10)
rmse <- sqrt(mean(cv$residual^2))
rmse


