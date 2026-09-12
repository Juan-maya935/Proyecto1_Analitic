#Instalación de librerias
install.packages("ggplot2")
install.packages("geostan")
install.packages("gridExtra")

install.packages("")

#Cargar librerias
library(ggplot2)
library(geostan)
library(gridExtra)
library(openxlsx)
library(dplyr)


#Directorio de trabajo



df_contamina <- read.xlsx("df_contamina.xlsx")
#df_contamina <- df_contamina %>% sample_n(120)


#Análisis Exploratorio de Datos espaciales

# Histograma de PM2.5
ggplot(df_contamina, aes(x = PM2.5)) +
  geom_histogram(binwidth = 2, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribución de PM2.5", x = "PM2.5 (μg/m³)", y = "Frecuencia")

#Prueba de normalidad
shapiro.test(df_contamina$PM2.5)


# Boxplot de PM2.5
ggplot(df_contamina, aes(y = PM2.5)) +
  geom_boxplot(fill = "orange", color = "black") +
  theme_minimal() +
  labs(title = "Boxplot de PM2.5", y = "PM2.5 (μg/m³)")

# Mapa de calor de PM2.5
ggplot(df_contamina, aes(x = x, y = y, fill = PM2.5)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Mapa de Concentración de PM2.5", fill = "PM2.5 (μg/m³)")




# Construcción de la matriz de pesos espaciales
n <- nrow(df_contamina)
W <- matrix(0, n, n)

# Definir vecinos basados en distancia euclidiana
calc_dist <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

# Definir umbral de vecindad basado en la densidad de los datos
umbral <- quantile(dist(df_contamina[, c("x", "y")]), 0.05)  # Considerando el 5% de las distancias más cortas

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j && calc_dist(df_contamina$x[i], df_contamina$y[i],
                            df_contamina$x[j], df_contamina$y[j]) < umbral) {
      W[i, j] <- 1
    }
  }
}

class(W)

# Normalización de la matriz de pesos
row_sums <- rowSums(W)
W <- matrix(ifelse(row_sums == 0, 0, W / row_sums), nrow = n, ncol = n)


# Cálculo del Índice de Moran
z_mean <- mean(df_contamina$PM2.5)
num <- 0
denom <- sum((df_contamina$PM2.5 - z_mean)^2)

for (i in 1:n) {
  for (j in 1:n) {
    num <- num + W[i, j] * (df_contamina$PM2.5[i] - z_mean) * (df_contamina$PM2.5[j] - z_mean)
  }
}


# Índice de Moran
z_mean <- mean(df_contamina$PM2.5)
num <- 0
denom <- sum((df_contamina$PM2.5 - z_mean)^2)

for (i in 1:n) {
  for (j in 1:n) {
    num <- num + W[i, j] * (df_contamina$PM2.5[i] - z_mean) * (df_contamina$PM2.5[j] - z_mean)
  }
}


moran_I <- (n / sum(W)) * (num / denom)



#Interpretación:
#Valor cercano a cero: Un Índice de Moran cercano a cero indica que no hay autocorrelación
#espacial significativa. Esto significa que los datos están distribuidos de manera aleatoria
#en el espacio, sin un patrón claro de agrupamiento o dispersión.

#Valor positivo: Si el Índice de Moran fuera positivo, indicaría autocorrelación espacial positiva
#es decir, que valores similares tienden a agruparse en el espacio.

#Valor negativo: Un valor negativo indicaría autocorrelación espacial negativa, lo que significa
#que valores diferentes tienden a estar cerca unos de otros.


# Cálculo del Índice de Geary
num_geary <- 0
for (i in 1:n) {
  for (j in 1:n) {
    num_geary <- num_geary + W[i, j] * (df_contamina$PM2.5[i] - df_contamina$PM2.5[j])^2
  }
}

guary_C <- ((n - 1) / (2 * sum(W))) * (num_geary / denom)


#Interpretación del Índice de Geary:
#Valor cercano a 1:indica que no hay autocorrelación espacial significativa.
#Esto sugiere que los datos están distribuidos de manera aleatoria en el espacio,
#sin un patrón claro de agrupamiento o dispersión.

#Valor menor que 1: Un valor menor que 1 indica autocorrelación espacial positiva.
#Esto significa que valores similares tienden a agruparse en el espacio.

#Valor mayor que 1:   Un valor mayor que 1 sugiere autocorrelación espacial negativa.
#Esto implica que valores diferentes tienden a estar cerca unos de otros.

#El valor de 1.386772, indica una ligera tendencia hacia la autocorrelación espacial negativa,
#pero cercana a 1, lo que sugiere que la autocorrelación no es fuerte y los datos podrían
#considerarse relativamente aleatorios en su distribución espacial.


# Cálculo de semivariograma
max_dist <- max(dist(df_contamina[, c("x", "y")]))
bins <- seq(1, max_dist, length.out = 20)
semivars <- rep(0, length(bins))

for (h in 1:(length(bins)-1)) {
  num <- 0
  count <- 0
  
  for (i in 1:n) {
    for (j in (i+1):n) {
      d <- calc_dist(df_contamina$x[i], df_contamina$y[i],
                     df_contamina$x[j], df_contamina$y[j])
      if (!is.na(d) && !is.na(bins[h]) && !is.na(bins[h+1]) && d >= bins[h] && d < bins[h+1]) {
        num <- num + (df_contamina$PM2.5[i] - df_contamina$PM2.5[j])^2
        count <- count + 1
      }
    }
  }
  
  semivars[h] <- ifelse(count > 0, num / (2 * count), NA)
}

# Gráfica del semivariograma
semivariograma_df <- data.frame(h = bins, gamma = semivars)

ggplot(semivariograma_df, aes(x = h, y = gamma)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(title = "Semivariograma Experimental", x = "Distancia (h)", y = "Semivarianza")

#Intepretación:
#El semivariograma muestra una estructura espacial en los datos, donde la semivarianza
#aumenta rápidamente con la distancia, indicando que los puntos cercanos están más correlacionados.
#La semivarianza se estabiliza alrededor de 150 (meseta), sugiriendo que la variabilidad máxima
#se alcanza a una distancia de aproximadamente 10 a 15 unidades (alcance).
#Más allá de esta distancia, los datos pierden correlación espacial.
#La ausencia de un efecto pepita significativo indica poca variabilidad a distancias muy cortas,
#lo que refleja una distribución espacial relativamente suave y definida.



