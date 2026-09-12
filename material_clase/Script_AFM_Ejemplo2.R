# =====================================================================
# ANALISIS FACTORIAL MULTIPLE (AFM / MFA) IMPLEMENTADO DESDE CERO EN R
# =====================================================================

set.seed(42)
options(digits = 4)

#   X = [ X1 | X2 | X3 ],   Xk en R^(I x Jk)

I <- 18
materias <- c("Matematicas", "Lengua", "Ciencias", "Ingles")
periodos <- c("Periodo_1", "Periodo_2", "Periodo_3")
K <- length(periodos)

estudiantes <- sprintf("Est_%02d", 1:I)

habilidad_base <- rnorm(I, mean = 3.5, sd = 0.6)
tendencia <- rnorm(I, mean = 0.15, sd = 0.10)

grupos <- list()  # lista de data.frames, uno por periodo (I x 4)
for (k in seq_along(periodos)) {
  per <- periodos[k]
  cols <- list()
  for (mat in materias) {
    ruido <- rnorm(I, 0, 0.35)
    efecto_materia <- rnorm(1, 0, 0.25)
    nota <- habilidad_base + (k - 1) * tendencia + efecto_materia + ruido
    nota <- pmin(pmax(nota, 0), 5)
    cols[[paste0(mat, "_", per)]] <- nota
  }
  grupos[[per]] <- as.data.frame(cols, row.names = estudiantes)
}

X <- do.call(cbind, grupos)
colnames(X) <- unlist(lapply(grupos, colnames))  # quita el prefijo de cbind
rownames(X) <- estudiantes
J <- ncol(X)


# Pesos de los individuos: D = (1/I) * Identidad (uniforme, ACP clasico)
D <- diag(1 / I, I, I)


estandarizar <- function(Z) {
  # Centra y reduce (media 0, varianza poblacional 1) cada columna
  Zm <- scale(Z, center = TRUE, scale = FALSE)
  sdev <- apply(Z, 2, function(col) sqrt(mean((col - mean(col))^2)))
  sweep(Zm, 2, sdev, "/")
}

acp_manual <- function(Z, D) {
  # ACP clasico via el triplete (Z, M = I, D).
  # Resuelve   Z' D Z u = lambda u   (matriz simetrica -> eigen)
  Z <- as.matrix(Z)
  C <- t(Z) %*% D %*% Z                 # covarianza ponderada (M = I)
  ee <- eigen(C, symmetric = TRUE)      # ya ordena de mayor a menor
  valores <- ee$values
  vectores <- ee$vectors
  # Convencion de signo: la mayor carga en valor absoluto es positiva
  for (s in seq_len(ncol(vectores))) {
    j_max <- which.max(abs(vectores[, s]))
    if (vectores[j_max, s] < 0) vectores[, s] <- -vectores[, s]
  }
  componentes <- Z %*% vectores
  list(valores = valores, vectores = vectores, componentes = componentes)
}


Z_grupos <- list()      
lambda1_grupos <- list() # primer valor propio de cada analisis separado

for (per in periodos) {
  Zk <- estandarizar(grupos[[per]])
  Z_grupos[[per]] <- Zk
  res_k <- acp_manual(Zk, D)
  lambda1_grupos[[per]] <- res_k$valores[1]
  cat(sprintf("  %s: lambda1^(k) = %.4f   (valores propios: %s)\n",
              per, res_k$valores[1],
              paste(round(res_k$valores, 3), collapse = ", ")))
}


#Ponderación

alpha <- list()
for (per in periodos) {
  alpha[[per]] <- 1 / lambda1_grupos[[per]]
  cat(sprintf("  alpha_%s = 1/%.4f = %.4f  -> factor sqrt(alpha) = %.4f\n",
              per, lambda1_grupos[[per]], alpha[[per]], sqrt(alpha[[per]])))
}

# Reescalamos cada tabla: Yk = Zk / sqrt(lambda1^(k))
Y_grupos <- list()
for (per in periodos) {
  Y_grupos[[per]] <- Z_grupos[[per]] / sqrt(lambda1_grupos[[per]])
}
Y <- do.call(cbind, Y_grupos)


#Siempre debemos verificar
for (per in periodos) {
  res_v <- acp_manual(Y_grupos[[per]], D)
  cat(sprintf("  %s: primer valor propio tras ponderar = %.6f\n",
              per, res_v$valores[1]))
}


#Análisis global ponderado
res_global <- acp_manual(Y, D)
valores_globales <- res_global$valores
vectores_globales <- res_global$vectores
Fmat <- res_global$componentes

inercia_total <- sum(valores_globales)
pct_var <- 100 * valores_globales / inercia_total
n_ejes <- 6


tabla_valores <- data.frame(
  Eje = paste0("Dim", 1:n_ejes),
  Valor_propio = round(valores_globales[1:n_ejes], 3),
  Pct_inercia = round(pct_var[1:n_ejes], 2),
  Pct_acumulado = round(cumsum(pct_var[1:n_ejes]), 2)
)
print(tabla_valores, row.names = FALSE)

q <- 2  # con 2 dimensiones se acumula casi el 90%

F_df <- as.data.frame(Fmat[, 1:q])
colnames(F_df) <- paste0("Dim", 1:q)
rownames(F_df) <- estudiantes
print(round(head(F_df), 3))


# Individuos

# F^k_s(i) = K * Yk(i,:) %*% u_s|k   (bloque de cargas del grupo k)
bloques_col <- list()
inicio <- 1
for (per in periodos) {
  n_j <- ncol(Y_grupos[[per]])
  bloques_col[[per]] <- inicio:(inicio + n_j - 1)
  inicio <- inicio + n_j
}

F_parciales <- list()
for (per in periodos) {
  u_bloque <- vectores_globales[bloques_col[[per]], 1:q, drop = FALSE]
  Yk_mat <- as.matrix(Y_grupos[[per]])
  F_parciales[[per]] <- K * (Yk_mat %*% u_bloque)
}

# Propiedad de baricentro: F_s(i) = (1/K) * sum_k F^k_s(i)
promedio_parciales <- Reduce(`+`, F_parciales) / K
diferencia_max <- max(abs(promedio_parciales - Fmat[, 1:q]))


for (per in periodos) {
  cat(sprintf("  %s: [%.3f, %.3f]\n", per,
              F_parciales[[per]][1, 1], F_parciales[[per]][1, 2]))
}


#Representación variables y grupos

# correlacion entre variables ORIGINALES
#     (estandarizadas, sin reponderar) y los ejes globales F_s
Z_todas <- do.call(cbind, Z_grupos)
correlaciones <- matrix(NA, nrow = ncol(Z_todas), ncol = q,
                        dimnames = list(colnames(Z_todas), paste0("Dim", 1:q)))
for (var in colnames(Z_todas)) {
  for (s in 1:q) {
    correlaciones[var, s] <- cor(Z_todas[, var], Fmat[, s])
  }
}

corrplot(correlaciones, method = "number")

# Contribuciones de cada variable a cada eje

contribuciones <- matrix(NA, nrow = ncol(Y), ncol = q,
                         dimnames = list(colnames(Y), paste0("Dim", 1:q)))
for (s in 1:q) {
  contribuciones[, s] <- 100 * (vectores_globales[, s]^2)
}


#Coeficiente Lg
Lg <- matrix(NA, nrow = K, ncol = q,
             dimnames = list(periodos, paste0("Dim", 1:q)))
for (per in periodos) {
  Lg[per, ] <- colSums(contribuciones[bloques_col[[per]], , drop = FALSE])
}

#Coeficiente RV entre pares de grupos (similitud estructural)
coef_RV <- function(Zk, Zk2) {
  Zk <- as.matrix(Zk); Zk2 <- as.matrix(Zk2)
  Wk <- Zk %*% t(Zk)
  Wk2 <- Zk2 %*% t(Zk2)
  num <- sum(diag(Wk %*% Wk2))
  den <- sqrt(sum(diag(Wk %*% Wk)) * sum(diag(Wk2 %*% Wk2)))
  num / den
}

RV <- matrix(NA, nrow = K, ncol = K, dimnames = list(periodos, periodos))
for (p1 in periodos) {
  for (p2 in periodos) {
    RV[p1, p2] <- coef_RV(Z_grupos[[p1]], Z_grupos[[p2]])
  }
}

corrplot(RV, method = "number")

# Visualización

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))


barplot(pct_var[1:n_ejes], names.arg = paste0("D", 1:n_ejes),
        col = "#2b6cb0", main = "Paso 3: Valores propios (analisis global)",
        ylab = "% de inercia explicada", xlab = "Dimension")
lines(x = seq(0.7, by = 1.2, length.out = n_ejes), y = pct_var[1:n_ejes],
      type = "o", pch = 16)


colores <- c(Periodo_1 = "#c53030", Periodo_2 = "#2f855a", Periodo_3 = "#2b6cb0")
plot(Fmat[, 1], Fmat[, 2], col = "gray50", pch = 16,
     xlab = "Dim1", ylab = "Dim2",
     main = "Paso 4: Baricentro (Est_01)")
abline(h = 0, v = 0, col = "lightgray")
idx <- 1
for (per in periodos) {
  points(F_parciales[[per]][idx, 1], F_parciales[[per]][idx, 2],
         col = colores[per], pch = 17, cex = 1.6)
  segments(Fmat[idx, 1], Fmat[idx, 2],
           F_parciales[[per]][idx, 1], F_parciales[[per]][idx, 2],
           col = colores[per])
}
points(Fmat[idx, 1], Fmat[idx, 2], pch = 8, cex = 1.8, lwd = 2)
legend("topleft", legend = c("Global", periodos), pch = c(8, 17, 17, 17),
       col = c("black", colores), cex = 0.7, bty = "n")


plot(0, 0, type = "n", xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1),
     xlab = "Dim1", ylab = "Dim2", asp = 1,
     main = "Paso 5: Circulo de correlacion")
symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "gray")
abline(h = 0, v = 0, col = "lightgray")
for (per in periodos) {
  vars_per <- grep(per, rownames(correlaciones), value = TRUE)
  for (v in vars_per) {
    arrows(0, 0, correlaciones[v, 1], correlaciones[v, 2],
           col = colores[per], length = 0.08)
    text(correlaciones[v, 1] * 1.1, correlaciones[v, 2] * 1.1,
         labels = strsplit(v, "_")[[1]][1], col = colores[per], cex = 0.6)
  }
}


plot(Lg[, 1], Lg[, 2], col = colores[periodos], pch = 16, cex = 2.2,
     xlab = "Dim1 (Lg)", ylab = "Dim2 (Lg)",
     xlim = c(min(Lg[,1]) - 3, max(Lg[,1]) + 3),
     ylim = c(min(Lg[,2]) - 3, max(Lg[,2]) + 3),
     main = "Paso 5: Representacion de los grupos (Lg)")
text(Lg[, 1], Lg[, 2], labels = periodos, pos = 4, cex = 0.8)
abline(h = 0, v = 0, col = "lightgray")


