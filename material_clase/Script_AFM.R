# =========================================================
# Tema:AFM
# Curso: Analítica de Datos
# Prof: Johann A. Ospina
# =========================================================

# Conjunto de datos


texto_datos <- '
Vino	Label	Soil	Odor.Intensity.before.shaking	Aroma.quality.before.shaking	Fruity.before.shaking	Flower.before.shaking	Spice.before.shaking	Visual.intensity	Nuance	Surface.feeling	Odor.Intensity	Quality.of.odour	Fruity	Flower	Spice	Plante	Phenolic	Aroma.intensity	Aroma.persistency	Aroma.quality	Attack.intensity	Acidity	Astringency	Alcohol	Balance	Smooth	Bitterness	Intensity	Harmony	Overall.quality	Typical
2EL	Saumur	Env1	3.074	3	2.714	2.28	1.96	4.321	4	3.269	3.407	3.308	2.885	2.32	1.84	2	1.65	3.259	2.963	3.2	2.963	2.107	2.429	2.5	3.25	2.731	1.926	2.857	3.143	3.393	3.25
1CHA	Saumur	Env1	2.964	2.821	2.375	2.28	1.68	3.222	3	2.808	3.37	3	2.56	2.44	1.739	2	1.381	2.962	2.808	2.926	3.036	2.107	2.179	2.654	2.926	2.5	1.926	2.893	2.964	3.214	3.036
1FON	Bourgueuil	Env1	2.857	2.929	2.56	1.96	2.077	3.536	3.393	3	3.25	2.929	2.769	2.192	2.25	1.75	1.25	3.077	2.8	3.077	3.222	2.179	2.25	2.643	3.321	2.679	2	3.074	3.143	3.536	3.179
1VAU	Chinon	Env2	2.808	2.593	2.417	1.913	2.16	2.893	2.786	2.538	3.16	2.88	2.391	2.083	2.167	2.304	1.476	2.542	2.583	2.478	2.704	3.179	2.185	2.5	2.333	1.68	1.963	2.462	2.038	2.464	2.25
1DAM	Saumur	Reference	3.607	3.429	3.154	2.154	2.04	4.393	4.036	3.385	3.536	3.36	3.16	2.231	2.148	1.762	1.6	3.615	3.296	3.462	3.464	2.571	2.536	2.786	3.464	3.036	2.071	3.643	3.643	3.741	3.444
2BOU	Bourgueuil	Reference	2.857	3.111	2.577	2.04	2.077	4.464	4.259	3.407	3.179	3.385	2.8	2.24	2.148	1.75	1.476	3.214	3.148	3.321	3.286	2.393	2.643	2.857	3.286	2.857	2.179	3.464	3.5	3.643	3.393
1BOI	Bourgueuil	Reference	3.214	3.222	2.962	2.115	2.04	4.143	3.929	3.25	3.429	3.5	3.038	2.2	2.385	1.826	1.476	3.25	3.222	3.385	3.393	2.607	2.607	2.778	3.464	2.857	1.929	3.643	3.556	3.714	3.357
3EL	Saumur	Env1	3.12	2.852	2.5	2.2	2.185	4.214	3.857	3.077	3.654	3.077	2.52	2.32	2.444	2.08	1.905	3.28	3.16	2.962	3.25	2.179	2.63	2.778	3.179	2.786	2	3.321	3.296	3.393	3.071
DOM1	Chinon	Env1	2.857	2.815	2.808	1.923	2.074	4.037	3.893	3.28	3.357	3.346	3	2.04	2.125	1.875	1.524	3.148	2.893	3.308	3.286	2.286	2.407	2.741	3.143	2.821	1.964	3.148	3.286	3.2	3.5
1TUR	Saumur	Env2	2.893	3	2.571	1.846	1.68	3.704	3.407	3.111	3.222	3.259	2.926	2.04	2.042	2	1.773	3.077	2.704	2.778	2.893	2.357	2.25	2.704	3.214	2.5	2.185	2.857	2.963	3.179	2.964
4EL	Saumur	Env2	3.25	3.286	2.714	1.926	1.962	3.857	3.643	3.259	3.607	3.385	2.889	2.115	2.16	1.955	1.571	3.286	3.036	3.222	3.321	2.429	2.571	2.893	3.192	2.857	2.214	3.357	3.071	3.571	3.5
PER1	Saumur	Env2	3.393	3.179	2.769	2.038	1.92	4.714	4.5	3.321	3.481	3.385	2.962	2	2.2	2.042	1.545	3.321	3.071	3.143	3.357	2.429	2.607	2.821	3.107	2.889	2.037	3.25	3.393	3.148	3.556
2DAM	Saumur	Reference	3.179	3.286	2.778	2.231	1.76	4.222	4.071	3.462	3.481	3.423	2.963	2.269	2.154	1.957	1.571	3.481	3.259	3.269	3.393	2.286	2.5	2.821	3.5	3.286	2	3.407	3.643	3.571	3.929
1POY	Saumur	Reference	3.071	3.107	2.731	2.12	1.8	4.714	4.536	3.429	3.357	3.444	2.885	2.12	2.346	1.826	1.55	3.269	3.08	3.192	3.519	2.111	2.536	2.778	3.444	3.231	2.071	3.667	3.786	3.929	3.481
1ING	Bourgueuil	Env1	3.107	3.143	2.846	2.185	1.962	4.071	3.893	3.462	3.357	3.37	2.846	2.24	2.28	1.75	1.524	3.333	3.037	3.37	3.185	2.286	2.643	2.929	3.286	2.821	2.107	3.321	3.296	3.643	3.296
1BEN	Bourgueuil	Reference	2.929	3.179	2.852	2	2.037	3.889	3.429	3.143	3.286	3.308	3.115	2.269	2	1.917	1.4	3.04	2.96	3.2	3.393	2.393	2.357	2.704	3.321	3	2	3.214	3.214	3.75	3.571
2BEA	Chinon	Reference	3.036	3.179	3.037	2.231	1.667	3.786	3.607	3.357	3.444	3.5	3.185	2.16	2.24	1.913	1.75	3.52	3.296	3.462	3.071	2.571	2.321	2.929	3.333	2.821	2.143	3.321	3.25	3.536	3.269
1ROC	Chinon	Env2	3.071	2.926	2.741	2	1.88	3.679	3.393	3.192	3.37	3.36	2.963	2.308	1.917	2	1.429	3.25	2.92	2.88	3.071	2.393	2.321	2.821	3.143	2.607	2.143	3.037	3.074	3.464	3.444
2ING	Bourgueuil	Env1	2.643	2.786	2.536	1.889	1.808	2.607	2.536	2.444	2.889	2.8	2.5	1.962	2.111	2.08	1.318	2.68	2.308	2.556	2.179	2.25	1.964	2.25	2.464	1.821	1.679	2.179	2.107	2.37	2.321
T1	Saumur	Env4	3.696	3.192	2.833	1.826	2.385	4.321	4	3.333	3.737	3.08	2.833	1.773	2.44	2.292	1.571	3.437	2.958	2.6	2.963	2.407	2.643	2.963	2.571	2.071	2.222	3.037	2.741	2.643	2.571
T2	Saumur	Env4	3.708	2.926	2.52	2.04	2.667	4.321	4.107	3.259	3.727	2.885	2.6	2.083	2.609	2.174	1.65	3.095	3.136	2.545	3.333	2.571	2.667	2.704	2.769	2.308	2.667	3.333	3	2.852	2.75
'

# Cargar datos
vino <- read.table(text = texto_datos, header = TRUE, sep = "\t", row.names = 1)


# Grupos
lista_grupos <- list(
  Olor_antes = as.matrix(vino[, c("Odor.Intensity.before.shaking", "Aroma.quality.before.shaking",
                                  "Fruity.before.shaking", "Flower.before.shaking", "Spice.before.shaking")]),
  Visual     = as.matrix(vino[, c("Visual.intensity", "Nuance", "Surface.feeling")]),
  Olor_despues = as.matrix(vino[, c("Odor.Intensity", "Quality.of.odour", "Fruity", "Flower", "Spice",
                                    "Plante", "Phenolic", "Aroma.intensity", "Aroma.persistency", "Aroma.quality")]),
  Gusto      = as.matrix(vino[, c("Attack.intensity", "Acidity", "Astringency", "Alcohol", "Balance",
                                  "Smooth", "Bitterness", "Intensity", "Harmony")])
)

I <- nrow(vino)             # numero de individuos (21 vinos)
K <- length(lista_grupos)   # numero de grupos activos (4)


# Estandarizar grupos
estandarizar <- function(X) scale(X, center = TRUE, scale = TRUE)

# ACP con descomposción en valores singugales
acp_svd <- function(X) {
  n <- nrow(X)
  Z <- X / sqrt(n - 1)
  desc <- svd(Z)
  list(valores = desc$d^2, ejes = desc$v, U = desc$u, d = desc$d)
}

#Análisis por grupos
lambda1 <- numeric(K)
grupos_std <- vector("list", K)
for (k in seq_len(K)) {
  Xk <- estandarizar(lista_grupos[[k]])
  grupos_std[[k]] <- Xk
  lambda1[k] <- acp_svd(Xk)$valores[1]
}
names(lambda1) <- names(lista_grupos)


# Ponderación
alpha <- 1 / lambda1
grupos_pond <- Map(function(Xk, a) Xk * sqrt(a), grupos_std, alpha)

#Análisis general ponderado
X_global <- do.call(cbind, grupos_pond)
res_global <- acp_svd(X_global)

q <- 2  # seleccionar las dimensiones como en ACP
valores_propios <- res_global$valores
inercia_pct <- 100 * valores_propios / sum(valores_propios)


U_ejes <- res_global$ejes[, 1:q, drop = FALSE]
F_global <- X_global %*% U_ejes
colnames(F_global) <- paste0("Eje", 1:q)
rownames(F_global) <- rownames(vino)



#Individuos po grupo

indices_grupo <- list(); inicio <- 1
for (k in seq_len(K)) {
  fin <- inicio + ncol(grupos_pond[[k]]) - 1
  indices_grupo[[k]] <- inicio:fin
  inicio <- fin + 1
}

F_parciales <- vector("list", K)
for (k in seq_len(K)) {
  bloque_u <- U_ejes[indices_grupo[[k]], , drop = FALSE]
  F_parciales[[k]] <- K * (grupos_pond[[k]] %*% bloque_u)
  colnames(F_parciales[[k]]) <- paste0("Eje", 1:q)
  rownames(F_parciales[[k]]) <- rownames(vino)
}
names(F_parciales) <- names(lista_grupos)


#Cos2 p contri de individuos
suma_cuadrados_i <- rowSums(F_global^2)
cos2_individuos <- F_global^2 / suma_cuadrados_i
p_i <- 1 / (I - 1)


# Variables y contribuciones
X_std_todas <- do.call(cbind, grupos_std)
correlaciones <- cor(X_std_todas, F_global)

grupo_de_variable <- rep(names(lista_grupos), sapply(grupos_std, ncol))
alpha_por_variable <- rep(alpha, sapply(grupos_std, ncol))
ctr_variables <- sweep(correlaciones^2 * alpha_por_variable, 2, valores_propios[1:q], "/") * 100


Lg <- aggregate(correlaciones^2 * alpha_por_variable,
                by = list(grupo = grupo_de_variable), FUN = sum)


RV <- function(Xk, Xkp) {
  Wk  <- Xk  %*% t(Xk)
  Wkp <- Xkp %*% t(Xkp)
  sum(diag(Wk %*% Wkp)) / sqrt(sum(diag(Wk %*% Wk)) * sum(diag(Wkp %*% Wkp)))
}
RV_matriz <- matrix(NA, K, K, dimnames = list(names(lista_grupos), names(lista_grupos)))
for (k in seq_len(K)) for (kp in seq_len(K)) RV_matriz[k, kp] <- RV(grupos_std[[k]], grupos_std[[kp]])


# Gráf de indiv
colores <- c(Saumur = "firebrick", Bourgueuil = "steelblue", Chinon = "darkgreen")

plot(F_global, type = "n", xlab = paste0("Eje 1 (", round(inercia_pct[1],1), "%)"),
     ylab = paste0("Eje 2 (", round(inercia_pct[2],1), "%)"),
     main = "AFM - Vinos (nube global de individuos)")
abline(h = 0, v = 0, lty = 2, col = "grey70")
text(F_global, labels = rownames(F_global), col = colores[vino$Label], cex = 0.8, font = 2)
legend("topleft", legend = names(colores), col = colores, pch = 19, bty = "n", cex = 0.8)

