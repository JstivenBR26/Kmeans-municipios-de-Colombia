# =============================================================================
# Simulación K-Means · municipios SINTÉTICOS de Colombia
# Réplica en R de las 7 etapas de la presentación (KMeans_Presentacion.html)
#
# AVISO: los datos son inventados (generados con generar_datos.py, semilla 2026).
#        No son datos reales de Colombia ni deben leerse como hallazgos empíricos.
#
# Uso: abrir R en la carpeta del repositorio y ejecutar source("R/simulacion_kmeans.R")
#      Las figuras quedan en la carpeta "figuras/".
# Paquetes: solo "cluster" (silueta; viene con R). El K-medias difuso se programa en R base.
# =============================================================================

if (!requireNamespace("cluster", quietly = TRUE)) install.packages("cluster", repos = "https://cloud.r-project.org")
library(cluster)

dir.create("figuras", showWarnings = FALSE)
pal <- c("#6ea8ff", "#ffb703", "#3ddc97", "#ff6b9d", "#b388ff")

# ---- Etapa 1: el problema ----------------------------------------------------
# Agrupar 1.100 municipios (sintéticos) según su perfil agropecuario/rural,
# para pensar políticas diferenciadas.
datos    <- read.csv("data/municipios_sinteticos.csv", stringsAsFactors = FALSE)
contorno <- read.csv("data/poly_colombia.csv")   # silueta esquemática, no cartografía oficial
n <- nrow(datos)

# ---- Etapa 2: variables ------------------------------------------------------
vars <- c("pct_va_agro", "pct_area_agro", "va_agro_pc", "pob_rural",
          "rendimiento", "conectividad", "nbi")
print(summary(datos[vars]))
print(round(cor(datos[vars]), 2))               # p. ej., conectividad vs nbi ~ -0,81

# ---- Etapa 3: preparación (estandarización z-score) --------------------------
# Sin estandarizar, la variable con mayor escala (pob_rural) dominaría la distancia euclidiana.
Z <- scale(datos[vars])                          # media 0, desviación muestral 1

# Proyección a 2 componentes principales: SOLO para visualizar
pca <- prcomp(Z)
cat("Varianza explicada PC1-PC2:", round(summary(pca)$importance[2, 1:2], 3), "\n")

mapa <- function(col, main, cex = 0.6) {
  plot(contorno$lon, contorno$lat, type = "n", asp = 1, xlab = "lon", ylab = "lat", main = main)
  polygon(contorno$lon, contorno$lat, col = "grey92", border = "grey50")
  points(datos$lon, datos$lat, pch = 19, cex = cex, col = col)
}

# ---- Etapa 4: K-Means con K = 4 ----------------------------------------------
set.seed(2026)
km4 <- kmeans(Z, centers = 4, nstart = 25, iter.max = 100)   # 25 inicios aleatorios
datos$grupo <- km4$cluster
print(table(datos$grupo))

png("figuras/01_kmeans_k4.png", 1400, 640, res = 130)
par(mfrow = c(1, 2))
plot(pca$x[, 1:2], col = pal[km4$cluster], pch = 19, cex = 0.6,
     main = "K-Means (K = 4) · plano PC1-PC2")
mapa(pal[km4$cluster], "Mapa esquemático (ubicación ilustrativa)")
dev.off()

# ---- Etapa 5: ¿cuántos grupos? -----------------------------------------------
# (a) Coeficiente de fusión: nivel de unión de un método jerárquico auxiliar (Ward).
#     El libro no da una fórmula cerrada; es la altura de fusión que reporta el software.
hc  <- hclust(dist(Z), method = "ward.D2")
fus <- sapply(1:12, function(k) hc$height[n - k])   # altura al pasar de k+1 a k grupos
fus_norm <- fus / max(fus)
print(round(fus_norm, 3))
salto <- fus_norm[3:1] - fus_norm[4:2]                 # CF(K-1) - CF(K) para K = 4, 3, 2
names(salto) <- paste0("K=", 4:2)
print(round(salto, 3))

# (b) Inercia (codo) y (c) silueta media
Ks <- 2:8
inercia <- sapply(Ks, function(k) { set.seed(2026); kmeans(Z, k, nstart = 25, iter.max = 100)$tot.withinss })
D <- dist(Z)
silue <- sapply(Ks, function(k) { set.seed(2026)
  mean(silhouette(kmeans(Z, k, nstart = 25, iter.max = 100)$cluster, D)[, 3]) })

png("figuras/02_elegir_k.png", 1500, 480, res = 130)
par(mfrow = c(1, 3))
plot(1:12, fus_norm, type = "b", pch = 19, xlab = "Número de grupos", ylab = "Coef. de fusión (normalizado)",
     main = "Coeficiente de fusión (Ward)")
plot(Ks, inercia, type = "b", pch = 19, xlab = "K", ylab = "Inercia (WCSS)", main = "Codo")
plot(Ks, silue, type = "b", pch = 19, xlab = "K", ylab = "Silueta media", main = "Silueta")
dev.off()
# Ningún criterio decide solo: K = 4 se escoge además por interpretabilidad.

# ---- Etapa 6: interpretación de los grupos -----------------------------------
perfil <- aggregate(datos[vars], by = list(grupo = datos$grupo), FUN = mean)
perfil$n <- as.vector(table(datos$grupo))
print(round(perfil, 2))   # centroides en unidades originales; el nombre de cada grupo lo pone el analista

# ---- Etapa 7: el problema de los conjuntos difusos ---------------------------
# Fuzzy C-Means en R base (sin paquetes extra). Criterio: min sum f_ik^m * d^2(x_i, mu_k).
# Centros: mu_k = sum f_ik^m x_i / sum f_ik^m. Actualización de pertenencias con exponente m = 2
# (forma estándar de Bezdek; el libro transcribe el criterio y los centros, no esta actualización).
fcm <- function(X, C0, m = 2, iter = 200, tol = 1e-6) {
  C <- C0
  for (it in seq_len(iter)) {
    d2 <- sapply(seq_len(nrow(C)), function(k) rowSums((X - matrix(C[k, ], nrow(X), ncol(X), byrow = TRUE))^2))
    d2[d2 < 1e-12] <- 1e-12
    U <- 1 / (d2^(1 / (m - 1)) * rowSums(d2^(-1 / (m - 1))))
    W <- U^m
    Cn <- t(W) %*% X / colSums(W)
    if (max(abs(Cn - C)) < tol) { C <- Cn; break }
    C <- Cn
  }
  list(membership = U, centers = C, cluster = max.col(U))
}
fz <- fcm(Z, km4$centers, m = 2)                      # se inicia en los centros de K-Means
U  <- fz$membership
maxU <- apply(U, 1, max)
cat("Municipios con pertenencia máxima < 0,6:", sum(maxU < 0.6), "de", n, "\n")

frontera <- order(apply(U, 1, function(u) { s <- sort(u, decreasing = TRUE); s[1] - s[2] }))[1:5]
print(data.frame(id = datos$id[frontera], round(U[frontera, ], 2)))

png("figuras/03_difuso.png", 1400, 640, res = 130)
par(mfrow = c(1, 2))
mapa(pal[fz$cluster], "Fuzzy C-Means: grupo de mayor pertenencia")
mapa(ifelse(maxU < 0.6, "black", "grey70"), "En negro: pertenencia máxima < 0,6 (frontera)", cex = 0.7)
dev.off()

cat("Listo. Figuras en la carpeta 'figuras/'.\n")

