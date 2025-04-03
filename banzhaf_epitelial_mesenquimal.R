
# ============================================== #
#    Banzhaf para Células Epiteliales y Mesenquimales
# ============================================== #

# Paquetes
library(Matrix)
library(dplyr)

# Asumimos que los datos ya están leídos:
# counts, genes, barcodes, gene_names, etc.
# Y que tienes cargada la función estimacion_banzhaf_bitmask()

# ------------------------- #
# Definir genes EMT
# ------------------------- #
genes_emt <- c("Cdh1", "Epcam", "Ocln", "Cldn1", "Dsg2", "Krt18", "Krt19",
               "Vim", "Fn1", "Cdh2", "Mmp2", "Mmp9", "S100a4",
               "Snai1", "Snai2", "Zeb1", "Zeb2", "Twist1", "Twist2", "Foxc2", "Tcf4")

# ------------------------- #
# Filtrar por tipo celular usando barcodes$V20
# ------------------------- #
epiteliales_idx <- which(barcodes$V20 == "Epithelial Cells")
mesenquimales_idx <- which(barcodes$V20 == "Connective Tissue Progenitors")

# ------------------------- #
# Obtener conjuntos activos
# ------------------------- #
library(Matrix)

counts_epi <- as(counts[, epiteliales_idx], "CsparseMatrix")
active_sets_epi <- get_active_gene_sets_sparse(counts_epi, genes$V3)

counts_mes <- as(counts[, mesenquimales_idx], "CsparseMatrix")
active_sets_mes <- get_active_gene_sets_sparse(counts_mes, genes$V3)

# ------------------------- #
# Filtrar genes EMT presentes en cada grupo
# ------------------------- #
genes_epi <- intersect(genes_emt, unique(unlist(active_sets_epi)))
genes_mes <- intersect(genes_emt, unique(unlist(active_sets_mes)))

# ------------------------- #
# Parámetros
# ------------------------- #
n_samples <- 500
subset_size <- 3

# ------------------------- #
# Calcular valores de Banzhaf
# ------------------------- #
banzhaf_epi <- estimacion_banzhaf_bitmask(
  maximal_simplices = active_sets_epi,
  genes = genes_epi,
  n_samples = n_samples,
  subset_size = subset_size
)

banzhaf_mes <- estimacion_banzhaf_bitmask(
  maximal_simplices = active_sets_mes,
  genes = genes_mes,
  n_samples = n_samples,
  subset_size = subset_size
)

# ------------------------- #
# Armar tabla comparativa
# ------------------------- #
emt_epi <- data.frame(
  gene = names(banzhaf_epi),
  banzhaf_epi = as.numeric(banzhaf_epi)
)

emt_mes <- data.frame(
  gene = names(banzhaf_mes),
  banzhaf_mes = as.numeric(banzhaf_mes)
)

emt_comparacion_tipo <- merge(emt_epi, emt_mes, by = "gene", all = TRUE)
emt_comparacion_tipo$diferencia <- emt_comparacion_tipo$banzhaf_mes - emt_comparacion_tipo$banzhaf_epi
emt_comparacion_tipo <- emt_comparacion_tipo[order(-emt_comparacion_tipo$diferencia), ]

# Guardar CSV
write.csv(emt_comparacion_tipo, "banzhaf_epitelial_vs_mesenquimal.csv", row.names = FALSE)

# Mostrar
print(emt_comparacion_tipo)


library(ggplot2)

ggplot(emt_comparacion_tipo, aes(x = banzhaf_epi, y = banzhaf_mes, label = gene)) +
  geom_point(color = "steelblue", size = 3.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_text(nudge_y = 0.015, size = 4, check_overlap = TRUE) +
  labs(
    title = "Comparación entre tipos celulares",
    x = "Índice de riesgo en células epiteliales",
    y = "Índice de riesgo en células mesenquimales"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )
ggsave("banzhaf_epi_vs_mesenquimal.png", width = 10, height = 7, dpi = 300, bg = "white")


