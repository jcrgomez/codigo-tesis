# ===================== #
#       PAQUETES        #
# ===================== #
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(future.apply)

# ===================== #
#     PARÁMETROS        #
# ===================== #
setwd("/Users/protor/Documents/Doctorado/celulas/datos")  # Ajusta esta ruta
n_samples <- 500
subset_size <- 3
min_cel_genes <- 50
max_cel_genes <- 500
min_gene_freq <- 500
block_size <- 200
output_dir <- "banzhaf_resultados"
dir.create(output_dir, showWarnings = FALSE)

# ===================== #
#   LECTURA DE DATOS    #
# ===================== #
counts <- readMM("GSE166692_sciSpace_count_matrix.mtx")
genes <- read.delim("GSE166692_sciSpace_gene_metadata.tsv", header = FALSE)
barcodes <- read.delim("GSE166692_sciSpace_cell_metadata.tsv", header = FALSE)
gene_names <- genes$V3

# ===================== #
#  FILTRADO CELULAR     #
# ===================== #
mt_genes <- grep("^mt-", gene_names, value = TRUE)
mt_rows <- which(gene_names %in% mt_genes)
total_counts <- Matrix::colSums(counts)
mt_counts <- Matrix::colSums(counts[mt_rows, ])
percent_mito <- mt_counts / total_counts * 100

healthy_cells <- which(percent_mito <= 10)
mito_cells    <- which(percent_mito > 10)
counts_healthy <- as(counts[, healthy_cells], "CsparseMatrix")
counts_mito    <- as(counts[, mito_cells], "CsparseMatrix")

# ===================== #
#  CONJUNTOS ACTIVOS    #
# ===================== #
get_active_gene_sets_sparse <- function(counts_matrix, gene_names) {
  active_sets <- vector("list", ncol(counts_matrix))
  for (j in seq_len(ncol(counts_matrix))) {
    col_start <- counts_matrix@p[j] + 1
    col_end   <- counts_matrix@p[j + 1]
    gene_indices <- counts_matrix@i[col_start:col_end] + 1
    active_sets[[j]] <- gene_names[gene_indices]
  }
  return(active_sets)
}

active_sets_healthy <- get_active_gene_sets_sparse(counts_healthy, genes$V3)
active_sets_mito    <- get_active_gene_sets_sparse(counts_mito, genes$V3)

# ===================== #
#    FILTRADO FINAL     #
# ===================== #
longitudes <- lengths(active_sets_healthy)
celdas_filtradas <- which(longitudes >= min_cel_genes & longitudes <= max_cel_genes)
active_sets_filtrado <- active_sets_healthy[celdas_filtradas]

gene_frecuencias <- table(unlist(active_sets_filtrado))
genes_filtrados <- names(gene_frecuencias[gene_frecuencias >= min_gene_freq])

# ===================== #
#     FUNCIONES CORE    #
# ===================== #
asignar_indices_genes <- function(genes) setNames(seq_along(genes), genes)

set_to_bitmask <- function(S, gen_indices, n_genes) {
  mask <- rep(0L, n_genes)
  indices <- gen_indices[S]
  mask[indices] <- 1L
  return(mask)
}

build_bitmask_index <- function(maximal_simplices, gen_indices, n_genes) {
  lapply(maximal_simplices, set_to_bitmask, gen_indices = gen_indices, n_genes = n_genes)
}

is_bitmask_subset <- function(mask_S, bitmask_index) {
  for (mask_T in bitmask_index) {
    if (length(mask_T) != length(mask_S)) next
    if (all(bitwAnd(mask_T, mask_S) == mask_S)) return(TRUE)
  }
  return(FALSE)
}

v_riesgo_bitmask <- function(S, bitmask_index, gen_indices, n_genes) {
  mask_S <- set_to_bitmask(S, gen_indices, n_genes)
  if (is_bitmask_subset(mask_S, bitmask_index)) return(0) else return(1)
}

estimacion_banzhaf_bitmask <- function(maximal_simplices, genes, 
                                       n_samples = 10000, 
                                       subset_size = 3, 
                                       seed = 123) {
  set.seed(seed)
  
  n_genes <- length(genes)
  gen_indices <- setNames(seq_along(genes), genes)
  bitmask_index <- build_bitmask_index(maximal_simplices, gen_indices, n_genes)
  
  banzhaf_scores <- numeric(n_genes)
  names(banzhaf_scores) <- genes
  
  for (i in seq_along(genes)) {
    gen <- genes[i]
    universo <- setdiff(genes, gen)
    universo_size <- length(universo)
    
    # Adaptar subset_size si hay pocos genes
    k <- min(subset_size, universo_size)
    usar_replace <- universo_size < subset_size
    
    contribuciones <- numeric(n_samples)
    
    for (j in seq_len(n_samples)) {
      if (universo_size == 0) {
        contribuciones[j] <- 0
        next
      }
      S <- sample(universo, k, replace = usar_replace)
      v1 <- v_riesgo_bitmask(S, bitmask_index, gen_indices, n_genes)
      v2 <- v_riesgo_bitmask(c(S, gen), bitmask_index, gen_indices, n_genes)
      contribuciones[j] <- v2 - v1
    }
    
    banzhaf_scores[i] <- mean(contribuciones)
  }
  
  return(sort(banzhaf_scores, decreasing = TRUE))
}



# ===================== #
#   ESTIMACIÓN BLOQUE   #
# ===================== #
procesar_bloque_banzhaf <- function(genes_bloque, maximal_simplices, 
                                    n_samples, subset_size, block_id) {
  n_genes <- length(genes_bloque)
  gen_indices <- asignar_indices_genes(genes_bloque)
  bitmask_index <- build_bitmask_index(maximal_simplices, gen_indices, n_genes)
  
  resultados <- future_lapply(genes_bloque, function(gen) {
    contribuciones <- numeric(n_samples)
    for (j in seq_len(n_samples)) {
      S <- sample(setdiff(genes_bloque, gen), subset_size)
      v1 <- v_riesgo_bitmask(S, bitmask_index, gen_indices, n_genes)
      v2 <- v_riesgo_bitmask(c(S, gen), bitmask_index, gen_indices, n_genes)
      contribuciones[j] <- v2 - v1
    }
    mean(contribuciones)
  })
  
  names(resultados) <- genes_bloque
  saveRDS(resultados, file = file.path(output_dir, paste0("banzhaf_bloque_", block_id, ".rds")))
  return(invisible(NULL))
}

# ===================== #
#   LANZAMIENTO BLOQUES #
# ===================== #
plan(multisession)  # Usa todos los núcleos disponibles
bloques <- split(genes_filtrados, ceiling(seq_along(genes_filtrados) / block_size))

for (i in seq_along(bloques)) {
  file_out <- file.path(output_dir, paste0("banzhaf_bloque_", i, ".rds"))
  if (!file.exists(file_out)) {
    cat("Procesando bloque", i, "/", length(bloques), "...\n")
    procesar_bloque_banzhaf(bloques[[i]], active_sets_filtrado, 
                            n_samples, subset_size, i)
  } else {
    cat("Bloque", i, "ya existe. Saltando.\n")
  }
}

# ===================== #
#    CSV Y VISUALIZACIÓN #
# ===================== #

# 1. Combinar bloques
archivos <- list.files(output_dir, full.names = TRUE)
resultados_todos <- do.call(c, lapply(archivos, readRDS))
banzhaf_scores <- unlist(resultados_todos)

# 2. Calcular frecuencia para los genes en resultados
gene_frecuencias_filtrado <- table(unlist(active_sets_filtrado))
frecuencia_genes <- gene_frecuencias_filtrado[names(banzhaf_scores)]

# 3. DataFrame final
df_banzhaf <- data.frame(
  gene = names(banzhaf_scores),
  banzhaf = as.numeric(banzhaf_scores),
  frecuencia = as.integer(frecuencia_genes)
)
df_banzhaf <- df_banzhaf[order(-df_banzhaf$banzhaf), ]

# 4. Guardar CSV
write.csv(df_banzhaf, file = file.path(output_dir, "banzhaf_final.csv"), row.names = FALSE)

# 5. Visualizar top-20
top20 <- head(df_banzhaf, 20)

ggplot(top20, aes(x = reorder(gene, banzhaf), y = banzhaf)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 genes por valor de Banzhaf (riesgo)",
       x = "Gen", y = "Valor de Banzhaf") +
  theme_minimal(base_size = 14)

#genes asociados con la EMT

genes_emt <- c("Cdh1", "Epcam", "Ocln", "Cldn1", "Dsg2", "Krt18", "Krt19",
               "Vim", "Fn1", "Cdh2", "Mmp2", "Mmp9", "S100a4",
               "Snai1", "Snai2", "Zeb1", "Zeb2", "Twist1", "Twist2", "Foxc2", "Tcf4")

emt_detectados <- df_banzhaf[df_banzhaf$gene %in% genes_emt, ]
emt_detectados <- emt_detectados[order(-emt_detectados$banzhaf), ]
print(emt_detectados)


#Ahora realizamos el mismo estudio pero para celulas mitocondriales

# ========================== #
# 1. Lista de genes EMT      #
# ========================== #

genes_emt <- c("Cdh1", "Epcam", "Ocln", "Cldn1", "Dsg2", "Krt18", "Krt19",
               "Vim", "Fn1", "Cdh2", "Mmp2", "Mmp9", "S100a4",
               "Snai1", "Snai2", "Zeb1", "Zeb2", "Twist1", "Twist2", "Foxc2", "Tcf4")

# ========================== #
# 2. Procesar células mitocondriales #
# ========================== #

# Filtrar por número de genes por célula
longitudes_mito <- lengths(active_sets_mito)
celdas_mito_filtradas <- which(longitudes_mito >= min_cel_genes & longitudes_mito <= max_cel_genes)
active_sets_mito_filtrado <- active_sets_mito[celdas_mito_filtradas]

# Filtrar genes por frecuencia
gene_frecuencias_mito <- table(unlist(active_sets_mito_filtrado))
genes_filtrados_mito <- names(gene_frecuencias_mito[gene_frecuencias_mito >= min_gene_freq])

# ========================== #
# 3. Calcular valor de Banzhaf (mito) #
# ========================== #

genes_emt_mito <- intersect(genes_emt, genes_filtrados_mito)
genes_emt_mito <- intersect(genes_emt, names(gene_frecuencias_mito))


banzhaf_mito <- estimacion_banzhaf_bitmask(
  maximal_simplices = active_sets_mito_filtrado,
  genes = genes_emt_mito,
  n_samples = n_samples,
  subset_size = subset_size
)

# ========================== #
# 4. Armar tabla comparativa #
# ========================== #

# Extraer datos para genes EMT en saludables
emt_saludables <- df_banzhaf[df_banzhaf$gene %in% genes_emt, ]
emt_saludables <- emt_saludables[order(emt_saludables$gene), ]

# Crear data frame para mito
frecuencias_mito <- gene_frecuencias_mito[names(banzhaf_mito)]
emt_mito <- data.frame(
  gene = names(banzhaf_mito),
  banzhaf_mito = as.numeric(banzhaf_mito),
  frecuencia_mito = as.integer(frecuencias_mito)
)
emt_mito <- emt_mito[order(emt_mito$gene), ]

# Combinar
emt_comparacion <- merge(emt_saludables, emt_mito, by = "gene", all = TRUE)
colnames(emt_comparacion) <- c("gene", "banzhaf_saludable", "frecuencia_saludable",
                               "banzhaf_mito", "frecuencia_mito")

# Mostrar tabla ordenada por diferencia de Banzhaf
emt_comparacion$diferencia <- emt_comparacion$banzhaf_mito - emt_comparacion$banzhaf_saludable
emt_comparacion <- emt_comparacion[order(-emt_comparacion$diferencia), ]
print(emt_comparacion)

#mapa de calor
emt_comparacion$tipo <- NA
emt_comparacion$tipo[emt_comparacion$gene %in% c("Cdh1", "Epcam", "Ocln", "Cldn1", "Dsg2", "Krt18", "Krt19")] <- "Epitelial"
emt_comparacion$tipo[emt_comparacion$gene %in% c("Vim", "Fn1", "Cdh2", "Mmp2", "Mmp9", "S100a4")] <- "Mesenquimal"
emt_comparacion$tipo[emt_comparacion$gene %in% c("Snai1", "Snai2", "Zeb1", "Zeb2", "Twist1", "Twist2", "Foxc2", "Tcf4")] <- "Regulador"

library(ggplot2)

ggplot(emt_comparacion, aes(x = banzhaf_saludable, y = banzhaf_mito, color = tipo, label = gene)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_abline(slope = 1, linetype = "dashed", color = "gray") +
  geom_text(nudge_y = 0.015, size = 3.5, check_overlap = TRUE) +
  labs(title = "Comparación de valores de Banzhaf en genes EMT",
       x = "Valor de Banzhaf (células saludables)",
       y = "Valor de Banzhaf (células mitocondriales)",
       color = "Tipo EMT") +
  theme_minimal(base_size = 14)

ggsave("banzhaf_comparacion_emt.png", width = 10, height = 7)

emt_comparacion$presente_en <- ifelse(is.na(emt_comparacion$banzhaf_saludable), 
                                      "Solo mitocondriales", "Ambos")

ggplot(emt_comparacion, aes(x = banzhaf_saludable, y = banzhaf_mito, 
                            color = tipo, shape = presente_en, label = gene)) +
  geom_point(size = 4, alpha = 0.9, na.rm = FALSE) +
  geom_abline(slope = 1, linetype = "dashed", color = "gray") +
  geom_text(nudge_y = 0.015, size = 3.5, check_overlap = TRUE, na.rm = FALSE) +
  scale_shape_manual(values = c("Ambos" = 16, "Solo mitocondriales" = 17)) +
  labs(title = "Comparación de valores de Banzhaf en genes EMT",
       x = "Valor de Banzhaf (células saludables)",
       y = "Valor de Banzhaf (células mitocondriales)",
       color = "Tipo EMT", shape = "Presencia del gen") +
  theme_minimal(base_size = 14)

#indice de riesgo en todas las celulas
# 1. Juntar las células saludables y mitocondriales
active_sets_todas <- c(active_sets_healthy, active_sets_mito)

# 2. Tomar una muestra aleatoria de N células
set.seed(123)  # reproducibilidad
n_muestra <- 10000  # ajusta según capacidad computacional
muestra_indices <- sample(seq_along(active_sets_todas), n_muestra)
active_sets_muestra <- active_sets_todas[muestra_indices]

# 3. Definir genes EMT
genes_emt <- c("Cdh1", "Epcam", "Ocln", "Cldn1", "Dsg2", "Krt18", "Krt19",
               "Vim", "Fn1", "Cdh2", "Mmp2", "Mmp9", "S100a4",
               "Snai1", "Snai2", "Zeb1", "Zeb2", "Twist1", "Twist2", "Foxc2", "Tcf4")

# 4. Filtrar genes EMT presentes en la muestra
genes_emt_muestra <- intersect(genes_emt, unique(unlist(active_sets_muestra)))

# 5. Calcular el valor de Banzhaf solo sobre estos genes
banzhaf_emt_muestra <- estimacion_banzhaf_bitmask(
  maximal_simplices = active_sets_muestra,
  genes = genes_emt_muestra,
  n_samples = 500,       # ajustable
  subset_size = 3
)

# 6. Mostrar resultados
banzhaf_emt_muestra <- sort(banzhaf_emt_muestra, decreasing = TRUE)
print(banzhaf_emt_muestra)
