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
setwd("/Users/Documents/...")  # Ajusta esta ruta
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
