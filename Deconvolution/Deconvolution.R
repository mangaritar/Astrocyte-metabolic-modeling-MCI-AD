# ─────────────────────────────────────────────────────────────
# 1. Install and load required packages
# ─────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment", "scater", "org.Hs.eg.db", "AnnotationDbi"))

install.packages(c("reticulate", "pheatmap", "ggplot2", "reshape2", "dplyr", "openxlsx", "readxl"))

library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)
library(openxlsx)
library(readxl)
library(CDSeq)

# ─────────────────────────────────────────────────────────────
# 2. Load input data (ASTROCYTES)
# ─────────────────────────────────────────────────────────────
ref_matrix <- read.csv("hippocampus_reference_matrix.csv", row.names = 1, check.names = FALSE)
pseudo_counts <- read.csv("pseudo_counts_geneSymbol_GSE28146_FIXED.csv", sep = ";", row.names = 1, check.names = FALSE)

# ─────────────────────────────────────────────────────────────
# 3. Map ENSEMBL to gene symbols
# ─────────────────────────────────────────────────────────────
ensembl_ids <- rownames(ref_matrix)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

ref_matrix$GeneSymbol <- gene_symbols
ref_matrix <- ref_matrix[!is.na(ref_matrix$GeneSymbol) & !duplicated(ref_matrix$GeneSymbol), ]
rownames(ref_matrix) <- ref_matrix$GeneSymbol
ref_matrix$GeneSymbol <- NULL

# ─────────────────────────────────────────────────────────────
# 4. Match common genes
# ─────────────────────────────────────────────────────────────
common_genes <- intersect(rownames(ref_matrix), rownames(pseudo_counts))
ref_common <- ref_matrix[common_genes, ]
pseudo_common <- pseudo_counts[common_genes, ]

# ─────────────────────────────────────────────────────────────
# 5. Define sample conditions
# ─────────────────────────────────────────────────────────────
conditions <- c(rep("Control", 8), rep("Incipient", 7), rep("Moderate", 8), rep("Severe", 7))
sample_ids <- colnames(pseudo_counts)

stopifnot(length(sample_ids) == length(conditions))

sample_condition <- data.frame(Condition = conditions)
rownames(sample_condition) <- sample_ids

# ─────────────────────────────────────────────────────────────
# 6. Run CDSeq per condition for ASTROCYTES
# ─────────────────────────────────────────────────────────────
pseudo_counts_t <- t(pseudo_counts)
pseudo_df <- as.data.frame(pseudo_counts_t)
pseudo_df$Condition <- sample_condition[rownames(pseudo_df), "Condition"]

cdseq_results_astro <- list()

for (cond in unique(sample_condition$Condition)) {
  message("Running CDSeq for condition: ", cond)
  
  subset_data <- subset(pseudo_df, Condition == cond)
  subset_data$Condition <- NULL
  expr_matrix <- t(subset_data)
  
  cdseq_results_astro[[cond]] <- CDSeq(
    bulk_data = expr_matrix,
    cell_type_number = 6,
    mcmc_iterations = 500
  )
}

save(cdseq_results_astro, file = "cdseq_results_astro.RData")

# ─────────────────────────────────────────────────────────────
# 7. Export estimated GEPs
# ─────────────────────────────────────────────────────────────
dir.create("GEP_astrocytes", showWarnings = FALSE)

for (cond in names(cdseq_results_astro)) {
  write.csv(cdseq_results_astro[[cond]]$estGEP,
            file = file.path("GEP_astrocytes", paste0("GEP_astro_", cond, ".csv")),
            row.names = TRUE)
}

# ─────────────────────────────────────────────────────────────
# 8. Select most astrocyte-like cell type per condition
# ─────────────────────────────────────────────────────────────
astro_profile <- ref_common[, "astrocyte", drop = FALSE]
astro_profile <- astro_profile[order(rownames(astro_profile)), , drop = FALSE]

best_geps <- list()

for (cond in names(cdseq_results_astro)) {
  est_gep <- cdseq_results_astro[[cond]]$estGEP
  est_gep <- est_gep[order(rownames(est_gep)), ]
  
  shared_genes <- intersect(rownames(astro_profile), rownames(est_gep))
  
  correlations <- apply(est_gep[shared_genes, ], 2, function(x) {
    cor(x, astro_profile[shared_genes, "astrocyte"], method = "pearson")
  })
  
  best_index <- which.max(correlations)
  best_geps[[cond]] <- est_gep[, best_index, drop = FALSE]
  colnames(best_geps[[cond]]) <- cond
}

gep_best_astro <- do.call(cbind, best_geps)

# ─────────────────────────────────────────────────────────────
# 9. Load astrocyte marker genes (from your file)
# ─────────────────────────────────────────────────────────────
markers_df <- read_excel("astrocytesmarkers (1).xlsx")
marker_genes <- unique(markers_df$Gene)

# ─────────────────────────────────────────────────────────────
# 10. Subset for marker genes and generate heatmap
# ─────────────────────────────────────────────────────────────
gep_filtered <- gep_best_astro[rownames(gep_best_astro) %in% marker_genes, ]
gep_scaled <- t(scale(t(gep_filtered)))  # z-score

pheatmap(gep_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Astrocyte Marker Genes (CDSeq GEPs)",
         fontsize_row = 7,
         fontsize_col = 12)
