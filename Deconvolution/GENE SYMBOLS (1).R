if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")

# Instalar si no lo tienes
if (!require("hgu133plus2.db")) BiocManager::install("hgu133plus2.db")

# Cargar paquete
library(hgu133plus2.db)
library(dplyr)
library(readr)

getwd()

pseudo_counts <- read.csv("pseudo_counts_GSE28146.csv", row.names = 1, check.names = FALSE)

# Mapear los IDs de sonda a símbolos de genes
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = rownames(pseudo_counts),
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Añadir la columna GeneSymbol a la matriz
pseudo_counts$GeneSymbol <- gene_symbols

# Eliminar filas sin asignación de símbolo
pseudo_counts_clean <- pseudo_counts[!is.na(pseudo_counts$GeneSymbol), ]


library(dplyr)

# Rehacer la agrupación desde pseudo_counts_clean (que aún tiene GeneSymbol)
pseudo_counts_gene <- pseudo_counts_clean %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

# Verifica que los nombres de genes están presentes
head(pseudo_counts_gene$GeneSymbol)

# Asignar correctamente los nombres de fila
pseudo_counts_matrix <- as.data.frame(pseudo_counts_gene)
rownames(pseudo_counts_matrix) <- pseudo_counts_matrix$GeneSymbol
pseudo_counts_matrix <- pseudo_counts_matrix[, -1]

# Guardar nuevamente
write.csv(pseudo_counts_matrix, "pseudo_counts_geneSymbol_GSE28146_FIXED.csv")
