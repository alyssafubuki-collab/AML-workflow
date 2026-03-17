library(Seurat)

cat("Loading reference...\n")
ref <- readRDS("/data/chellal/aml_project/AML_ref_seurat_fixed2.rds")

cat("Original dimensions (genes x cells):\n")
print(dim(ref))

# ======================
# Tirage aléatoire 20k cellules
# ======================
set.seed(123)

all_cells <- colnames(ref)
sampled_cells <- sample(all_cells, size = 20000, replace = FALSE)

ref_20k <- subset(ref, cells = sampled_cells)

cat("Subset dimensions:\n")
print(dim(ref_20k))

# Sauvegarde
saveRDS(ref_20k, file = "/data/chellal/aml_project/AML_ref_fixed2_20k_sampled.rds")

cat("DONE\n")
