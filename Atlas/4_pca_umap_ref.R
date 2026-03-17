library(Seurat)
library(future)

options(future.globals.maxSize = 20 * 1024^3)  # 20 GB

cat("Loading raw reference...\n")
ref <- readRDS("/data/chellal/aml_project/AML_ref_fixed2_20k_sampled.rds")

cat("Dimensions:\n")
print(dim(ref))

# ======================
# Normalisation
# ======================
cat("Normalizing...\n")
ref <- NormalizeData(ref)

# ======================
# Variable features
# ======================
cat("Finding variable features...\n")
ref <- FindVariableFeatures(ref, nfeatures = 3000)

# ======================
# Scaling
# ======================
cat("Scaling data...\n")
ref <- ScaleData(ref)

# ======================
# PCA
# ======================
cat("Running PCA...\n")
ref <- RunPCA(ref, npcs = 30)

# ======================
# Neighbors + Clusters
# ======================
cat("Finding neighbors and clusters...\n")
ref <- FindNeighbors(ref, dims = 1:20)
ref <- FindClusters(ref, resolution = 0.5)

# ======================
# UMAP
# ======================
cat("Running UMAP...\n")
ref <- RunUMAP(ref, dims = 1:20)

# ======================
# Save plots
# ======================
cat("Saving UMAP plots...\n")

pdf("UMAP_clusters.pdf", width = 8, height = 6)
print(DimPlot(ref, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
dev.off()

pdf("UMAP_celltypes.pdf", width = 8, height = 6)
print(DimPlot(ref, reduction = "umap", group.by = "Cell.Type"))
dev.off()

# Sauvegarder objet avec embeddings
saveRDS(ref, file = "AML_downsampled_ref_fixed2_with_umap.rds")

cat("DONE\n")
