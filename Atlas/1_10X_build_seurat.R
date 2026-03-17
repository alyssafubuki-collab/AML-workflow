library(Seurat)

counts <- Read10X("matrix_files")
head(rownames(counts))


# ----------------------------
# 1. Lire la matrice 10X
# ----------------------------
inputMatrixDir <- "/data/chellal/aml_project/matrix_files"

if(!dir.exists(inputMatrixDir)){
  stop(paste("ERREUR : Le dossier n'existe pas :", inputMatrixDir))
}

cat("Reading 10X matrix...\n")
counts <- Read10X(
  data.dir = inputMatrixDir,
  gene.column = 1,        # utiliser 2 si features.tsv a gene symbol
  unique.features = FALSE,
  strip.suffix = FALSE
)

cat("First gene names:\n")
print(head(rownames(counts)))
head(rownames(counts))
str(counts)
class(counts)
dim(counts)
list.files("/data/chellal/aml_project/matrix_files")


# Vérifier doublons
dup_genes <- sum(duplicated(rownames(counts)))
cat("Nombre de gènes dupliqués :", dup_genes, "\n")
if(dup_genes > 0){
  rownames(counts) <- make.unique(rownames(counts))
  cat("Doublons corrigés avec make.unique()\n")
}

# ----------------------------
# 2. Créer Seurat object
# ----------------------------
cat("Creating Seurat object...\n")
seurat_ref <- CreateSeuratObject(
  counts = counts,
  project = "AML_ref",
  min.cells = 0,
  min.features = 0
)
cat("Object dimensions (cells x genes):\n")
print(dim(seurat_ref))

# ----------------------------
# 3. Ajouter les metadata
# ----------------------------
meta_file <- "/data/chellal/aml_project/metadata.csv"

if(!file.exists(meta_file)){
  stop(paste("ERREUR : Le fichier metadata n'existe pas :", meta_file))
}

cat("Reading metadata...\n")
meta <- read.csv(meta_file, row.names = 1)

# Vérifier que les barcodes correspondent
if(!all(rownames(meta) %in% colnames(seurat_ref))){
  warning("Attention : certains barcodes du metadata ne sont pas dans Seurat object !")
}

cat("Adding metadata to Seurat object...\n")
seurat_ref <- AddMetaData(seurat_ref, metadata = meta)

# ----------------------------
# 4. Sauvegarder
# ----------------------------
cat("Saving Seurat object as AML_ref_seurat_fixed2.rds\n")
saveRDS(seurat_ref, file = "AML_ref_seurat_fixed2.rds")

cat("DONE\n")
