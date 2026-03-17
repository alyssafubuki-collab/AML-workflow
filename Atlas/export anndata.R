import scanpy as sc
from scipy import io
import os

# CrÃ©er dossier sortie
os.makedirs("matrix_files", exist_ok=True)

# Lire fichier h5ad (CHEMIN CLUSTER)
adata = sc.read_h5ad("/data/chellal/AML_scAtlas2_Sept2025.h5ad")

print(adata)

# Barcodes
with open('matrix_files/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')

# Features
with open('matrix_files/features.tsv', 'w') as f:
    for item in ['\t'.join([x, x, 'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')

# Matrix
io.mmwrite('matrix_files/matrix', adata.X)

# Compression
os.system("gzip matrix_files/*")

# Metadata
adata.obs.to_csv('metadata.csv')

zcat matrix_files/matrix.mtx.gz | head -3
