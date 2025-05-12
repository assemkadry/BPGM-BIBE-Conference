# ---------------------------------------
# Single-Cell RNA-Seq Analysis
# Author: Assem K. Elsherif
# Data Source: GSE168652
# Date September 25, 2023
# ---------------------------------------

# Load Seurat library for single-cell RNA-seq analysis
library(Seurat)

# -------------------------
# 1. Load and Create Seurat Objects
# -------------------------
# Loop through each sample directory ("normal", "tumor")
# Read data using Read10X and create Seurat objects with minimum filtering thresholds
for (file in c("normal", "tumor")){
  seurat_data <- Read10X(data.dir = paste0("/Users/assemkadry/Desktop/scRNA_Seq_Analysis/GSE168652_RAW/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.cells = 3,
                                   min.features = 200,
                                   project = file)
  assign(file, seurat_obj)
}

# -------------------------
# 2. Merge Seurat Objects
# -------------------------
# Merge normal and tumor Seurat objects with cell ID prefixes
merged_seurat <- merge(x = normal, 
                       y = tumor, 
                       add.cell.id = c("normal", "tumor"))

# Check dimensions: 22100 genes across 24498 cells
merged_seurat

# Save merged object as RDS
saveRDS(merged_seurat, "/Users/assemkadry/Desktop/scRNA_Seq_Analysis/merged_seurat.rds")

# Load merged Seurat object from RDS (optional if reloading)
merged_seurat <- readRDS("/Users/assemkadry/Desktop/scRNA_Seq_Analysis/merged_seurat.rds")

# -------------------------
# 3. Quality Control
# -------------------------
# Visualize metadata
View(merged_seurat@meta.data)

# Calculate mitochondrial gene percentage using "^MT-" prefix
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio

# Visual QC plots: violin plots and feature scatter plots
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)

# Check relationships between features
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# -------------------------
# 4. Cell Filtering
# -------------------------
# Keep cells with 200–2500 features and <5% mitochondrial reads
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mitoRatio < 5)

# -------------------------
# 5. Normalization and HVGs
# -------------------------
# Normalize gene expression
merged_seurat <- NormalizeData(merged_seurat)

# Identify 2000 most variable genes using VST method
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

# Visualize top 10 highly variable genes
top10 <- head(VariableFeatures(merged_seurat), 10)
plot1 <- VariableFeaturePlot(merged_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# -------------------------
# 6. Scaling and PCA
# -------------------------
# Scale all genes
all.genes <- rownames(merged_seurat)
merged_seurat <- ScaleData(merged_seurat, features = all.genes)

# Run PCA using variable features
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))

# Examine PCA output
print(merged_seurat[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(merged_seurat, dims = 1, cells = 500, balanced = TRUE)

# Determine number of PCs to use
ElbowPlot(merged_seurat)

# -------------------------
# 7. Clustering and Dimensionality Reduction
# -------------------------
# Build kNN graph and find clusters
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

# Run tSNE for visualization
merged_seurat <- RunTSNE(merged_seurat,dims = 1:30)

# Plot clusters and sample identities
DimPlot(merged_seurat, reduction = "tsne",label = TRUE, repel = TRUE,label.size = 6)
DimPlot(merged_seurat, group.by ="orig.ident" , reduction = "tsne", repel = TRUE,label.size = 6)
DimPlot(merged_seurat, split.by ="orig.ident" , reduction = "tsne",label = TRUE, repel = TRUE,label.size = 6)

# -------------------------
# 8. Marker Gene Analysis
# -------------------------
# Find marker genes for all clusters
merged_seurat.markers <- FindAllMarkers(merged_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Show top 2 markers per cluster
merged_seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Find markers for cluster 0 only
cluster0.markers <- FindMarkers(merged_seurat, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)

# Visualize expression of BPGM gene
FeaturePlot(merged_seurat, features = c("BPGM"), min.cutoff = 'q10')
VlnPlot(merged_seurat, features = c("BPGM"))

# -------------------------
# 9. Manual Cell Type Annotation
# -------------------------
# Marker references (commented for clarity)
# Endothelial: VWF, EMCN, FLT1, EGFL7, PECAM1
# Cancer: CDH1, EPCAM, CDKN2A
# Lymphocytes: CD27, PRF1
# Macrophages: CD163, FCGR2A
# Fibroblasts: APOD, LUM, COL1A2
# Smooth Muscle: PRLR, SFRP4, ACTG2
# Endometrial Stromal: PLN, SUSD2, RGS5

# View pre-annotation clusters
DimPlot(merged_seurat, reduction = 'tsne', label = TRUE, repel = TRUE)

# Visualize marker expression
FeaturePlot(merged_seurat, features = c("put your markers here"), min.cutoff = 'q10')

# Rename cluster identities manually
merged_seurat <- RenameIdents(merged_seurat, `1` = 'Unknown')
merged_seurat <- RenameIdents(merged_seurat, `2` = 'Smooth muscle cells')
merged_seurat <- RenameIdents(merged_seurat, `3` = 'Lymphocytes')
merged_seurat <- RenameIdents(merged_seurat, `4` = 'Endothelail cells')
merged_seurat <- RenameIdents(merged_seurat, `5` = 'Fibroblasts')
merged_seurat <- RenameIdents(merged_seurat, `6` = 'Endometrial stromal cells')
merged_seurat <- RenameIdents(merged_seurat, `7` = 'Unknown')
merged_seurat <- RenameIdents(merged_seurat, `8` = 'Macrophage')
merged_seurat <- RenameIdents(merged_seurat, `0` = 'Cancer cells')

# View annotated clusters
DimPlot(merged_seurat, reduction = 'tsne', label = TRUE,repel = TRUE,label.size = 7)
DimPlot(merged_seurat, split.by ="orig.ident" , reduction = "tsne",label = TRUE, repel = TRUE,label.size = 7)

# Replot expression of key genes after annotation
FeaturePlot(merged_seurat, features = c("BPGM"), min.cutoff = 'q10')
VlnPlot(merged_seurat, features = c("BPGM"))

# Show cell count per cluster
table(merged_seurat$seurat_clusters)

# -------------------------
# 10. BPGM Group DE Analysis in Cancer Cells
# -------------------------
# Subset only "Cancer cells"
cancer_cells <- subset(merged_seurat, ident = "Cancer cells")

# Extract raw counts of BPGM
BPGM_counts <- cancer_cells@assays[["RNA"]]@counts["BPGM", ]
table(BPGM_counts)

# Remove zero-count cells and summarize
BPGM_counts <- BPGM_counts[BPGM_counts > 0]
median(BPGM_counts)
summary(BPGM_counts)

# Create high/low BPGM groups based on count cutoff
cancer_cells$BPGM_group <- ifelse(cancer_cells@assays[["RNA"]]@counts["BPGM", ] >= 1, "BPGM High", "BPGM Low")

# Set group identities
cancer_cells <- SetIdent(cancer_cells, value = "BPGM High", cells = which(cancer_cells$BPGM_group == "BPGM High"))
cancer_cells <- SetIdent(cancer_cells, value = "BPGM Low", cells = which(cancer_cells$BPGM_group == "BPGM Low"))

# Compare gene expression between BPGM groups
BPGM_markers <- FindMarkers(cancer_cells, ident.1 = "BPGM High", ident.2 = "BPGM Low")
top_BPGM_markers <- head(BPGM_markers, n = 10)
print(top_BPGM_markers)

# Add column for absolute log2FC
BPGM_markers$abs_avg_log2FC <- abs(BPGM_markers$avg_log2FC)

# View full table of DEGs
view(BPGM_markers)

# Visualize BPGM, SUN1, EGLN3 expression in merged and cancer-only data
VlnPlot(merged_seurat, features = c("BPGM", "SUN1", "EGLN3"))
VlnPlot(cancer_cells, features = c("BPGM", "SUN1", "EGLN3"))

# Export DEGs to CSV
write.csv(as.matrix(BPGM_markers), file="BPGM_markers.csv", quote = F, row.names = T)