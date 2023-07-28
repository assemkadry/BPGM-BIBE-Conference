#Create a Seurat object for each sample
for (file in c("normal", "tumor")){
  seurat_data <- Read10X(data.dir = paste0("/Users/assemkadry/Desktop/scRNA_Seq_Analysis/GSE168652_RAW/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.cells = 3,
                                   min.features = 200,
                                   project = file)
  assign(file, seurat_obj)
}

# Create a merged Seurat object
merged_seurat <- merge(x = normal, 
                       y = tumor, 
                       add.cell.id = c("normal", "tumor"))

merged_seurat
# 22100 features across 24498 samples

#Save RDS file
saveRDS(merged_seurat, "/Users/assemkadry/Desktop/scRNA_Seq_Analysis/merged_seurat.rds")

#Load RDS file
merged_seurat <- readRDS("/Users/assemkadry/Desktop/scRNA_Seq_Analysis/merged_seurat.rds")

# 1. QC -------
View(merged_seurat@meta.data)
# table(merged_seurat$orig.ident)

# MT%
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio
View(merged_seurat@meta.data)

# Visualize QC metrics as a violin plot
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# 2. Filtering -----------------
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mitoRatio < 5)
# table(merged_seurat$orig.ident)

# 3. Normalize data ----------
merged_seurat <- NormalizeData(merged_seurat)
#str(merged_seurat)

# 4. Identify highly variable features --------------
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(merged_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 5. Scaling -------------
all.genes <- rownames(merged_seurat)
merged_seurat <- ScaleData(merged_seurat, features = all.genes)
#str(merged_seurat)

#merged_seurat <- ScaleData(merged_seurat, vars.to.regress = "mitoRatio")

# 6. Perform Linear dimensionality reduction --------------
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))

# visualize PCA results
print(merged_seurat[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(merged_seurat, dims = 1, cells = 500, balanced = TRUE)

# determine dimensionality of the data
ElbowPlot(merged_seurat)

# 7. Clustering ------------
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(merged_seurat), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
#merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)
merged_seurat <- RunTSNE(merged_seurat,dims = 1:30)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(merged_seurat, reduction = "tsne",label = TRUE, repel = TRUE,label.size = 6)
DimPlot(merged_seurat, group.by ="orig.ident" , reduction = "tsne", repel = TRUE,label.size = 6)
DimPlot(merged_seurat, split.by ="orig.ident" , reduction = "tsne",label = TRUE, repel = TRUE,label.size = 6)

# 8. Cluster biomarkers -------
#find All markers

merged_seurat.markers <- FindAllMarkers(merged_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
merged_seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# find all markers of cluster 0
cluster0.markers <- FindMarkers(merged_seurat, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)

# let's visualize the main feature
FeaturePlot(merged_seurat, features = c("BPGM"), min.cutoff = 'q10')
VlnPlot(merged_seurat, features = c("BPGM"))

# 9. Assigning cell type identity to clusters (Manual) ---------

#Used Markers:

#Endothelial: from molecular marker (VWF, EMCN, FLT1) or #Endothelial cells: (EGFL7, EMCN, PECAM1)
#Cancer cells: (CDH1, EPCAM, CDKN2A), 
#Lymphocytes: (CD27, PRF1)
#Macrophages: (CD163, FCGR2A) 
#Fibroblasts: (APOD,LUM,COL1A2)
#Smooth muscle cells: (PRLR,SFRP4,ACTG2) 
#Endometrial stromal cells: (PLN,SUSD2,RGS5)

#see the plot before annotation
DimPlot(merged_seurat, reduction = 'tsne', label = TRUE, repel = TRUE)
#for test the markers
FeaturePlot(merged_seurat, features = c("put your markers here"), min.cutoff = 'q10')

merged_seurat <- RenameIdents(merged_seurat, `1` = 'Unknown')
merged_seurat <- RenameIdents(merged_seurat, `2` = 'Smooth muscle cells')
merged_seurat <- RenameIdents(merged_seurat, `3` = 'Lymphocytes')
merged_seurat <- RenameIdents(merged_seurat, `4` = 'Endothelail cells')
merged_seurat <- RenameIdents(merged_seurat, `5` = 'Fibroblasts')
merged_seurat <- RenameIdents(merged_seurat, `6` = 'Endometrial stromal cells')
merged_seurat <- RenameIdents(merged_seurat, `7` = 'Unknown')
merged_seurat <- RenameIdents(merged_seurat, `8` = 'Macrophage')
merged_seurat <- RenameIdents(merged_seurat, `0` = 'Cancer cells')

#visualize the result after annotation
DimPlot(merged_seurat, reduction = 'tsne', label = TRUE,repel = TRUE,label.size = 7)
DimPlot(merged_seurat, split.by ="orig.ident" , reduction = "tsne",
        label = TRUE, repel = TRUE,label.size = 7)


# Re-visualize the main feature
FeaturePlot(merged_seurat, features = c("BPGM"), min.cutoff = 'q10')
VlnPlot(merged_seurat, features = c("BPGM"))

#counts of cells in each cluster 
table(merged_seurat$seurat_clusters)

# 10. Differential expression testing ----

# Subset cancer cells from the merged Seurat object
cancer_cells <- subset(merged_seurat, ident = "Cancer cells")
# Extract the BPGM counts from the subsetted Seurat object
BPGM_counts <- cancer_cells@assays[["RNA"]]@counts["BPGM", ]
table(BPGM_counts)
#Remove cells with zero counts
BPGM_counts <- BPGM_counts[BPGM_counts > 0]
median(BPGM_counts)
summary(BPGM_counts)

# Create BPGM binary groups based on cutoff
cancer_cells$BPGM_group <- ifelse(cancer_cells@assays[["RNA"]]@counts["BPGM", ] >= 1, "BPGM High", "BPGM Low")

# Assign cluster identities to BPGM high and low expression groups
cancer_cells <- SetIdent(cancer_cells, value = "BPGM High", cells = which(cancer_cells$BPGM_group == "BPGM High"))
cancer_cells <- SetIdent(cancer_cells, value = "BPGM Low", cells = which(cancer_cells$BPGM_group == "BPGM Low"))

# Perform DEG analysis between BPGM High and Low groups
BPGM_markers <- FindMarkers(cancer_cells, ident.1 = "BPGM High", ident.2 = "BPGM Low")

# View the top DEGs
top_BPGM_markers <- head(BPGM_markers, n = 10)

# Print the top DEGs
print(top_BPGM_markers)

# Retrieve the column names (gene names) of the BPGM_markers data frame
colnames(BPGM_markers)
# Calculate the absolute average log2 fold change for each gene and add it as a new column
BPGM_markers$abs_avg_log2FC <- abs(BPGM_markers$avg_log2FC)
# View the BPGM_markers data frame
view(BPGM_markers)
# Generate violin plots to visualize the expression of genes BPGM, SUN1, and EGLN3 in the cancer_cells dataset
VlnPlot(cancer_cells, features = c("BPGM","SUN1", "EGLN3"))

#export the DEGs
write.csv(as.matrix(BPGM_markers),file="BPGM_markers.csv",quote = F,row.names = T)
