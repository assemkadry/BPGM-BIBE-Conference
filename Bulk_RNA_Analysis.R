# ---------------------------------------
# Bulk RNA-Seq Analysis of TCGA Data
# Author: Assem K. Elsherif
# Date September 25, 2023
# ---------------------------------------

# ------------------------
# Required Libraries
# ------------------------

# For expression processing and annotation
library(tidyverse)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# For DE analysis
library(DESeq2)

# For visualization
library(EnhancedVolcano)
library(pheatmap)
library(cowplot)

# For enrichment analysis
library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(ggplot2)
library(forcats)

# ------------------------
# 1. Load the data
# ------------------------

# Set the path to the mRNA-Seq data directory
data.path="/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/TCGA/gdc_download_20220928_184027.813651"

# Get a list of all .tsv files in the data directory (recursively)
files <- list.files(path=data.path, recursive=TRUE, pattern="tsv")

# Read the first file as initialization
file = files[1]
file.id = strsplit(file, "/")[[1]][1]  # extract file ID from folder structure

# Open the .tsv file and skip first 6 lines (header info)
tsv.con <- file(file.path(data.path, files[1]), "r")
temp <- read.table(tsv.con, sep = "\t", header = FALSE, skip = 6)

# Initialize expression matrix with gene and count columns
mrna.exp = temp
mrna.exp <- mrna.exp[, c(1, 4)]
rownames(mrna.exp) = mrna.exp[, 1]     # use gene ID as row names
mrna.exp = mrna.exp[-1]                # remove first column after storing row names
colnames(mrna.exp) = c(file.id)        # assign column name based on file ID

# Read and append remaining files
for (i in 2:length(files)) {
  file = files[i]
  file.id = strsplit(file, "/")[[1]][1]
  
  tsv.con <- file(file.path(data.path, files[i]), "r")
  temp <- read.table(tsv.con, sep = "\t", header = FALSE, skip = 6)
  
  temp <- temp[, c(1, 4)]
  rownames(temp) = temp[, 1]
  temp = temp[-1]
  colnames(temp) = c(file.id)
  
  mrna.exp = cbind(mrna.exp, temp)  # append column-wise
}

# View the resulting expression matrix
View(mrna.exp)

# Close the last opened file connection
close(tsv.con)

# ------------------------
# 2. Mapping Ensembl ID to Gene Symbol
# ------------------------

# Extract base Ensembl ID (remove version suffix)
ensemble.id = sapply(rownames(mrna.exp), function(x) strsplit(as.character(x), "\\.")[[1]][1])
View(ensemble.id)

# Add Ensembl ID as a new column
mrna.exp = cbind(ensemble.id, mrna.exp)

# Map Ensembl IDs to gene symbols using org.Hs.eg.db
mapper <- mapIds(org.Hs.eg.db, keys=ensemble.id, column="SYMBOL", keytype="ENSEMBL", multiVals="first") 
mapper.df = as.data.frame(mapper)
mapper.df <- cbind(ensemble.id, mapper.df)
names(mapper.df) = c("ensemble.id", "symbol")

# Merge gene symbol info with expression data
mrna.exp2 = merge(mrna.exp, mapper.df, by="ensemble.id", all.x=TRUE)

# Remove Ensembl column
mrna.exp2 = mrna.exp2[-1]

# Remove genes without symbols
mrna.exp2 = mrna.exp2[!is.na(mrna.exp2$symbol), ]

# ------------------------
# 3. Collapse duplicate gene symbols
# ------------------------

# Check number of duplicated gene symbols
x = duplicated(mrna.exp2$symbol)  
sum(x)

# Prepare for aggregation
mrna.exp.data = mrna.exp2[-dim(mrna.exp2)[2]]  # remove last column (symbol)
mrna.exp.data = apply(mrna.exp.data, 2, as.numeric)

# Aggregate expression by gene symbol (mean)
mrna.exp.data.agg = aggregate(mrna.exp.data, list(mrna.exp2$symbol), FUN=mean)

# Set rownames to gene symbols
rownames(mrna.exp.data.agg) = mrna.exp.data.agg$Group.1
mrna.exp.data.agg = mrna.exp.data.agg[-1]

file.ids = colnames(mrna.exp.data.agg)

# ------------------------
# 4. Load phenotype data (sample sheet)
# ------------------------

# Read TCGA sample sheet
pheno <- read_delim("/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/TCGA/gdc_sample_sheet.2022-09-28.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
View(pheno)
table(pheno$`Sample Type`)

# Replace spaces in column names with dots
pheno.names = names(pheno)
names(pheno) = sapply(pheno.names, function(x) gsub(" ", ".", x))
table(pheno$Sample.Type)

# ------------------------
# 5. Match expression with phenotype data
# ------------------------

# Match expression file IDs with sample sheet
file.ids.pheno = pheno$File.ID
index.files = match(file.ids, file.ids.pheno)

# Rename expression matrix columns using sample IDs
names(mrna.exp.data.agg) = pheno$Sample.ID[index.files]

# ------------------------
# 6. Extract tumor samples only
# ------------------------

tumor.samples_ID = pheno[pheno$Sample.Type %in% c("Primary Tumor", "Metastatic"), ]$Sample.ID
tumor.exp = mrna.exp.data.agg[, names(mrna.exp.data.agg) %in% tumor.samples_ID]

pheno.sub = pheno[pheno$Sample.ID %in% tumor.samples_ID, c("Sample.ID", "Sample.Type")]

# Convert expression values to integers
exp.sub = apply(tumor.exp, 2, as.integer)
rownames(exp.sub) = rownames(tumor.exp)

# ------------------------
# 7. Create BPGM expression groups
# ------------------------

# Log-transform expression
log.exp.sub <- log2(exp.sub + 1)

# Extract BPGM expression and calculate median
BPGM.exp <- log.exp.sub["BPGM", ]
BPGM.median <- median(BPGM.exp)

# Assign samples to High/Low BPGM groups
BPGM_group <- ifelse(BPGM.exp > BPGM.median, "BPGM_High", "BPGM_Low")
pheno.sub$BPGM <- as.factor(BPGM_group)
table(pheno.sub$BPGM)

# Compare group-wise BPGM expression means
BPGM_high_exp = log.exp.sub[, pheno.sub$BPGM == "BPGM_High"]
BPGM_high_mean = mean(BPGM_high_exp["BPGM", ])

BPGM_low_exp = log.exp.sub[, pheno.sub$BPGM == "BPGM_Low"]
BPGM_low_mean = mean(BPGM_low_exp["BPGM", ])

# ------------------------
# 8. DE Analysis using DESeq2
# ------------------------

# Create DESeqDataSet object and run DE analysis
cond1 = "BPGM_High"
cond2 = "BPGM_Low"
dds = DESeqDataSetFromMatrix(countData = exp.sub, colData = pheno.sub, design = ~ BPGM)
dds.run = DESeq(dds)

# Get DE results with contrast
res = results(dds.run)
res = results(dds.run, contrast = c("BPGM", cond1, cond2))

# Clean up results and filter DEGs
res = res[complete.cases(res), ]
summary(res)
res.df = as.data.frame(res)
res.degs = res.df[res.df$padj < 0.05 & abs(res.df$log2FoldChange) > 1, ]

# Extract expression data for DEGs
exp.degs = exp.sub[rownames(exp.sub) %in% rownames(res.degs), ]

# ------------------------
# 9. DEG Utilities: Up/Down Regulated Genes
# ------------------------

# Function to extract upregulated genes
get_upregulated <- function(df) {
  key <- intersect(rownames(df)[which(df$log2FoldChange >= 1)], rownames(df)[which(df$padj <= 0.05)])
  results <- as.data.frame(df[rownames(df) %in% key, ])
  return(results)
}

# Function to extract downregulated genes
get_downregulated <- function(df) {
  key <- intersect(rownames(df)[which(df$log2FoldChange <= -1)], rownames(df)[which(df$padj <= 0.05)])
  results <- as.data.frame(df[rownames(df) %in% key, ])
  return(results)
}

# Save DEG results to CSV
up_reg <- get_upregulated(res.df)
down_reg <- get_downregulated(res.df)

write.csv(as.matrix(res.degs), file = "res.degs.csv", quote = FALSE, row.names = TRUE)
write.csv(as.matrix(up_reg), file = "up_reg.csv", quote = FALSE, row.names = TRUE)
write.csv(as.matrix(down_reg), file = "down_reg.csv", quote = FALSE, row.names = TRUE)

# ------------------------
# 10. Visualization
# ------------------------

# 10.1 Volcano Plot
# ------------------------
# Visualize DEGs using EnhancedVolcano
EnhancedVolcano(res, lab = rownames(res),
                title = "",
                subtitle = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-6, 6),
                legendLabels = c('Not sig.', 'LFC', 'p-value', 'p-value & LFC'),
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                pCutoff = 0.05,
                FCcutoff = 1,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                labFace = 'bold',
                labSize = 5.0,
                parseLabels = TRUE,
                selectLab = c('BPGM'))  # highlight BPGM
                

# 10.2 Heatmap
# ------------------------
# Prepare annotation for heatmap columns
colnames = colnames(exp.degs) 
case.vector = rep("BPGM_high_exp", length(colnames(BPGM_high_exp)))
ctrl.vector = rep("BPGM_low_exp", length(colnames(BPGM_low_exp)))
Sample = c(case.vector, ctrl.vector)

# Create sample annotation dataframe
annotation = as.data.frame(Sample)
rownames(annotation) = colnames

# Define annotation colors
Sample = c("lightgreen", "navy")
names(Sample) = c("BPGM_high_exp", "BPGM_low_exp")
ann_colors = list(Sample = Sample)

# Scale expression data (gene-wise)
m2 = scale(t(exp.degs), center = TRUE, scale = TRUE)
m2 = t(m2)

# Plot heatmap of DEGs using pheatmap
pheatmap(m2, annotation = annotation, annotation_colors = ann_colors,
         fontsize_row = 5, fontsize_col = 9)

# Save R environment (optional)
save.image()

# ------------------------
# 11. Gene Set Enrichment Analysis with ClusterProfiler
# ------------------------

# 11.1 GO Enrichment (gseGO)
# ------------------------

# Set the desired organism database
organism = "org.Hs.eg.db"

# Prepare input: read DE result
df = read.csv("res.degs.csv", header = TRUE)
names(df)[names(df) == "X"] <- "GeneSymbol"

# Extract log2FoldChange vector with gene symbols as names
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$GeneSymbol
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

# Show available ID types for the selected organism
keytypes(org.Hs.eg.db)

# Run GO enrichment using fgsea method
set.seed(50)
gse <- gseGO(geneList = gene_list, 
             ont = "ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism,
             by = "fgsea",
             eps = 1e-10,
             pAdjustMethod = "BH")

# Preview GO enrichment output
head(gse)

# Save GO descriptions to file
require(DOSE)
go <- as.data.frame(gse@result[["Description"]])
write.csv(as.matrix(go), file = "go.csv", quote = FALSE, row.names = TRUE)

# Plot top GO terms as dotplot
go_plot <- dotplot(gse, showCategory = 10, orderBy = "setSize", split = ".sign", 
                   label_format = function(x) stringr::str_wrap(x, width = 80)) + 
                   facet_grid(. ~ .sign) + ggtitle("GO")
plot(go_plot)

# Choose specific GO terms for focused plot
selected_indices <- c(1,2,3,4,5,15,21,22,23,28,35,49,54,58,61,83,87,88)
selected_go <- gse$Description[selected_indices]

# Plot selected GO terms
go_plot_desired <- dotplot(gse, showCategory = selected_go, orderBy = "setSize", split = ".sign", 
                           label_format = function(x) stringr::str_wrap(x, width = 80)) + 
                           facet_grid(. ~ .sign) + ggtitle("GO")
plot(go_plot_desired)

# 11.2 KEGG Enrichment (gseKEGG)
# ------------------------

# Convert SYMBOL → ENTREZID
df = read.csv("res.degs.csv", header = TRUE)
names(df)[names(df) == "X"] <- "GeneSymbol"
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$GeneSymbol

ids = bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
colnames(dedup_ids) = c("GeneSymbol", "EntrezID")

df2 = merge(df, dedup_ids, by = "GeneSymbol")
kegg_gene_list = df2$log2FoldChange
names(kegg_gene_list) = df2$EntrezID
kegg_gene_list = na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"

# Run KEGG enrichment
set.seed(50)
kegg = gseKEGG(geneList = kegg_gene_list,
               organism = kegg_organism,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               by = "fgsea",
               eps = 1e-10,
               keyType = "ncbi-geneid")

# Save KEGG term descriptions
head(kegg)
pathways_kegg <- as.data.frame(kegg@result[["Description"]])
write.csv(as.matrix(pathways_kegg), file = "pathways_kegg.csv", quote = FALSE, row.names = TRUE)

# Custom theme for KEGG plot
mytheme <- theme(
  panel.border = element_blank(),
  panel.grid.major = element_line(linetype = 'dotted', colour = '#808080'),
  panel.grid.major.y = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line.x = element_line()
)

# Plot KEGG pathways as NES bar chart
kegg_plot <- ggplot(kegg, showCategory = 25, aes(NES, fct_reorder(Description, NES), fill = pvalue)) +
  geom_col() +
  geom_segment(mapping = aes(x = -2.1, xend = ifelse(sign(NES) > 0, 0, NES), yend = Description)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours = c('#b3eebe', "#46bac2", '#371ea3'), guide = guide_colorbar(reverse = TRUE)) +
  theme_dose(12) + mytheme +
  xlab("Normalized Enrichment Score") + ylab(NULL) +
  ggtitle("KEGG Enriched Pathways")
plot(kegg_plot)

# 11.3 MSigDb GSEA
# ------------------------

# Prepare input
df = read.csv("res.degs.csv", header = TRUE)
names(df)[names(df) == "X"] <- "GeneSymbol"
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$GeneSymbol

ids = bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
colnames(dedup_ids) = c("GeneSymbol", "EntrezID")

df3 = merge(df, dedup_ids, by = "GeneSymbol")
MSigDb_gene_list = df3$log2FoldChange
names(MSigDb_gene_list) = df3$EntrezID
MSigDb_gene_list = na.omit(MSigDb_gene_list)
MSigDb_gene_list = sort(MSigDb_gene_list, decreasing = TRUE)

# Load MSigDb C2 gene sets
C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C2_t2g)

# Run GSEA using C2 collection
set.seed(50)
msig <- GSEA(geneList = MSigDb_gene_list, TERM2GENE = C2_t2g,
             pvalueCutoff = 0.05,
             eps = 1e-10,
             by = "fgsea",
             pAdjustMethod = "none")

# Plot selected gene sets
gsea_plot <- gseaplot2(msig, geneSetID = c(42, 77, 95))
print(gsea_plot)

# Combine all enrichment plots in a single figure
cowplot::plot_grid(go_plot_desired, kegg_plot, print(gsea_plot), labels = c("A", "B", "C"))