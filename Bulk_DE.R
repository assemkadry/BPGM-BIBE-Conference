# 1.Load the data ----
#Set the path to the mRNA-Seq data directory
data.path="/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/TCGA/gdc_download_20220928_184027.813651"
# Get a list of files in the data directory
files <- list.files(path=data.path,recursive=T, pattern = "tsv")

# read the first file for the first time
file=files[1]
file.id=strsplit(file,"/")[[1]][1]

#read tsv the file
tsv.con <- file(file.path(data.path,files[1]),"r")
temp <- read.table(tsv.con, sep = "\t", header = F, skip = 6)

#create a storing object mrna.exp to save the whole read counts of each file read in an iteration
mrna.exp=temp
mrna.exp <- mrna.exp[, c(1, 4)]
rownames(mrna.exp)=mrna.exp[,1]
mrna.exp=mrna.exp[-1]
colnames(mrna.exp)=c(file.id)

#Read and process the remaining files

for(i in 2: length(files))
{
  #refer to the next file (note that we start from index 2, because we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  
  # read the next file  
  tsv.con <- file(file.path(data.path,files[i]),"r")
  temp <- read.table(tsv.con, sep = "\t", header = F, skip = 6)
  
  # remove the first column, bec we had it already
  
  temp <- temp[, c(1, 4)]
  rownames(temp)=temp[,1]
  temp=temp[-1]
  colnames(temp)=c(file.id)
  mrna.exp=cbind(mrna.exp,temp)
}
# View the mrna.exp data in a table view
View(mrna.exp)
# Close the file connection
close(tsv.con)

##mapping of ensembel.id to gene symbol
# Prepare the ensemble id to be in the same format as the one in the database
ensemble.id=sapply(rownames(mrna.exp), function(x) strsplit(as.character(x),"\\.")[[1]][1])
View(ensemble.id)
# Add ensemble.id as a column to mrna.exp
mrna.exp=cbind(ensemble.id,mrna.exp)
# Map ensemble.id to gene symbol using the org.Hs.eg.db database
mapper<- mapIds(org.Hs.eg.db, keys=ensemble.id, column="SYMBOL",keytype="ENSEMBL", multiVals="first") 
mapper.df=as.data.frame(mapper)
mapper.df <- cbind(ensemble.id, mapper.df)
names(mapper.df)=c("ensemble.id","symbol")
# Merge mrna.exp and mapper.df based on ensemble.id
mrna.exp2=merge(mrna.exp,mapper.df,by="ensemble.id",all.x=T) 

# Drop the first column (ensemble.id)
mrna.exp2=mrna.exp2[-1]
# Remove rows with missing gene symbols
mrna.exp2=mrna.exp2[ ! is.na(mrna.exp2$symbol),]

# Check for duplicated gene symbols
x=duplicated(mrna.exp2$symbol)  
sum(x)

# Duplicate gene symbols may be due to different transcripts, so we need to aggregate them
mrna.exp.data=mrna.exp2[-dim(mrna.exp2)[2]]
mrna.exp.data=apply(mrna.exp.data,2, as.numeric)

# Remove duplication by aggregation
mrna.exp.data.agg= aggregate(mrna.exp.data, list(mrna.exp2$symbol),FUN=mean)
# Set gene symbols as row names
rownames(mrna.exp.data.agg)=mrna.exp.data.agg$Group.1
mrna.exp.data.agg=mrna.exp.data.agg[-1]
file.ids=colnames(mrna.exp.data.agg)

#save(mrna.exp.data.agg, file="/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/TCGA/mrna.exp.data.agg.RData")
#load("/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/TCGA/mrna.exp.data.agg.RData")

# Load the mRNA sample sheets
pheno <- read_delim("/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/TCGA/gdc_sample_sheet.2022-09-28.tsv", "\t", escape_double = FALSE, trim_ws = TRUE,show_col_types = FALSE)
View(pheno)
table(pheno$`Sample Type`)

# Rename column names by replacing spaces with dots
pheno.names=names(pheno)
names(pheno)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno$Sample.Type)

#we will rename the columns of our expression data with the sample ids columns of the pheno file
# Match file ids to file.ids.pheno

file.ids.pheno=pheno$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mrna.exp.data.agg)=pheno$Sample.ID[index.files]

# Get tumor samples ID
tumor.samples_ID= pheno[ pheno$Sample.Type %in% c("Primary Tumor","Metastatic"),]$Sample.ID

# Retrieve expression data and pheno data for tumor samples
tumor.exp=mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% tumor.samples_ID]
# Retrieve pheno data for tumor samples
pheno.sub=pheno[pheno$Sample.ID %in% c(tumor.samples_ID), c("Sample.ID", "Sample.Type")]
# Convert expression data to integer and set row names
exp.sub=apply (tumor.exp, 2,as.integer)
rownames(exp.sub)=rownames(tumor.exp)

##Relative gene expression
#Normalize the expression data using log transformation method.
log.exp.sub <- log2(exp.sub + 1)
#Calculate the expression of BPGM for each sample.
BPGM.exp <- log.exp.sub["BPGM", ]
#Calculate the median value of BPGM expression across all samples.
BPGM.median <- median(BPGM.exp)
#Assign patients with BPGM expression values above the median as BPGM_High and patients with BPGM expression values below the median as BPGM_Low.
BPGM_group <- ifelse(BPGM.exp > BPGM.median, "BPGM_High", "BPGM_Low")
pheno.sub$BPGM <- BPGM_group
#print table showing the number of samples in each group
table(pheno.sub$BPGM)
pheno.sub$BPGM <- as.factor(pheno.sub$BPGM)

# subset expression data for samples classified as "BPGM_High"
BPGM_high_exp = log.exp.sub[, pheno.sub$BPGM == "BPGM_High"]
# calculate median expression value of BPGM for "BPGM_High" samples
BPGM_high_mean = mean(BPGM_high_exp["BPGM", ])

# subset expression data for samples classified as "BPGM_Low"
BPGM_low_exp = log.exp.sub[, pheno.sub$BPGM == "BPGM_Low"]
# calculate median expression value of BPGM for "BPGM_Low" samples
BPGM_low_mean = mean(BPGM_low_exp["BPGM", ])

#2. DE Analysis -----
#DE analysis using DeSeq2
# Define conditions for comparison
cond1="BPGM_High"
cond2="BPGM_Low"
# Create DESeqDataSet object
dds = DESeqDataSetFromMatrix( countData = exp.sub , colData = pheno.sub , design = ~ BPGM)
# Run DESeq analysis
dds.run = DESeq(dds)

#direct results or specifying the contrast (to make a res object based on two specific BPGM high and low groups)
res=results(dds.run)
# Obtain DE results for specific contrast (BPGM high vs BPGM low)
res=results(dds.run, contrast = c("BPGM",cond1 ,cond2))

# Remove rows with null values
res=res[complete.cases(res), ]
# Summary of DE results
summary(res)
# Convert DE results to a data frame
res.df=as.data.frame(res)
# Select differentially expressed genes based on significance and fold change thresholds
res.degs=res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange)>1,]

#expression of these degs for plots
exp.degs= exp.sub[rownames(exp.sub) %in% rownames(res.degs), ]

#Function to split results
get_upregulated <- function(df){
  key <- intersect(rownames(df)[which(df$log2FoldChange>=1)],
                   rownames(df)[which(df$padj<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

get_downregulated <- function(df){
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-1)],
                   rownames(df)[which(df$padj<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

up_reg <- get_upregulated(res.df)
down_reg <- get_downregulated(res.df)

#export the DEGs
write.csv(as.matrix(res.degs),file="res.degs.csv",quote = F,row.names = T)
write.csv(as.matrix(up_reg),file="up_reg.csv",quote = F,row.names = T)
write.csv(as.matrix(down_reg),file="down_reg.csv",quote = F,row.names = T)


#3. Visualization -----
#3.1 Volcano Plot ----

EnhancedVolcano(res,lab = rownames(res),
                title = "",
                subtitle = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-6, 6),
                legendLabels=c('Not sig.','LFC','p-value','p-value & LFC'),
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize =4.0,
                pCutoff = 0.05,
                FCcutoff = 1,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                labFace = 'bold',
                labSize = 5.0,
                parseLabels = TRUE,
                selectLab = c('BPGM'))


#3.3. Heatmap ------
# Get the expression profiles of the differentially expressed genes (DEGs) only
colnames=colnames(exp.degs) 
case.vector=rep("BPGM_high_exp", length(colnames(BPGM_high_exp)))
ctrl.vector=rep("BPGM_low_exp", length(colnames(BPGM_low_exp)))
Sample=c(case.vector, ctrl.vector) 
# Create an annotation data frame for sample labels
annotation=as.data.frame(Sample)
rownames(annotation)=colnames
annotation

# Specify colors for the sample labels
Sample = c("lightgreen", "navy")
names(Sample) = c("BPGM_high_exp", "BPGM_low_exp")
ann_colors = list(Sample = Sample)
# Scale the expression data
m2=scale(t(exp.degs),center=T,scale=T)

m2=t(m2)
# Generate the heatmap using pheatmap
pheatmap(m2, annotation = annotation, annotation_colors = ann_colors, fontsize_row = 5,fontsize_col = 9)
# Save the current workspace image
save.image()


#4. Gene Set Enrichment Analysis with ClusterProfiler ----
#4.1. Gene Set Enrichment usig GO ----

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"

#Prepare Input
# reading in data from desq2
df = read.csv("res.degs.csv", header=TRUE)
names(df)[names(df) == "X"] <- "GeneSymbol"

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$GeneSymbol
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Identify available key types in the specified organism database
keytypes(org.Hs.eg.db)

# for reproducibility
set.seed(50) 
# Perform Gene Set Enrichment Analysis using GO
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism,
             by = "fgsea",
             eps = 1e-10,
             pAdjustMethod = "BH")

# Display the first few rows of the Gene Set Enrichment Analysis results
head(gse)
# Load the DOSE package
require(DOSE)
# Extract the gene set descriptions and save them to a CSV file
go <- as.data.frame(gse@result[["Description"]])
write.csv(as.matrix(go),file="go.csv",quote = F,row.names = T)

# Generate a dotplot of the Gene Set Enrichment Analysis results, showing the top N categories and splitting by sign
go_plot <- dotplot(gse, showCategory=10, orderBy = "setSize",split=".sign",label_format = function(x) stringr::str_wrap(x, width=80)) + facet_grid(.~.sign) + ggtitle("GO")
plot(go_plot)

#slecet the desired GO terms to be plotted
selected_indices <-c(1,2,3,4,5,15,21,22,23,28,35,49,54,58,61,83,87,88)
selected_go <- gse$Description[selected_indices]

# Generate a dotplot of the Gene Set Enrichment Analysis results, showing the desired categories and splitting by sign
go_plot_desired <- dotplot(gse, showCategory=selected_go,orderBy = "setSize",split=".sign",label_format = function(x) stringr::str_wrap(x, width=80)) + facet_grid(.~.sign) + ggtitle("GO")

# Display the plot
plot(go_plot_desired)

#4.2. Gene Set Enrichment usig KEGG ----

#Prepare Input
#no need to repeat the following steps if done before (from reading: name the vetor)
# reading in data from deseq2
df = read.csv("res.degs.csv", header=TRUE)
names(df)[names(df) == "X"] <- "GeneSymbol"

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$GeneSymbol

# Convert gene IDs for gseKEGG function
ids = bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# remove duplicate IDS (here we use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has the respective entrez IDs for the gene symbols.
colnames(dedup_ids) = c("GeneSymbol", "EntrezID")
df2 = merge(df, dedup_ids, by = "GeneSymbol")

# Create a vector of the gene universe
kegg_gene_list = df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) = df2$EntrezID

# omit any NA values 
kegg_gene_list = na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
# Set the desired KEGG organism
kegg_organism = "hsa"
# for reproducibility
set.seed(50) 
# Perform Gene Set Enrichment Analysis using KEGG
kegg = gseKEGG(geneList     = kegg_gene_list,
              organism     = kegg_organism,
              minGSSize    = 3,
              maxGSSize    = 800,
              pvalueCutoff = 0.05,
              pAdjustMethod = "none",
              by = "fgsea",
              eps = 1e-10,
              keyType = "ncbi-geneid")

# Display the first few rows of the KEGG Gene Set Enrichment Analysis results
head(kegg)

# Extract the pathway descriptions and save them to a CSV file
pathways_kegg <- as.data.frame(kegg@result[["Description"]])
write.csv(as.matrix(pathways_kegg),file="pathways_kegg.csv",quote = F,row.names = T)

# Define custom theme for the plot
mytheme <- theme(
  panel.border = element_blank(),  # Remove panel border
  panel.grid.major = element_line(linetype = 'dotted', colour = '#808080'),  # Set dotted grid lines
  panel.grid.major.y = element_blank(),  # Remove major grid lines on the y-axis
  panel.grid.minor = element_blank(),  # Remove minor grid lines
  axis.line.x = element_line()  # Add x-axis line
)

# Create a KEGG plot using ggplot
kegg_plot <- ggplot(kegg, showCategory = 25, aes(NES, fct_reorder(Description, NES), fill = pvalue)) +
  geom_col() +  # Add columns based on the data
  geom_segment(mapping = aes(x = -2.1, xend = ifelse(sign(NES) > 0, 0, NES), yend = Description)) +
  scale_x_continuous(expand = c(0, 0)) +  # Set x-axis limits
  scale_fill_gradientn(colours = c('#b3eebe', "#46bac2", '#371ea3'), guide = guide_colorbar(reverse = TRUE)) +  # Set fill gradient colors
  theme_dose(12) +  # Apply a specific theme variant called "theme_dose" with size 12
  mytheme +  # Apply the custom theme defined earlier
  xlab("Normalized Enrichment Score") +  # Set x-axis label
  ylab(NULL) +  # Remove y-axis label
  ggtitle("KEGG Enriched Pathways")  # Set plot title

# Display the plot
plot(kegg_plot)


#4.3. Gene Set Enrichment usig GSEA (MSigDb) ----

#Prepare Input
#no need to repeat the following steps if done before (from reading: name the vetor)
# reading in data from deseq2
df = read.csv("res.degs.csv", header=TRUE)
names(df)[names(df) == "X"] <- "GeneSymbol"

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$GeneSymbol

# Convert gene IDs for ENTREZID
ids = bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# remove duplicate IDS (here we use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has the respective entrez IDs for the gene symbols.
colnames(dedup_ids) = c("GeneSymbol", "EntrezID")
df3 = merge(df, dedup_ids, by = "GeneSymbol")

# Create a vector of the gene universe
MSigDb_gene_list = df3$log2FoldChange

# Name vector with ENTREZ ids
names(MSigDb_gene_list) = df3$EntrezID

# omit any NA values 
MSigDb_gene_list = na.omit(MSigDb_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
MSigDb_gene_list = sort(MSigDb_gene_list, decreasing = TRUE)

#Retrieve specific collection
C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C2_t2g)

# for reproducibility
set.seed(50)
# Perform Gene Set Enrichment Analysis using GSEA
msig <- GSEA(geneList = MSigDb_gene_list,TERM2GENE = C2_t2g,
            pvalueCutoff = 0.05,
            eps = 1e-10,
            by = "fgsea",
            pAdjustMethod = "none")
# Display the first few rows of the GSEA Gene Set Enrichment Analysis results
head(msig)

# Generate a GSEA plot of the GSEA Gene Set Enrichment Analysis results for the shared (with kegg) N gene sets 
gsea_plot <- gseaplot2(msig, geneSetID = c(42,77,95))
print(gsea_plot)

#merge all Gene Set Enrichment Plot
cowplot::plot_grid(go_plot_desired, kegg_plot,print(gsea_plot), labels=c("A", "B","C"))

# 4.4. Gene Set Enrichment using WikiPathway ----
#Prepare Input
#no need to repeat the following steps if done before (from reading: name the vetor)
# reading in data from deseq2
df = read.csv("res.degs.csv", header=TRUE)
names(df)[names(df) == "X"] <- "GeneSymbol"

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$GeneSymbol

# Convert gene IDs for ENTREZID
ids = bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# remove duplicate IDS (here we use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has the respective entrez IDs for the gene symbols.
colnames(dedup_ids) = c("GeneSymbol", "EntrezID")
df4 = merge(df, dedup_ids, by = "GeneSymbol")

# Create a vector of the gene universe
WP_gene_list = df4$log2FoldChange

# Name vector with ENTREZ ids
names(WP_gene_list) = df4$EntrezID

# omit any NA values 
WP_gene_list = na.omit(WP_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
WP_gene_list = sort(WP_gene_list, decreasing = TRUE)

get_wp_organisms()

# for reproducibility
set.seed(50)

WP <- gseWP(geneList =WP_gene_list,
            pvalueCutoff = 0.05,
            organism = "Homo sapiens",
            pAdjustMethod = "none")


# 4.5. Gene Set Enrichment using Reactome ----
#Prepare Input
#no need to repeat the following steps if done before (from reading: name the vetor)
# reading in data from deseq2
df = read.csv("res.degs.csv", header=TRUE)
names(df)[names(df) == "X"] <- "GeneSymbol"

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$GeneSymbol

# Convert gene IDs for ENTREZID
ids = bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# remove duplicate IDS (here we use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has the respective entrez IDs for the gene symbols.
colnames(dedup_ids) = c("GeneSymbol", "EntrezID")
df5 = merge(df, dedup_ids, by = "GeneSymbol")

# Create a vector of the gene universe
Reactome_gene_list = df5$log2FoldChange

# Name vector with ENTREZ ids
names(Reactome_gene_list) = df5$EntrezID

# omit any NA values 
Reactome_gene_list = na.omit(Reactome_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
Reactome_gene_list = sort(Reactome_gene_list, decreasing = TRUE)

# for reproducibility
set.seed(50)

Reactome <- gsePathway(geneList =Reactome_gene_list,
            pvalueCutoff = 0.05,
            minGSSize = 3, 
            maxGSSize = 800,
            organism = "human",
            eps = 1e-10,
            by = "fgsea",
            pAdjustMethod = "none")
