# --------------------------------------------
# Bulk RNA-Seq Analysis of TCGA Data (Survival analysis and wilcoxon Test)
# Author: Assem K. Elsherif
# Date September 25, 2023
# --------------------------------------------

# ------------------------
# Required Libraries
# ------------------------


# Load required libraries
library(TCGAbiolinks)         # For querying TCGA data
library(survminer)            # For Kaplan-Meier survival plotting
library(survival)             # For survival analysis models
library(SummarizedExperiment) # For working with TCGA gene expression format
library(tidyverse)            # For data manipulation
library(DESeq2)               # For normalization and filtering
library(ggplot2)
library(data.table)

# --------------------------------------------
# 1. Download and Prepare Clinical Data
# --------------------------------------------

# Retrieve clinical data for TCGA-CESC patients
clinical_CESC <- GDCquery_clinic("TCGA-CESC")

# Check which survival-related columns are present
any(colnames(clinical_CESC) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_CESC) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_CESC[, c(9, 38, 44)]

# Summary of vital status (Alive/Dead)
table(clinical_CESC$vital_status)

# Create logical column for deceased (TRUE/FALSE)
clinical_CESC$deceased <- ifelse(clinical_CESC$vital_status == "Alive", FALSE, TRUE)
table(clinical_CESC$deceased)

# Create a unified survival time column
clinical_CESC$overall_survival <- ifelse(clinical_CESC$vital_status == "Alive",
                                         clinical_CESC$days_to_last_follow_up,
                                         clinical_CESC$days_to_death)

# --------------------------------------------
# 2. Download and Prepare Gene Expression Data
# --------------------------------------------

# Build query for downloading RNA-Seq counts for TCGA-CESC tumors
query_CESC_all = GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Metastatic"),
  access = "open"
)

# Get metadata for matched samples
output_CESC <- getResults(query_CESC_all)

# Download data (only if not done before)
query_CESC <- query_CESC_all
GDCdownload(query_CESC)

# Prepare count matrix from downloaded data
tcga_CESC_data <- GDCprepare(query_CESC, summarizedExperiment = TRUE)
CESC_matrix <- assay(tcga_CESC_data, "unstranded")
CESC_matrix[1:10, 1:10]

# Extract gene annotations and sample annotations
gene_metadata <- as.data.frame(rowData(tcga_CESC_data))
coldata <- as.data.frame(colData(tcga_CESC_data))

# --------------------------------------------
# 3. VST Normalization of Gene Expression
# --------------------------------------------

# Convert raw counts to DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = CESC_matrix,
                              colData = coldata,
                              design = ~ 1)

# Filter out genes with low expression
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Apply variance-stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)
CESC_matrix_vst <- assay(vsd)
CESC_matrix_vst[1:10, 1:10]

# -------------------------------
# 4. Extract BPGM Expression and Classify Samples
# -------------------------------

# Convert expression matrix to long format and filter for BPGM
CESC_BPGM <- CESC_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "BPGM")

# Compute median expression for BPGM
median_value <- median(CESC_BPGM$counts)

# Classify samples based on high vs low BPGM expression
CESC_BPGM$strata <- ifelse(CESC_BPGM$counts >= median_value, "HIGH", "LOW")

# Match sample IDs with clinical info
CESC_BPGM <- CESC_BPGM %>%
  mutate(case_id = gsub('-01.*', '', case_id)) %>%
  left_join(clinical_CESC, by = c("case_id" = "submitter_id"))

# -------------------------------
# 5. Perform Survival Analysis Based on BPGM Level
# -------------------------------

# Fit Kaplan-Meier survival model by BPGM expression strata
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = CESC_BPGM)

# Show model summary
fit

# Plot survival curves with p-value and risk table
ggsurvplot(fit,
           data = CESC_BPGM,
           pval = TRUE,
           risk.table = TRUE)

# Log-rank test
fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = CESC_BPGM)

# -------------------------------
# 6. Perform Wilcoxon Rank-Sum Test on BPGM Expression
# -------------------------------

#####################################################################
#
#       Null Hypothesis: No difference in medians between groups
#  Alternate Hypothesis: Significant median difference exists
#               Alpha = 0.05
#
#####################################################################

# Load BPGM expression data from Cancer vs Normal samples
expression_CESC_BPGM <- read.csv(file = "/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/Expression Data/expression_BPGM.csv", header = TRUE, sep = ",")
print(expression_CESC_BPGM)

# Test for normality in Cancer and Normal samples
shapiro.test(expression_CESC_BPGM$Cancer)  # p < 0.05 → not normal
shapiro.test(expression_CESC_BPGM$Normal)  # p > 0.05 → possibly normal

# Perform Mann-Whitney U test (non-parametric)
wilcox.test(expression_CESC_BPGM$Cancer, expression_CESC_BPGM$Normal, paired = FALSE)
wilcox.test(expression_CESC_BPGM$Normal, expression_CESC_BPGM$Cancer, paired = FALSE)

# Interpretation: p < 0.05 → reject null → BPGM expression differs significantly

# -------------------------------
# 7. Plot BPGM Expression in Tumor vs. Normal Samples
# -------------------------------

library(ggplot2)
library(data.table)

# Read BPGM expression data for TCGA and GTEx combined
expression_CESC_BPGM <- fread("/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/Expression Data/TCGA_GTEx(endocervix).txt")

# Select and rename columns
data <- expression_CESC_BPGM[, .(Sample, BPGM_expression_value)]
colnames(data) <- c("Sample", "Expression")

# Relabel sample groups as Tumor vs Normal
data$Sample <- factor(data$Sample, levels = c("CESC", "normal"),
                      labels = c("Tumor", "Normal"))

# Create violin plot of BPGM expression
BPGM_plot <- ggplot(data, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_violin(trim = FALSE) +                                  # show full distribution
  geom_jitter(width = 0.2, shape = 21) +                        # add sample points
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "pointrange", shape = 17, color = "white") +  # show mean ± SD
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = "", y = "BPGM Expression - TPM", fill = "Condition")

# Display the plot
plot(BPGM_plot)

# Save the plot to file
ggsave("violin_plot.png", BPGM_plot, width = 10, height = 7)