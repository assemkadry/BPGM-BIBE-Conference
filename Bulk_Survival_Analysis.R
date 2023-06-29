library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

#1. Survival Analysis ----

# getting clinical data for TCGA-CESC cohort -------------------
clinical_CESC <- GDCquery_clinic("TCGA-CESC")
any(colnames(clinical_CESC) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_CESC) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_CESC[,c(9,38,44)]

# looking at some variables associated with survival 
table(clinical_CESC$vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# change certain values the way they are encoded
clinical_CESC$deceased <- ifelse(clinical_CESC$vital_status == "Alive", FALSE, TRUE)
table(clinical_CESC$deceased)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_CESC$overall_survival <- ifelse(clinical_CESC$vital_status == "Alive",
                                         clinical_CESC$days_to_last_follow_up,
                                         clinical_CESC$days_to_death)

# get gene expression data -----------

# build a query to get gene expression data for entire cohort
query_CESC_all = GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor","Metastatic"),
  access = "open")

output_CESC <- getResults(query_CESC_all)


# # get gene expression data from ALL tumors 
query_CESC <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Metastatic"),
  access = "open")

# download data, do not run if you downloaded it before
GDCdownload(query_CESC)

# get counts
tcga_CESC_data <- GDCprepare(query_CESC, summarizedExperiment = TRUE)
CESC_matrix <- assay(tcga_CESC_data, "unstranded")
CESC_matrix[1:10,1:10]


# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_CESC_data))
coldata <- as.data.frame(colData(tcga_CESC_data))


# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = CESC_matrix,
                              colData = coldata,
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# vst 
vsd <- vst(dds, blind=FALSE)
CESC_matrix_vst <- assay(vsd)
CESC_matrix_vst[1:10,1:10]


# Get data for BPGM  gene and add gene metadata information to it -------------
CESC_BPGM  <- CESC_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "BPGM")

# get median value
median_value <- median(CESC_BPGM$counts)

# denote which cases have higher or lower expression than median count
CESC_BPGM$strata <- ifelse(CESC_BPGM$counts >= median_value, "HIGH", "LOW")

# Add clinical information to CESC_BPGM 
CESC_BPGM <- CESC_BPGM %>%
  mutate(case_id = gsub('-01.*', '', case_id)) %>%
  left_join(clinical_CESC, by = c("case_id" = "submitter_id"))

# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = CESC_BPGM)

fit
ggsurvplot(fit,
           data = CESC_BPGM ,
           pval = T,
           risk.table = T)

fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = CESC_BPGM)

#2. Wilcoxon Test ----

#####################################################################
#
#       Null Hypothesis: The median difference between pairs of observations is zero
#  Alternate Hypothesis: The median difference between pairs of observations is not zero
#
#    Alpha = 0.05
#
#####################################################################

expression_CESC_BPGM <- read.csv(file = "/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/Expression Data/expression_BPGM.csv", header = TRUE, sep = ",")
print(expression_CESC_BPGM)
#
# Check for Normality
shapiro.test(expression_CESC_BPGM$Cancer)     # p < 0.05, data are not normal
shapiro.test(expression_CESC_BPGM$Normal)   # p > 0.05, data are normal

# Perform Mann-Whitney Test using "wilcox.test()" function
wilcox.test(expression_CESC_BPGM$Cancer, expression_CESC_BPGM$Normal, paired =  FALSE)
wilcox.test(expression_CESC_BPGM$Normal, expression_CESC_BPGM$Cancer, paired =  FALSE)

#   Report: W = 0.9734, p = 6.498e-09
#   Decision: Reject Null Hypothesis
# Conclusion: There is a significant difference in Expression of BPGM between Cancer and Normal sample
#

#3. Expression Plot ----

library(ggplot2)
library(data.table)

# Reading in the data with fread from data.table
expression_CESC_BPGM <- fread("/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis/Expression Data/expression_CESC_BPGM.txt")

# Selecting columns we need
data <- expression_CESC_BPGM[, .(Sample, BPGM_expression_value)]

# Renaming columns
colnames(data) <- c("Sample", "Expression")

# Create a factor variable to differentiate cancer and normal samples
data$Sample <- factor(data$Sample, levels = c("CESC", "normal"),
                      labels = c("Tumor", "Normal"))

# Plot
BPGM_plot <- ggplot(data, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, shape = 21) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "pointrange", shape = 17, color = "white") +
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

# Visualize the resulted plot
plot(BPGM_plot)


# Save the plot as an image
ggsave("violin_plot.png", BPGM_plot, width = 10, height = 7)



