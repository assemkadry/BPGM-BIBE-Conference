library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

#single data frame
setwd("/Users/assemkadry/Desktop/Bulk_RNA_Seq_Analysis//Expression Data")
# set the directory where the files are located
dir <- "~/Desktop/Bulk_RNA_Seq_Analysis//Expression Data"

# get a list of all files in the directory that match the pattern "expression_*_TFAP2A.txt"
file_list <- list.files(dir, pattern = "expression_.*_TFAP2A.txt")

# create an empty data frame to store the data
expression_df <- data.frame()

# loop over each file in the list
for (file in file_list) {
  
  # extract the cancer type from the file name
  cancer <- gsub("expression_(\\w+)_TFAP2A.txt", "\\1", file)
  
  # read in the data for the current cancer type
  expression_data <- read_delim(file.path(dir, file), delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  
  # create a new column to indicate the sample group
  expression_data$Group <- ifelse(grepl(cancer, expression_data$Sample), paste(cancer, "_T", sep = ""), paste(cancer, "_N", sep = ""))
  
  # append the data to the expression_df data frame
  expression_df <- rbind(expression_df, expression_data)
}

expression_df <- expression_df %>% select(-Sample)
# Replace values in the "TFAP2A_expression_value" column with the calculated log2 values
expression_df$TFAP2A_expression_value <- log2(expression_df$TFAP2A_expression_value + 1)

# view the data
head(expression_df)

# summarize the data by group and count the number of samples in each group
sample_count <- expression_df %>%
  group_by(Group) %>%
  summarize(count = n()) %>%
  ungroup()

# create labels with sample name and count
x_labels <- paste(sample_count$Group, " (n=", sample_count$count, ")", sep="")

#z_labels <- sapply(file_list, function(file) {
#  cancer <- gsub("expression_(\\w+)_TFAP2A.txt", "\\1", file)
#  return(cancer)
#})

ggplot(expression_df, aes(x = Group, y = TFAP2A_expression_value, fill = ifelse(grepl("N\\s*$", Group), "green", "red"))) +
  geom_boxplot() +
  scale_fill_identity() +
  ggtitle("TFAP2A Expression") +
  labs(x = "", y = "Transcript Per Million (TPM)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5, size = 14)) +
  scale_x_discrete(labels = x_labels)




