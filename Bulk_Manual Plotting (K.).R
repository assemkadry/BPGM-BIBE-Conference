# rm(list = ls())

library(rstudioapi)
library (ggplot2)

# change the working directory to the script path
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))

# Create a list of the txt files in the working directory
file_list <- list.files(pattern = "\\.txt$")

# Create a list to store the extracted values for each label
Data <- list()

# Create a data frame to hold all the values in coordinates for plotting
Dataframe <- data.frame()

# Create empty lists to evaluate the means' lines coordinates
x_vals <- list()
y_vals <- list()

# Create empty lists to evaluate x_axis labels
point_value <- NULL
point_label <- NULL

# Create empty lists to evaluate secondary x_axis labels
sec_point_value <- NULL
sec_point_label <- NULL

# Create variables for the values x-axis range
start <- 0
end <- 2

# Loop through each file name in the file list
for (file_df in file_list) 
{
  # Read the file using its file name
  file_df <- read.table(file_df)
  
  # Extract the label from the first column of the data frame (as long as it's not "normal" or "Sample")
  label <- unique(file_df[file_df[, 1] != "normal" & file_df[, 1] != "Sample", 1])[1]
  
  # Create variables to store the values for each label
  tumor_values <- NULL
  normal_values <- NULL
  
  # Loop through each row of the dataset starting from the second row
  for (i in 2:(nrow(file_df)-1)) 
  {
    # Check the label of the current row
    if (file_df[2:nrow(file_df), 1][i] == "normal") 
    {
      # If the label is "normal", add the value from the third column to the normal_values variable
      normal_values <- c(normal_values, file_df[2:nrow(file_df), 3][i])
      
    } 
    else if (file_df[2:nrow(file_df), 1][i] == label) 
    {
      # If the label matches the label extracted from the first column, add the value from the third column to the tumor_values variable
      tumor_values <- c(tumor_values, file_df[2:nrow(file_df), 3][i])
    }
  }
  
  # Transform the values by log2 and sort them in ascending order
  tumor_values <- sort(log2(as.numeric(tumor_values) + 1))
  normal_values <- sort(log2(as.numeric(normal_values) + 1))
  
  # Assign x-axis values with proper spacing
  tumor_values_x_axis <- seq(start, end, by = 2 / (length(tumor_values) - 1))
  normal_values_x_axis <- seq(start + 2.5, end + 2.5, by = 2 / (length(normal_values) - 1))

  # Calculate the mean of the scaled values
  tumor_values_mean <- mean(tumor_values)
  normal_values_mean <- mean(normal_values)

  # Add the extracted values for each label to the Data
  Data[[label]][["tumor"]] <- list(values_scaled = tumor_values, values_x_axis = tumor_values_x_axis, values_mean = tumor_values_mean)
  Data[[label]][["normal"]] <- list(values_scaled = normal_values, values_x_axis = normal_values_x_axis, values_mean = normal_values_mean)
  
  # Combine the scaled values and x-axis values into data frames
  tumor_df <- data.frame(values_scaled = tumor_values, values_x_axis = tumor_values_x_axis, label = label, type = "tumor")
  normal_df <- data.frame(values_scaled = normal_values, values_x_axis = normal_values_x_axis, label = label, type = "normal")
  
  # Append the data frames to the normal and tumor data frames
  Dataframe <- rbind(Dataframe, tumor_df, normal_df)
  
  # Append the new lines' coordinates
  x_vals <- c(x_vals, (start+end)/2, (start+end+5)/2)
  y_vals <- c(y_vals, tumor_values_mean, normal_values_mean)

  # Append the new x_axis labels
  point_value <- c(point_value, start, start + 2.5)
  point_label <- c(point_label, sprintf("T (n=%s)", length(tumor_values)), sprintf("N (n=%s)", length(normal_values)))

  # Append the new secondary x_axis labels
  sec_point_value <- c(sec_point_value, end + 0.5)
  sec_point_label <- c(sec_point_label, label)
  
  # Assign the new range for the new file
  start <- start + 6
  end <- end + 6
}

# Create the scatterplot and store it
scatterplot <- ggplot(Dataframe, aes(x = values_x_axis, y = values_scaled, color = type)) +
  geom_point() +
  scale_color_manual(values = c("green", "red")) +
  theme_bw()
  
# Loop on the ranges that has the mean coordinates
lines_list <- lapply(seq_along(x_vals), function(i) {
  geom_segment(aes(x = x_vals[[i]] - 0.5, xend = x_vals[[i]] + 0.5, y = y_vals[[i]], yend = y_vals[[i]]), linewidth = 1, color = "black")
})

# Add the lines to the scatterplot and modify the labels
scatterplot <- scatterplot + lines_list + 
  ylab("Transcripts Per Million (TPM)
       log2 (TPM + 1)") + 
  scale_x_continuous(breaks = point_value, labels = point_label, name = NULL, expand = c(0, 0), 
                     sec.axis = sec_axis(trans = ~., breaks = sec_point_value, labels = sec_point_label)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x.top = element_text(angle = 45, hjust = 0))
  
# Display the final plot
scatterplot