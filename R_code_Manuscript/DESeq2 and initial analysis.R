# Load necessary libraries
library(dplyr)

# Set the folder path
folder_path <- "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Count Files"  # Replace with your folder path

# List all .txt files in the folder
file_list <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)

# Read all files and combine them column-wise
combined_table <- file_list %>%
  lapply(read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  bind_cols()

# View the combined table
print(combined_table)

# Optionally, save the combined table as a file
write.table(combined_table, "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/combined_table.txt", sep = "\t", row.names = FALSE)




# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)

# Load your data (replace with your file paths)
counts <- read.delim("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/combined_table_clean1.txt", row.names = 1)
colData <- read.delim("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/colData_2.txt") # Sample info: treatment, time, replicate

# Ensure colData is correctly formatted as factors
colData$treatment <- factor(colData$treatment, levels = c("Control", "CFPR", "PYR", "Combined"))
colData$time <- factor(colData$time, levels = c("1", "4", "8", "24", "48", "72"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ time + treatment)

# Perform variance-stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Extract PCA for all samples
pca_data <- plotPCA(vsd, intgroup = c("time", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# General PCA Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment, shape = time)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of RNAseq Data", 
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  theme(plot.title = element_text(hjust = 0.5))

# PCA Analysis for Individual Time Points
time_points <- levels(colData$time)
pca_list <- list()

for (t in time_points) {
  # Subset VST data for the specific time point
  subset_vsd <- vsd[, colData$time == t]
  
  # Perform PCA
  pca_data_time <- plotPCA(subset_vsd, intgroup = "treatment", returnData = TRUE)
  percentVar_time <- round(100 * attr(pca_data_time, "percentVar"))
  
  # Add sample names to the PCA data
  pca_data_time$sample <- rownames(colData(subset_vsd))
  
  # Define custom colors for treatments
  custom_colors <- c("Control" = "#999999", "CFPR" = "#E69F00", "PYR" = "#009E73", "Combined" = "#CC79A7")
  
  # Plot PCA for the current time point with sample names
  p <- ggplot(pca_data_time, aes(x = PC1, y = PC2, color = treatment)) +
    geom_point(size = 3) +
    geom_text(aes(label = sample), hjust = 0.5, vjust = -0.5, size = 3) +  # Add sample names
    theme_minimal() +
    labs(title = paste("PCA at Time Point:", t),
         x = paste0("PC1: ", percentVar_time[1], "% variance"),
         y = paste0("PC2: ", percentVar_time[2], "% variance")) +
    scale_color_manual(values = custom_colors) +  # Use custom colors
    theme(plot.title = element_text(hjust = 0.5))+
    theme_light()
  
  # Save plot to list
  pca_list[[t]] <- p
}

# Save PCA plots for each time point
for (t in time_points) {
  ggsave(filename = paste0("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/PCA_TimePoint_", t, ".png"), plot = pca_list[[t]], width = 8, height = 6)
}

# Display PCA plots in RStudio viewer
pca_list

# Load necessary libraries
library(DESeq2)
library(dplyr)

# Load data
counts <- read.table("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/combined_table_clean1.txt", header = TRUE, row.names = 1)
colData <- read.table("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/colData_2.txt", header = TRUE, row.names = 1)

    # Specify comparison parameters
time_point <- "72"   # Specify the time point
treatment <- "Combined"   # Specify the treatment to compare
p_adj_threshold <- 0.05   # Adjusted p-value threshold for significance
log2_fc_threshold <- 0    # Log2 fold change threshold for up-/down-regulation

# Subset colData to include the selected time point and relevant treatment
filtered_colData <- colData %>%
  filter(time == time_point & (treatment == "Control" | treatment == !!treatment))

# Ensure rownames of colData match column names in counts
filtered_counts <- counts[, rownames(filtered_colData)]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = filtered_colData,
  design = ~ treatment
)

# Perform DESeq2 analysis
dds <- DESeq(dds)

# Get results for the treatment vs control comparison
res <- results(dds, contrast = c("treatment", treatment, "Control"))

# Filter significant results
sig_genes <- res[which(res$padj < p_adj_threshold & !is.na(res$padj)), ]
upregulated <- sig_genes[which(sig_genes$log2FoldChange > log2_fc_threshold), ]
downregulated <- sig_genes[which(sig_genes$log2FoldChange < -log2_fc_threshold), ]

# Summary of results
total_sig <- nrow(sig_genes)   # Total significant genes
num_upregulated <- nrow(upregulated)  # Upregulated genes
num_downregulated <- nrow(downregulated)  # Downregulated genes

# Print summary
cat("Summary for", treatment, "vs Control at", time_point, ":\n")
cat("Total significant genes:", total_sig, "\n")
cat("Upregulated genes:", num_upregulated, "\n")
cat("Downregulated genes:", num_downregulated, "\n")

# Save results to files
write.csv(as.data.frame(res), file = paste0('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/',treatment, "_vs_Control_", time_point, "_full_results.csv"))
#write.csv(as.data.frame(sig_genes), file = paste0('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/',treatment, "_vs_Control_", time_point, "_significant_results.csv"))

# Optional: Visualize results
plotMA(res, main = paste0("MA Plot: ", treatment, " vs Control at ", time_point))






# Optional: Visualize results
plotMA(res, main = paste0("MA Plot: ", treatment, " vs Control at ", time_point))


# Create DESeq2 dataset with time component = more complex model
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ treatment + time + treatment:time)

# Relevel the Condition factor to ensure "Control" is the baseline
dds$Condition <- relevel(dds$treatment, ref = "Control")

# Run DESeq2 normalization and dispersion estimation
dds <- DESeq(dds)

# Identify rows where betaConv is explicitly FALSE (non-converged genes)
non_converged_genes <- rownames(dds)[mcols(dds)$betaConv == FALSE]

# Print the non-converged genes
cat("Non-converged genes:\n", non_converged_genes)

# Initialize list to store results
results_list <- list()

# Retain only genes with a defined betaConv (TRUE or FALSE)
dds <- dds[!is.na(mcols(dds)$betaConv), ]

# Optionally, exclude only the NA and FALSE cases
dds <- dds[mcols(dds)$betaConv == TRUE, ]

# Initialize list to store results
results_list <- list()

# Perform differential expression analysis for each treatment at each time point
conditions <- c("CFPR", "PYR", "Combined")
time_points <- c("1", "4", "8", "24", "48", "72")

# Initialize a data frame to store the summary of DEGs
deg_summary <- data.frame(
  Condition = character(),
  Time = character(),
  Total_DEGs = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  stringsAsFactors = FALSE
)

for (cond in conditions) {
  for (time in time_points) {
    if (time == "1") {
      # For time = 1, compare treatments to Control using main effect of Condition
      res <- results(dds, contrast = c("treatment", cond, "Control"))
    } else {
      # For other time points, use the interaction term
      interaction_name <- paste0("treatment", cond, ".time", time)
      
      # Check if the interaction term exists
      if (interaction_name %in% resultsNames(dds)) {
        res <- results(dds, name = interaction_name)
      } else {
        message("Interaction term not found for ", cond, " at time ", time)
        next
      }
    }
    
    # Save results (optional)
    file_name <- paste0("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/DESeq2_", cond, "_vs_Control_Time_", time, ".csv")
    write.csv(as.data.frame(res), file = file_name, row.names = TRUE)
    
    # Count DEGs
    total_degs <- sum(res$padj < 0.05, na.rm = TRUE)
    upregulated <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
    downregulated <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
    
    # Append to summary table
    deg_summary <- rbind(deg_summary, data.frame(
      Condition = cond,
      Time = time,
      Total_DEGs = total_degs,
      Upregulated = upregulated,
      Downregulated = downregulated
    ))
  }
}

write.csv(deg_summary, file = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/DEG_Summary_Table.csv", row.names = FALSE)


# Function to read in DESeq2 results and filter for overexpressed genes (padj < 0.05 and log2FoldChange > 0)
read_deseq_results <- function(file_path) {
  res <- read.csv(file_path, row.names = 1)  # Read the CSV and set the row names as gene IDs
  res <- res[!is.na(res$padj), ]  # Remove rows with NA in the padj column
  res_overexpressed <- rownames(res[res$padj < 0.05 & res$log2FoldChange > 0, ])  # Filter for overexpressed genes
  return(res_overexpressed)
}

# Specify the file paths for each treatment and time point
files_CFPR <- list(
  "1" = "DESeq2_CFPR_vs_Control_Time_1.csv",
  "4" = "DESeq2_CFPR_vs_Control_Time_4.csv",
  "8" = "DESeq2_CFPR_vs_Control_Time_8.csv",
  "24" = "DESeq2_CFPR_vs_Control_Time_24.csv",
  "48" = "DESeq2_CFPR_vs_Control_Time_48.csv",
  "72" = "DESeq2_CFPR_vs_Control_Time_72.csv"
)

# Read in the overexpressed genes for CFPR at each time point
overexpressed_CFPR <- lapply(files_CFPR, read_deseq_results)

# Remove empty sets from the list before plotting
overexpressed_CFPR <- overexpressed_CFPR[sapply(overexpressed_CFPR, length) > 0]




# Specify the file paths for each time point across treatments
files_timepoints <- list(
  "CFPR" = list(
    "1" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/CFPR_vs_Control_1_full_results.csv",
    "4" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/CFPR_vs_Control_4_full_results.csv",
    "8" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/CFPR_vs_Control_8_full_results.csv",
    "24" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/CFPR_vs_Control_24_full_results.csv",
    "48" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/CFPR_vs_Control_48_full_results.csv",
    "72" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/CFPR_vs_Control_72_full_results.csv"
  ),
  "PYR" = list(
    "1" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/PYR_vs_Control_1_full_results.csv",
    "4" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/PYR_vs_Control_4_full_results.csv",
    "8" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/PYR_vs_Control_8_full_results.csv",
    "24" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/PYR_vs_Control_24_full_results.csv",
    "48" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/PYR_vs_Control_48_full_results.csv",
    "72" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/PYR_vs_Control_72_full_results.csv"
  ),
  "Combined" = list(
    "1" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Combined_vs_Control_1_full_results.csv",
    "4" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Combined_vs_Control_4_full_results.csv",
    "8" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Combined_vs_Control_8_full_results.csv",
    "24" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Combined_vs_Control_24_full_results.csv",
    "48" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Combined_vs_Control_48_full_results.csv",
    "72" = "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Combined_vs_Control_72_full_results.csv"
  )
)

# Read in the overexpressed genes for each time point across treatments
overexpressed_timepoints <- lapply(files_timepoints, function(treatment_files) {
  lapply(treatment_files, read_deseq_results)
})

overexpressed_timepoints = overexpressed_timepoints[sapply(overexpressed_timepoints, length) > 0]

library(VennDiagram)
# Create the Venn Diagram for Time Point 1
venn.plot_time1 <- venn.diagram(
  x = list(
    CFPR = overexpressed_timepoints$CFPR$`72`,
    PYR = overexpressed_timepoints$PYR$`72`,
    Combined = overexpressed_timepoints$Combined$`72`
  ),
  category.names = c("CFPR", "IG1", "IG2"),
  filename = NULL,
  output = TRUE,
  col = "black",
  fill = c("#0072B2", "#E69F00", "#CC79A7"),
  alpha = 0.5,
  cex = 1,
  main="72 hour",
  main.fontface = "bold",
  main.fontfamily = "sans",
  fontface = "bold",
  fontfamily = "sans",
  main.cex = 1.5,
  label.col = "black",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
  
)

grid::grid.newpage()
# Plot the Venn diagram for Time Point 1
grid.draw(venn.plot_time1)

# Function to write overlapping genes to file
write_overlap_genes <- function(overlap_set, file_name) {
  write.csv(overlap_set, file = file_name, row.names = FALSE)
}

overlaps=c()
# Example for Time Point 1 - CFPR, PYR, and Combined
overlap_time1 <- Reduce(intersect, list(
  overexpressed_timepoints$CFPR$`72`,
  overexpressed_timepoints$PYR$`72`,
  overexpressed_timepoints$Combined$`72`
))

overlaps = cbind(unname(overlap_time1),rep('Three-way',length(overlap_time1)))

overlap_time1 <- Reduce(intersect, list(
  overexpressed_timepoints$CFPR$`72`,
  overexpressed_timepoints$PYR$`72`
))

overlaps = rbind(overlaps,cbind(unname(overlap_time1),rep('CFPR-PYR',length(overlap_time1))))

overlap_time1 <- Reduce(intersect, list(
  overexpressed_timepoints$CFPR$`72`,
  overexpressed_timepoints$Combined$`72`
))

overlaps = rbind(overlaps,cbind(unname(overlap_time1),rep('CFPR-Combined',length(overlap_time1))))

overlap_time1 <- Reduce(intersect, list(
  overexpressed_timepoints$Combined$`72`,
  overexpressed_timepoints$PYR$`72`
))

overlaps = rbind(overlaps,cbind(unname(overlap_time1),rep('Combined-PYR',length(overlap_time1))))

colnames(overlaps) = c('ID','Type')

# Write the overlapping genes to a CSV file
write_overlap_genes(overlaps, "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/overlap_time72_genes.csv")

