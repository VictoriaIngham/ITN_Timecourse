library(ggplot2)
library(dplyr)
library(tidyr)
library(dplyr)
library(tibble)
library(readr)
library(DESeq2)

# Define the folder containing your comparison files
folder_path <- "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/ImpulseDE"

# Define the list of comparison names corresponding to your files
comparisons <- c("CFPR_vs_Con", "IG2_vs_Con", "IG1_vs_Con")

# Initialize an empty list to store the data frames
impulseDE2_results <- list()

# Loop through each comparison file and read it
for (comparison in comparisons) {
  # Construct the file path for the current comparison
  file_path <- file.path(folder_path, paste0(comparison, ".txt"))
  
  # Read the file (assuming it's tab-delimited and contains 'Gene.ID' and 'p_value' columns)
  temp_data <- read.delim(file_path)
  
  # Add a column indicating the comparison
  temp_data <- temp_data %>%
    mutate(Comparison = comparison)
  
  # Append the data frame to the list
  impulseDE2_results[[comparison]] <- temp_data
}

# Combine all the results into one data frame
combined_results <- bind_rows(impulseDE2_results)


# Define significance threshold
pval_threshold <- 0.05

# Subset for sigmoidal genes (Monotonous, not transient)
sigmoidal <- combined_results %>%
  filter(isMonotonous == TRUE, isTransient == FALSE, padj < pval_threshold)

# # Upregulated: Higher expression in case than control
# Sigmoidal_up <- sigmoidal %>%
#   filter(sigmoidTOconst_padj < pval_threshold)
# 
# # Downregulated: Lower expression in case than control
# Sigmoidal_down <- sigmoidal %>%
#   filter(impulseTOsigmoid_padj < pval_threshold)

# Subset for pulse genes (Transient, not Monotonous)
pulse <- combined_results %>%
  filter(isTransient == TRUE, isMonotonous == FALSE, padj < pval_threshold)

# # Upregulated pulses
# Impulse_up <- pulse %>%
#   filter(sigmoidTOconst_padj < pval_threshold)
# 
# # Downregulated pulses
# Impulse_down <- pulse %>%
#   filter(impulseTOsigmoid_padj < pval_threshold)

# # Write out the results
# write_tsv(Sigmoidal_up, file.path(folder_path, "Sigmoidal_upvsCon.txt"))
# write_tsv(Sigmoidal_down, file.path(folder_path, "Sigmoidal_downvsCon.txt"))
# write_tsv(Impulse_up, file.path(folder_path, "Impulse_upvsCon.txt"))
# write_tsv(Impulse_down, file.path(folder_path, "Impulse_downvsCon.txt"))







# Function to find genes appearing in multiple comparisons
find_common_genes <- function(df) {
  df %>%
    group_by(Gene) %>%
    summarise(Count = n(), Comparisons = paste(unique(Comparison), collapse = ", ")) %>%
    filter(Count > 1)  # Keep only genes appearing in multiple comparisons
}
# 
# # Identify overlapping genes for each category
# Sigmoidal_up_common <- find_common_genes(Sigmoidal_up)
# Sigmoidal_down_common <- find_common_genes(Sigmoidal_down)
# Impulse_up_common <- find_common_genes(Impulse_up)
# Impulse_down_common <- find_common_genes(Impulse_down)

##ignore the up and down
Sigmoidal_common = find_common_genes(sigmoidal)
Impulse_common = find_common_genes(pulse)

# # Write out the results
# write_tsv(Sigmoidal_up_common, file.path(folder_path, "Sigmoidal_up_common_vsCon.txt"))
# write_tsv(Sigmoidal_down_common, file.path(folder_path, "Sigmoidal_down_commonvsCon.txt"))
# write_tsv(Impulse_up_common, file.path(folder_path, "Impulse_up_commonvsCon.txt"))
# write_tsv(Impulse_down_common, file.path(folder_path, "Impulse_down_commonvsCon.txt"))
# 
# 
# write_tsv(Impulse_common, file.path(folder_path, "Impulse_all_commonvsCon.txt"))
# write_tsv(Sigmoidal_common, file.path(folder_path, "Sigmoidal_all_commonvsCon.txt"))

# Load your data (replace with your file paths)
counts <- read.delim("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/combined_table_clean1.txt", row.names = 1)
colData <- read.delim("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/colData_2.txt") # Sample info: treatment, time, replicate

# Ensure colData is correctly formatted as factors
colData$treatment <- factor(colData$treatment, levels = c("Control", "CFPR", "PYR", "Combined"))
colData$time <- factor(colData$time, levels = c("1", "4", "8", "24", "48", "72"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ time + treatment)

# Perform variance-stabilizing transformation (VST)
vst_counts <- assay(vst(dds, blind = FALSE))

# Convert normalized counts to long format for ggplot
tidy_counts <- as.data.frame(vst_counts) %>%
  rownames_to_column(var = "Geneid") %>%
  pivot_longer(cols = -Geneid, names_to = "colData", values_to = "Count") %>%
  left_join(colData, by = "colData") 

goi_table = read.delim('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/ImpulseDE/Impulse_GOI_list.txt',header=T)

# Filter common gene tables for ATP synthase-related genes
filter_for_goi <- function(df) {
  df %>% filter(Gene %in% goi_table$Gene.ID)
}

# sigmoidal_up_goi <- filter_for_goi(Sigmoidal_up_common)
# sigmoidal_down_goi <- filter_for_goi(Sigmoidal_down_common)
# impulse_up_goi <- filter_for_goi(Impulse_up_common)
# impulse_down_goi <- filter_for_goi(Impulse_down_common)

##again ignore up and down
sigmoidal_goi <- filter_for_goi(Sigmoidal_common)
impulse_goi <- filter_for_goi(Impulse_common)


custom_colors <- c("Control" = "#4a4747", "CFPR" = "#E69F00", "PYR" = "#009E73", "Combined" = "#CC79A7")


# Function to plot expression curves for a gene set
plot_expression <- function(gene_list, title) {
  plot_data <- tidy_counts %>% filter(Geneid %in% gene_list$Gene)

  if (nrow(plot_data) == 0) {
    print(paste("No genes found for", title))
    return(NULL)
  }

  plot_data <- plot_data %>%
    mutate(TimeNumeric = as.numeric(as.character(time)))

  # Ensure there is a 'Comparison' column that provides the correct info for each gene
  plot_data <- plot_data %>%
    left_join(gene_list %>%
                select(Gene, Comparisons), by = c("Geneid" = "Gene"))

  # Create a label column specifically for each gene's comparison
  plot_data$comparison_label <- plot_data$Comparisons

  # Now, let's create the plot
  ggplot(plot_data, aes(x = TimeNumeric, y = Count, color = treatment, group = interaction(Geneid, treatment))) +
    geom_smooth(aes(fill = treatment), method = "loess", se = TRUE, alpha = 0.2) +
    facet_wrap(~ Geneid, scales = "free_y") +  # Facet by Geneid
    scale_x_continuous(breaks = c(1, 4, 8, 24, 48, 72), labels = c("1", "4", "8", "24", "48", "72")) +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    labs(
      title = paste("Expression Trends:", title),
      x = "Hours Post Exposure",
      y = "Normalized Expression"
    ) +
    theme_minimal() +
    # Add comparison text for each facet (gene) in the top-right corner
    geom_text(aes(x = 72, y = max(Count), label = comparison_label),
              data = plot_data, inherit.aes = FALSE,
              hjust = 1, vjust = 1, size = 3, color = "black", fontface = "italic")
}

##Ignore up and down
plot_expression(sigmoidal_goi, "Sigmoidal")
plot_expression(impulse_goi, "Impulse")



###For saving
plot_expression <- function(gene_list, title) {
  plot_data <- tidy_counts %>% filter(Geneid %in% gene_list$Gene)
  
  if (nrow(plot_data) == 0) {
    print(paste("No genes found for", title))
    return(NULL)
  }
  
  plot_data <- plot_data %>%
    mutate(TimeNumeric = as.numeric(as.character(time))) %>%
    left_join(gene_list %>% select(Gene, Comparisons), by = c("Geneid" = "Gene")) %>%
    mutate(comparison_label = Comparisons)
  
  # Create the plot object
  p <- ggplot(plot_data, aes(x = TimeNumeric, y = Count, color = treatment, group = interaction(Geneid, treatment))) +
    geom_smooth(aes(fill = treatment), method = "loess", se = TRUE, alpha = 0.2) +
    facet_wrap(~ Geneid, scales = "free_y") +
    scale_x_continuous(breaks = c(1, 4, 8, 24, 48, 72), labels = c("1", "4", "8", "24", "48", "72")) +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    labs(
      title = paste("Expression Trends:", title),
      x = "Hours Post Exposure",
      y = "Normalized Expression"
    ) +
    theme_minimal() +
    geom_text(aes(x = 72, y = max(Count), label = comparison_label), 
              data = plot_data, inherit.aes = FALSE,
              hjust = 1, vjust = 1, size = 3, color = "black", fontface = "italic")
  
  return(p)  # 👈 Return the plot object
}

# Save sigmoidal plot
p1 <- plot_expression(sigmoidal_goi, "Sigmoidal")
if (!is.null(p1)) {
  ggsave("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Figures/Impulse DE figures/Sigmoidal_all.pdf", p1, width = 10, height = 6, dpi = 300,device="pdf")
}

# Save impulse plot
p2 <- plot_expression(impulse_goi, "Impulse")
if (!is.null(p2)) {
  ggsave("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Figures/Impulse DE figures/Impulse_all.pdf", p2, width = 10, height = 6, dpi = 300,device="pdf")
}

# # Generate plots for genes appearing in multiple comparisons
# plot_expression(sigmoidal_up_goi, "Sigmoidal Upregulated")
# plot_expression(sigmoidal_down_goi, "Sigmoidal Downregulated")
# plot_expression(impulse_up_goi, "Impulse Upregulated")
# plot_expression(impulse_down_goi, "Impulse Downregulated")
