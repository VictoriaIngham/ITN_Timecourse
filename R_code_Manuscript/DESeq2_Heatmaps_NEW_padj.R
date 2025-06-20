# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(gplots)  # For heatmap.2 clustering

# Assuming goi_table and rnaseq_data are already loaded
goi_table = read.delim('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Gene lists of interest/list_66_percent_in_one_treatment.txt', header = T)
rnaseq_data = read.csv('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Analysis of all data/final_combined_results_fixed.csv', header = T)

# Merge the genes of interest with the RNAseq data
merged_data <- inner_join(goi_table, rnaseq_data, by = c("Gene.ID" = "ID"))

# Reshape log2FoldChange and padj data
logfc_data <- melt(merged_data, id.vars = c("Gene.ID", "GeneName.x"),
                   measure.vars = grep("log2FoldChange", names(merged_data), value = TRUE),
                   variable.name = "Treatment_Time", value.name = "log2FoldChange")

padj_data <- melt(merged_data, id.vars = c("Gene.ID", "GeneName.x"),
                  measure.vars = grep("^padj", names(merged_data), value = TRUE),
                  variable.name = "Padj_Treatment_Time", value.name = "padj")

# Ensure matching order
stopifnot(all(logfc_data$Gene.ID == padj_data$Gene.ID))

# Merge the two
melted_data <- logfc_data
melted_data$padj <- padj_data$padj

# Extract Treatment and Time as before
melted_data <- melted_data %>%
  mutate(
    Treatment = str_extract(Treatment_Time, "(?<=log2FoldChange\\.)[A-Za-z]+"),
    Time = str_extract(Treatment_Time, "\\d+(?=$)")
  )

melted_data$Time <- factor(melted_data$Time, levels = c("1", "4", "8", "24", "48", "72"))


# Remove weak/non-significant changes based on adj.pval threshold
melted_data$log2FoldChange[melted_data$padj > 0.05 | is.na(melted_data$padj)] <- 0

# Remove genes with no differential expression
filtered_data <- melted_data %>%
  group_by(Gene.ID) %>%
  filter(any(log2FoldChange != 0)) %>%
  ungroup()

# Cap fold changes for color scale
filtered_data$log2FoldChange <- pmax(pmin(filtered_data$log2FoldChange, 2.5), -2.5)

# Reorder gene names based on clustering
heatmap_data <- dcast(filtered_data, GeneName.x ~ Treatment + Time, value.var = "log2FoldChange", fun.aggregate = mean)
heatmap_matrix <- as.matrix(heatmap_data[, -1])
row.names(heatmap_matrix) <- heatmap_data$GeneName.x
heatmap_matrix[is.na(heatmap_matrix)] <- 0
row_dendrogram <- hclust(dist(heatmap_matrix), method = "complete")
ordered_rows <- order.dendrogram(as.dendrogram(row_dendrogram))
filtered_data$GeneName.x <- factor(filtered_data$GeneName.x, 
                                   levels = row.names(heatmap_matrix)[ordered_rows],
                                   ordered = TRUE)

# Plot
p <- ggplot(filtered_data, aes(x = Time, y = factor(GeneName.x, levels = levels(filtered_data$GeneName.x)), fill = log2FoldChange)) +
  geom_tile(color = "lightgrey", size = 0.25) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limit = c(-2.5, 2.5), name = "log2 Fold Change") +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(~ Treatment, scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 4),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15),
        strip.background = element_rect(fill = NA, color = NA),
        panel.spacing = unit(0.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "Gene Expression Heatmap", x = "Time (hrs)", y = "Gene Name")

# Save
ggsave('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Figures/DEseq2 Heatmaps/genes_in_66_percent_out.png', plot = p, width = 10, height = 8, dpi = 300)
