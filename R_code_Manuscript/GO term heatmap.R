library(clusterProfiler)
library(org.Ag.eg.db)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(circlize)

# Read in your data
data <- read.csv('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/GO Terms.csv')  # Replace with the path to your input CSV

# Perform GO Slim on the entire dataset
gene_list <- unique(data$GO_ID)  # Use GO IDs from your data

# Perform enrichment analysis (no need to specify the ontology if analyzing all)
ego <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Ag.eg.db,
  keyType       = "GOALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  minGSSize     = 1,  # Reduce minimum gene set size
  maxGSSize     = 1000,  # Increase maximum gene set size
  readable      = TRUE
)

# Simplify the GO terms
ego_simplified <- simplify(ego, cutoff = 0.7, by = "p.adjust", measure = "Wang")

# Extract simplified GO terms
simplified_go_terms <- ego_simplified@result$ID

# Filter the original data to retain only simplified GO terms
data_filtered <- data %>% filter(GO_ID %in% simplified_go_terms)

# Prepare the data for the heatmap: Reshape it into long format and include NA for missing values
heatmap_data <- data_filtered %>%
  select(GO_ID,Description, padj, fold, Time.point, Treatment) %>%
  pivot_longer(cols = c(padj, fold), names_to = "Metric", values_to = "Value") %>%
  unite("Time_Treatment", Time.point, Treatment, sep = "_") %>%
  pivot_wider(names_from = c(Time_Treatment, Metric), values_from = Value) %>%
  arrange(GO_ID)

# Now, we need to restructure it such that there are three rows per GO_ID for each treatment
expanded_data <- data_filtered %>%
  filter(Time.point %in% c("1hr", "4hr", "8hr", "24hr", "48hr", "72hr")) %>%
  complete(GO_ID, Time.point = c("1hr", "4hr", "8hr", "24hr", "48hr", "72hr"), Treatment = c("PYR", "Dual", "CFPR")) %>%
  arrange(GO_ID, Time.point, Treatment)

# Pivot the data so each treatment has its own column for each time point
heatmap_data_pivoted <- expanded_data %>%
  select(GO_ID, Time.point, Treatment, padj) %>%
  pivot_wider(names_from = c(Time.point, Treatment), values_from = padj) %>%
  arrange(GO_ID)  # Ensure it's sorted by GO_ID

# Now, we have the correct shape: rows for GO terms, columns for time point/treatment combos

# Create a matrix for the heatmap (e.g., p-values)
p_value_matrix = heatmap_data_pivoted[,-1]


time_levels <- c("1hr", "4hr", "8hr", "24hr", "48hr", "72hr")
                 
p_value_matrix = p_value_matrix[, order(factor(sub("_.*", "", colnames(p_value_matrix)), levels = time_levels))]     

p_value_matrix = as.matrix(p_value_matrix)



descriptions = heatmap_data[which(heatmap_data$GO_ID %in% simplified_go_terms == T),]

descriptions = descriptions[order(descriptions$GO_ID),]

rownames(p_value_matrix)  = descriptions$Description

p_value_matrix = p_value_matrix[order(rownames(p_value_matrix)),]

#write.table(p_value_matrix,'/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Reduced terms for heatmap.txt',sep='\t',row.names=T)

##Manually remove GO terms that aren't of interest in Excel


# Define column split to create the visual separation every 3 columns (for 3 treatments)
column_split <- rep(1:6, each = 3)  # 6 time points, 3 treatments per time point

p_value_matrix <- read.delim("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Reduced terms for heatmap.txt", header = TRUE, sep = "\t")

rownames(p_value_matrix) = p_value_matrix[,1]
p_value_matrix = p_value_matrix[,-1]
p_value_matrix = as.matrix(p_value_matrix)

# Heatmap customization (e.g., p-values)
Heatmap(
  p_value_matrix,
  name = "P-Value",
  na_col = "cornsilk",
  col = colorRamp2(c(0,0.0000001,0.000001, 0.0001,0.001 ), c("darkblue","blue3","blue","cadetblue","lightblue")),
  border = 'black',
  show_heatmap_legend = T,
  row_split = rownames(p_value_matrix),
  column_split = column_split,  # Split columns every 3 (treatments per time point)
  show_row_names = TRUE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_title = "GO terms for up-regulated transcripts",
  row_title = "GO Terms",
  heatmap_width = unit(20,"cm"),
  left_annotation = NULL,  # Add annotations if needed
  row_gap = unit(0.5, "mm"),  # Add space between rows (thicker line effect)
  column_gap = unit(1, "mm")  # Space between groups of columns (this is for the visual separation)
)

####Plot genes of interest
# Read in the gene lists with GO Terms

GO_terms = read.delim('/Users/vickyingham/Seafile/Wasim & Vicky/GO_term_download.txt',header=T)


library(DESeq2)

# Load RNAseq count data
counts <- read.delim('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/combined_table_clean1.txt', row.names = 1,header=T)  # Replace with the actual file path

# Load sample metadata
coldata <- read.delim('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/colData_2.txt',header=T)  # Replace with the actual metadata file path


# Create a DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ treatment)

# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)



# Define the GO term of interest
go_term_of_interest <- "GO:0006119"  # Replace with your GO term

filtered_genes <- GO_terms %>%
  filter(
    (!is.na(GO_Term) & grepl(go_term_of_interest, GO_Term))
  ) %>%
  select(Gene.ID) %>%
  distinct()

##Give a list of genes of interest
filtered_genes = read.delim('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Gene lists of interest/ATPases.txt',header=T)


# Filter normalized counts to keep only the genes of interest
filtered_genes = as.matrix(filtered_genes$Gene.ID)
counts_data <- normalized_counts[rownames(normalized_counts) %in% filtered_genes, ]

counts_data = as.data.frame(counts_data)

counts_data$Gene.ID = rownames(counts_data)

# Extract metadata from column headers
counts_long <- counts_data %>%
  pivot_longer(cols = -Gene.ID, names_to = "Sample", values_to = "Counts") %>%
  mutate(
    Treatment = sub("_.*", "", Sample),
    Timepoint = sub(".*_(\\d+hr).*", "\\1", Sample),
    Replicate = sub(".*_(R\\d+)$", "\\1", Sample)
  ) %>%
  filter(!is.na(Counts))  # Remove rows with missing counts

# Ensure time points are ordered correctly
time_order <- c("1hr", "4hr", "8hr", "24hr", "48hr", "72hr")
counts_long$Timepoint <- factor(counts_long$Timepoint, levels = time_order)




library(ggplot2)
library(tidyr)
library(dplyr)


# ggplot(counts_long, aes(x = Timepoint, y = Counts, color = Treatment, group = interaction(Treatment, Replicate))) +
#   geom_smooth(se = FALSE, method = "loess", aes(group = Treatment)) +  # Smoothed lines for each treatment
#   geom_point(size = 1, alpha = 0.6) +  # Points for each replicate
#   facet_wrap(~Gene.ID, scales = "free_y") +  # One plot per gene
#   theme_minimal() +
#   labs(
#     title = "Smoothed Expression Patterns by Treatment and Timepoint",
#     x = "Time Point",
#     y = "Normalized Expression (log2)",
#     color = "Treatment"
#   ) +
#   scale_color_manual(
#     values = c("IG1" = "#1A85FF", "IG2" = "#E66100", "CFPR" = "#5D3A9B", "Cont" = "#D41159")
#   ) +  # Add custom colors for each treatment
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "bottom"
#   )



###Try to account for differences in time steps
# Convert timepoints to numeric values for continuous spacing
counts_long <- counts_long %>%
  mutate(
    Timepoint = recode(Timepoint,
                       "1hr" = 1, 
                       "4hr" = 4, 
                       "8hr" = 8, 
                       "24hr" = 24, 
                       "48hr" = 48, 
                       "72hr" = 72
    )
  )

# Plot with continuous time
ggplot(counts_long, aes(x = Timepoint, y = Counts, color = Treatment, group = interaction(Treatment, Replicate))) +
  geom_smooth(se = TRUE, method = "loess", aes(group = Treatment, fill = Treatment), level = 0.9, alpha = 0.1) +  # Smoothed lines with 95% CI and lighter alpha
  geom_point(size = 1, alpha = 0.6) +  # Points for each replicate
  facet_wrap(~Gene.ID, scales = "free_y") +  # One plot per gene
  theme_minimal() +
  labs(
    title = "ATPase",
    x = "Time Point (hrs)",
    y = "Normalized Expression (log2)",
    color = "Treatment",
    fill = "Treatment"
  ) +
  scale_color_manual(
    values = c("IG1" = "#1A85FF", "IG2" = "#E66100", "CFPR" = "#D41159","Cont" = "black")
  ) +
  scale_fill_manual(  # Match fill color of confidence interval to line color
    values = c("IG1" = "#1A85FF", "IG2" = "#E66100", "CFPR" = "#D41159","Cont" = "black")
  ) +
  scale_x_continuous(
    breaks = c(1, 4, 8, 24, 48, 72),  # Specify major breaks for clarity
    labels = c("1hr", "4hr", "8hr", "24hr", "48hr", "72hr")  # Keep the timepoint labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )