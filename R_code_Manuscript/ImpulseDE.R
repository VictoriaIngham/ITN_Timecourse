# Load your data (replace with your file paths)
counts <- read.delim("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/combined_table_clean1.txt", row.names = 1)
colData <- read.delim("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/colData_2.txt") # Sample info: treatment, time, replicate

library(ImpulseDE2)
library(ggplot2)

# Prepare metadata (ensure columns: Sample, Condition, Timepoint)
annotation_table <- data.frame(
  Sample = colData$colData,
  Condition = colData$treatment,
  Time = colData$time,
  Replicate = gsub(".*_R(\\d+)", "rep\\1", colData$colData) # Extract replicate number
)

annotation_table <- annotation_table[order(annotation_table$Time), ]

count_matrix <- as.matrix(counts)

# Set a threshold for minimum gene counts
min_total_counts <- 870 # Adjust this threshold as needed

# Calculate row sums (total counts per gene)
gene_sums <- rowSums(count_matrix)

# Identify genes to keep
keep_genes <- gene_sums >= min_total_counts

# Filter the count matrix
count_matrix <- count_matrix[keep_genes, ]


# CFPR vs Control
subset_samples_cfpr_control <- annotation_table$Sample[annotation_table$Condition %in% c("CFPR", "Control")]
subset_counts_cfpr_control <- count_matrix[, subset_samples_cfpr_control]
subset_annotation_cfpr_control <- annotation_table[annotation_table$Sample %in% subset_samples_cfpr_control, ]
subset_annotation_cfpr_control$Condition <- ifelse(subset_annotation_cfpr_control$Condition == "Control", "control", "case")

# CFPR vs PYR
subset_samples_cfpr_pyr <- annotation_table$Sample[annotation_table$Condition %in% c("CFPR", "PYR")]
subset_counts_cfpr_pyr <- count_matrix[, subset_samples_cfpr_pyr]
subset_annotation_cfpr_pyr <- annotation_table[annotation_table$Sample %in% subset_samples_cfpr_pyr, ]
subset_annotation_cfpr_pyr$Condition <- ifelse(subset_annotation_cfpr_pyr$Condition == "PYR", "control", "case")

# CFPR vs Combined
subset_samples_cfpr_combined <- annotation_table$Sample[annotation_table$Condition %in% c("CFPR", "Combined")]
subset_counts_cfpr_combined <- count_matrix[, subset_samples_cfpr_combined]
subset_annotation_cfpr_combined <- annotation_table[annotation_table$Sample %in% subset_samples_cfpr_combined, ]
subset_annotation_cfpr_combined$Condition <- ifelse(subset_annotation_cfpr_combined$Condition == "Combined", "control", "case")

# Control vs PYR
subset_samples_control_pyr <- annotation_table$Sample[annotation_table$Condition %in% c("Control", "PYR")]
subset_counts_control_pyr <- count_matrix[, subset_samples_control_pyr]
subset_annotation_control_pyr <- annotation_table[annotation_table$Sample %in% subset_samples_control_pyr, ]
subset_annotation_control_pyr$Condition <- ifelse(subset_annotation_control_pyr$Condition == "Control", "control", "case")

# Control vs Combined
subset_samples_control_combined <- annotation_table$Sample[annotation_table$Condition %in% c("Control", "Combined")]
subset_counts_control_combined <- count_matrix[, subset_samples_control_combined]
subset_annotation_control_combined <- annotation_table[annotation_table$Sample %in% subset_samples_control_combined, ]
subset_annotation_control_combined$Condition <- ifelse(subset_annotation_control_combined$Condition == "Control", "control", "case")

# PYR vs Combined
subset_samples_pyr_combined <- annotation_table$Sample[annotation_table$Condition %in% c("PYR", "Combined")]
subset_counts_pyr_combined <- count_matrix[, subset_samples_pyr_combined]
subset_annotation_pyr_combined <- annotation_table[annotation_table$Sample %in% subset_samples_pyr_combined, ]
subset_annotation_pyr_combined$Condition <- ifelse(subset_annotation_pyr_combined$Condition == "Combined", "control", "case")

# Run ImpulseDE2 manually for each comparison and print the results
impulse_results_cfpr_control <- runImpulseDE2(matCountData = subset_counts_cfpr_control,
                                               dfAnnotation = subset_annotation_cfpr_control,boolCaseCtrl = T, boolIdentifyTransients = T,scaNProc = 4)
library(dplyr)


# Subset the results to find genes significantly different between case and control
significant_genes <- impulse_results_cfpr_control$dfImpulseDE2Results %>%
  filter(padj < 0.05)  # Adjust the p-value threshold as needed

write.table(significant_genes,"/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/ImpulseDE/CFPR vs Con.txt",sep='\t',row.names=F)





impulse_results_cfpr_pyr <- runImpulseDE2(matCountData = subset_counts_cfpr_pyr,
                                           dfAnnotation = subset_annotation_cfpr_pyr,boolCaseCtrl = T, boolIdentifyTransients = T,scaNProc = 4)
significant_genes <- impulse_results_cfpr_pyr$dfImpulseDE2Results %>% filter(padj < 0.05)  # Adjust the p-value threshold as needed


write.table(significant_genes,"/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/ImpulseDE/CFPR vs PYR.txt",sep='\t',row.names=F)

impulse_results_cfpr_combined <- runImpulseDE2(matCountData = subset_counts_cfpr_combined,
                                                dfAnnotation = subset_annotation_cfpr_combined,boolCaseCtrl = T, boolIdentifyTransients = T,scaNProc = 4)
significant_genes <- impulse_results_cfpr_combined$dfImpulseDE2Results %>%
  filter(padj < 0.05)  # Adjust the p-value threshold as needed


write.table(significant_genes,"/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/ImpulseDE/CFPR vs IG2.txt",sep='\t',row.names=F)


impulse_results_control_pyr <- runImpulseDE2(matCountData = subset_counts_control_pyr,
                                               dfAnnotation = subset_annotation_control_pyr,boolCaseCtrl = T, boolIdentifyTransients = T,scaNProc = 4)
significant_genes <- impulse_results_control_pyr$dfImpulseDE2Results %>%
  filter(padj < 0.05)  # Adjust the p-value threshold as needed

write.table(significant_genes,"/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/ImpulseDE/PYR vs Con.txt",sep='\t',row.names=F)


impulse_results_control_combined <- runImpulseDE2(matCountData = subset_counts_control_combined,
                                                    dfAnnotation = subset_annotation_control_combined,boolCaseCtrl = T, boolIdentifyTransients = T,scaNProc = 4)
significant_genes <- impulse_results_control_combined$dfImpulseDE2Results %>%
  filter(padj < 0.05)  # Adjust the p-value threshold as needed

write.table(significant_genes,"/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/ImpulseDE/Combined vs Con.txt",sep='\t',row.names=F)


impulse_results_pyr_combined <- runImpulseDE2(matCountData = subset_counts_pyr_combined,
                                                dfAnnotation = subset_annotation_pyr_combined,boolCaseCtrl = T, boolIdentifyTransients = T,scaNProc = 4)
significant_genes <- impulse_results_pyr_combined$dfImpulseDE2Results %>%
  filter(padj < 0.05)  # Adjust the p-value threshold as needed


write.table(significant_genes,"/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/ImpulseDE/PYR vs combined.txt",sep='\t',row.names=F)
