# Load required packages
library(dplyr)
library(readr)
library(stringr)


# For a .csv file:
df <- read_csv('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Analysis of all data/final_combined_results_fixed.csv')

# For a tab-delimited file (e.g. from Excel or RNAseq output):
# df <- read_tsv("your_file.txt")

# --- Step 2: Define treatments and time points ---
treatments <- c("CFPR", "PYR", "Combined")
timepoints <- c("1", "4", "8", "24", "48", "72")

# --- Step 3: Count significant padj values (< 0.05) ---
count_significant <- function(df, treatment) {
  # Build column names like: padj.CFPR_vs_Control_1, etc.
  padj_cols <- paste0("padj.", treatment, "_vs_Control_", timepoints)
  
  # Safely subset these columns (check they exist)
  padj_cols <- padj_cols[padj_cols %in% colnames(df)]
  
  # Count how many time points have padj < 0.05 per gene
  sig_counts <- apply(df[, padj_cols], 1, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE))
  return(sig_counts)
}

# --- Step 4: Apply to all treatments ---
df$CFPR_sig <- count_significant(df, "CFPR")
df$PYR_sig <- count_significant(df, "PYR")
df$Combined_sig <- count_significant(df, "Combined")

# --- Step 5: Filter genes with ≥2 significant timepoints in any treatment ---
sig_genes <- df %>%
  filter(CFPR_sig >= 4 | PYR_sig >= 4 | Combined_sig >= 4)

# --- Optional: View or save the result ---
print(paste("Number of genes meeting criteria:", nrow(sig_genes)))
# View(sig_genes)  # Uncomment in RStudio
write_csv(sig_genes, "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Gene lists of interest/significant_genes_in_66_percent_in_one_treatment.csv")  # Save result


###TFs in 25%

TFs = read.delim('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Gene lists of interest/Unique_TFs.txt',header=T)

merged = merge(sig_genes,TFs, by.x='ID',by.y='Gene.ID')
write_csv(merged, "/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Gene lists of interest/TFs_in_25_percent.csv")  # Save result
