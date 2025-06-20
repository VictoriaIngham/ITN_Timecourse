df = read.delim('/Users/vickyingham/Dropbox/IG1, IG2 and CFPR paper/Figures/OxPhos_in_50_percent.txt',header=T)

library(tidyverse)



## Data Preparation and Transformation
# 1. Pivot the data from wide to long format
#    Select columns that contain 'log2FoldChange' and pivot them.
df_long <- df %>%
  pivot_longer(
    cols = starts_with("log2FoldChange"),
    names_to = "Measurement",
    values_to = "log2FoldChange" # Keep this as log2FoldChange
  ) %>%
  # No need for RealFoldChange calculation here, as we are returning to log2FoldChange
  # 2. Extract Condition and Time from the 'Measurement' column
  #    This step uses regular expressions to parse the column names.
  mutate(
    Condition = str_extract(Measurement, "CFPR|PYR|Combined"),
    Time = as.numeric(str_extract(Measurement, "\\d+(?=hr|$|$)"))
  ) %>%
  # 3. Clean up the 'Condition' column values (remove '_vs_Control' and ensure correct names)
  #    Also, replace any NA in 'Condition' if the pattern extraction failed for some reason.
  mutate(
    Condition = case_when(
      grepl("CFPR", Measurement) ~ "CFPR",
      grepl("PYR", Measurement) ~ "PYR",
      grepl("Combined", Measurement) ~ "Combined",
      TRUE ~ NA_character_
    )
  ) %>%
  # 4. Ensure 'Time' is numeric and 'Complex' and 'Condition' are factors for plotting
  mutate(
    Time = as.numeric(Time),
    Complex = as.factor(Complex),
    Condition = as.factor(Condition)
  ) %>%
  # 5. Remove rows where 'Condition' or 'Time' might be NA if extraction failed
  filter(!is.na(Condition) & !is.na(Time))

# Define the custom colors for each condition
colors <- c("CFPR" = "#E69F00", "PYR" = "#009E73", "Combined" = "#CC79A7")

# ---
## Plotting Each Complex on a Separate Graph
# Get unique complex names
unique_complexes <- unique(df_long$Complex)

# Loop through each complex and create a separate plot
for (comp in unique_complexes) {
  # Filter data for the current complex
  df_complex <- df_long %>%
    filter(Complex == comp)
  
  # Create the plot for the current complex
  p <- ggplot(df_complex, aes(x = Time, y = log2FoldChange, color = Condition)) +
    geom_smooth(aes(fill = Condition), se = TRUE, method = "loess", alpha = 0.2) +
    scale_color_manual(values = colors) + # Apply custom colors for lines
    scale_fill_manual(values = colors) +   # Apply custom colors for fill (confidence intervals)
    labs(
      title = paste0("Smoothed log2 Fold Change Over Time for ", comp), # Dynamic title
      x = "Time (hrs)",
      y = "log2 Fold Change",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_minimal() + # Use a minimal theme for a clean look
    theme(
      plot.title = element_text(hjust = 0.5), # Center plot title
      legend.position = "bottom", # Position legend at the bottom
      # Increase size of axis tick labels
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  # Print the plot (this will display it in your R environment, or save it if directed)
  print(p)
  
  # Optional: Save each plot to a file
  # Uncomment and modify the lines below to save plots
   ggsave(filename = paste0("/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/log2FoldChange_plot_", gsub(" ", "_", comp), ".png"),
          plot = p, width = 8, height = 6, units = "in", dpi = 300)
}