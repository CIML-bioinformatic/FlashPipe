## ####################################################
## This script aims to create a plot for the UMI RNA Count, for each well and each plate.
## ####################################################

## @knitr compute_RNA_UMI

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Pre-processing of ERCC RNA data from zUMIs
# 2. Creation of general graph of ERCC UMI RNA data
# 3. Creation of plate-by-plate graph of ERCC RNA data for UMIs
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# ######################################################
# ## 1. Pre-processing of ERCC RNA data from zUMIs #####
# ######################################################

# Create an empty dataframe to store the UMI Count and feature count
all_summary_RNA_df_all_plate <- data.frame(stringsAsFactors = FALSE)

# Initialize the list that will store all the dataframes for each plate to be used for the figures.
rna_count_list = list()

# Iterate over each plate in the list
for (plate in levels(PLATES_LIST)) {
  plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate)
  
  # Read the RNA_Count file (required for UMI plot)
  RNA_Count = RNA_count_df_list[[ plate]]
  
  summary_RNA_df = data.frame( RNA.UMI_count = colSums(RNA_Count),
                               RNA.UMI_count_log = log10(colSums(RNA_Count)),
                               WellID = names(RNA_Count))
  summary_RNA_df[ , COLUMN_HEADER_PLATE_NAME] = plate
  
  
  # Add the UMI data to the global dataframe all_summary_RNA_df_all_plate
  all_summary_RNA_df_all_plate <- rbind(all_summary_RNA_df_all_plate, summary_RNA_df)
  
  # Save the data in a list, to use it after (for the report)
  rna_count_list[[plate]] = summary_RNA_df
}

# ##########################################################
# ## 2. Creation of general graph of ERCC UMI RNA data #####
# ##########################################################

# Order the plate names by the order provided by the user
all_summary_RNA_df_all_plate[ , COLUMN_HEADER_PLATE_NAME] = factor(all_summary_RNA_df_all_plate[ , COLUMN_HEADER_PLATE_NAME], levels = levels(PLATES_LIST))

# Dot plot for UMI_Count
dotplot_all_RNA_umi_count_log <- ggplot(all_summary_RNA_df_all_plate, aes(x = .data[[COLUMN_HEADER_PLATE_NAME]], y = RNA.UMI_count_log)) + 
  geom_violin(aes(fill = .data[[COLUMN_HEADER_PLATE_NAME]])) + 
  geom_jitter(width = 0.2) +
  stat_summary(fun = "mean", geom = "crossbar", colour = "red", linewidth = 0.5, width = 0.08, aes(linetype = "Mean", group = 1)) + 
  stat_summary(fun = "median", geom = "crossbar", colour = "orange", linewidth = 0.5, width = 0.08, aes(linetype = "Median", group = 1)) +
  scale_y_continuous() +  
  scale_linetype_manual(values = c("Mean" = "solid", "Median" = "solid")) + 
  labs(linetype = "Statistics") + 
  xlab("Plate name") + ylab(" Log RNA UMI count") +
  ggtitle("Violin Plot of the log of RNA UMI count for all plates") + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(dotplot_all_RNA_umi_count_log)

# ######################################################################
# ## 3. Creation of plate-by-plate graph of ERCC RNA data for UMIs #####
# ######################################################################

# Compute the minimum and maximum values of ERCC to set a common range on figures
rna_umi_min = min( all_summary_RNA_df_all_plate$RNA.UMI_count_log)
rna_umi_max = max( all_summary_RNA_df_all_plate$RNA.UMI_count_log)

# With the log, certain values can be -Inf or +Inf. Numerically this makes sense, but biologically it doesn't.
# Therefore, for the limit, if these values are present, they will be set directly to 0 or 1 to avoid errors.
if (rna_umi_min == -Inf) {
  rna_umi_min = 0
}
if (rna_umi_max == Inf) {
  rna_umi_max = 1
}


# Print the RNA UMI per plate in plate plots
for (plate_name in levels(PLATES_LIST)) {
  
  cat("\n \n")
  cat("\n#### ", plate_name, "\n\n")
  
  # Plot the UMI Count in log10
  plates_plot_feature_UMI_Count_log <- plate_plot(
    data = rna_count_list[[plate_name]],
    position = WellID,
    value = RNA.UMI_count_log,
    plate_size = 96,
    plate_type = "round",
    limits = c(rna_umi_min, rna_umi_max),
    title = paste("Plate", plate_name, " : Log of RNA UMI count")
  )
  plates_plot_feature_UMI_Count_log = plates_plot_feature_UMI_Count_log + labs(fill = "log RNA UMI count")
  print( plates_plot_feature_UMI_Count_log)
  
  # Save the plot to file
  plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name, '05_computeRNAUMI')
  # Create the output directory if it doesn't exist
  if (!dir.exists(plate_output_directory)) {
    dir.create(plate_output_directory, recursive = TRUE)
  }
  
  plot_filename_UMI_Count_log <- file.path(plate_output_directory, paste0("PlatePlot_UMI_Count_log_", plate_name, ".png"))
  ggsave(plot_filename_UMI_Count_log, plot = plates_plot_feature_UMI_Count_log)
} 
