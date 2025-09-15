## ####################################################
## This script aims to create a plot for the Feature 
## Number for the RNA UMI, for each well and each plate.
## ####################################################

## @knitr compute_Feature_Number

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Pre-processing and feature number calculation of ERCC RNA data from zUMIs
# 2. Creation of general graph of ERCC UMI RNA data
# 3. Creation of Upset Plot of ERCC UMI RNA data
# 4. Creation of plate-by-plate graph of ERCC RNA data for UMIs
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# #####################################################################################
# ## 1. Pre-processing and feature number calculation of ERCC RNA data from zUMIs #####
# #####################################################################################

### The aim is to analyze all the features present in the data
### to produce 3 figures, a general graph for all plates
### a figure for each of the plates, and an Upset Plot (which groups the
### information linked to the genes present between the plates).

# Initialize variables used to store data (especially for figures)
all_RNA_feature_count <- data.frame( stringsAsFactors = FALSE)
all_features_set = c()
RNA_feature_count_list = list()

# For loop to detect the total number of features for each plate (then stock them in the variable for graphics)
for (plate_name in levels(PLATES_LIST)) {
  plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name)
  
  # Read the RNA_Feature_Count file (required for feature count plot)
  RNA_count_df = RNA_count_df_list[[ plate_name]]
  
  # Get the number of features and store it in a new column.
  RNA_feature_count = data.frame( WellID = names( RNA_count_df), feature_count = rep( 0, ncol(RNA_count_df)))
  RNA_feature_count$feature_count <- apply( RNA_count_df, 2, function(x) sum(x > 0))
  RNA_feature_count_list[[ plate_name]] = RNA_feature_count
  
  # Accumulate the genes expressed in all the plates  
  all_features_set = unique( append( all_features_set, row.names( RNA_count_df)))
  
  # Add the feature count data to the global dataframe all_UMI_Count
  all_RNA_feature_count <- rbind(all_RNA_feature_count, data.frame(Plate_Name = rep(plate_name, nrow(RNA_feature_count)),
                                                                  WellID = RNA_feature_count$WellID,    
                                                                  feature_count = RNA_feature_count$feature_count,
                                                                  plate_feature_count = rep( nrow(RNA_count_df), nrow(RNA_feature_count)))
                                )
  
}

# ##########################################################
# ## 2. Creation of general graph of ERCC UMI RNA data #####
# ##########################################################

# Separates in the dataframe all features_count by levels, i.e. each plate name
all_RNA_feature_count[, COLUMN_HEADER_PLATE_NAME] = factor(all_RNA_feature_count[,COLUMN_HEADER_PLATE_NAME], levels = levels(PLATES_LIST))

# Dot plot for feature count
plot_all_RNA_feature_count <- ggplot(all_RNA_feature_count, aes(x = .data[[COLUMN_HEADER_PLATE_NAME]], y = feature_count)) + 
  geom_violin(aes(fill = .data[[COLUMN_HEADER_PLATE_NAME]])) + 
  geom_jitter(width = 0.2) +
  stat_summary(fun = "mean", geom = "crossbar", colour = "red", linewidth = 0.5, width = 0.08, aes(linetype = "Mean", group = 1)) + 
  stat_summary(fun = "median", geom = "crossbar", colour = "orange", linewidth = 0.5, width = 0.08, aes(linetype = "Median", group = 1)) +
  scale_y_continuous() +  
  scale_linetype_manual(values = c("Mean" = "solid", "Median" = "solid")) + 
  labs(linetype = "Statistics") + 
  xlab("Plate name") + ylab("RNA feature count") +
  ggtitle("Violin Plot of RNA feature count for all plates") + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot_all_RNA_feature_count)

# #######################################################
# ## 3. Creation of Upset Plot of UMI RNA data ##########
# #######################################################

if( length( PLATES_LIST) > 1){
  cat("\n \n")
  cat("<b>Number total of gene through all the plates : ", length( all_features_set), "</b>.\n \n", sep = '')
  
  # Get a list of all the gene names across all my plates.
  all_genes <- unique(unlist(lapply(levels(PLATES_LIST), function(plate) {
    rownames(RNA_count_df_list[[plate]])
  })))
  
  # Create a dataframe containing the names of all genes with the names of all plates in columns.
  df_Up_SetPlot <- data.frame(matrix(0, nrow = length(all_genes), ncol = length(PLATES_LIST)))
  rownames(df_Up_SetPlot) <- all_genes
  colnames(df_Up_SetPlot) <- levels(PLATES_LIST)
  
  # We replace the 0's in the dataframe with a “1” if the genes are present in the plate.
  for (plate in levels(PLATES_LIST)) {
    RNA_Feature_Count <- RNA_count_df_list[[plate]]
    gene_names <- rownames(RNA_Feature_Count)
    df_Up_SetPlot[gene_names, plate] <- 1
  }
  
  upset_plot <- upset(df_Up_SetPlot, 
                      sets = colnames(df_Up_SetPlot),  # Plates as sets
                      sets.bar.color = "#56B4E9", 
                      order.by = "freq", 
                      main.bar.color = "#009E73", 
                      matrix.color = "#D55E00", 
                      keep.order = TRUE) # Keep the order of the plates
  print(upset_plot)
  
  cat("\n \n")
}

# ######################################################################
# ## 4. Creation of plate-by-plate graph of ERCC RNA data for UMIs #####
# ######################################################################

# Compute the minimum and maximum values of ERCC to set a common range on figures
rna_feature_count_min = min( all_RNA_feature_count$feature_count)
rna_feature_count_max = max( all_RNA_feature_count$feature_count)

# Print the ERCC UMI per plate in plate plots
for (plate_name in levels(PLATES_LIST)) {
  cat("\n \n")
  cat("\n#### ", plate_name, "\n\n")
  plot_plate_feature_Count <- plate_plot(
    data = RNA_feature_count_list[[ plate_name]],
    position = WellID,
    value = feature_count,
    plate_size = 96,
    plate_type = "round",
    limits = c(rna_feature_count_min, rna_feature_count_max),
    title = paste("Plate", plate_name, " : RNA feature count")
  )
  plot_plate_feature_Count = plot_plate_feature_Count + labs(fill = "RNA feature count")
  print( plot_plate_feature_Count)
  
  # Save the plot to file
  plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name, '06_computeFeatureNumber')
  # Create the output directory if it doesn't exist
  if (!dir.exists(plate_output_directory)) {
    dir.create(plate_output_directory, recursive = TRUE)
  }
  path_output_plot_feature_count <- file.path(plate_output_directory, paste0("PlatePlot_feature_Count_", plate_name, ".png"))
  ggsave(path_output_plot_feature_count, plot = plot_plate_feature_Count)
} 
