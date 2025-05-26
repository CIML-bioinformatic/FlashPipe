## ####################################################
## This script aims to create a plot for the Feature Number for the RNA UMI, for each well and each plate.
## ####################################################

## @knitr compute_Feature_Number

# TODO CHANGE COMMENT : Iterate over each plate in the list
all_RNA_feature_count <- data.frame( stringsAsFactors = FALSE)
all_features_set = c()
RNA_feature_count_list = list()
for (plate_name in levels(PLATES_LIST)) {
  plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name)
  
  # Read the RNA_Feature_Count file (required for feature count plot)
  RNA_count_df = RNA_count_df_list[[ plate_name]]
  
  # Get the number of features and store it in a new column.
  RNA_feature_count = data.frame( WellID = names( RNA_count_df), feature_count = rep( 0, ncol(RNA_count_df)))
  RNA_feature_count$feature_count <- apply( RNA_count_df, 2, function(x) sum(x > 0))
  RNA_feature_count_list[[ plate_name]] = RNA_feature_count
  
  # Accumulate the genes expresed in all the plates  
  all_features_set = unique( append( all_features_set, row.names( RNA_count_df)))
  
  # Add the feature count data to the global dataframe all_UMI_Count
  all_RNA_feature_count <- rbind(all_RNA_feature_count, data.frame(Plate_Name = rep(plate_name, nrow(RNA_feature_count)),
                                                                  WellID = RNA_feature_count$WellID,    
                                                                  feature_count = RNA_feature_count$feature_count,
                                                                  plate_feature_count = rep( nrow(RNA_count_df), nrow(RNA_feature_count)))
                                  )
  
}

## ....................................................
## dotplot for all the plate with feature count
## ....................................................

all_RNA_feature_count[, PLATE_NAME] = factor(all_RNA_feature_count[,PLATE_NAME], levels = levels(PLATES_LIST))

# Dot plot for feature count
plot_all_RNA_feature_count <- ggplot(all_RNA_feature_count, aes(x = .data[[PLATE_NAME]], y = feature_count)) + 
  geom_violin(aes(fill = .data[[PLATE_NAME]])) + 
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

## ....................................................
## upset plot for the feature of all the plate
## ....................................................

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
                    sets = colnames(df_Up_SetPlot),  # Les plaques comme ensembles
                    sets.bar.color = "#56B4E9", # Couleur des barres
                    order.by = "freq", # Ordre par fréquence des intersections
                    main.bar.color = "#009E73", # Couleur des barres principales
                    matrix.color = "#D55E00", # Couleur de la matrice de l'upset
                    keep.order = TRUE) # Garder l'ordre des plaques

# Save the graph in file
# upset_plot_filename <- file.path(PATH_ANALYSIS_OUTPUT, "UpSet_Plot_Gene_Presence.png")
# ggsave(upset_plot_filename, plot = upset_plot)
print(upset_plot)

cat("\n \n")

## ....................................................
## Plotplate for each plate 
## ....................................................

# Compute the minimum and maximum values of ERCC to set a common range on figures
rna_feature_count_min = min( all_RNA_feature_count$feature_count)
rna_feature_count_max = max( all_RNA_feature_count$feature_count)

# Print the ERCC UMI per plate in plate plots
for (plate_name in levels(PLATES_LIST)) {
  
  cat("\n \n")
  cat("\n#### ", plate_name, "\n\n")
  
  # Plot the feature count
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
  plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name)
  path_output_plot_feature_count <- file.path(plate_output_directory, paste0("PlatePlot_feature_Count_", plate_name, ".png"))
  ggsave(path_output_plot_feature_count, plot = plot_plate_feature_Count)
} 
